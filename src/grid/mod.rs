use std::{io::{BufReader, Read}, marker::PhantomData, os::raw::c_void};

use mpi::{ffi::{MPI_Recv, MPI_Send}, raw::AsRaw, topology::SimpleCommunicator, traits::{Communicator, Equivalence}};

use crate::basetree::BaseTree;


pub mod refine;
pub mod transfer;
pub mod iterate;
pub mod corners;
pub mod cell;
pub mod nodes;


#[derive(Clone, Debug, Default)]
pub(crate) struct CellData<T> {
    pub data: T,
    pub local_id: u32,
    pub global_id: u32,
    pub owner_rank: u32,
}


pub struct Grid<T> {

    communicator: SimpleCommunicator,

    base_tree: BaseTree,

    connectivity: *mut p4est_sys::p4est_connectivity,

    grid: *mut p4est_sys::p4est,
    
    ghosts: *mut p4est_sys::p4est_ghost,

    ghost_data: Vec<CellData<T>>,
    
    pt: PhantomData<T>,

    global_id_offset: usize,
}


impl<T> Grid<T> {

    pub fn new_unitsquare(communicator: SimpleCommunicator) -> Self where T: Sized {

        let comm = communicator.as_raw();

        let tree = BaseTree::new_unitsquare();

        let connectivity = unsafe { tree.build_connectivity() };

        let s = size_of::<CellData<T>>();

        let grid = unsafe {
            p4est_sys::p4est_new(comm.0 as *mut c_void as p4est_sys::MPI_Comm, connectivity, s, None, std::ptr::null_mut())
        };

        Self {
            base_tree: tree,
            communicator,
            connectivity,
            grid,
            ghosts: std::ptr::null_mut(),
            ghost_data: vec![],
            pt: PhantomData,
            global_id_offset: 0,
        }.with_updated_ids()
    }


    pub fn from_su2<F>(f: F, communicator: SimpleCommunicator) -> Result<Self, std::io::Error> where F: Read {

        let reader = BufReader::new(f);
        let tree = BaseTree::from_su2(reader)?;

        let connectivity = unsafe { tree.build_connectivity() };

        let s = size_of::<CellData<T>>();

        let comm = communicator.as_raw();

        let grid = unsafe {
            p4est_sys::p4est_new(comm.0 as *mut c_void as p4est_sys::MPI_Comm, connectivity, s, None, std::ptr::null_mut())
        };

        Ok(Self {
            base_tree: tree,
            communicator,
            connectivity,
            grid,
            ghosts: std::ptr::null_mut(),
            ghost_data: vec![],
            pt: PhantomData,
            global_id_offset: 0,
        }.with_updated_ids())
    }

    pub fn local_len(&self) -> usize {
        unsafe { (*self.grid).local_num_quadrants as usize }
    }
    pub fn global_len(&self) -> usize {
        unsafe { (*self.grid).global_num_quadrants as usize }
    }

    pub fn len_with_ghosts(&self) -> usize {
        self.local_len() + self.ghost_data.len()
    }

    pub fn partition(&mut self) where T: Clone + Default {
        if self.communicator.size() <= 1 {
            return;
        }

        unsafe {
            p4est_sys::p4est_partition(self.grid, 1, None);
        }
        self.update_ids();
        self.update_ghosts();
    }

    pub fn update_ghosts(&mut self) where T: Clone + Default {
        if self.communicator.size() <= 1 {
            return;
        }

        if self.ghosts != std::ptr::null_mut() {
            unsafe {
                p4est_sys::p4est_ghost_destroy(self.ghosts);
                self.ghosts = std::ptr::null_mut();
            }
        }

        unsafe {
            self.ghosts = p4est_sys::p4est_ghost_new(self.grid, p4est_sys::p4est_connect_type_t_P4EST_CONNECT_FACE);
        
            self.ghost_data.resize((*self.ghosts).ghosts.elem_count, CellData::<T>::default());
            
            p4est_sys::p4est_ghost_exchange_data(self.grid, self.ghosts, self.ghost_data.as_mut_ptr() as *mut c_void);
        }

        self.update_ghost_ids();

    }

    pub fn exchange_ghost_data(&mut self) {
        if self.ghosts != std::ptr::null_mut() {
            unsafe {
                p4est_sys::p4est_ghost_exchange_data(self.grid, self.ghosts, self.ghost_data.as_mut_ptr() as *mut c_void);
            }
        }
        self.update_ghost_ids();
    }

    pub(crate) fn update_ghost_ids(&mut self) {
        let local_len = self.local_len();

        for i in 0..self.ghost_data.len() {
            self.ghost_data[i].local_id = ( local_len + i ) as u32;
        }
    }

    pub fn boundary_len(&self) -> usize {
        self.base_tree.boundary_len()
    }
    pub fn boundary_name(&self, boundary: u16) -> &str {
        self.base_tree.boundary_name(boundary)
    }

}



impl<T> Drop for Grid<T> {
    fn drop(&mut self) {
        unsafe {
            p4est_sys::p4est_destroy(self.grid);
        }
        unsafe {
            p4est_sys::p4est_connectivity_destroy(self.connectivity);
        }
        unsafe {
            if self.ghosts != std::ptr::null_mut() {
                p4est_sys::p4est_ghost_destroy(self.ghosts);
            }
        }
    }
}




impl<T> Grid<T> {


    pub(crate) fn with_updated_ids(mut self) -> Self {
        self.update_ids();
        self
    }

    pub(crate) fn update_ids(&mut self) {

        let rank = self.communicator.rank();
        let size = self.communicator.size();

        let raw_comm = self.communicator.as_raw();

        let local_len = self.local_len();

        let mut self_offset: usize = 0;
        let mut user_data = [0_usize; 2];
        user_data[0] = self_offset;
        user_data[1] = rank as usize;
        
        if rank > 0 {
            unsafe {
                MPI_Recv(&mut self_offset as *mut usize as *mut c_void, 1, usize::equivalent_datatype().as_raw(), rank - 1, 0, raw_comm, std::ptr::null_mut());
            }
        }

        let next_offset = self_offset + local_len;

        if (rank + 1) < size {
            unsafe {
                MPI_Send(&next_offset as *const usize as *mut c_void, 1, usize::equivalent_datatype().as_raw(), rank + 1, 0, raw_comm);
            }
        }

        // we have the local offset for global id evaluation
        self.global_id_offset = self_offset;
        user_data[0] = self.global_id_offset;


        unsafe {
            #[cfg(feature = "2d")]
            p4est_sys::p4est_iterate(
                self.grid,
                self.ghosts,
                user_data.as_mut_ptr() as *mut c_void,
                Some(iter_volume_update_ids_fn::<T>),
                None,
                None,
            );

            #[cfg(feature = "3d")]
            p4est_sys::p4est_iterate(
                self.grid,
                self.ghosts,
                user_data.as_mut_ptr() as *mut c_void,
                Some(iter_volume_update_ids_fn::<T>),
                None,
                None,
                None,
            );
        }

    }

}






extern "C" fn iter_volume_update_ids_fn<T>(info: *mut p4est_sys::p4est_iter_volume_info, data: *mut std::ffi::c_void) {

    unsafe {
        

        let global_offset = *(data as *const usize);
        let rank = *(data as *const usize).offset(1);

        let grid = (*info).p4est;

        let tarrayelemsize = (*(*grid).trees).elem_size as i32;
        let tree = (*(*grid).trees).array.offset((tarrayelemsize * (*info).treeid) as isize) as *mut p4est_sys::p4est_tree;
        let local_offset = (*tree).quadrants_offset;

        let local_id = (*info).quadid + local_offset;

        let cell_data: *mut CellData<T> = (*(*info).quad).p.user_data as *mut c_void as *mut CellData<T>;
        let cell_data = &mut (*cell_data);

        cell_data.local_id = local_id as u32;
        cell_data.global_id = (local_id as u32) + (global_offset as u32);
        cell_data.owner_rank = rank as u32;

    }

}




impl<T> Clone for Grid<T> where T: Clone {
    fn clone(&self) -> Self {
        let p4est = unsafe {p4est_sys::p4est_copy(self.grid, 1)};
        let ghosts = if self.ghosts == std::ptr::null_mut() {
            std::ptr::null_mut()
        } else {
            unsafe {p4est_sys::p4est_ghost_new(p4est, p4est_sys::p4est_connect_type_t_P4EST_CONNECT_FACE)}
        };
        let ghost_data = self.ghost_data.clone();
        Self {
            base_tree: self.base_tree.clone(),
            communicator: self.communicator.duplicate(),
            connectivity: unsafe {p4est_sys::p4est_connectivity_copy(self.connectivity, 1)},
            grid: p4est,
            ghosts,
            ghost_data,
            global_id_offset: self.global_id_offset,
            pt: PhantomData,
        }
    }
}


