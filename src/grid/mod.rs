use std::{io::{BufReader, Read}, marker::PhantomData, os::raw::c_void};

use mpi::{ffi::{MPI_Recv, MPI_Send}, raw::AsRaw, topology::SimpleCommunicator, traits::{Communicator, Equivalence}};

use crate::basetree::BaseTree;


pub mod refine;
pub mod iterate;
pub mod corners;
pub mod cell;


#[derive(Clone, Debug)]
pub(crate) struct CellData<T> {
    pub data: T,
    pub local_id: u32,
    pub global_id: u32,
}


pub struct Grid<T> {

    communicator: SimpleCommunicator,

    base_tree: BaseTree,

    connectivity: *mut p4est_sys::p4est_connectivity,

    grid: *mut p4est_sys::p4est,
    
    ghosts: *mut p4est_sys::p4est_ghost,
    
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

    pub fn partition(&mut self) {
        unsafe {
            p4est_sys::p4est_partition(self.grid, 1, None);
        }
        self.update_ids();
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


        unsafe {
            #[cfg(feature = "2d")]
            p4est_sys::p4est_iterate(
                self.grid,
                self.ghosts,
                (&mut self_offset as *mut usize) as *mut c_void,
                Some(iter_volume_update_ids_fn::<T>),
                None,
                None,
            );

            #[cfg(feature = "3d")]
            p4est_sys::p4est_iterate(
                self.grid,
                self.ghosts,
                (&mut self_offset as *mut usize) as *mut c_void,
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

        let grid = (*info).p4est;

        let tarrayelemsize = (*(*grid).trees).elem_size as i32;
        let tree = (*(*grid).trees).array.offset((tarrayelemsize * (*info).treeid) as isize) as *mut p4est_sys::p4est_tree;
        let local_offset = (*tree).quadrants_offset;

        let local_id = (*info).quadid + local_offset;

        let cell_data: *mut CellData<T> = (*(*info).quad).p.user_data as *mut c_void as *mut CellData<T>;
        let cell_data = &mut (*cell_data);

        cell_data.local_id = local_id as u32;
        cell_data.global_id = (local_id as u32) + (global_offset as u32);

    }

}

