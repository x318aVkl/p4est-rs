
use std::ffi::c_void;

use p4est_sys::consts::CELL_CORNERS;

use crate::{basetree::BaseTree, grid::{CellData, Grid, corners::cell_corners}};

use crate::grid::cell::Cell;



extern "C" fn refine_uniform_fn(_grid: *mut p4est_sys::p4est, _treeid: i32, _quad: *mut p4est_sys::p4est_quadrant) -> i32 {
    1
}


static mut USER_REFINE_FN: Option<*mut std::ffi::c_void> = None;
static mut USER_COARSEN_FN: Option<*mut std::ffi::c_void> = None;
static mut USER_BASETREE: Option<*const BaseTree> = None;

extern "C" fn refine_fn<'a, F, T>(_grid: *mut p4est_sys::p4est, treeid: i32, quad: *mut p4est_sys::p4est_quadrant) -> i32 where F: FnMut(Cell<'a, T>) -> bool, T: 'a {
    unsafe {

        let cell_data: &CellData<T> = &*((*quad).p.user_data as *mut c_void as *const CellData<T>);

        let id = cell_data.local_id;
        let tree = &*USER_BASETREE.unwrap();
        let (corners, corners_int) = cell_corners(tree, treeid, quad);
        
        let cell = Cell {
            data: &cell_data.data, 
            local_id: id as usize, 
            global_id: id as usize,
            level: (*quad).level as u8,
            is_ghost: false,
            owner_rank: cell_data.owner_rank,
            corners,
            corners_int,
            raw_quad: quad,
            tree_id: treeid,
        };

        let f = USER_REFINE_FN.unwrap() as *mut F;

        if (*f)(cell) {1} else {0}
    }
}


extern "C" fn coarsen_fn<'a, F, T>(_grid: *mut p4est_sys::p4est, treeid: i32, quads: *mut *mut p4est_sys::p4est_quadrant) -> i32 where F: FnMut([Cell<'a, T>; CELL_CORNERS]) -> bool, T: 'a {
    unsafe {

        let cells: [Cell<'a, T>; CELL_CORNERS] = std::array::from_fn(|i| {

            let quad = *quads.offset(i as isize); 
            
            let cell_data: &CellData<T> = &*((*quad).p.user_data as *mut c_void as *const CellData<T>);

            let id = cell_data.local_id;
            let tree = &*USER_BASETREE.unwrap();
            let (corners, corners_int) = cell_corners(tree, treeid, quad);
                
            let cell = Cell {
                data: &cell_data.data, 
                local_id: id as usize, 
                global_id: id as usize,
                level: (*quad).level as u8,
                is_ghost: false,
                owner_rank: cell_data.owner_rank,
                corners,
                corners_int,
                raw_quad: quad,
                tree_id: treeid,
            };

            cell
        });

        let f = USER_COARSEN_FN.unwrap() as *mut F;

        if (*f)(cells) {1} else {0}
    }
}




impl<T> Grid<T> {

    pub fn refine_uniform(&mut self) where T: Clone + Default {
        unsafe {
            p4est_sys::p4est_refine(self.grid, 0, Some(refine_uniform_fn), None);
        }
        self.update_ids();
        self.update_ghosts();
    }

    pub fn refine<'a, F>(&'a mut self, f: F) where F: Fn(Cell<'a, T>) -> bool,T: Clone + Default {
        unsafe {
            USER_REFINE_FN = Some(&f as *const F as *mut c_void);
            USER_BASETREE = Some(&self.base_tree as *const BaseTree);

            p4est_sys::p4est_refine(
                self.grid,
                0,
                Some(refine_fn::<'a, F, T>),
                None,
            );

            USER_REFINE_FN = None;
            USER_BASETREE = None;
        }
        self.update_ids();
        self.update_ghosts();
    }

    pub fn coarsen<'a, F>(&'a mut self, f: F) where F: Fn([Cell<'a, T>; CELL_CORNERS]) -> bool, T: Clone + Default {
        unsafe {
            USER_COARSEN_FN = Some(&f as *const F as *mut c_void);
            USER_BASETREE = Some(&self.base_tree as *const BaseTree);

            p4est_sys::p4est_coarsen(
                self.grid,
                0,
                Some(coarsen_fn::<'a, F, T>),
                None,
            );

            USER_COARSEN_FN = None;
            USER_BASETREE = None;
        }
        self.update_ids();
        self.update_ghosts();
    }

    pub fn balance(&mut self) where T: Clone + Default {
        unsafe {
            p4est_sys::p4est_balance(self.grid, p4est_sys::p4est_connect_type_t_P4EST_CONNECT_FACE, None);
        }
        self.update_ids();
        self.update_ghosts();
    }

}



