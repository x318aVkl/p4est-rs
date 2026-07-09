use std::ffi::c_void;

use crate::{basetree::BaseTree, grid::{CellData, Grid, corners::cell_corners}};


use crate::grid::cell::Cell;




impl<T> Grid<T> {
    pub fn map_cells<'a, F>(&'a self, f: F) where F: FnMut(Cell<'a, T>) {
        unsafe {
            USER_VOLUME_FN = Some(&f as *const F as *mut c_void);

            let user_data = &self.base_tree as *const BaseTree as *mut c_void;

            #[cfg(feature = "2d")]
            p4est_sys::p4est_iterate(
                self.grid,
                self.ghosts,
                user_data,
                Some(iter_volume_fn::<F, T>),
                None,
                None,
            );

            #[cfg(feature = "3d")]
            p4est_sys::p4est_iterate(
                self.grid,
                self.ghosts,
                user_data,
                Some(iter_volume_fn::<F, T>),
                None,
                None,
                None,
            );

            USER_VOLUME_FN = None;
        }
    }
}



static mut USER_VOLUME_FN: Option<*mut std::ffi::c_void> = None;

extern "C" fn iter_volume_fn<'a, F, T>(info: *mut p4est_sys::p4est_iter_volume_info, data: *mut std::ffi::c_void) where F: FnMut(Cell<'a, T>), T: 'a + Sized {

    unsafe {

        let cell_data: *const CellData<T> = (*(*info).quad).p.user_data as *mut c_void as *const CellData<T>;
        let cell_data = &(*cell_data);

        let base_tree = &*(data as *const BaseTree);

        let id = cell_data.local_id;
        let gid = cell_data.global_id;
        let corners = cell_corners(base_tree, (*info).treeid, (*info).quad);
        
        let cell = Cell {
            data: &cell_data.data, 
            local_id: id as usize, 
            global_id: gid as usize,
            level: (*(*info).quad).level as u8,
            corners,
        };

        let f = USER_VOLUME_FN.unwrap() as *mut F;

        (*f)(cell);
    }

}

