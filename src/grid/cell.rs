
use std::os::raw::c_void;

use p4est_sys::consts::{CELL_CORNERS, DIM, FACE_CORNERS};

use crate::{basetree::BaseTree, grid::{CellData, corners::cell_corners}};



#[derive(Debug)]
pub struct Cell<'a, T> { 
    pub data: &'a T,
    pub local_id: usize,
    pub global_id: usize,
    pub level: u8,
    pub is_ghost: bool,
    pub owner_rank: u32,
    pub(crate) corners: [[f64; DIM]; CELL_CORNERS],
    pub(crate) corners_int: [[u32; DIM]; CELL_CORNERS],
    pub(crate) raw_quad: *const p4est_sys::p4est_quadrant,
    pub(crate) tree_id: i32,
}



impl<'a, T> Cell<'a, T> {

    pub fn corner(&self, i: usize) -> [f64; DIM] {
        self.corners[i]
    }

    pub fn corner_len(&self) -> usize {
        self.corners.len()
    }

    pub(crate) unsafe fn from_quad(
        base_tree: &BaseTree,
        treeid: i32,
        quad: *const p4est_sys::p4est_quadrant,
    ) -> Self {
        unsafe {
            
            let cell_data: *const CellData<T> = (*quad).p.user_data as *mut c_void as *const CellData<T>;
            let cell_data: &CellData<T> = &(*cell_data);

            let id = cell_data.local_id;
            let gid = cell_data.global_id;
            let (corners, corners_int) = cell_corners(base_tree, treeid, quad as *mut p4est_sys::p4est_quadrant);
            
            Self {
                data: &cell_data.data, 
                local_id: id as usize, 
                global_id: gid as usize,
                level: (*quad).level as u8,
                is_ghost: false,
                owner_rank: cell_data.owner_rank,
                corners,
                corners_int,
                raw_quad: quad,
                tree_id: treeid,
            }
        }
    }
}



#[derive(Debug)]
pub struct Face<'a, T> {
    pub cell0: Cell<'a, T>,
    pub cell1: Option<Cell<'a, T>>,
    pub id: usize,
    pub(crate) face_id: i8,
    pub(crate) face_id_side: u8,
    pub(crate) corners: [[f64; DIM]; FACE_CORNERS],
    pub(crate) corners_int: [[u32; DIM]; FACE_CORNERS],
    pub(crate) subface_id: Option<u8>,
    pub boundary: Option<u16>,
}


impl<'a, T> Face<'a, T> {
    pub fn corner(&self, i: usize) -> [f64; DIM] {
        self.corners[i]
    }
}


