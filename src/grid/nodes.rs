use p4est_sys::consts::{CELL_CORNERS, FACE_CORNERS};

use crate::grid::{Grid, cell::{Cell, Face}, corners::face_corner_ids};



pub struct NodeNumbering {
    nodes: *mut p4est_sys::p4est_nodes,
}

impl NodeNumbering {

    pub fn new<T>(grid: &Grid<T>) -> Self {
        Self {
            nodes: unsafe {p4est_sys::p4est_nodes_new(grid.grid, grid.ghosts)}
        }
    }


    pub fn cell_nodes<'a, T>(&self, cell: &Cell<'a, T>) -> [usize; CELL_CORNERS] {
        let cell_id = cell.local_id;

        let mut out = [0; CELL_CORNERS];

        for i in 0..CELL_CORNERS {
            let nidx = unsafe { *(*self.nodes).local_nodes.offset((CELL_CORNERS * cell_id + i) as isize) };
            out[i] = nidx as usize;
        }

        out
    }

    pub fn face_nodes<'a, T>(&self, face: &Face<'a, T>) -> [usize; FACE_CORNERS] {
        let face_id = face.face_id;
        let cell = if let Some(c1) = &face.cell1 {
            if c1.level > face.cell0.level {
                c1
            } else {
                &face.cell0
            }
        } else {
            &face.cell0
        };
        let fc = face_corner_ids(face_id);
        let mut out = [0; FACE_CORNERS];
        let cnodes = self.cell_nodes(cell);
        for i in 0..FACE_CORNERS {
            out[i] = cnodes[fc[i] as usize];
        }
        out
    }

}

impl Drop for NodeNumbering {
    fn drop(&mut self) {
        unsafe {
            p4est_sys::p4est_nodes_destroy(self.nodes);
        }
    }
}


