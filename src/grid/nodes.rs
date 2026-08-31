use std::collections::{HashMap, HashSet};

use p4est_sys::consts::{CELL_CORNERS, DIM, FACE_CORNERS};

use crate::grid::{Grid, cell::{Cell, Face}, corners::face_corner_ids};



pub struct NodeNumbering {

    // computes a unique id 
    node_map: HashMap<[u32; 4], u32>,
    node_positions: Vec<[f64; DIM]>,

    cell_nodes: Vec<usize>,
    cell_nodes_starts: Vec<usize>,

}

impl NodeNumbering {

    pub fn new<T>(grid: &Grid<T>) -> Self {
        Self {
            node_map: HashMap::new(),
            node_positions: vec![],
            cell_nodes: vec![],
            cell_nodes_starts: vec![0],
        }.compute(&grid)
    }

    pub fn cell_nodes<'a, T>(&self, cell: &Cell<'a, T>) -> &[usize] {
        &self.cell_nodes[self.cell_nodes_starts[cell.local_id]..self.cell_nodes_starts[cell.local_id + 1]]
    }

    pub fn face_nodes<'a, T>(&self, face: &Face<'a, T>) -> [usize; FACE_CORNERS] {
        
        let face_corners = face.corners_int;
        let tree = if face.face_id_side == 0 {face.cell0.tree_id} else {face.cell1.as_ref().unwrap().tree_id};

        let mut out = [usize::MAX; FACE_CORNERS];
        for i in 0..FACE_CORNERS {
            let mut hash = [u32::MAX; 4];
            for j in 0..DIM {
                hash[j] = face_corners[i][j] as u32;
            }
            hash[3] = tree as u32;

            match self.node_map.get(&hash) {
                Some(n) => {
                    out[i] = *n as usize;
                },
                None => {
                    panic!("Node not found in map {:?}", hash);
                }
            }
        }

        out
    }

}


impl NodeNumbering {
    fn compute<T>(mut self, grid: &Grid<T>) -> Self {

        let mut cell_nodes = vec![HashSet::new(); grid.len_with_ghosts()];

        let mut node_id = 0;
        grid.map_faces(|face| {

            let face_corners = face.corners_int;
            let tree = if face.face_id_side == 0 {face.cell0.tree_id} else {face.cell1.as_ref().unwrap().tree_id};

            let c0 = face.cell0.local_id;
            let c1 = match face.cell1.as_ref() {
                Some(c1) => Some(c1.local_id),
                None => None,
            };

            for i in 0..FACE_CORNERS {

                let mut hash = [u32::MAX; 4];
                for j in 0..DIM {
                    hash[j] = face_corners[i][j] as u32;
                }
                hash[3] = tree as u32;

                cell_nodes[c0].insert(node_id);
                if let Some(c1) = c1 {
                    cell_nodes[c1].insert(node_id);
                }

                self.node_map.insert(hash, node_id);
                self.node_positions.push(face.corners[i]);
                node_id += 1;
            }

        });

        for nodes in cell_nodes {
            for n in nodes {
                self.cell_nodes.push(n as usize);
            }
            self.cell_nodes_starts.push(self.cell_nodes.len());
        }

        self
    }
}
