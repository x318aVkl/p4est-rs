use std::collections::HashMap;

use p4est_sys::consts::{CELL_CORNERS, FACE_CORNERS};

use crate::grid::{Grid, cell::{Cell, Face}, corners::face_corner_ids};



pub struct NodeNumbering {
    nodes: *mut p4est_sys::p4est_lnodes,

    // computes the element node local ids, includes the ghost elements
    all_element_node_local_ids: Vec<usize>,

    // map from [n0, X, n1], (n0, n1, usize::MAX) -> X, where X is non local, hanging
    // also maps from [n0, n1, n2] -> F, where F is the center node of a face
    hanging_map: HashMap<[usize; 3], usize>,
}

impl NodeNumbering {

    pub fn new<T>(grid: &Grid<T>) -> Self {
        Self {
            nodes: unsafe {p4est_sys::p4est_lnodes_new(grid.grid, grid.ghosts, 1)},
            all_element_node_local_ids: vec![usize::MAX; grid.len_with_ghosts() * CELL_CORNERS],
            hanging_map: HashMap::new(),
        }.compute(&grid)
    }

    pub fn cell_nodes<'a, T>(&self, cell: &Cell<'a, T>) -> [usize; CELL_CORNERS] {
        let cell_id = cell.local_id;

        let mut out = [usize::MAX; CELL_CORNERS];

        let n_local_cells = unsafe {(*self.nodes).num_local_elements} as usize;
        if cell_id >= n_local_cells {
            return out;
        }

        for i in 0..CELL_CORNERS {
            let nidx = unsafe { *(*self.nodes).element_nodes.offset((CELL_CORNERS * cell_id + i) as isize) };
            out[i] = nidx as usize;
        }

        for i in 0..CELL_CORNERS {
            out[i] = self.all_element_node_local_ids[cell_id * CELL_CORNERS + i];
        }

        out
    }

    pub fn face_nodes<'a, T>(&self, face: &Face<'a, T>) -> [usize; FACE_CORNERS] {
        let face_id = face.face_id;
        let cell = if face.face_id_side == 0 {
            &face.cell0
        } else {
            face.cell1.as_ref().unwrap()
        };

        if cell.is_ghost {
            // actually use the other cell
            let other_cell = if face.face_id_side == 0 {
                face.cell1.as_ref().unwrap()
            } else {
                &face.cell0
            };
            let other_face_id = opposite_id(face_id);
            let fc = face_corner_ids(other_face_id);
            let mut out = [0; FACE_CORNERS];
            let cnodes = self.cell_nodes(other_cell);
            for i in 0..FACE_CORNERS {
                out[i] = cnodes[fc[i] as usize];
            }
            if other_cell.level < cell.level {
                // the reference cell is coarser
                // this is not good, we need the small face nodes
                let subout = subface_corners_local(
                    face.subface_id.unwrap(), 
                    out, 
                    &self.hanging_map,
                    None,
                );
                out = subout;
            }
            return out;
        }
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
            p4est_sys::p4est_lnodes_destroy(self.nodes);
        }
    }
}




impl NodeNumbering {
    fn compute<T>(mut self, grid: &Grid<T>) -> Self {

        grid.map_cells(|cell| {
            let id = cell.local_id;

            for i in 0..CELL_CORNERS {
                let nidx = unsafe { *(*self.nodes).element_nodes.offset((CELL_CORNERS * id + i) as isize) };
                self.all_element_node_local_ids[CELL_CORNERS *id + i] = nidx as usize;
            }
        });

        // correct the hanging nodes in certain faces
        let dummy_map = HashMap::new();
        let mut currglobid = unsafe { (*self.nodes).num_local_nodes } as usize;
        grid.map_faces(|face| {
            let c_ref = if face.face_id_side == 0 {
                &face.cell0
            } else {
                face.cell1.as_ref().unwrap()
            };

            if c_ref.is_ghost {
                // we must do something about that, reference cell is a ghost
                let c_other = if face.face_id_side == 0 {
                    face.cell1.as_ref().unwrap()
                } else {
                    &face.cell0
                };

                if c_other.is_ghost {
                    panic!("erhmmm both cells are ghosts, cant find face local nodes then");
                }

                if c_ref.level < c_other.level {
                    // must handle hanging nodes

                    let other_face_id = opposite_id(face.face_id);
                    let fc = face_corner_ids(other_face_id);
                    let mut out = [0; FACE_CORNERS];
                    let cnodes = self.cell_nodes(c_other);
                    for i in 0..FACE_CORNERS {
                        out[i] = cnodes[fc[i] as usize];
                    }

                    let _ = subface_corners_local(
                        face.subface_id.unwrap(), 
                        out, 
                        &dummy_map,
                        Some((&mut currglobid, &mut self.hanging_map)),
                    );
                }
            }
        });


        self
    }
}


fn opposite_id(side: i8) -> i8 {
    match side {
        0 => 1,
        1 => 0,
        2 => 3,
        3 => 2,
        4 => 5,
        5 => 4,
        _ => panic!("invalid face side {}", side)
    }
}

// 10 means the hanging node in the center of the edge
#[cfg(feature = "2d")]
fn subface_corners(subface: u8) -> &'static [u8] {
    match subface {
        0 => &[0, 8],
        1 => &[8, 1],
        _ => panic!("error, subface invalid {}", subface),
    }
}
//   2 ---13--- 3
//   |          |
//  10    20    11
//   |          |
//   0 ---12--- 1
#[cfg(feature = "3d")]
fn subface_corners(subface: u8) -> &'static [u8] {
    match subface {
        0 => &[0, 12, 10, 20],
        1 => &[12, 1, 20, 11],
        2 => &[10, 20, 2, 13],
        3 => &[20, 11, 13, 3],
        _ => panic!("error, subface invalid {}", subface),
    }
}

enum CornerDecode {
    Direct(u8),
    FetchFromMap([usize; 3]),
}

fn subface_corners_local(
    subface: u8, 
    parent_face_local: [usize; FACE_CORNERS], 
    hanging_map: &HashMap<[usize; 3], usize>,
    mut insert_map_builder: Option<(&mut usize, &mut HashMap<[usize; 3], usize>)>,
) -> [usize; FACE_CORNERS] {
    let mut out = [usize::MAX; FACE_CORNERS];

    let corners = subface_corners(subface);

    for (i, c) in corners.iter().enumerate() {
        let c = *c;
        let c_op;
        #[cfg(feature = "2d")]
        {
            c_op = if c < 4 {
                CornerDecode::Direct(c)
            } else if c == 8 {
                CornerDecode::FetchFromMap([parent_face_local[0], parent_face_local[1], usize::MAX])
            } else {
                panic!("unrecognized corner decode command {}", c)
            };
        }
        #[cfg(feature = "3d")]
        {
            c_op = if c < 4 {
                CornerDecode::Direct(c)
            } else if c == 8 {
                CornerDecode::FetchFromMap([parent_face_local[0], parent_face_local[1], usize::MAX])
            } else if c == 10 {
                CornerDecode::FetchFromMap([parent_face_local[0], parent_face_local[2], usize::MAX])
            } else if c == 11 {
                CornerDecode::FetchFromMap([parent_face_local[1], parent_face_local[3], usize::MAX])
            } else if c == 12 {
                CornerDecode::FetchFromMap([parent_face_local[0], parent_face_local[1], usize::MAX])
            } else if c == 13 {
                CornerDecode::FetchFromMap([parent_face_local[2], parent_face_local[3], usize::MAX])
            } else if c == 20 {
                let mut pclone = parent_face_local;
                pclone.sort();
                CornerDecode::FetchFromMap([pclone[0], pclone[1], pclone[2]])
            } else {
                panic!("unrecognized corner decode command {}", c)
            };
        }

        match c_op {
            CornerDecode::Direct(id) => {
                out[i] = parent_face_local[id as usize];
            },
            CornerDecode::FetchFromMap(mut hash) => {
                hash.sort();

                let id = if let Some((currid, map)) = insert_map_builder.as_mut() {
                    map.entry(hash).or_insert_with(|| {
                        **currid += 1;
                        **currid - 1
                    })
                } else {
                    hanging_map.get(&hash).expect("found hash in hanging map")
                };
                out[i] = *id;
            }
        }
    }

    out
}