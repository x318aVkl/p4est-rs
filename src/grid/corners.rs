use p4est_sys::consts::{CELL_CORNERS, DIM, FACE_CORNERS};

use crate::basetree::BaseTree;



pub(crate) unsafe fn qcoord_to_vertex(tree: &BaseTree, treeid: i32, node: [i32; DIM]) -> [f64; DIM] {

    let treeid = treeid;

    let mut vxyz = [0.0; DIM];

    let root_len = (2_u32).pow(p4est_sys::P4EST_MAXLEVEL);

    for i in 0..DIM {
        vxyz[i] = (node[i] as f64) / (root_len as f64);
    }

    let element = tree.element(treeid as usize);

    element.tree_relative_position_to_global_position(vxyz)
}


pub(crate) unsafe fn collect_quadrant_corner_nodes(quad: *mut p4est_sys::p4est_quadrant, nodes: &mut [[i32; DIM]; CELL_CORNERS]) {
    for i in 0..CELL_CORNERS {
        unsafe {
            let mut node = *quad;
            p4est_sys::p4est_quadrant_corner_node(quad, i as i32, (&mut node) as *mut p4est_sys::p4est_quadrant);

            let mut node_coords = [-1_i32; DIM];
            node_coords[0] = node.x;
            node_coords[1] = node.y;

            #[cfg(feature = "3d")]
            { node_coords[2] = node.z; }

            nodes[i] = node_coords;
        }
    }
}
pub(crate) unsafe fn cell_corners(tree: &BaseTree, treeid: i32, quad: *mut p4est_sys::p4est_quadrant) -> ([[f64; DIM]; CELL_CORNERS], [[u32; DIM]; CELL_CORNERS]) {
    let mut corners = [[0.0; DIM]; CELL_CORNERS];

    let mut corner_nodes = [[0_i32; DIM]; CELL_CORNERS];
    unsafe { collect_quadrant_corner_nodes(quad, &mut corner_nodes); }

    for i in 0..CELL_CORNERS {
        unsafe { corners[i] = qcoord_to_vertex(tree, treeid, corner_nodes[i]); }
    }

    let mut c = [[0_u32; DIM]; CELL_CORNERS];
    for i in 0..CELL_CORNERS {
        for j in 0..DIM {
            c[i][j] = corner_nodes[i][j] as u32;
        }
    }

    (corners, c)
}

//  2 ----- 3
//  |       |
//  |       |
//  0 ----- 1
#[cfg(feature = "2d")]
pub(crate) fn face_corner_ids(face: i8) -> [i32; FACE_CORNERS] {
    match face {

        0 => [2, 0],
        1 => [1, 3],
        2 => [0, 1],
        3 => [3, 2],

        _ => panic!("invalid face id {}", face)
    }
}

#[cfg(feature = "3d")]
pub(crate) fn face_corner_ids(face: i8) -> [i32; FACE_CORNERS] {
    match face {

        0 => [0, 4, 6, 2],
        1 => [1, 3, 7, 5],
        2 => [0, 1, 5, 4],
        3 => [3, 2, 6, 7],
        4 => [0, 2, 3, 1],
        5 => [4, 5, 7, 6],

        _ => panic!("invalid face id {}", face)
    }
}

pub(crate) fn face_corners_from_cell<T>(face: i8, cell_corners: [[T; DIM]; CELL_CORNERS]) -> [[T; DIM]; FACE_CORNERS] where T: Default + Copy {
    let mut out = [[T::default(); DIM]; FACE_CORNERS];

    for (k, i) in face_corner_ids(face).into_iter().enumerate() {
        out[k] = cell_corners[i as usize];
    }

    out
}
