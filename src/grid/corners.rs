use p4est_sys::consts::{CELL_CORNERS, DIM};

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
pub(crate) unsafe fn cell_corners(tree: &BaseTree, treeid: i32, quad: *mut p4est_sys::p4est_quadrant) -> [[f64; DIM]; CELL_CORNERS] {
    let mut corners = [[0.0; DIM]; CELL_CORNERS];

    let mut corner_nodes = [[0_i32; DIM]; CELL_CORNERS];
    unsafe { collect_quadrant_corner_nodes(quad, &mut corner_nodes); }

    for i in 0..CELL_CORNERS {
        unsafe { corners[i] = qcoord_to_vertex(tree, treeid, corner_nodes[i]); }
    }

    corners
}