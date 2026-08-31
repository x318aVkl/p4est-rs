use std::ffi::c_void;

use p4est_sys::consts::FACE_CORNERS;

use crate::{basetree::BaseTree, grid::{CellData, Grid, cell::Face, corners::{cell_corners, face_corners_from_cell}}};


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


    pub fn map_faces<'a, F>(&'a self, f: F) where F: FnMut(Face<'a, T>) {
        unsafe {
            USER_FACE_FN = Some(&f as *const F as *mut c_void);

            let face_info = FaceIterInfo {
                geometry: &self.base_tree as *const BaseTree,
                ghost_data: self.ghost_data.as_ptr() as *mut CellData<T>,
                local_size: self.local_len(),
                current_face_id: 0,
            };

            let user_data = &face_info as *const FaceIterInfo<T> as *mut c_void;

            #[cfg(feature = "2d")]
            p4est_sys::p4est_iterate(
                self.grid,
                self.ghosts,
                user_data,
                None,
                Some(iter_face_fn::<'a, F, T>),
                None,
            );

            #[cfg(feature = "3d")]
            p4est_sys::p4est_iterate(
                self.grid,
                self.ghosts,
                user_data,
                None,
                Some(iter_face_fn::<'a, F, T>),
                None,
                None,
            );

            USER_FACE_FN = None;
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
        let (corners, corners_int) = cell_corners(base_tree, (*info).treeid, (*info).quad);
        
        let cell = Cell {
            data: &cell_data.data, 
            local_id: id as usize, 
            global_id: gid as usize,
            is_ghost: false,
            owner_rank: cell_data.owner_rank,
            level: (*(*info).quad).level as u8,
            corners,
            corners_int,
            raw_quad: (*info).quad,
            tree_id: (*info).treeid,
        };

        let f = USER_VOLUME_FN.unwrap() as *mut F;

        (*f)(cell);
    }

}



struct FaceIterInfo<T> {
    geometry: *const BaseTree,
    ghost_data: *mut CellData<T>,
    local_size: usize,
    current_face_id: usize,
}


static mut USER_FACE_FN: Option<*mut std::ffi::c_void> = None;

extern "C" fn iter_face_fn<'a, F, T>(info: *mut p4est_sys::p4est_iter_face_info, data: *mut std::ffi::c_void) where F: FnMut(Face<'a, T>), T: 'a + Sized {

    unsafe {

        let sides = (*info).sides;
        let side0: *mut p4est_sys::p4est_iter_face_side = sides.array.offset((sides.elem_size * 0) as isize) as *mut p4est_sys::p4est_iter_face_side;
        
        let mut quads_0: [Option<(*const p4est_sys::p4est_quadrant, usize, i32, bool)>; FACE_CORNERS] = [None; FACE_CORNERS];
        let mut side0len = 1;

        let ghost_info_ptr = data as *mut FaceIterInfo<T>;

        let geom = & *((*ghost_info_ptr).geometry);

        let ghost_data = (*ghost_info_ptr).ghost_data;

        let local_len = (*ghost_info_ptr).local_size as u32;

        let face_id = &mut (*ghost_info_ptr).current_face_id;

        if (*side0).is_hanging == 0 {
            let qid = (*side0).is.full.quadid as usize;
            quads_0[0] = Some(((*side0).is.full.quad, qid, (*side0).treeid, (*side0).is.full.is_ghost == 1));
        } else {
            // it is hanging, N_CHILD quads
            for i in 0..FACE_CORNERS {
                let qid = (*side0).is.hanging.quadid[i] as usize;
                quads_0[i] = Some(((*side0).is.hanging.quad[i], qid, (*side0).treeid, (*side0).is.hanging.is_ghost[i] == 1));
            }
            side0len = FACE_CORNERS;
        }

        if sides.elem_count == 1 {
            // this is a boundary!

            for i in 0..side0len {
                let (quad_0, q0id, s0tid, is_ghost_0) = quads_0[i].unwrap();

                let (corners0, corners_int_0) = cell_corners(geom,  (*side0).treeid, quad_0 as *mut p4est_sys::p4est_quadrant);

                let cell_data = if is_ghost_0 {
                        &mut *(ghost_data.offset(q0id as isize) as *mut CellData<T>) as &mut CellData<T>
                    } else {
                        &mut *((*quad_0).p.user_data as *mut CellData<T>)
                    };
                
                let face_corners = face_corners_from_cell((*side0).face, corners0);
                let face_corners_int = face_corners_from_cell((*side0).face, corners_int_0);

                let face = Face {
                    cell0: Cell {
                        data: &(*cell_data).data,
                        local_id: (*cell_data).local_id as usize,
                        global_id: (*cell_data).global_id as usize,
                        level: (*quad_0).level as u8,
                        is_ghost: is_ghost_0,
                        owner_rank: (*cell_data).owner_rank,
                        corners: corners0,
                        corners_int: corners_int_0,
                        raw_quad: quad_0,
                        tree_id: s0tid,
                    },
                    cell1: None,
                    id: *face_id,
                    face_id: (*side0).face,
                    face_id_side: 0,
                    corners: face_corners,
                    corners_int: face_corners_int,
                    subface_id: if side0len > 1 {Some(i as u8)} else {None},
                };

                let f = USER_FACE_FN.unwrap() as *mut F;

                (*f)(face);

                *face_id += 1;
            }
            
        } else if sides.elem_count == 2 {

            // this is an internal face, two sides
            let side1: *mut p4est_sys::p4est_iter_face_side = sides.array.offset((sides.elem_size * 1) as isize) as *mut p4est_sys::p4est_iter_face_side;

            let mut quads_1: [Option<(*const p4est_sys::p4est_quadrant, usize, i32, bool)>; FACE_CORNERS] = [None; FACE_CORNERS];
            let mut side1len = 1;

            if (*side1).is_hanging == 0 {
                quads_1[0] = Some(((*side1).is.full.quad, (*side1).is.full.quadid as usize, (*side1).treeid, (*side1).is.full.is_ghost == 1));
            } else {
                // it is hanging, N_CHILD quads
                for i in 0..FACE_CORNERS {
                    quads_1[i] = Some(((*side1).is.hanging.quad[i], (*side1).is.hanging.quadid[i] as usize, (*side1).treeid, (*side1).is.hanging.is_ghost[i] == 1));
                }
                side1len = FACE_CORNERS;
            }

            for i in 0..side0len {
                for j in 0..side1len {

                    let (quad_0, q0id, s0tid, is_ghost_0) = quads_0[i].unwrap();
                    let (quad_1, q1id, s1tid, is_ghost_1) = quads_1[j].unwrap();

                    let (corners0, corners_int_0) = cell_corners(geom,  (*side0).treeid, quad_0 as *mut p4est_sys::p4est_quadrant);
                    let (corners1, corners_int_1) = cell_corners(geom,  (*side1).treeid, quad_1 as *mut p4est_sys::p4est_quadrant);

                    let quad_0_data = if is_ghost_0 {
                        &mut *(ghost_data.offset(q0id as isize) as *mut CellData<T>) as &mut CellData<T>
                    } else {
                        &mut *((*quad_0).p.user_data as *mut CellData<T>)
                    };

                    let quad_1_data = if is_ghost_1 {
                        &mut *(ghost_data.offset(q1id as isize) as *mut CellData<T>) as &mut CellData<T>
                    } else {
                        &mut *((*quad_1).p.user_data as *mut CellData<T>)
                    };

                    let ((face_corners, lface_id, lface_id_side), face_corners_int) = if side0len >= side1len {
                        ((face_corners_from_cell((*side0).face, corners0), (*side0).face, 0), face_corners_from_cell((*side0).face, corners_int_0))
                    } else {
                        ((face_corners_from_cell((*side1).face, corners1), (*side1).face, 1), face_corners_from_cell((*side1).face, corners_int_1))
                    };

                    // edit, allow both side 0 and side 1 to be ghosts
                    // this allows collection of all ghost data, required for lnodes
                    //if (quad_0_data.local_id < local_len) || (quad_1_data.local_id < local_len) {

                        let face = if quad_0_data.local_id < quad_1_data.local_id {
                            Face {
                                cell0: Cell {
                                    data: &quad_0_data.data,
                                    local_id: quad_0_data.local_id as usize,
                                    global_id: quad_0_data.global_id as usize,
                                    level: (*quad_0).level as u8,
                                    is_ghost: is_ghost_0,
                                    owner_rank: quad_0_data.owner_rank,
                                    corners: corners0,
                                    corners_int: corners_int_0,
                                    raw_quad: quad_0,
                                    tree_id: s0tid,
                                },
                                cell1: Some(Cell {
                                    data: &quad_1_data.data,
                                    local_id: quad_1_data.local_id as usize,
                                    global_id: quad_1_data.global_id as usize,
                                    level: (*quad_1).level as u8,
                                    is_ghost: is_ghost_1,
                                    owner_rank: quad_1_data.owner_rank,
                                    corners: corners1,
                                    corners_int: corners_int_1,
                                    raw_quad: quad_1,
                                    tree_id: s1tid,
                                }),
                                id: *face_id,
                                face_id: lface_id,
                                face_id_side: lface_id_side,
                                corners: face_corners,
                                corners_int: face_corners_int,
                                subface_id: if (side0len > 1) || (side1len > 1) {Some(side0.max(side1) as u8)} else {None},
                            }
                        } else {
                            Face {
                                cell0: Cell {
                                    data: &quad_1_data.data,
                                    local_id: quad_1_data.local_id as usize,
                                    global_id: quad_1_data.global_id as usize,
                                    level: (*quad_1).level as u8,
                                    is_ghost: is_ghost_1,
                                    owner_rank: quad_1_data.owner_rank,
                                    corners: corners1,
                                    corners_int: corners_int_1,
                                    raw_quad: quad_1,
                                    tree_id: s1tid,
                                },
                                cell1: Some(Cell {
                                    data: &quad_0_data.data,
                                    local_id: quad_0_data.local_id as usize,
                                    global_id: quad_0_data.global_id as usize,
                                    level: (*quad_0).level as u8,
                                    is_ghost: is_ghost_0,
                                    owner_rank: quad_0_data.owner_rank,
                                    corners: corners0,
                                    corners_int: corners_int_0,
                                    raw_quad: quad_0,
                                    tree_id: s0tid,
                                }),
                                id: *face_id,
                                face_id: lface_id,
                                face_id_side: if lface_id_side == 0 {1} else {0},   // FLIP IT, its bad but anyway
                                corners: face_corners,
                                corners_int: face_corners_int,
                                subface_id: if (side0len > 1) || (side1len > 1) {Some(side0.max(side1) as u8)} else {None},
                            }
                        };

                        let f = USER_FACE_FN.unwrap() as *mut F;

                        (*f)(face);

                        *face_id += 1;

                   //}

                }
            }

        } else {
            panic!("face has not 1 or two sides?? wth thats an issue, better call saul!");
        }

    }

}

