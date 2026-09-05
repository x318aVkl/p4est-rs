use std::os::raw::c_void;

use mpi::{raw::AsRaw, traits::Communicator};
use p4est_sys::consts::CELL_CORNERS;

use crate::grid::{CellData, Grid, cell::Cell};



pub fn transfer_data_custom_partition<'a, V, T>(
    old_mesh: &'a Grid<T>,
    old_data: &[V],
    new_mesh: &'a Grid<T>,
    new_data: &mut [V],
) -> Result<(), i32> {

    // transfer data after partitioning
    unsafe {
        p4est_sys::p4est_transfer_fixed(
            (*new_mesh.grid).global_first_quadrant,
            (*old_mesh.grid).global_first_quadrant,
            new_mesh.communicator.as_raw().0 as *mut c_void as p4est_sys::MPI_Comm,
            0,
            new_data.as_mut_ptr() as *mut c_void,
            old_data.as_ptr() as *const c_void,
            size_of::<V>(),
        );
    }

    Ok(())
}


pub fn transfer_data_custom_adapt<'a, FC, FR, V, T>(
    old_mesh: &'a Grid<T>,
    old_data: &[V],
    new_mesh: &'a Grid<T>,
    new_data: &mut [V],
    mut coarsening_function: FC,
    mut refining_function: FR,
) -> Result<(), i32>
where 
FC: FnMut([(Cell<'a, T>, &V); CELL_CORNERS], (Cell<'a, T>, &mut V)),
FR: FnMut((Cell<'a, T>, &V), [(Cell<'a, T>, &mut V); CELL_CORNERS]),
V: Clone + Default,
{

    // transfer data after adaptation
    // iterate over both mesh at the same time

unsafe {

    for tt in (*new_mesh.grid).first_local_tree..=(*new_mesh.grid).last_local_tree {

        let tarrayelemsize = (*(*(*old_mesh).grid).trees).elem_size as i32;
        let ptree = (*(*(*old_mesh).grid).trees).array.offset((tarrayelemsize * tt) as isize) as *mut p4est_sys::p4est_tree;

        let mut pquad = (*ptree).quadrants.array as *mut p4est_sys::p4est_quadrant;
        
        let tarrayelemsize = (*(*(*new_mesh).grid).trees).elem_size as i32;
        let ntree = (*(*(*new_mesh).grid).trees).array.offset((tarrayelemsize * tt) as isize) as *mut p4est_sys::p4est_tree;
        
        let mut nquad = (*ntree).quadrants.array as *mut p4est_sys::p4est_quadrant;

        let mut ptraverse: usize = 0;
        let ptreesize = (*ptree).quadrants.elem_count;

        let mut ntraverse: usize = 0;
        let ntreesize = (*ntree).quadrants.elem_count;

        // simultaneous loop over old and new quadrants
        loop {

            //println!("{} : {} {}", new_mesh.communicator().rank(), (*pquad).level , (*nquad).level);

            if (*pquad).level == (*nquad).level {
                // old and new quadrants are the same size, copy value
                let nid = (*((*nquad).p.user_data as *mut CellData<T>)).local_id as usize;
                let pid = (*((*pquad).p.user_data as *mut CellData<T>)).local_id as usize;
                new_data[nid] = old_data[pid].clone();

                nquad = nquad.add(1);
                ntraverse += 1;

                pquad = pquad.add(1);
                ptraverse += 1;
            } else if ((*pquad).level + 1) == (*nquad).level {
                // new quadrant are refined from the old one, use refine function to determine the new quadrant data
    
                #[allow(invalid_reference_casting)]
                #[cfg(feature = "2d")]
                {
                    refining_function(
                        (
                            Cell::from_quad(&old_mesh.base_tree, tt, pquad),
                            &old_data[(*((*pquad).p.user_data as *mut CellData<T>)).local_id as usize],
                        ),
                        [
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(0)),
                                &mut *(&new_data[(*((*nquad.offset(0)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(1)),
                                &mut *(&new_data[(*((*nquad.offset(1)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(2)),
                                &mut *(&new_data[(*((*nquad.offset(2)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(3)),
                                &mut *(&new_data[(*((*nquad.offset(3)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                        ]
                    );
                }
                #[allow(invalid_reference_casting)]
                #[cfg(feature = "3d")]
                {
                    refining_function(
                        (
                            Cell::from_quad(&old_mesh.base_tree, tt, pquad),
                            &old_data[(*((*pquad).p.user_data as *mut CellData<T>)).local_id as usize],
                        ),
                        [
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(0)),
                                &mut *(&new_data[(*((*nquad.offset(0)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(1)),
                                &mut *(&new_data[(*((*nquad.offset(1)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(2)),
                                &mut *(&new_data[(*((*nquad.offset(2)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(3)),
                                &mut *(&new_data[(*((*nquad.offset(3)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(4)),
                                &mut *(&new_data[(*((*nquad.offset(4)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(5)),
                                &mut *(&new_data[(*((*nquad.offset(5)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(6)),
                                &mut *(&new_data[(*((*nquad.offset(6)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                            (
                                Cell::from_quad(&new_mesh.base_tree, tt, nquad.offset(7)),
                                &mut *(&new_data[(*((*nquad.offset(7)).p.user_data as *mut CellData<T>)).local_id as usize] as *const V as *mut V),
                            ),
                        ]
                    );
                }

                nquad = nquad.add(CELL_CORNERS);
                ntraverse += CELL_CORNERS;

                pquad = pquad.add(1);
                ptraverse += 1;
            } else if (*pquad).level == ((*nquad).level + 1) {
                // new quadrant is coarsened from the old one, use coarsening function

                #[allow(invalid_reference_casting)]
                #[cfg(feature = "2d")]
                {
                    coarsening_function(
                        [
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(0)),
                                &old_data[(*((*pquad.offset(0)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(1)),
                                &old_data[(*((*pquad.offset(1)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(2)),
                                &old_data[(*((*pquad.offset(2)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(3)),
                                &old_data[(*((*pquad.offset(3)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                        ],
                        (
                            Cell::from_quad(&new_mesh.base_tree, tt, nquad),
                            &mut new_data[(*((*nquad).p.user_data as *mut CellData<T>)).local_id as usize],
                        ),
                    );
                }
                #[allow(invalid_reference_casting)]
                #[cfg(feature = "3d")]
                {
                    coarsening_function(
                        [
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(0)),
                                &old_data[(*((*pquad.offset(0)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(1)),
                                &old_data[(*((*pquad.offset(1)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(2)),
                                &old_data[(*((*pquad.offset(2)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(3)),
                                &old_data[(*((*pquad.offset(3)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(4)),
                                &old_data[(*((*pquad.offset(4)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(5)),
                                &old_data[(*((*pquad.offset(5)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(6)),
                                &old_data[(*((*pquad.offset(6)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                            (
                                Cell::from_quad(&old_mesh.base_tree, tt, pquad.offset(7)),
                                &old_data[(*((*pquad.offset(7)).p.user_data as *mut CellData<T>)).local_id as usize],
                            ),
                        ],
                        (
                            Cell::from_quad(&new_mesh.base_tree, tt, nquad),
                            &mut new_data[(*((*nquad).p.user_data as *mut CellData<T>)).local_id as usize],
                        ),
                    );
                }

                nquad = nquad.add(1);
                ntraverse += 1;

                pquad = pquad.add(CELL_CORNERS);
                ptraverse += CELL_CORNERS;
            } else {
                // panic, this should not happen
                panic!("Error in rank {}, new and old quadrant levels {}, {} are invalid", new_mesh.communicator.rank(), (*nquad).level, (*pquad).level);
            }

            if ptraverse >= ptreesize {
                break;
            }
            if ntraverse >= ntreesize {
                break;
            }
            
        }
    }
    
}



    Ok(())
}