use std::{collections::{HashMap, HashSet}, io::BufRead};

use p4est_sys::consts::{CELL_CORNERS, CELL_FACES, DIM};





// A base tree that describes potentially high order elements
// all elements are quads or hexahedra
pub struct BaseTree {
    nodes: Vec<[f64; DIM]>,
    elements: Vec<usize>,
    element_starts: Vec<usize>,
}


#[derive(Debug)]
pub struct Element<'a> {
    all_nodes: &'a [[f64; DIM]],
    elem_nodes: &'a [usize],
}


impl std::fmt::Debug for BaseTree {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\nTriangulation with:\n- {} nodes\n- {} elements\n- order {}\n", self.nodes.len(), self.element_starts.len() - 1, self.element(0).order())
    }
}


impl<'a> Element<'a> {
    pub fn order(&self) -> u8 {
        let elen = self.elem_nodes.len();
        let elen = elen as f64;
        let blen = elen.powf(1.0 / (DIM as f64));
        (blen as u8) - 1
    }

    pub fn corners(&self) -> [usize; CELL_CORNERS] {
        let mut corners = [0; CELL_CORNERS];

        corners[0] = self.elem_nodes[0];
        corners[1] = self.elem_nodes[1];
        corners[2] = self.elem_nodes[3];
        corners[3] = self.elem_nodes[2];

        corners
    }


    pub fn tree_relative_position_to_global_position(
        &self,
        norm_pos: [f64; DIM],
    ) -> [f64; DIM] {
        // norm pos input is between zero and one

        // convert input to between -1 and +1
        for _ in 0..DIM {
            //norm_pos[i] = 2.0 * norm_pos[i] - 1.0;
        }

        // lagrange basis evaluations
        let mut out = [0.0; DIM];
        
        for i in 0..self.elem_nodes.len() {
            // add the contribution from this node
            let li = basis(norm_pos, i);

            let node = self.all_nodes[self.elem_nodes[i]];
            for k in 0..DIM {
                out[k] += node[k] * li;
            }
        }

        out
    }
}


impl BaseTree {

    pub(crate) unsafe fn build_connectivity(&self) -> *mut p4est_sys::p4est_connectivity {
        unsafe {

            let mut unique_corners = HashSet::<usize>::new();
            for elem in 0..self.element_len() {
                for c in self.element(elem).corners() {
                    unique_corners.insert(c);
                }
            }
            let unique_corners = unique_corners.into_iter().collect::<Vec<_>>();
            let mut global_to_corner = HashMap::new();

            for (i, c) in unique_corners.iter().enumerate() {
                global_to_corner.insert(*c, i);
            }

            #[cfg(feature = "2d")]
            let conn = p4est_sys::p4est_connectivity_new(
                unique_corners.len() as i32, 
                self.element_len() as i32, 
                0,
                0,
            );
            #[cfg(feature = "3d")]
            let conn = p4est_sys::p4est_connectivity_new(
                unique_corners.len() as i32, 
                self.element_len() as i32, 
                0,
                0,
                0,
                0,
            );

            // fill conn->vertices
            // and conn->tree_to_vertex
            for (i, c) in unique_corners.iter().enumerate() {
                let c = *c;

                let node = self.nodes[c];

                *(*conn).vertices.offset((3 * i + 0) as isize) = node[0];
                *(*conn).vertices.offset((3 * i + 1) as isize) = node[1];
            }

            for elem in 0..self.element_len() {

                let corners = self.element(elem).corners();
                for n in 0..CELL_CORNERS {

                    let local_node = *global_to_corner.get(&corners[n]).unwrap();

                    *(*conn).tree_to_vertex.offset((CELL_CORNERS * elem + n) as isize) = local_node as i32;
                }
            }

            
            for tree in 0..((*conn).num_trees) as usize {
                for face in 0..CELL_FACES {
                    *((*conn).tree_to_tree.offset((CELL_FACES * tree + face) as isize)) = tree as i32;
                    *((*conn).tree_to_face.offset((CELL_FACES * tree + face) as isize)) = face as i8;
                }
            }


            // check the connectivity is valid
            let ecode = p4est_sys::p4est_connectivity_is_valid(conn);
            if ecode != 1 {
                panic!("error, connectivity is not valid, code {}", ecode);
            }


            // complete the connecvitiy
            p4est_sys::p4est_connectivity_complete(conn);

            conn
        }
    }

}


impl BaseTree {

    pub fn node(&self, index: usize) -> &[f64; DIM] {
        &self.nodes[index]
    }

    pub fn element<'a>(&'a self, index: usize) -> Element<'a> {
        let ls = self.element_starts[index];
        let le = self.element_starts[index + 1];
        Element { all_nodes: &self.nodes, elem_nodes: &self.elements[ls..le] }
    }

    pub fn element_len(&self) -> usize {
        self.element_starts.len() - 1
    }


    pub fn from_su2<F>(f: F) -> Result<Self, std::io::Error> where F: BufRead {

        let mut nodes = vec![];
        let mut elements = vec![];
        let mut element_starts = vec![0];

        let mut section = "".to_string();
        for line in f.lines() {
            let line = line?;

            let line = line.trim();

            if line.len() == 0 {
                continue;
            }

            // figure out the section were in
            if line.contains("=") {
                let mut ls = line.split("= ");
                section = ls.nth(0).unwrap().to_string();
                let id = ls.nth(0).unwrap().parse::<usize>().unwrap();

                if section == "NDIME" {
                    assert_eq!(id, DIM);
                }
            } else {


                if section == "NELEM" {
                    // read an element
                    let mut ls = line.split(" ");
                    let elem_id = ls.nth(0).unwrap().parse::<u8>().unwrap();

                    #[cfg(feature = "2d")]
                    { assert_eq!(elem_id, 9); }


                    #[cfg(feature = "3d")]
                    { assert_eq!(elem_id, 12); }

                    for value_str in ls {
                        elements.push(
                            value_str.parse().unwrap()
                        );
                    }
                    elements.pop();
                    element_starts.push(elements.len());
                } else if section == "NPOIN" {
                    // read a node
                    
                    let ls = line.split(" ");

                    let mut node = [0.0; DIM];

                    for (k, v) in ls.enumerate() {
                        if k >= DIM {break;}
                        node[k] = v.parse().unwrap();
                    }

                    nodes.push(node);

                }

            }


        }


        Ok(Self { nodes, elements, element_starts })
    }

}


impl BaseTree {

    #[cfg(feature = "2d")]
    pub fn new_unitsquare() -> Self {
        Self { nodes: vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ], 
        elements: vec![
            0, 1, 2, 3,
        ], 
        element_starts: vec![0, 4] }
    }

    #[cfg(feature = "3d")]
    pub fn new_unitsquare() -> Self {
        Self { nodes: vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 1.0],
            [0.0, 1.0, 1.0],
            [1.0, 1.0, 1.0],
        ], 
        elements: vec![
            0, 1, 2, 3, 4, 5, 6, 7
        ], 
        element_starts: vec![0, 8] }
    }

}





// basis function for tensor element
#[cfg(feature = "2d")]
pub fn basis(
    point: [f64; DIM],
    i: usize
) -> f64 {
    let x = point[0];
    let y = point[1];
    if i == 0 {
        (1.0 - x) * (1.0 - y)
    } else if i == 1 {
        x * (1.0 - y)
    } else if i == 2 {
        x * y
    } else if i == 3 {
        (1.0 - x) * y
    } else {
        panic!()
    }
}

#[cfg(feature = "3d")]
pub fn basis(
    point: [f64; DIM],
    i: usize
) -> f64 {
    let x = point[0];
    let y = point[1];
    let z = point[2];
    if i == 0 {
        (1.0 - x) * (1.0 - y) * (1.0 - z)
    } else if i == 1 {
        x * (1.0 - y) * (1.0 - z)
    } else if i == 2 {
        x * y  * (1.0 - z)
    } else if i == 3 {
        (1.0 - x) * y  * (1.0 - z)
    } else if i == 4 {
        (1.0 - x) * (1.0 - y) * z
    } else if i == 5 {
        x * (1.0 - y) * z
    } else if i == 6 {
        x * y  * z
    } else if i == 7 {
        (1.0 - x) * y  * z
    } else {
        panic!()
    }
}





