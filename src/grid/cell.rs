
use p4est_sys::consts::{CELL_CORNERS, DIM};



#[derive(Debug)]
pub struct Cell<'a, T> { 
    pub data: Option<&'a T>,
    pub local_id: usize,
    pub global_id: usize,
    pub level: u8,
    pub(super) corners: [[f64; DIM]; CELL_CORNERS],
}



impl<'a, T> Cell<'a, T> {

    pub fn corner(&self, i: usize) -> [f64; DIM] {
        self.corners[i]
    }

}

