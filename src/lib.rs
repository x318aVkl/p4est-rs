

pub mod grid;
pub mod basetree;


pub mod consts {
    pub use p4est_sys::consts::DIM;
}

pub mod env {
    use std::os::raw::c_void;
    use mpi::{raw::AsRaw, topology::SimpleCommunicator};
    use p4est_sys::sc_init;


    pub fn initialize(mpi_comm: &SimpleCommunicator) {
        unsafe {
            sc_init(mpi_comm.as_raw().0 as *mut c_void as p4est_sys::MPI_Comm, 1, 1, None, p4est_sys::SC_LP_ERROR as i32);
        }
        unsafe {
            p4est_sys::p4est_init(None, p4est_sys::SC_LP_ERROR as i32);
        }
    }

}
