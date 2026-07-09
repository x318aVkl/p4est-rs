use p4est_sys::consts::DIM;





pub fn reorder_raw_vtk_element(
    element: Vec<usize>
) -> Vec<usize> {
    // gmsh outputs the nodes in a bullit order
    // we need to reorder nodes, annoying

    let npoints = element.len();

    let order = (((npoints as f64).powf(1.0 / (DIM as f64))).round() as usize) - 1;

    let mut result = element.clone();

    if order == 1 {
        // simple reordering
        result[0] = element[0];
        result[1] = element[1];
        result[2] = element[3];
        result[3] = element[2];

        #[cfg(feature = "3d")]
        { 
            result[4] = element[4];
            result[5] = element[5];
            result[6] = element[7];
            result[7] = element[6];
        }

    } else if order == 2 {
        #[cfg(feature = "3d")]
        {
            result[0] = element[0];
            result[1] = element[8];
            result[2] = element[1];
            result[3] = element[11];
            result[4] = element[24];
            result[5] = element[9];
            result[6] = element[3];
            result[7] = element[10];
            result[8] = element[2];
            result[9] = element[16];
            result[10] = element[22];
            result[11] = element[17];
            result[12] = element[20];
            result[13] = element[26];
            result[14] = element[21];
            result[15] = element[19];
            result[16] = element[23];
            result[17] = element[18];
            result[18] = element[4];
            result[19] = element[12];
            result[20] = element[5];
            result[21] = element[15];
            result[22] = element[25];
            result[23] = element[13];
            result[24] = element[7];
            result[25] = element[14];
            result[26] = element[6];
        }

    } else if order == 3 {

        #[cfg(feature = "3d")]
        {
            result[0] = element[0];
            result[1] = element[8];
            result[2] = element[9];
            result[3] = element[1];
            result[4] = element[10];
            result[5] = element[32];
            result[6] = element[35];
            result[7] = element[14];
            result[8] = element[11];
            result[9] = element[33];
            result[10] = element[34];
            result[11] = element[15];
            result[12] = element[3];
            result[13] = element[19];
            result[14] = element[18];
            result[15] = element[2];
            result[16] = element[12];
            result[17] = element[36];
            result[18] = element[37];
            result[19] = element[16];
            result[20] = element[40];
            result[21] = element[56];
            result[22] = element[57];
            result[23] = element[44];
            result[24] = element[43];
            result[25] = element[59];
            result[26] = element[58];
            result[27] = element[45];
            result[28] = element[22];
            result[29] = element[49];
            result[30] = element[48];
            result[31] = element[20];
            result[32] = element[13];
            result[33] = element[39];
            result[34] = element[38];
            result[35] = element[17];
            result[36] = element[41];
            result[37] = element[60];
            result[38] = element[61];
            result[39] = element[47];
            result[40] = element[42];
            result[41] = element[63];
            result[42] = element[62];
            result[43] = element[46];
            result[44] = element[23];
            result[45] = element[50];
            result[46] = element[51];
            result[47] = element[21];
            result[48] = element[4];
            result[49] = element[24];
            result[50] = element[25];
            result[51] = element[5];
            result[52] = element[26];
            result[53] = element[52];
            result[54] = element[53];
            result[55] = element[28];
            result[56] = element[27];
            result[57] = element[55];
            result[58] = element[54];
            result[59] = element[29];
            result[60] = element[7];
            result[61] = element[31];
            result[62] = element[30];
            result[63] = element[6];
        }

    } else {
        panic!("Unsupported tree element order {}, must be >= 1 and <= 3", order);
    }

    //println!("nodes = {:?}", result);

    result
}


