

pub use basis::*;


pub mod basis {
    use p4est_sys::consts::DIM;

    const Q_POINTS: &[&[f64]] = &[
        &[
            -1.0, 1.0,
        ],
        &[
            -1.0, 0.0, 1.0,
        ],
        &[
            -1.0, -1.0/3.0, 1.0/3.0, 1.0,
        ],
    ];

    const W_CST: &[&[f64]] = &[
        &[-0.5, 0.5],
        &[0.5, -1.0, 0.5],
        &[-0.5625, 1.6875, -1.6875, 0.5625],
    ];

    pub fn basis(
        point: [f64; DIM],
        i: usize,
        order: usize,
    ) -> f64 {
        
        let mut prod = 1.0;
        
        let mut k = i;
        //println!("{}", i);
        for dimi in 0..DIM {
            let kdim = k % (order + 1);
            //println!("   {}", kdim);

            prod *= basis_1d(point[dimi], kdim, order);

            k /= order + 1;
        }

        prod
    }

    pub fn basis_1d(
        x: f64,
        i: usize,
        order: usize,
    ) -> f64 {

        let q_points = Q_POINTS[order - 1];

        let dx = x - q_points[i];
        if dx.abs() < 1e-14 {
            return 1.0;
        }

        let mut prod = 1.0;
        for j in 0..q_points.len() {
            prod *= x - q_points[j];
        }

        if prod.abs() < 1e-14 {
            return 0.0;
        }

        prod * W_CST[order - 1][i] / dx
        
    }

}



