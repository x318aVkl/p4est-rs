use std::io::{Write, BufWriter};

use mpi::traits::Communicator;
use p4est::grid::{Grid, nodes::NodeNumbering};
use p4est_sys::consts::CELL_CORNERS;

fn main() -> Result<(), std::io::Error> {

    let universe = mpi::initialize().unwrap();
    let world = universe.world();

    p4est::env::initialize(&world);

    //let mut grid = Grid::<()>::new_unitsquare(world.duplicate());
    let mut grid = Grid::<()>::from_su2(std::fs::File::open("test.su2").unwrap(), world.duplicate()).unwrap();

    for _ in 1..=3 {
        grid.refine_uniform();
        grid.partition();
    }
    grid.refine(|cell| {
        let mut c = 0.;
        for i in 0..CELL_CORNERS {
            c += cell.corner(i)[0];
        }
        c /= CELL_CORNERS as f64;
        c < 0.0
    });
    grid.partition();

    let nodes = NodeNumbering::new(&grid);

    //println!("Grid len = {} / {}", grid.local_len(), grid.global_len());

    //grid.map_cells(|cell| {
    //    let cell_nodes = nodes.cell_nodes(&cell);
    //    println!("rank {} cell {} {}, nodes: {:?}", world.rank(), cell.local_id, cell.global_id, cell_nodes);
    //});

    let file = std::fs::File::create(format!("data_{}.py", world.rank()).as_str()).unwrap();
    let mut writer = BufWriter::new(file);

    write!(writer, "data = [\n").unwrap();

    grid.map_faces(|face| {
        let face_nodes = nodes.face_nodes(&face);

        write!(writer, "    [").unwrap();
        for i in 0..face_nodes.len() {
            let v = face.corner(i);
            write!(writer, "{:?}", v).unwrap();
            if (i+1) != face_nodes.len() {
                write!(writer, ", ").unwrap();
            }
        }
        write!(writer, "],\n").unwrap();
    });

    writeln!(writer, "]").unwrap();

    Ok(())
}
