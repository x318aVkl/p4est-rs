use mpi::traits::Communicator;
use p4est::grid::{Grid, nodes::NodeNumbering};

fn main() {

    let universe = mpi::initialize().unwrap();
    let world = universe.world();

    p4est::env::initialize(&world);

    let mut grid = Grid::<()>::new_unitsquare(world.duplicate());
    //let mut grid = Grid::<()>::from_su2(std::fs::File::open("test.su2").unwrap(), world.duplicate()).unwrap();

    for _ in 1..=2 {
        grid.refine_uniform();
        grid.partition();
    }
    grid.refine(|cell| {
        cell.corner(0)[0] < 0.5
    });
    grid.partition();

    let nodes = NodeNumbering::new(&grid);

    println!("Grid len = {} / {}", grid.local_len(), grid.global_len());

    grid.map_cells(|cell| {
        let cell_nodes = nodes.cell_nodes(&cell);
        println!("rank {} cell {} {}, nodes: {:?}", world.rank(), cell.local_id, cell.global_id, cell_nodes);
    });

    println!("");

    grid.map_faces(|face| {
        let face_nodes = nodes.face_nodes(&face);
        println!("rank {} face {} : {:?} {:?}, nodes {:?}", world.rank(), face.id, (face.cell0.local_id, face.cell0.is_ghost), match face.cell1 {Some(c) => (c.local_id as i32, c.is_ghost), None => (-1_i32, false)}, face_nodes);
    });

}
