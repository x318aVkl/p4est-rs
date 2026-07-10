use mpi::traits::Communicator;
use p4est::grid::Grid;

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

    println!("Grid len = {} / {}", grid.local_len(), grid.global_len());

    grid.map_cells(|cell| {
        println!("rank {} cell {} {} {:?}", world.rank(), cell.local_id, cell.global_id, cell.corner(0));
    });

    println!("");

    grid.map_faces(|face| {
        println!("rank {} face {} : {} {}", world.rank(), face.id, face.cell0.local_id, match face.cell1 {Some(c) => c.local_id as i32, None => -1});
    });

}
