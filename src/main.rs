mod config;
mod grid;

fn main() {
    let (mut rl, thread) = raylib::init()
        .size(config::SIZE, config::SIZE)
        .title("fluidrs")
        .build();

    let sample_grid = grid::Grid::new(config::SIZE, config::SCALE);
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        sample_grid.draw(&mut d);
    }
}
