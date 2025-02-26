mod config;
mod grid;

fn main() {
    let (mut rl, thread) = raylib::init()
        .size(
            config::SIZE.try_into().unwrap(),
            config::SIZE.try_into().unwrap(),
        )
        .title("fluidrs")
        .build();

    let simulator = grid::Simulator::new(0.0, 0.0, 0.0);
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        simulator.draw(&mut d);
    }
}
