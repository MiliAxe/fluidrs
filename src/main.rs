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

    let mut simulator = grid::Simulator::new(5.0, 2.0, 10.0);
    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        simulator.step();
        simulator.draw(&mut d);
    }
}
