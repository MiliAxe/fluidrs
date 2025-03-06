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

    let mut simulator = grid::Simulator::new(config::VISC, config::DIFF, config::DT);
    while !rl.window_should_close() {
        simulator.update_mouse_density(&rl);
        simulator.step();
        let mut d = rl.begin_drawing(&thread);
        simulator.draw(&mut d);
    }
}
