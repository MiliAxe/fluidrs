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
    simulator.run(&mut rl, thread);
}
