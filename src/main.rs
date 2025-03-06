mod config;
mod grid;

fn main() {
    let mut simulator = grid::Simulator::new(config::VISC, config::DIFF, config::DT);
    simulator.run();
}
