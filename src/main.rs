mod config;
mod grid;
mod simulator;

use simulator::Simulator;

fn main() {
    let mut simulator = Simulator::new(config::VISC, config::DIFF, config::DT);
    simulator.run();
}
