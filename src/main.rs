mod app;
mod config;
mod grid;
mod simulator;

use app::App;

fn main() {
    let mut app = App::new();
    app.run();
}
