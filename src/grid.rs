use rand::{seq::IndexedRandom, Rng};
use raylib::{ffi::GetApplicationDirectory, prelude::*};
use std::ops;

struct Velocity {
    x: f32,
    y: f32,
}

pub struct Cell {
    color: raylib::prelude::Color,
    velocity: Velocity,
    velocity_0: Velocity,
    scratch: f32,
    density: f32,
}

pub struct Grid {
    cells: Vec<Vec<Cell>>,
    scale: usize,
    size: usize,
}

pub struct Simulator {
    grid: Grid,
    viscosity: f32,
    diffusion: f32,
    deltatime: f32,
}

impl ops::Add<Velocity> for Velocity {
    type Output = Velocity;

    fn add(self, rhs: Velocity) -> Self::Output {
        Velocity {
            x: self.x + rhs.x,
            y: self.x + rhs.x,
        }
    }
}

// TODO: we should somehow implement the diffuse, project and advect functions. I don't know
//       what is the best way to do it in rust. and where i should do it.

// Diffuse function basically calculates a constant, and passes a lot of parameters to a linear
// solving function...So we basically need to implement linear solve first?

// In order to implement lin solve, I will need to implement set_bnd. which sets the boundries basically
// but a huge problem here that we have is the fact that the code that mike has written has a major conflict
// with my design. and that is the fact that set_bnd in mike ash's design is using a parameter called 'b'
// this b makes the function have multiple purposes which is not really aligned to the uncle bob's rule
// of "doing one thing and doing it well". So I think the best thing to do here is to first understand
// what the function is doing and the reformatted to my rust liking.

enum ActionType {
    VelX,
    VelY,
    Density,
}

impl Grid {
    fn add_density(&mut self, x: usize, y: usize, amount: f32) {
        self.cells[y][x].density += amount;
    }

    fn add_velocity(&mut self, x: usize, y: usize, amount: Velocity) {
        let old_velocity = &mut self.cells[y][x].velocity;
        old_velocity.x += amount.x;
        old_velocity.y += amount.y;
    }

    fn set_bnd(&mut self, action_type: ActionType) {
        // I am totally aware of the fact that this function is super
        // unreadable. This is merely a copy paste of mike ash's code
        // for my own implementation
        for x in 1..self.size - 1 {
            match action_type {
                ActionType::VelX => {
                    self.cells[0][x].velocity.x = self.cells[1][x].velocity.x;
                    self.cells[self.size - 1][x].velocity.x =
                        self.cells[self.size - 2][x].velocity.x;
                }
                ActionType::VelY => {
                    self.cells[0][x].velocity.y = self.cells[1][x].velocity.y;
                    self.cells[self.size - 1][x].velocity.y =
                        -self.cells[self.size - 2][x].velocity.y;
                }
                ActionType::Density => {
                    self.cells[0][x].density = self.cells[1][x].density;
                    self.cells[self.size - 1][x].density = self.cells[self.size - 2][x].density;
                }
            }
        }

        for y in 1..self.size - 1 {
            match action_type {
                ActionType::VelX => {
                    self.cells[y][0].velocity.x = self.cells[y][1].velocity.x;
                    self.cells[y][self.size - 1].velocity.x =
                        -self.cells[y][self.size - 2].velocity.x;
                }
                ActionType::VelY => {
                    self.cells[y][0].velocity.y = self.cells[y][1].velocity.y;
                    self.cells[y][self.size - 1].velocity.y =
                        self.cells[y][self.size - 2].velocity.y;
                }
                ActionType::Density => {
                    self.cells[y][0].density = self.cells[y][1].density;
                    self.cells[y][self.size - 1].density = self.cells[y][self.size - 2].density;
                }
            }
        }

        match action_type {
            ActionType::VelX => {
                self.cells[0][0].velocity.x =
                    0.5 * self.cells[0][1].velocity.x + 0.5 * self.cells[1][0].velocity.x;
                self.cells[0][self.size - 1].velocity.x = 0.5
                    * self.cells[1][self.size - 1].velocity.x
                    + 0.5 * self.cells[0][self.size - 2].velocity.x;
                self.cells[self.size - 1][0].velocity.x = 0.5
                    * self.cells[self.size - 2][0].velocity.x
                    + 0.5 * self.cells[self.size - 1][1].velocity.x;
                self.cells[self.size - 1][self.size - 1].velocity.x = 0.5
                    * self.cells[self.size - 1][self.size - 2].velocity.x
                    + 0.5 * self.cells[self.size - 2][self.size - 1].velocity.x;
            }
            ActionType::VelY => {
                self.cells[0][0].velocity.y =
                    0.5 * self.cells[0][1].velocity.y + 0.5 * self.cells[1][0].velocity.y;
                self.cells[0][self.size - 1].velocity.y = 0.5
                    * self.cells[1][self.size - 1].velocity.y
                    + 0.5 * self.cells[0][self.size - 2].velocity.y;
                self.cells[self.size - 1][0].velocity.y = 0.5
                    * self.cells[self.size - 2][0].velocity.y
                    + 0.5 * self.cells[self.size - 1][1].velocity.y;
                self.cells[self.size - 1][self.size - 1].velocity.y = 0.5
                    * self.cells[self.size - 1][self.size - 2].velocity.y
                    + 0.5 * self.cells[self.size - 2][self.size - 1].velocity.y;
            }
            ActionType::Density => {
                self.cells[0][0].density =
                    0.5 * self.cells[0][1].density + 0.5 * self.cells[1][0].density;
                self.cells[0][self.size - 1].density = 0.5 * self.cells[1][self.size - 1].density
                    + 0.5 * self.cells[0][self.size - 2].density;
                self.cells[self.size - 1][0].density = 0.5 * self.cells[self.size - 2][0].density
                    + 0.5 * self.cells[self.size - 1][1].density;
                self.cells[self.size - 1][self.size - 1].density = 0.5
                    * self.cells[self.size - 1][self.size - 2].density
                    + 0.5 * self.cells[self.size - 2][self.size - 1].density;
            }
        }
    }

    // Used for solving an equation linearly. I don't know what is happening here really
    fn lin_solve(&mut self, action_type: ActionType, a: f32, c: f32) {
        let c_recip = 1.0 / c;
        for k in 0..super::config::ITER {
            for y in 1..self.size - 1 {
                for x in 1..self.size - 1 {}
            }
        }
    }
}

impl Simulator {
    pub fn new(viscosity: f32, diffusion: f32, deltatime: f32) -> Simulator {
        Simulator {
            grid: Grid::new(),
            viscosity,
            diffusion,
            deltatime,
        }
    }

    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        self.grid.draw(d);
    }
}

impl Grid {
    pub fn new() -> Grid {
        let size = super::config::SIZE;
        let scale = super::config::SCALE;
        let grid_count = size / scale;
        let mut grid: Vec<Vec<Cell>> = Vec::new();

        for _i in 0..grid_count {
            let mut row = Vec::new();
            for _j in 0..grid_count {
                let mut rng = rand::rng();
                let r: u8 = rng.random();
                let g: u8 = rng.random();
                let b: u8 = rng.random();

                row.push(Cell {
                    color: Color { r, g, b, a: 255 },
                    velocity: Velocity { x: 0.0, y: 0.0 },
                    velocity_0: Velocity { x: 0.0, y: 0.0 },
                    scratch: 0.0,
                    density: 0.0,
                });
            }
            grid.push(row);
        }

        Grid {
            cells: grid,
            scale,
            size,
        }
    }

    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        for (row_idx, row) in self.cells.iter().enumerate() {
            for (col_idx, value) in row.iter().enumerate() {
                let y = (row_idx as i32) * (self.scale as i32);
                let x = (col_idx as i32) * (self.scale as i32);
                d.draw_rectangle(x, y, self.scale as i32, self.scale as i32, value.color);
            }
        }
    }
}
