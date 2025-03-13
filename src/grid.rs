use raylib::prelude::*;

pub struct Velocity {
    pub x: f32,
    pub y: f32,
}

pub struct Cell {
    pub color: raylib::prelude::Color,
    pub velocity_0: Velocity,
    pub velocity: Velocity,
    pub density: f32,
    pub density_0: f32,
}

pub struct Grid {
    pub cells: Vec<Vec<Cell>>,
    pub scale: usize,
    pub size: usize,
}

pub enum ActionType {
    VelX,
    VelY,
    Density,
}

impl Grid {
    pub fn add_velocity(&mut self, x: usize, y: usize, amount: Vector2) {
        if x >= self.size || y >= self.size {
            return;
        }
        let old_velocity = &mut self.cells[y][x].velocity;
        old_velocity.x += amount.x;
        old_velocity.y += amount.y;
    }

    pub fn add_density(&mut self, x: usize, y: usize, amount: f32) {
        if x >= self.size || y >= self.size {
            return;
        }
        self.cells[y][x].density += amount;
    }

    pub fn update_cell_colors(&mut self) {
        for row in self.cells.iter_mut() {
            for cell in row.iter_mut() {
                let color_intensity = cell.density.clamp(0.0, 255.0) as u8;

                cell.color = Color {
                    r: color_intensity,
                    g: color_intensity,
                    b: color_intensity,
                    a: 255,
                };
            }
        }
    }

    pub fn set_bnd_generic<FGet, FSet>(&mut self, action_type: ActionType, get: FGet, set: FSet)
    where
        FGet: Fn(&Cell) -> f32,
        FSet: Fn(&mut Cell, f32),
    {
        // I am totally aware of the fact that this function is super
        // unreadable. This is merely a copy paste of mike ash's code
        // for my own implementation
        for x in 1..(self.size - 1) {
            let left_value = get(&self.cells[1][x]);
            let right_value = get(&self.cells[self.size - 2][x]);
            match action_type {
                ActionType::VelX => {
                    set(&mut self.cells[0][x], left_value);
                    set(&mut self.cells[self.size - 1][x], right_value);
                }
                ActionType::VelY => {
                    set(&mut self.cells[0][x], -left_value);
                    set(&mut self.cells[self.size - 1][x], -right_value);
                }
                ActionType::Density => {
                    set(&mut self.cells[0][x], left_value);
                    set(&mut self.cells[self.size - 1][x], right_value);
                }
            }
        }

        for y in 1..(self.size - 1) {
            let top_value = get(&self.cells[y][1]);
            let bottom_value = get(&self.cells[y][self.size - 2]);
            match action_type {
                ActionType::VelX => {
                    set(&mut self.cells[y][0], -top_value);
                    set(&mut self.cells[y][self.size - 1], -bottom_value);
                }
                ActionType::VelY => {
                    set(&mut self.cells[y][0], top_value);
                    set(&mut self.cells[y][self.size - 1], bottom_value);
                }
                ActionType::Density => {
                    set(&mut self.cells[y][0], top_value);
                    set(&mut self.cells[y][self.size - 1], bottom_value);
                }
            }
        }

        let top_left = 0.5 * (get(&mut self.cells[0][1]) + get(&mut self.cells[1][0]));
        let top_right =
            0.5 * (get(&mut self.cells[1][self.size - 1]) + get(&mut self.cells[0][self.size - 2]));
        let bottom_left =
            0.5 * (get(&mut self.cells[self.size - 2][0]) + get(&mut self.cells[self.size - 1][1]));
        let bottom_right = 0.5
            * (get(&mut self.cells[self.size - 1][self.size - 2])
                + get(&mut self.cells[self.size - 2][self.size - 1]));

        set(&mut self.cells[0][0], top_left);
        set(&mut self.cells[0][self.size - 1], top_right);
        set(&mut self.cells[self.size - 1][0], bottom_left);
        set(&mut self.cells[self.size - 1][self.size - 1], bottom_right);
    }

    pub fn lin_solve_generic<FSet, FGet, FPrevGet>(
        &mut self,
        action_type: ActionType,
        get: FGet,
        prev_get: FPrevGet,
        set: FSet,
        a: f32,
        c: f32,
    ) where
        FSet: Fn(&mut Cell, f32),
        FGet: Fn(&Cell) -> f32,
        FPrevGet: Fn(&Cell) -> f32,
    {
        let c_recip = 1.0 / c;
        for _k in 0..super::config::ITER {
            for y in 1..(self.size - 1) {
                for x in 1..(self.size - 1) {
                    let previous_element = prev_get(&self.cells[y][x]);

                    let left_neighbor = get(&self.cells[y][x - 1]);
                    let right_neighbor = get(&self.cells[y][x + 1]);
                    let bottom_neighbor = get(&self.cells[y - 1][x]);
                    let top_neighbor = get(&self.cells[y + 1][x]);

                    set(
                        &mut self.cells[y][x],
                        (previous_element
                            + a * (left_neighbor
                                + right_neighbor
                                + bottom_neighbor
                                + top_neighbor))
                            * c_recip,
                    );
                }
            }
        }

        self.set_bnd_generic(action_type, get, set);
    }
    pub fn new() -> Grid {
        let size = super::config::SIZE;
        let scale = super::config::SCALE;
        let grid_count = size / scale;
        let mut grid: Vec<Vec<Cell>> = Vec::new();

        for _i in 0..grid_count {
            let mut row = Vec::new();
            for _j in 0..grid_count {
                row.push(Cell {
                    color: raylib::prelude::Color::BLACK,
                    velocity: Velocity { x: 0.0, y: 0.0 },
                    velocity_0: Velocity { x: 0.0, y: 0.0 },
                    density: 0.0,
                    density_0: 0.0,
                });
            }
            grid.push(row);
        }

        Grid {
            cells: grid,
            scale,
            size: grid_count,
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
