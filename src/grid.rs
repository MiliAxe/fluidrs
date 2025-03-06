use raylib::prelude::*;

struct Velocity {
    x: f32,
    y: f32,
}

pub struct Cell {
    color: raylib::prelude::Color,
    velocity_0: Velocity,
    velocity: Velocity,
    density: f32,
    density_0: f32,
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

enum ActionType {
    VelX,
    VelY,
    Density,
}

impl Grid {
    fn add_velocity(&mut self, x: usize, y: usize, amount: Velocity) {
        let old_velocity = &mut self.cells[y][x].velocity;
        old_velocity.x += amount.x;
        old_velocity.y += amount.y;
    }

    fn add_density(&mut self, x: usize, y: usize, amount: f32) {
        self.cells[y][x].density += amount;
    }

    fn update_cell_colors(&mut self) {
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

    fn set_bnd_generic<FGet, FSet>(&mut self, action_type: ActionType, get: FGet, set: FSet)
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

    fn lin_solve_generic<FSet, FGet, FPrevGet>(
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
}

enum ProjectAction {
    ToPrevious,
    ToCurrent,
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

    fn diffuse(&mut self, action_type: ActionType) {
        let spread_factor = match action_type {
            ActionType::VelX | ActionType::VelY => self.viscosity,
            ActionType::Density => self.diffusion,
        };

        let a = self.deltatime
            * spread_factor
            * (self.grid.size as f32 - 2.0)
            * (self.grid.size as f32 - 2.0);

        let get = match action_type {
            ActionType::VelX => |cell: &Cell| -> f32 { cell.velocity_0.x },
            ActionType::VelY => |cell: &Cell| -> f32 { cell.velocity_0.y },
            ActionType::Density => |cell: &Cell| -> f32 { cell.density_0 },
        };

        let prev_get = match action_type {
            ActionType::VelX => |cell: &Cell| -> f32 { cell.velocity.x },
            ActionType::VelY => |cell: &Cell| -> f32 { cell.velocity.y },
            ActionType::Density => |cell: &Cell| -> f32 { cell.density },
        };
        let set = match action_type {
            ActionType::VelX => |cell: &mut Cell, value: f32| cell.velocity_0.x = value,
            ActionType::VelY => |cell: &mut Cell, value: f32| cell.velocity_0.y = value,
            ActionType::Density => |cell: &mut Cell, value: f32| cell.density_0 = value,
        };

        self.grid
            .lin_solve_generic(action_type, get, prev_get, set, a, 1.0 + 6.0 * a);
    }

    fn project(&mut self, action_type: ProjectAction) {
        let get_x = match action_type {
            ProjectAction::ToPrevious => |cell: &Cell| cell.velocity_0.x,
            ProjectAction::ToCurrent => |cell: &Cell| cell.velocity.x,
        };
        let get_y = match action_type {
            ProjectAction::ToPrevious => |cell: &Cell| cell.velocity_0.y,
            ProjectAction::ToCurrent => |cell: &Cell| cell.velocity.y,
        };
        let set_x = match action_type {
            ProjectAction::ToPrevious => |cell: &mut Cell, value: f32| cell.velocity_0.x = value,
            ProjectAction::ToCurrent => |cell: &mut Cell, value: f32| cell.velocity.x = value,
        };
        let set_y = match action_type {
            ProjectAction::ToPrevious => |cell: &mut Cell, value: f32| cell.velocity_0.y = value,
            ProjectAction::ToCurrent => |cell: &mut Cell, value: f32| cell.velocity.y = value,
        };

        let get_div = match action_type {
            ProjectAction::ToPrevious => |cell: &Cell| cell.velocity.x,
            ProjectAction::ToCurrent => |cell: &Cell| cell.velocity_0.x,
        };
        let get_pressure = match action_type {
            ProjectAction::ToPrevious => |cell: &Cell| cell.velocity.y,
            ProjectAction::ToCurrent => |cell: &Cell| cell.velocity_0.y,
        };
        let set_pressure = match action_type {
            ProjectAction::ToPrevious => |cell: &mut Cell, value: f32| cell.velocity.y = value,
            ProjectAction::ToCurrent => |cell: &mut Cell, value: f32| cell.velocity_0.y = value,
        };
        let set_div = match action_type {
            ProjectAction::ToPrevious => |cell: &mut Cell, value: f32| cell.velocity.x = value,
            ProjectAction::ToCurrent => |cell: &mut Cell, value: f32| cell.velocity_0.x = value,
        };

        let size = self.grid.size as f32;

        for y in 1..self.grid.size - 1 {
            for x in 1..self.grid.size - 1 {
                let div_value = -0.5
                    * (get_x(&self.grid.cells[y][x + 1]) - get_x(&self.grid.cells[y][x - 1])
                        + get_y(&self.grid.cells[y + 1][x])
                        - get_y(&self.grid.cells[y - 1][x]))
                    / size;

                set_div(&mut self.grid.cells[y][x], div_value);
                set_pressure(&mut self.grid.cells[y][x], 0.0);
            }
        }

        self.grid
            .set_bnd_generic(ActionType::Density, get_div, set_div);
        self.grid
            .set_bnd_generic(ActionType::Density, get_pressure, set_pressure);
        self.grid.lin_solve_generic(
            ActionType::Density,
            get_pressure,
            get_div,
            set_pressure,
            1.0,
            6.0,
        );

        for y in 1..self.grid.size - 1 {
            for x in 1..self.grid.size - 1 {
                let current_x = get_x(&self.grid.cells[y][x]);
                let current_y = get_y(&self.grid.cells[y][x]);

                let x_difference = 0.5
                    * (get_pressure(&self.grid.cells[y][x + 1])
                        - get_pressure(&self.grid.cells[y][x - 1]))
                    * size;
                let y_difference = 0.5
                    * (get_pressure(&self.grid.cells[y + 1][x])
                        - get_pressure(&self.grid.cells[y - 1][x]))
                    * size;

                set_x(&mut self.grid.cells[y][x], current_x - x_difference);
                set_y(&mut self.grid.cells[y][x], current_y - y_difference);
            }
        }

        self.grid.set_bnd_generic(ActionType::VelX, get_x, set_x);
        self.grid.set_bnd_generic(ActionType::VelY, get_y, set_y);
    }

    fn advect(&mut self, action_type: ActionType) {
        let size = self.grid.size as f32;
        let dt = self.deltatime;

        let get_prev = match action_type {
            ActionType::VelX => |cell: &Cell| cell.velocity_0.x,
            ActionType::VelY => |cell: &Cell| cell.velocity_0.y,
            ActionType::Density => |cell: &Cell| cell.density_0,
        };
        let set_current = match action_type {
            ActionType::VelX => |cell: &mut Cell, value: f32| cell.velocity.x = value,
            ActionType::VelY => |cell: &mut Cell, value: f32| cell.velocity.y = value,
            ActionType::Density => |cell: &mut Cell, value: f32| cell.density = value,
        };
        let get_vel_x = match action_type {
            ActionType::VelX => |cell: &Cell| cell.velocity_0.x,
            ActionType::VelY => |cell: &Cell| cell.velocity_0.x,
            ActionType::Density => |cell: &Cell| cell.velocity.x,
        };
        let get_vel_y = match action_type {
            ActionType::VelX => |cell: &Cell| cell.velocity_0.y,
            ActionType::VelY => |cell: &Cell| cell.velocity_0.y,
            ActionType::Density => |cell: &Cell| cell.velocity.y,
        };

        let mut new_values = vec![vec![0.0; self.grid.size]; self.grid.size];

        (1..self.grid.size - 1).for_each(|y| {
            for x in 1..self.grid.size - 1 {
                let cell = &self.grid.cells[y][x];
                let vx = get_vel_x(cell);
                let vy = get_vel_y(cell);

                let mut x_prev = x as f32 - dt * vx * (size - 2.0);
                let mut y_prev = y as f32 - dt * vy * (size - 2.0);

                x_prev = x_prev.max(0.5).min(size - 1.5);
                y_prev = y_prev.max(0.5).min(size - 1.5);

                let x0 = x_prev.floor() as usize;
                let x1 = (x0 + 1).min(self.grid.size - 1);
                let y0 = y_prev.floor() as usize;
                let y1 = (y0 + 1).min(self.grid.size - 1);

                let sx = x_prev - x0 as f32;
                let sy = y_prev - y0 as f32;
                let s0 = 1.0 - sx;
                let t0 = 1.0 - sy;

                let value = s0 * t0 * get_prev(&self.grid.cells[y0][x0])
                    + sx * t0 * get_prev(&self.grid.cells[y0][x1])
                    + s0 * sy * get_prev(&self.grid.cells[y1][x0])
                    + sx * sy * get_prev(&self.grid.cells[y1][x1]);

                new_values[y][x] = value;
            }
        });

        (1..self.grid.size - 1).for_each(|y| {
            for x in 1..self.grid.size - 1 {
                set_current(&mut self.grid.cells[y][x], new_values[y][x]);
            }
        });

        let get_current = match action_type {
            ActionType::VelX => |cell: &Cell| cell.velocity.x,
            ActionType::VelY => |cell: &Cell| cell.velocity.y,
            ActionType::Density => |cell: &Cell| cell.density,
        };

        self.grid
            .set_bnd_generic(action_type, get_current, set_current);
    }

    pub fn step(&mut self) {
        self.diffuse(ActionType::VelX);
        self.diffuse(ActionType::VelY);

        self.project(ProjectAction::ToPrevious);

        self.advect(ActionType::VelX);
        self.advect(ActionType::VelY);

        self.project(ProjectAction::ToCurrent);

        self.diffuse(ActionType::Density);
        self.advect(ActionType::Density);

        self.grid.update_cell_colors();
    }


    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        self.grid.draw(d);
    }

    pub fn update_mouse_density(&mut self, rl: &RaylibHandle) {
        if rl.is_mouse_button_down(MouseButton::MOUSE_BUTTON_LEFT) {
            let mouse_pos = rl.get_mouse_position();

            let cell_x = (mouse_pos.x / self.grid.scale as f32) as usize;
            let cell_y = (mouse_pos.y / self.grid.scale as f32) as usize;

            // add density to a 10x10 square
            for i in 0..5 {
                for j in 0..5 {
                    self.grid.add_density(cell_x + i, cell_y + j, 1000.0);
                }
            }
        }
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
