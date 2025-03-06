use crate::grid::{ActionType, Cell, Grid, Velocity};
use imgui::Context;
use raylib::prelude::*;
use raylib_imgui_rs::Renderer;

enum ProjectAction {
    ToPrevious,
    ToCurrent,
}
pub struct Simulator {
    grid: Grid,
    viscosity: f32,
    diffusion: f32,
    deltatime: f32,
    raylib_handle: RaylibHandle,
    raylib_thread: RaylibThread,
    imgui_context: imgui::Context,
    imgui_renderer: raylib_imgui_rs::Renderer,
}

impl Simulator {
    pub fn new(viscosity: f32, diffusion: f32, deltatime: f32) -> Simulator {
        let (mut raylib_handle, raylib_thread) = raylib::init()
            .size(
                super::config::SIZE.try_into().unwrap(),
                super::config::SIZE.try_into().unwrap(),
            )
            .title("fluidrs")
            .build();

        let mut imgui_context = Context::create();

        let imgui_renderer =
            Renderer::create(&mut imgui_context, &mut raylib_handle, &raylib_thread);

        Simulator {
            grid: Grid::new(),
            viscosity,
            diffusion,
            deltatime,
            raylib_handle,
            raylib_thread,
            imgui_context,
            imgui_renderer,
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

    pub fn draw(&mut self) {
        let mut d = self.raylib_handle.begin_drawing(&self.raylib_thread);
        self.grid.draw(&mut d);

        self.imgui_renderer.render(&mut self.imgui_context, &mut d);
    }

    pub fn update_mouse_density(&mut self) {
        if self
            .raylib_handle
            .is_mouse_button_down(MouseButton::MOUSE_BUTTON_LEFT)
        {
            let mouse_pos = self.raylib_handle.get_mouse_position();

            let cell_x = (mouse_pos.x / self.grid.scale as f32) as usize;
            let cell_y = (mouse_pos.y / self.grid.scale as f32) as usize;

            // add density to a 10x10 square
            for i in 0..3 {
                for j in 0..3 {
                    self.grid.add_density(cell_x + i, cell_y + j, 100.0);
                    self.grid.add_density(cell_x - i, cell_y - j, 100.0);
                    self.grid.add_density(cell_x + i, cell_y - j, 100.0);
                    self.grid.add_density(cell_x - i, cell_y + j, 100.0);
                }
            }
        }
    }

    pub fn update_mouse_velocity(&mut self, last_mouse_pos: &mut Option<Vector2>) {
        if self
            .raylib_handle
            .is_mouse_button_down(MouseButton::MOUSE_BUTTON_LEFT)
        {
            if last_mouse_pos.is_none() {
                *last_mouse_pos = Some(self.raylib_handle.get_mouse_position());
            }

            let current_mouse_pos = self.raylib_handle.get_mouse_position();
            let velocity_vector = Vector2::new(
                (current_mouse_pos.x - last_mouse_pos.unwrap().x) * 5.0,
                (current_mouse_pos.y - last_mouse_pos.unwrap().y) * 5.0,
            );

            let cell_x = (current_mouse_pos.x / self.grid.scale as f32) as usize;
            let cell_y = (current_mouse_pos.y / self.grid.scale as f32) as usize;

            self.grid.add_velocity(cell_x, cell_y, velocity_vector);

            *last_mouse_pos = Some(current_mouse_pos);
        } else {
            *last_mouse_pos = Option::<Vector2>::None;
        }
    }

    pub fn run(&mut self) {
        let mut last_mouse_pos = Option::<Vector2>::None;

        while !self.raylib_handle.window_should_close() {
            self.update_mouse_density();
            self.update_mouse_velocity(&mut last_mouse_pos);

            self.step();
            self.draw();
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
