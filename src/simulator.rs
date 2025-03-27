use crate::grid::{ActionType, Cell, Grid};
use colorgrad::BasisGradient;
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
    brush_size: usize,
    brush_density: f32,
    brush_type: BrushType,
    gradient_index: usize,
    last_mouse_pos: Option<Vector2>,
}

#[derive(PartialEq, Eq)]
pub enum BrushType {
    Filled,
    Outline,
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
            brush_size: 1,
            brush_density: 100.0,
            gradient_index: 0,
            brush_type: BrushType::Filled,
            last_mouse_pos: Option::<Vector2>::None,
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

    fn draw_filled_circle(&mut self, cell_x: usize, cell_y: usize) {
        let radius = self.brush_size as i32; // Or calculate an actual circle radius

        for dy in -radius..=radius {
            for dx in -radius..=radius {
                // Only apply if within the circle of given radius. You can adjust the condition as needed.
                if dx * dx + dy * dy <= radius * radius {
                    let draw_x = (cell_x as i32 + dx).max(0) as usize;
                    let draw_y = (cell_y as i32 + dy).max(0) as usize;
                    self.grid.add_density(draw_x, draw_y, self.brush_density);
                }
            }
        }
    }

    fn draw_outline_circle(&mut self, cell_x: usize, cell_y: usize) {
        let mut x: i32 = 0;
        let mut y: i32 = -(self.brush_size as i32);

        let mut p: i32 = -(self.brush_size as i32);

        let calculcate_draw_pos =
            |cell_pos: usize, cell_dif: i32| (cell_pos as i32 + cell_dif).max(0) as usize;

        while x < -y {
            if p > 0 {
                y += 1;
                p += 2 * (x + y) + 1;
            } else {
                p += 2 * x + 1;
            }

            self.grid.add_density(
                calculcate_draw_pos(cell_x, x),
                calculcate_draw_pos(cell_y, y),
                self.brush_density,
            );
            self.grid.add_density(
                calculcate_draw_pos(cell_x, x),
                calculcate_draw_pos(cell_y, -y),
                self.brush_density,
            );
            self.grid.add_density(
                calculcate_draw_pos(cell_x, -x),
                calculcate_draw_pos(cell_y, y),
                self.brush_density,
            );
            self.grid.add_density(
                calculcate_draw_pos(cell_x, -x),
                calculcate_draw_pos(cell_y, -y),
                self.brush_density,
            );

            self.grid.add_density(
                calculcate_draw_pos(cell_x, y),
                calculcate_draw_pos(cell_y, x),
                self.brush_density,
            );
            self.grid.add_density(
                calculcate_draw_pos(cell_x, y),
                calculcate_draw_pos(cell_y, -x),
                self.brush_density,
            );
            self.grid.add_density(
                calculcate_draw_pos(cell_x, -y),
                calculcate_draw_pos(cell_y, x),
                self.brush_density,
            );
            self.grid.add_density(
                calculcate_draw_pos(cell_x, -y),
                calculcate_draw_pos(cell_y, -x),
                self.brush_density,
            );

            x += 1;
        }
    }

    pub fn update_mouse_density(&mut self) {
        if self
            .raylib_handle
            .is_mouse_button_down(MouseButton::MOUSE_BUTTON_LEFT)
        {
            let mouse_pos = self.raylib_handle.get_mouse_position();

            let cell_x = (mouse_pos.x / self.grid.scale as f32) as usize;
            let cell_y = (mouse_pos.y / self.grid.scale as f32) as usize;

            match self.brush_type {
                BrushType::Filled => self.draw_filled_circle(cell_x, cell_y),
                BrushType::Outline => self.draw_outline_circle(cell_x, cell_y),
            }
        }
    }

    fn available_gradients() -> Vec<(&'static str, BasisGradient)> {
        vec![
            ("Inferno", colorgrad::preset::inferno()),
            ("Viridis", colorgrad::preset::viridis()),
            ("Plasma", colorgrad::preset::plasma()),
            ("Magma", colorgrad::preset::magma()),
            ("Gray", colorgrad::preset::greys()),
            ("Blues", colorgrad::preset::blues()),
            ("Greens", colorgrad::preset::greens()),
            ("Reds", colorgrad::preset::reds()),

        ]
    } 

    fn handle_ui(&mut self) -> &mut imgui::Ui {
        self.imgui_renderer
            .update(&mut self.imgui_context, &mut self.raylib_handle);
        let ui = self.imgui_context.frame();

        let mut viscosity_slider_value = (self.viscosity * 10000.0) as i32;
        let mut diffusion_slider_value = (self.diffusion * 100000.0) as i32;

        let gradients = Simulator::available_gradients();
        let gradient_names: Vec<&str> = gradients.iter().map(|(name, _)| *name).collect();

        ui.window("Settings")
            .size([300.0, 100.0], imgui::Condition::FirstUseEver)
            .build(|| {
                ui.slider("Viscosity", 0, 100, &mut viscosity_slider_value);
                ui.slider("Diffusion", 0, 100, &mut diffusion_slider_value);
                ui.slider("Brush Size", 1, 100, &mut self.brush_size);
                ui.slider("Brush Density", 1.0, 1000.0, &mut self.brush_density);

                let mut current_brush = match self.brush_type {
                    BrushType::Filled => 0,
                    BrushType::Outline => 1,
                };
                if ui.combo(
                    "Brush Type",
                    &mut current_brush,
                    &["Filled Circle", "Circle Outline"],
                    |item| std::borrow::Cow::Borrowed(item),
                ) {
                    self.brush_type = match current_brush {
                        0 => BrushType::Filled,
                        1 => BrushType::Outline,
                        _ => panic!("Invalid brush type"),
                    };
                }

                if ui.combo("Gradient", &mut self.gradient_index, &gradient_names, |item| std::borrow::Cow::Borrowed(item)) {
                    self.grid.color_grad = gradients[self.gradient_index].1.clone();
                }});

        self.viscosity = viscosity_slider_value as f32 / 10000.0;
        self.diffusion = diffusion_slider_value as f32 / 100000.0;

        ui
    }

    pub fn update_mouse_velocity(&mut self) {
        if self
            .raylib_handle
            .is_mouse_button_down(MouseButton::MOUSE_BUTTON_RIGHT)
        {
            if self.last_mouse_pos.is_none() {
                self.last_mouse_pos = Some(self.raylib_handle.get_mouse_position());
            }

            let current_mouse_pos = self.raylib_handle.get_mouse_position();
            let velocity_vector = Vector2::new(
                (current_mouse_pos.x - self.last_mouse_pos.unwrap().x) * 2.0,
                (current_mouse_pos.y - self.last_mouse_pos.unwrap().y) * 2.0,
            );

            let cell_x = (current_mouse_pos.x / self.grid.scale as f32) as usize;
            let cell_y = (current_mouse_pos.y / self.grid.scale as f32) as usize;

            // self.grid.add_velocity(cell_x, cell_y, velocity_vector);
            for i in 0..self.brush_size / 2 {
                for j in 0..self.brush_size / 2 {
                    self.grid.add_velocity(
                        cell_x.saturating_add(i),
                        cell_y.saturating_add(j),
                        velocity_vector,
                    );
                    self.grid.add_velocity(
                        cell_x.saturating_sub(i),
                        cell_y.saturating_sub(j),
                        velocity_vector,
                    );
                    self.grid.add_velocity(
                        cell_x.saturating_add(i),
                        cell_y.saturating_sub(j),
                        velocity_vector,
                    );
                    self.grid.add_velocity(
                        cell_x.saturating_sub(i),
                        cell_y.saturating_add(j),
                        velocity_vector,
                    );
                }
            }

            self.last_mouse_pos = Some(current_mouse_pos);
        } else {
            self.last_mouse_pos = Option::<Vector2>::None;
        }
    }

    fn handle_mouse_simulator(&mut self) {
        self.update_mouse_density();
        self.update_mouse_velocity();
    }

    pub fn run(&mut self) {
        while !self.raylib_handle.window_should_close() {
            self.step();

            let want_capture_mouse = {
                let ui = self.handle_ui();
                ui.io().want_capture_mouse
            };

            if !want_capture_mouse {
                self.handle_mouse_simulator();
            }

            self.draw();
        }
    }
}
