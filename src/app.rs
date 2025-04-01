use crate::{grid::Grid, simulator::Simulator};
use colorgrad::Gradient;
use imgui::Context;
use raylib::prelude::*;
use raylib_imgui_rs::Renderer;
use std::sync::Arc;

#[derive(PartialEq, Eq)]
pub enum BrushType {
    Filled,
    Outline,
}

pub struct App {
    simulator: Simulator,
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

impl App {
    pub fn new() -> App {
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

        let simulator = Simulator::new(crate::config::VISC, crate::config::DIFF, crate::config::DT);

        App {
            simulator,
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

    fn draw_grid(grid: &Grid, d: &mut RaylibDrawHandle) {
        for (row_idx, row) in grid.cells.iter().enumerate() {
            for (col_idx, value) in row.iter().enumerate() {
                let y = (row_idx as i32) * (grid.scale as i32);
                let x = (col_idx as i32) * (grid.scale as i32);
                d.draw_rectangle(x, y, grid.scale as i32, grid.scale as i32, value.color);
            }
        }
    }

    pub fn draw(&mut self) {
        let mut d = self.raylib_handle.begin_drawing(&self.raylib_thread);
        // self.raylib_handle.draw(&self.raylib_thread, |mut d| {
        //     App::draw_grid(&self.simulator.grid, &mut d);
        //     self.imgui_renderer.render(&mut self.imgui_context, &mut d);
        // });
        App::draw_grid(&self.simulator.grid, &mut d);
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
                    self.simulator
                        .grid
                        .add_density(draw_x, draw_y, self.brush_density);
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

            self.simulator.grid.add_density(
                calculcate_draw_pos(cell_x, x),
                calculcate_draw_pos(cell_y, y),
                self.brush_density,
            );
            self.simulator.grid.add_density(
                calculcate_draw_pos(cell_x, x),
                calculcate_draw_pos(cell_y, -y),
                self.brush_density,
            );
            self.simulator.grid.add_density(
                calculcate_draw_pos(cell_x, -x),
                calculcate_draw_pos(cell_y, y),
                self.brush_density,
            );
            self.simulator.grid.add_density(
                calculcate_draw_pos(cell_x, -x),
                calculcate_draw_pos(cell_y, -y),
                self.brush_density,
            );

            self.simulator.grid.add_density(
                calculcate_draw_pos(cell_x, y),
                calculcate_draw_pos(cell_y, x),
                self.brush_density,
            );
            self.simulator.grid.add_density(
                calculcate_draw_pos(cell_x, y),
                calculcate_draw_pos(cell_y, -x),
                self.brush_density,
            );
            self.simulator.grid.add_density(
                calculcate_draw_pos(cell_x, -y),
                calculcate_draw_pos(cell_y, x),
                self.brush_density,
            );
            self.simulator.grid.add_density(
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

            let cell_x = (mouse_pos.x / self.simulator.grid.scale as f32) as usize;
            let cell_y = (mouse_pos.y / self.simulator.grid.scale as f32) as usize;

            match self.brush_type {
                BrushType::Filled => self.draw_filled_circle(cell_x, cell_y),
                BrushType::Outline => self.draw_outline_circle(cell_x, cell_y),
            }
        }
    }

    fn available_gradients() -> Vec<(&'static str, Arc<dyn Gradient>)> {
        vec![
            ("Inferno", Arc::new(colorgrad::preset::inferno())),
            ("Viridis", Arc::new(colorgrad::preset::viridis())),
            ("Plasma", Arc::new(colorgrad::preset::plasma())),
            ("Magma", Arc::new(colorgrad::preset::magma())),
            ("Gray", Arc::new(colorgrad::preset::greys())),
            ("Blues", Arc::new(colorgrad::preset::blues())),
            ("Greens", Arc::new(colorgrad::preset::greens())),
            ("Reds", Arc::new(colorgrad::preset::reds())),
            ("Spectral", Arc::new(colorgrad::preset::spectral())),
            ("Warm", Arc::new(colorgrad::preset::warm())),
            ("Cool", Arc::new(colorgrad::preset::cool())),
            ("Rainbow", Arc::new(colorgrad::preset::rainbow())),
            ("Turbo", Arc::new(colorgrad::preset::turbo())),
        ]
    }

    fn handle_ui(&mut self) -> &mut imgui::Ui {
        self.imgui_renderer
            .update(&mut self.imgui_context, &mut self.raylib_handle);
        let ui = self.imgui_context.frame();

        let mut viscosity_slider_value = (self.simulator.viscosity * 10000.0) as i32;
        let mut diffusion_slider_value = (self.simulator.diffusion * 100000.0) as i32;

        let gradients = App::available_gradients();
        let gradient_names: Vec<&str> = gradients.iter().map(|(name, _)| *name).collect();

        ui.window("Settings")
            .size([300.0, 100.0], imgui::Condition::FirstUseEver)
            .build(|| {
                ui.slider("Viscosity", 0, 100, &mut viscosity_slider_value);
                ui.slider("Diffusion", 0, 100, &mut diffusion_slider_value);
                ui.slider("Brush Size", 1, 100, &mut self.brush_size);
                ui.slider("Brush Density", 1.0, 255.0, &mut self.brush_density);

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

                if ui.combo(
                    "Gradient",
                    &mut self.gradient_index,
                    &gradient_names,
                    |item| std::borrow::Cow::Borrowed(item),
                ) {
                    self.simulator.grid.color_grad = gradients[self.gradient_index].1.clone();
                }
            });

        self.simulator.viscosity = viscosity_slider_value as f32 / 10000.0;
        self.simulator.diffusion = diffusion_slider_value as f32 / 100000.0;

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

            let cell_x = (current_mouse_pos.x / self.simulator.grid.scale as f32) as usize;
            let cell_y = (current_mouse_pos.y / self.simulator.grid.scale as f32) as usize;

            // self.grid.add_velocity(cell_x, cell_y, velocity_vector);
            for i in 0..self.brush_size / 2 {
                for j in 0..self.brush_size / 2 {
                    self.simulator.grid.add_velocity(
                        cell_x.saturating_add(i),
                        cell_y.saturating_add(j),
                        velocity_vector,
                    );
                    self.simulator.grid.add_velocity(
                        cell_x.saturating_sub(i),
                        cell_y.saturating_sub(j),
                        velocity_vector,
                    );
                    self.simulator.grid.add_velocity(
                        cell_x.saturating_add(i),
                        cell_y.saturating_sub(j),
                        velocity_vector,
                    );
                    self.simulator.grid.add_velocity(
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
            self.simulator.step();

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
