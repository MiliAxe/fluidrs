use raylib::prelude::*;
use std::sync::Arc;

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
    pub color_grad: Arc<dyn colorgrad::Gradient>,
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
                let color_intensity = cell.density.clamp(0.0, 255.0);
                let color_intensity_normalize = color_intensity / 255.0;
                let color = self.color_grad.at(color_intensity_normalize).to_rgba8();

                cell.color = raylib::prelude::Color {
                    r: color[0],
                    g: color[1],
                    b: color[2],
                    a: color[3],
                };
            }
        }
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
            color_grad: Arc::new(colorgrad::preset::inferno()),
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
