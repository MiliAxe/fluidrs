use rand::Rng;
use raylib::prelude::*;

pub struct Cell {
    color: raylib::prelude::Color,
}

pub struct Grid {
    cells: Vec<Vec<Cell>>,
    scale: i32,
}

impl Grid {
    pub fn new(size: i32, scale: i32) -> Grid {
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
                });
            }
            grid.push(row);
        }

        Grid { cells: grid, scale }
    }

    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        for (row_idx, row) in self.cells.iter().enumerate() {
            for (col_idx, value) in row.iter().enumerate() {
                let y = (row_idx as i32) * self.scale;
                let x = (col_idx as i32) * self.scale;
                d.draw_rectangle(x, y, self.scale, self.scale, value.color);
            }
        }
    }
}
