use rand::Rng;
use raylib::prelude::*;
use std::ops;

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

    // project is called twice in the step function. Once for projecting
    // v0 and once for projecting v. Because of that, we have two types of
    // actions. one is ToPrevious and one is ToCurrent. This isn't really
    // following clean coding principles because it really makes the implementation
    // tied to the project. But this at least makes the step function way more readable
    // which is what matters.
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
            ProjectAction::ToPrevious => |cell: &mut Cell, value: f32| cell.velocity_0.x = value,
            ProjectAction::ToCurrent => |cell: &mut Cell, value: f32| cell.velocity.x = value,
        };
        let set_div = match action_type {
            ProjectAction::ToPrevious => |cell: &mut Cell, value: f32| cell.velocity_0.y = value,
            ProjectAction::ToCurrent => |cell: &mut Cell, value: f32| cell.velocity.y = value,
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
                let current_y = get_x(&self.grid.cells[y][x]);

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

        for y in 1..self.grid.size - 1 {
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
        }

        for y in 1..self.grid.size - 1 {
            for x in 1..self.grid.size - 1 {
                set_current(&mut self.grid.cells[y][x], new_values[y][x]);
            }
        }

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
    }

    // Now it is time to write project... Here are what I should first consider and answer
    // Q: What are the inputs to project? what is given as input to project?
    // A: In the blog, we have Three velocity inputs and two other arrays of pressure
    //    field and divergence which are... temp storage? I should examine if I really
    //    need to have the temp spaces as struct fields or not. Before everything else.
    //    let's see how project is called in the step function.

    //    project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
    //    project(Vx, Vy, Vz, Vx0, Vy0, 4, N);

    // So first we call it with the 0 Velocities and pass the Vx and Vy as temp arrays?
    // Then we call it with the non-0 Velocities and pass the same velocities without 0
    // So a main question is, are the p and div parameters that are passed here used somewhere else?
    // Namely, is Vx and Vy used somewhere else after the first project call? If they are going to
    // be used, it is going to be in the advect function call. Vx and Vy are passed as the *d parameter.
    // The value of d is not really used! It is just used to store a value independent of the d value
    // in the advect function. The same thing goes for Vy. It just stores a new value. We also pass
    // it to set_bnd afterwards to set the boundries of it.

    // Let's see it after the second call. They are not used! Unless... We run the step function
    // again and it is used in diffuse? What happens there? The value is used here, Vx0 and Vy0
    // is passed as *x parameter to lin_solve and lin_solve uses that to update values.

    // I think I should start to develop a solid understanding of the input and output parameters.
    // and that will help me understand things better here. Let's start to analyze the step function

    // void FluidCubeStep(FluidCube *cube)
    // {
    //     int N          = cube->size;
    //     float visc     = cube->visc;
    //     float diff     = cube->diff;
    //     float dt       = cube->dt;
    //     float *Vx      = cube->Vx;
    //     float *Vy      = cube->Vy;
    //     float *Vz      = cube->Vz;
    //     float *Vx0     = cube->Vx0;
    //     float *Vy0     = cube->Vy0;
    //     float *Vz0     = cube->Vz0;
    //     float *s       = cube->s;
    //     float *density = cube->density;
    //
    //     diffuse(1, Vx0, Vx, visc, dt, 4, N);
    //     diffuse(2, Vy0, Vy, visc, dt, 4, N);
    //     diffuse(3, Vz0, Vz, visc, dt, 4, N);
    //
    //     project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
    //
    //     advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
    //     advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
    //     advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
    //
    //     project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
    //
    //     diffuse(0, s, density, diff, dt, 4, N);
    //     advect(0, density, s, Vx, Vy, Vz, dt, N);
    // }

    // Here we are first having 3 calls to diffuse. What is the purpose of them
    // here? we should thrive to understand what is going on with the inputs here.
    // The first parameter is b, which defines what kind of input we are giving it.
    // 1 for Vx, 2 for Vy and 3 for Vz. This is trivial but what matters and is important
    // is the next two parameters. We are passing Vx0 as x* and Vx as *x0! The rest isn't
    // important so we will focus on the second and third parameter for now.
    // The value of Vx0 is updated with a little bit of help of the value of Vx in the first
    // diffusion and so on... So after the diffusion, the arrays of Vx0 and Vy0 and Vz0 contain
    // the new diffused vectors. So to sum it up:

    // We use Vx to update Vx0 to a diffused array

    // This diffusion can ruin our rule of incompressible fluids so now we apply project.
    // Let's see what project does and what it updates. Project first sets all the values
    // of p parameter to 0! so this means that Vx just becomes all zero! Afterwards it sets
    // the values of divergence (parameter div) to Vy. So here too, Vy gets filled up with
    // new values. the old values are gone. Afterwards, we set the boundries for p and div,
    // which basically updates them. We also solve a linear equation and update p using the
    // help of div. So far, we are setting values to Vx and Vy (we also have to consider
    // the question of whether Vx and Vy are used later on!)
    // After this, we update the values of Vx0, Vy0 and Vz0. So far, we are just working with
    // Vx0, Vy0 and Vz0!

    // Now let's examine advect. This is quite simple to examine. the only place that we update
    // values is with the d parameter. Let's see what the d parameter is. Here, Vx, Vy and Vz
    // are our d parameters. They are getting new values! But let's see what other values are used
    // here. The V0 velocities are used here.

    // Now we get another project call but this time with Vx, Vy and Vz. it used Vx0 and Vy0 as
    // p and div.

    // Hold on... Is it clicking for me? Now Vx, Vy and Vz should contain the correct values. Now
    // after all of this hassle, we go for diffusing and then advecting our dye density. So what matters
    // at the end is Vx, Vy and Vz! But what happens to Vx0, Vy0 and Vz0? Are they used for the Next loop?
    // Of course! As I explained at the first of this thought process, We use the current value of V0 with
    // a little bit of help from V to update everything. So... They are used?? But what happened to them
    // before?

    // So turns out it doesn't really matter what initial value we give to lin solve! The whole purpose
    // of it is to iteratively make the result better. So I have kind of started to understand everything
    // now.

    // 1. First diffuse the dye using a linear solving function. Stores the new value in Vx0, Vy0 and Vz0
    //    using Vx, Vy and Vz.
    // 2. Fix the screwed Vx0, Vy0 and Vz0 using project. Using two other temp variables of Vx and Vy.
    // 3. Advect the Vx0 into Vx. So now we will have the updated velocity vectors in Vx, Vy and Vz.
    // 4. Fix the issues in Vx, Vy and Vz.
    // 5. Use the final velocities to move the dye around

    // 2. How does project use lin_solve and set_bnd? Is my current implementation enough
    //    for them?

    // Now this is the main question. We have to see the usages that are for lin_solve and
    // set_bnd

    // lin_solve:
    //      diffuse:
    //          x: v0, y: v
    //          x: s, y: density
    //      project:
    //
    //

    // set_bnd:
    //      lin_solve:
    //      project:
    //      advect:

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
                    color: raylib::prelude::Color { r, g, b, a: 255 },
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
