use livemod::LiveModHandle;
use macroquad::prelude::*;
use noise::{NoiseFn, Perlin};
use std::{
    cmp::Ordering,
    ops::{Deref, DerefMut},
    time::{Duration, Instant},
};

#[derive(Debug, Clone, Copy)]
struct Particle {
    pos: Vec2,
    birth_time: Instant,
    life_span: Duration,
}

#[derive(Debug, Clone, Copy)]
struct Vortex {
    pos: Vec2,
    range: f32,
    direction: i8,
    birth_time: Instant,
}

#[derive(Debug, Clone)]
struct Tri {
    v: Vec<Vec2>,
    c: Vec<Color>,
}

impl Tri {
    fn should_contain(&self, v: Vec2) -> bool {
        self.circ_center().distance(v) <= self.circ_radius()
    }

    fn circ_center(&self) -> Vec2 {
        let point_a = (self.v[0] + self.v[1]) / 2.;
        let point_b = (self.v[2] + self.v[1]) / 2.;
        let vec_a = (self.v[0] - self.v[1]).perp();
        let vec_b = (self.v[2] - self.v[1]).perp();

        let b = (point_a.x + (point_b.y - point_a.y) * vec_a.x / vec_a.y - point_b.x)
            / (vec_b.x - (vec_b.y * vec_a.x) / (vec_a.y));

        point_b + b * vec_b
    }
    fn circ_radius(&self) -> f32 {
        self.circ_center().distance(self.v[0])
    }
    fn raster(&self, depth: u16) {
        if depth == 1 {
            draw_triangle(
                self.v[0],
                self.v[1],
                self.v[2],
                //1.,
                Color::new(
                    self.c[0].r * 0.3333 + self.c[1].r * 0.3333 + self.c[2].r * 0.3333,
                    self.c[0].g * 0.3333 + self.c[1].g * 0.3333 + self.c[2].g * 0.3333,
                    self.c[0].b * 0.3333 + self.c[1].b * 0.3333 + self.c[2].b * 0.3333,
                    self.c[0].a * 0.3333 + self.c[1].a * 0.3333 + self.c[2].a * 0.3333,
                ),
            )
        } else {
            let v_ab = (self.v[0] + self.v[1]) / 2.;
            let c_ab = Color::from_vec((self.c[0].to_vec() + self.c[1].to_vec()) / 2.);

            let v_bc = (self.v[1] + self.v[2]) / 2.;
            let c_bc = Color::from_vec((self.c[1].to_vec() + self.c[2].to_vec()) / 2.);

            let v_ca = (self.v[2] + self.v[0]) / 2.;
            let c_ca = Color::from_vec((self.c[2].to_vec() + self.c[0].to_vec()) / 2.);

            Tri {
                v: vec![self.v[0], v_ab, v_ca],
                c: vec![self.c[0], c_ab, c_ca],
            }
            .raster(depth - 1);

            Tri {
                v: vec![self.v[1], v_ab, v_bc],
                c: vec![self.c[1], c_ab, c_bc],
            }
            .raster(depth - 1);

            Tri {
                v: vec![self.v[2], v_bc, v_ca],
                c: vec![self.c[2], c_bc, c_ca],
            }
            .raster(depth - 1);

            Tri {
                v: vec![v_ab, v_bc, v_ca],
                c: vec![c_ab, c_bc, c_ca],
            }
            .raster(depth - 1);
        }
    }
    fn divide(&self, v: Vec2, c: Color) -> Vec<Tri> {
        let mut ret = Vec::new();
        let mut color_point_angle: Vec<(Color, Vec2, f32)> = self
            .v
            .iter()
            .zip(self.c.iter())
            .map(|(e, i)| (*i, *e, (*e - v).angle_between(Vec2::new(1., 0.))))
            .collect();

        color_point_angle.sort_by(|(_, _, a), (_, _, b)| {
            if a < b {
                Ordering::Less
            } else {
                Ordering::Greater
            }
        });

        for i in 0..color_point_angle.len() {
            let (p1_c, p1_v, _) = color_point_angle[i];
            let (p2_c, p2_v, _) = color_point_angle[(i + 1) % color_point_angle.len()];

            ret.push(Tri {
                v: vec![v, p1_v, p2_v],
                c: vec![c, p1_c, p2_c],
            })
        }

        ret
    }

    fn min_side(&self) -> f32 {
        let (ret, _) = self
            .v
            .iter()
            .fold((0., self.v[self.v.len() - 1]), |(minimum, last), elm| {
                (f32::min(minimum, elm.distance(last)), *elm)
            });

        ret
    }

    fn perimeter(&self) -> f32 {
        let (ret, _) = self
            .v
            .iter()
            .fold((0., self.v[self.v.len() - 1]), |(sum, last), elm| {
                (sum + elm.distance(last), *elm)
            });
        ret
    }

    fn merge(&self, other: Tri) -> Tri {
        let mut ret = self.clone();
        ret.v.extend(&other.v);
        ret.c.extend(&other.c);
        let all_v = ret.v.clone();

        ret.v.sort_by(|a, b| {
            if a.x == b.x {
                if a.y == b.y {
                    Ordering::Equal
                } else if a.y > b.y {
                    Ordering::Greater
                } else {
                    Ordering::Less
                }
            } else {
                if a.x < b.x {
                    Ordering::Less
                } else {
                    Ordering::Greater
                }
            }
        });
        ret.v.dedup();

        ret.c = ret
            .v
            .iter()
            .map(|v| ret.c[all_v.iter().position(|e| v == e).unwrap_or(0)])
            .collect();

        ret
    }
}
fn teselate(sw: f32, sh: f32, positions: Vec<Vec2>, colors: Vec<Color>) -> Vec<Tri> {
    // println!("starting to calc tris");

    for p in &positions {
        draw_circle_lines(p.x, p.y, 10., 0.1, YELLOW)
    }

    let mut tris = Vec::new();
    tris.push(Tri {
        v: vec![
            Vec2::new(-1., -1.),
            Vec2::new(2. * sw, 0.),
            Vec2::new(0., 2. * sw),
        ],
        c: vec![WHITE, WHITE, WHITE],
    });

    for i in 0..positions.len() {
        // println!("start loop");

        let mut new_tris = Vec::new();
        let mut big_tri = Tri {
            v: vec![],
            c: vec![],
        };

        // println!("a, {}", tris.len());

        for t in tris {
            if t.should_contain(positions[i]) {
                big_tri = big_tri.merge(t);
            } else {
                new_tris.push(t);
            }
        }

        // println!("b,{}", new_tris.len());

        let division = big_tri.divide(positions[i], colors[i]);
        // println!("c,{}", division.len());

        new_tris.extend(division);
        // println!("d,{}", new_tris.len());

        tris = new_tris;
        // println!("c,{}", tris.len());

        // println!("{:?}", tris);
    }

    tris = tris
        .iter()
        .filter(|e| {
            for v in &e.v {
                if v.x > sw || v.x < 0. || v.y > sh || v.x < 0. {
                    return false;
                }
            }
            true
        })
        .filter(|e| e.perimeter() < 300.)
        .map(|e| e.clone())
        .collect();

    tris
}
#[macroquad::main("fireplace")]
async fn main() {
    let livemod = LiveModHandle::new_gui();

    let mut perlin_time_scale_a = livemod.create_variable("perlin_time_scale_a", 0_f64);
    let mut perlin_space_scale_a = livemod.create_variable("perlin_space_scale_a", 0_f64);
    let mut perlin_time_scale_b = livemod.create_variable("perlin_time_scale_b", 0_f64);
    let mut perlin_space_scale_b = livemod.create_variable("perlin_space_scale_b", 0_f64);
    let mut max_life_span = livemod.create_variable("max_life_span", 0_f64);
    let mut vortex_density = livemod.create_variable("vortex_density", 0_f32);
    let mut vortex_strength = livemod.create_variable("vortex_strength", 0_f32);
    {
        *(perlin_time_scale_a.lock_mut().deref_mut()) = 582.;
        *(perlin_space_scale_a.lock_mut().deref_mut()) = 243.;
        *(perlin_time_scale_b.lock_mut().deref_mut()) = 320.;
        *(perlin_space_scale_b.lock_mut().deref_mut()) = 100.;
        *(max_life_span.lock_mut().deref_mut()) = 2035.;
        *(vortex_density.lock_mut().deref_mut()) = 97.;
        *(vortex_strength.lock_mut().deref_mut()) = 3000.;
    }
    let mut fire_particles = Vec::new();
    let mut fire_vortexes = Vec::new();

    let mut last_loop_time = 0.;
    let start = Instant::now();

    let perlin = Perlin::new();

    loop {
        let loop_start = Instant::now();
        // println!("tick");

        // create particles
        let sw = screen_width();
        let sh = screen_height();

        let ms = start.elapsed().as_millis() as f64;

        //if rand::gen_range(0, 100) > 90 {
        for _ in 0..100 {
            let x = rand::gen_range(100., sw - 100.);
            fire_particles.push(Particle {
                pos: Vec2::new(x, sh),
                birth_time: Instant::now(),
                life_span: Duration::new(
                    0,
                    ((((1.
                        + perlin.get([
                            ms / (perlin_time_scale_b.lock().deref() + 0.1),
                            (x as f64) / (0.1 + perlin_space_scale_b.lock().deref()),
                        ]))
                        * 1000000.)
                        + ((1.
                            + perlin.get([
                                ms / (perlin_time_scale_a.lock().deref() + 0.1),
                                (x as f64) / (0.1 + perlin_space_scale_a.lock().deref()),
                            ]))
                            * 1000000.)
                        - 50000.)
                        * max_life_span.lock().deref()) as u32,
                ),
            });
            // println!("{:?}",fire_particles);
        }

        if rand::gen_range(0., 100.) > *vortex_density.lock().deref() {
            fire_vortexes.push(Vortex {
                pos: Vec2::new(rand::gen_range(0., sw), sh),
                range: rand::gen_range(100., 200.),
                direction: -1 + 2 * (rand::gen_range::<i16>(0, 2) as i8),
                birth_time: Instant::now(),
            });
        }

        // remove OOB

        fire_particles = fire_particles
            .iter()
            .filter(|e| e.birth_time.elapsed() < e.life_span)
            .filter(|e| e.pos.y >= 0.)
            .filter(|e| e.pos.y <= screen_height())
            .filter(|e| e.pos.x >= 0.)
            .filter(|e| e.pos.x <= screen_width())
            .map(|e| *e)
            .collect();

        fire_vortexes = fire_vortexes
            .iter()
            .filter(|e| e.pos.y >= -screen_height())
            .filter(|e| e.pos.y <= screen_height())
            .filter(|e| e.pos.x >= 0.)
            .filter(|e| e.pos.x <= screen_width())
            .map(|e| *e)
            .collect();

        // vortexes

        for v in &fire_vortexes {
            for p in &mut fire_particles {
                if v.pos.distance(p.pos) < v.range {
                    let influence = (v.range - v.pos.distance(p.pos))
                        * if v.birth_time.elapsed().as_millis() < 1000 {
                            v.birth_time.elapsed().as_millis() as f32 / 1000.
                        } else {
                            1.
                        };

                    p.pos -= v.pos;

                    let m = p.pos.length();
                    let a = -p.pos.angle_between(Vec2::new(1., 0.))
                        + (influence / vortex_strength.lock().deref())
                            * (v.direction as f32)
                            * last_loop_time
                            / 50.;

                    p.pos = v.pos + Vec2::new(a.cos() * m, a.sin() * m);
                }
            }
        }

        // atraction

        for v in &fire_vortexes {
            let other = fire_vortexes[rand::gen_range(0, fire_vortexes.len())];
        }

        // // viscosity

        // let max_d = 20.;

        // for a in 0..fire_particles.len() {
        //     for b in 0..fire_particles.len() {
        //         if a != b
        //             && max_d > fire_particles[a].pos.distance(fire_particles[b].pos)
        //             && rand::gen_range(0., 1.) > 0.5
        //         {
        //             let intensity = fire_particles[a].pos.distance(fire_particles[b].pos) / max_d;

        //             let speed_a = fire_particles[a].pos - fire_particles[a].old_pos;
        //             let speed_b = fire_particles[b].pos - fire_particles[b].old_pos;

        //             let speed_ab = (speed_a + speed_b) / 2.;

        //             let new_speed_a = intensity * speed_ab + (1. - intensity) * speed_a;
        //             let new_speed_b = intensity * speed_ab + (1. - intensity) * speed_b;

        //             fire_particles[a].old_pos = fire_particles[a].pos - new_speed_a;
        //             fire_particles[b].old_pos = fire_particles[b].pos - new_speed_b;
        //         }
        //     }
        // }

        // heat dynamics

        for p in &mut fire_particles {
            p.pos += Vec2::new(0., -0.2)
                * (1.
                    - (p.birth_time.elapsed().as_millis() as f32 / p.life_span.as_millis() as f32))
                * (last_loop_time + 1.)
        }

        // move particles

        // for p in &mut fire_particles {
        //     let new_old_pos = p.pos;
        //     p.pos += -p.old_pos + p.pos;
        //     p.old_pos = new_old_pos;
        // }

        for c in &mut fire_vortexes {
            c.pos += Vec2::new(0., -0.1) * last_loop_time;
        }

        // calculate tri mesh

        // let tris = teselate(
        //     sw,
        //     sh,
        //     fire_particles.iter().map(|e| e.pos).collect(),
        //     fire_particles
        //         .iter()
        //         .map(|e| {
        //             fire_color(
        //                 (e.birth_time.elapsed().as_millis()) as f32
        //                     / e.life_span.as_millis() as f32,
        //             )
        //         })
        //         .collect(),
        // );

        // draw

        //clear_background(Color::new(0., 0., 0., 0.5));

        // for t in &tris {
        //     t.raster(4)
        // }

        // draw particles

        for p in &fire_particles {
            let col = fire_color(
                p.birth_time.elapsed().as_millis() as f32 / p.life_span.as_millis() as f32,
            );
            // println!("{},{},{}: r{} g{} b{} a{}", p.pos.x, p.pos.y, 10.,col.r,col.g,col.b,col.a);
            draw_circle(p.pos.x, p.pos.y, 10., col);
        }

        // for p in &fire_vortexes {
        //     draw_circle(p.pos.x, p.pos.y, 10., YELLOW);
        // }

        last_loop_time = loop_start.elapsed().as_millis() as f32;
        next_frame().await;
    }
}

fn fire_color(i: f32) -> Color {
    let fade = Color::from_rgba(0, 0, 0, 0);

    let steps = vec![0., 0.1, 0.2, 0.5, 1.]; //vec![0., 0.5, 1.];
    let colors = vec![WHITE, YELLOW, ORANGE, RED, fade];

    let j = steps
        .iter()
        .position(|elm| elm >= &i)
        .unwrap_or(steps.len() - 1)
        .max(1);
    let first_weight = (i - steps[j - 1]) / (steps[j] - steps[j - 1]);
    let second_weight = 1. - first_weight;
    let first_index = j - 1;
    let second_index = j;

    // println!("{},{}", first_weight, second_weight);

    let ret = Color::new(
        (colors[first_index].r * second_weight) + (colors[second_index].r * first_weight),
        (colors[first_index].g * second_weight) + (colors[second_index].g * first_weight),
        (colors[first_index].b * second_weight) + (colors[second_index].b * first_weight),
        (colors[first_index].a * second_weight) + (colors[second_index].a * first_weight),
    );

    // println!("col: {:?}", ret);

    ret
}
