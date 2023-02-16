// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

use std::fs::File;
use std::io::Write;

use crate::simulation::System;
use crate::particle::Particle;

pub struct WangLandau {
    cv: Vec<f64>,
    histogram: Vec<u64>,
    bias: Vec<f64>,
    pub alpha: f64,
    wl_file: String,
    smoothness: f64,
    alpha_limit: f64,
}

impl WangLandau {

    pub fn new() -> WangLandau {

        let cv = vec![0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0];

        let cv_len = cv.len();

        WangLandau { cv,
                    histogram: vec![0u64; cv_len],
                    bias: vec![0.0; cv_len],
                    alpha: 1e-3,
                    wl_file: "wl.dat".to_string(),
                    smoothness: 0.0001,
                    alpha_limit: 1e-6,
                }
    }

    fn get_bin_index(&self, particles: &Vec<Particle>) -> Option<usize> {

        let center_x = (System::center(particles)[0] * 10.0).round() / 10.0;

        self.cv.iter().position(|&s| s == center_x)
    }

    pub fn energy(&self, particles: &Vec<Particle>) -> f64 {

        match self.get_bin_index(particles) {
            Some(x) => return self.bias[x],
            None => return f64::INFINITY,
        }
    }

    fn clear_histogram(&mut self) {

        let cv_len = self.cv.len();

        self.histogram = vec![0u64; cv_len];
    }

    pub fn update(&mut self, particles: &Vec<Particle>) {

        match self.get_bin_index(particles) {
            Some(x) => {
                self.histogram[x] += 1;
                self.bias[x] += self.alpha;
            },
            None => {
                eprintln!("Internal Multivalent error. System is not in CV range. WL calculation failed.");
            }   
        }

        // check histogram smoothness
        let max_h = *self.histogram.iter().max().unwrap();
        let min_h = *self.histogram.iter().min().unwrap();
        let smooth = (max_h as f64 / min_h as f64).log2();

        if smooth < self.smoothness {
            self.clear_histogram();
            self.alpha /= 2.0;
        }
    }

    pub fn check_alpha_limit(&self) -> bool {
        if self.alpha < self.alpha_limit {
            true
        } else {
            false
        }
    }

    pub fn write_free_energy(&self) {
        
        let mut output = match File::create(&self.wl_file) {
            Ok(x) => x,
            Err(_) => {
                eprintln!("\nError. File `{}` could not be created.", &self.wl_file);
                return
            }
        };

        if let Err(_) = write!(output, "# ALPHA: {}\n", self.alpha) {
            eprintln!("\nError. Could not write line into file `{}`.", &self.wl_file);
        }

        for i in 0..self.cv.len() {
            if let Err(_) = write!(output, "{} {} {}\n", self.cv[i], -self.bias[i], self.histogram[i]) {
                eprintln!("\nError. Could not write line into file `{}`.", &self.wl_file);
            }
        }
    }

}