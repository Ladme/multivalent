// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

use linreg;
use statistical;
use std::fs::File;
use std::io::{self, Write};

use crate::particle::Particle;
use crate::simulation::{Dimensionality, System};

pub struct Diffusion {
    /// vector of diffusion coefficients calculated from every block
    diffusion: Vec<f64>,
    /// vector of Mean Squared Deviations for every sweep
    msd: Vec<f64>,
    /// vector of particles containing initial positions of particles
    init_positions: Vec<Particle>,
    /// number of equilibration sweeps in each simulation repeat
    eq_sweeps: u32,
    /// frequency of MSD calculation
    msd_freq: u32,
    /// vector of sweeps at which MSD was calculated
    sweeps: Vec<u32>,
    /// dimensionality of the calculation
    dimensionality: Dimensionality,
    /// number of simulations in one block
    repeats_per_block: u32,
}

impl Diffusion {

    /// Generate a new structure for the calculation of diffusion coefficient.
    pub fn new(system: &System, init_positions: &Vec<Particle>) -> Diffusion {

        let len = system.prod_sweeps / system.msd_freq;
        let mut sweeps = Vec::with_capacity(len as usize);

        for i in 0..len {
            sweeps.push((i + 1) as u32 * system.msd_freq + system.eq_sweeps);
        }

        Diffusion { diffusion: Vec::with_capacity((system.repeats / system.diff_block) as usize),
                    msd:  vec![0.0; (system.prod_sweeps / system.msd_freq) as usize],
                    eq_sweeps: system.eq_sweeps,
                    msd_freq: system.msd_freq,
                    sweeps: sweeps,
                    init_positions: init_positions.clone(),
                    dimensionality: system.dimensionality,
                    repeats_per_block: system.diff_block,
                  }
    }

    /// Calculates diffusion coefficient from the provided MSD data and adds it to the diffusion vector.
    /// 
    /// ## Implementation
    /// Fits a line to the MSD curve. Slope of the line (divided by dimensionality) is the diffusion coefficient.
    pub fn calc_diffusion(&mut self) {
        

        let (slope, _): (f64, f64) = linreg::linear_regression(&self.sweeps, &self.msd).expect("\nError. Internal linalg error. Could not fit line through MSD data.");

        let diff = match self.dimensionality {
            Dimensionality::TWO => slope / 4.0,
            Dimensionality::ONE => slope / 2.0,
        };

        self.diffusion.push(diff);

    }

    /// Calculates average diffusion coefficient from the collected diffusion coefficients. Also provides standard deviation.
    pub fn get_average_diffusion(&self) -> (f64, f64) {

        if self.diffusion.len() == 1 {
            return (self.diffusion[0], 0.0);
        }

        let av_diff = statistical::mean(&self.diffusion);
        let std_diff = statistical::standard_deviation(&self.diffusion, Some(av_diff));

        (av_diff, std_diff)
    }

    /// Sets the MSD for all sweeps to 0.
    pub fn clear_msd(&mut self) {
        let len = self.msd.len();
        self.msd = vec![0.0; len];
    }

    /// Calculates MSD for the current configuration of particles.
    pub fn calc_msd(&mut self, particles: &Vec<Particle>, current_sweep: u32) {

        for (i, part) in particles.iter().enumerate() {
            let x = part.position[0] - self.init_positions[i].position[0];
            let y = part.position[1] - self.init_positions[i].position[1];

            self.msd[((current_sweep - self.eq_sweeps) / self.msd_freq - 1) as usize] += x * x + y * y;
        }
    }

    /// Normalizes MSD to be independent of a) the number of simulations it was obtained using and b) the number of particles.
    pub fn normalize_msd(&mut self) {
        for val in &mut self.msd {
            *val /= (self.repeats_per_block * self.init_positions.len() as u32) as f64;
        }
    }

    /// Get diffusion coefficient from the diffusion vector at target index.
    pub fn get_diff(&self, index: usize) -> f64 {
        self.diffusion[index]
    }

    pub fn write_msd(&self, filename: &str, block_id: u32) -> Result<(), io::Error> {

        // create MSD file
        let mut output = File::create(filename)?;

        assert_eq!(self.sweeps.len(), self.msd.len());

        write!(output, "$ label msd {}\n", block_id)?;
        // write msd data for each sweep
        for i in 0..self.sweeps.len() {
            write!(output, "{} {}\n", self.sweeps[i], self.msd[i])?;
        }

        Ok(())


    }

}