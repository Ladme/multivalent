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
    /// initial center of geometry of the system
    pub initial_center: [f64; 2],
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
    pub fn new(system: &System) -> Diffusion {

        let len = system.prod_sweeps / system.msd_freq;
        let mut sweeps = Vec::with_capacity(len as usize);

        for i in 0..len {
            sweeps.push((i + 1) as u32 * system.msd_freq);
        }

        Diffusion { diffusion: Vec::with_capacity((system.repeats / system.diff_block) as usize),
                    msd:  vec![0.0; (system.prod_sweeps / system.msd_freq) as usize],
                    msd_freq: system.msd_freq,
                    sweeps,
                    initial_center: [0.0, 0.0],
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

        let center = System::center(particles);

        let x = center[0] - self.initial_center[0];
        let y = center[1] - self.initial_center[1];

        self.msd[(current_sweep / self.msd_freq - 1) as usize] += x * x + y * y;
    }

    /// Normalizes MSD to be independent of the number of simulations it was obtained using.
    pub fn normalize_msd(&mut self) {
        for val in &mut self.msd {
            *val /= self.repeats_per_block as f64;
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


/*
*************************************
            UNIT TESTS
*************************************
*/


#[cfg(test)]
mod tests {

    use crate::parser::parse_input;
    use super::Diffusion;

    const INPUT_FILE: &str = "test_files/test_input_energy";

    #[test]
    fn test_clear_msd() {

        let mut system = parse_input(INPUT_FILE).expect("Could not find input file.");
        let mut diffusion = Some(Diffusion::new(&system));

        system.run_production( &mut diffusion, &mut None, 1);

        if let Some(unwrapped) = diffusion.as_mut() {

            for val in &unwrapped.msd {
                assert_ne!(*val, 0.0);
            }
    
            unwrapped.clear_msd();
    
            for val in &unwrapped.msd {
                assert_eq!(*val, 0.0);
            }

        }   
    }

    #[test]
    fn test_calc_msd() {

        let system = parse_input(INPUT_FILE).expect("Could not find input file.");

        let mut diffusion = Diffusion::new(&system);

        diffusion.initial_center = [0.0, 0.0];
        diffusion.calc_msd(&system.particles, 100);

        let result = diffusion.msd[0];
        let expected = 0.082;

        assert!((result - expected).abs() < 0.0001);
    }

    #[test]
    fn test_normalize_msd() {

        let mut system = parse_input(INPUT_FILE).expect("Could not find input file.");
        let mut diffusion = Some(Diffusion::new(&system));

        system.run_production( &mut diffusion, &mut None, 1);
        system.run_production(&mut diffusion, &mut None, 2);

        if let Some(unwrapped) = diffusion.as_mut() {

            let original = unwrapped.msd.clone();

            unwrapped.normalize_msd();

            for i in 0..unwrapped.msd.len() {
                assert!((original[i] - 2.0 * unwrapped.msd[i]).abs() < 0.0001);
            }
        }
    }

}