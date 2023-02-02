// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

use std::fmt;
use crate::particle::Particle;
use crate::simulation::System;

/// Bond is a harmonic potential of the following form:
/// 
/// `(1/2) * k * (x - x_0)^2`, 
/// 
/// where `k` is `force_const`, 
/// 
/// `x_0` is `eq_dist`,
/// 
/// and `x` is the distance between the bonded particles.
pub struct Bond {
    /// bonded particles
    pub particles: [usize; 2],
    /// force constant
    force_const: f64,
    /// equilibrium distance between the particles
    eq_dist: f64,
}

impl Bond {

    /// Creates a new bond.
    pub fn new(particles: [usize; 2], force_const: f64, eq_dist: f64) -> Bond {
        Bond { particles, force_const, eq_dist }
    }

    /// Calculates energy of a bond connecting two particles.
    pub fn energy(&self, particles: &Vec<Particle>) -> f64 {

        let part1 = &particles[self.particles[0]];
        let part2 = &particles[self.particles[1]];

        let x = System::distance(part1, part2) - self.eq_dist;

        0.5 * self.force_const * x * x
    }

}

/// Allows the usage of print* macros for Bond.
impl fmt::Display for Bond {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Bonded Particles: {}, {}\nForce Constant: {}\nEquilibrium Distance: {}\n",
                   self.particles[0], self.particles[1], self.force_const, self.eq_dist)
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

    const INPUT_FILE: &str = "test_files/test_input_energy";

    #[test]
    fn test_bond_energy() {

        let system = parse_input(INPUT_FILE).expect("Could not find input file.");
        let expected = [0.5, 0.0, 0.01, 0.000345];
        for (i, bond) in system.bonds.iter().enumerate() {
            assert!((bond.energy(&system.particles) - expected[i]).abs() < 0.0001);
        }

    }
}