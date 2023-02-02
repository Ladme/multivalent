// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

#[derive(Clone)]
pub struct MoveStatistics {
    /// number of accepted MC moves for every particle
    pub accepted: Vec<u64>,
    /// number of rejected MC moves for every particle
    pub rejected: Vec<u64>,
}

impl MoveStatistics {

    /// Prepares and prints statistics of translation moves in the simulation(s).
    pub fn report(&self) {

        let total_acc: u64 = self.accepted.iter().sum();
        let total_rej: u64 = self.rejected.iter().sum();
        let total = total_acc + total_rej;

        println!("\nMove Statistics:");
        println!(">> Total moves: {}", total);
        println!(">> Accepted moves: {} ({:.2} %) ", total_acc, (100.0 * total_acc as f64 / total as f64));
        println!(">> Rejected moves: {} ({:.2} %) ", total_rej, (100.0 * total_rej as f64 / total as f64));

        println!("\nMove Statistics for Individual Particles:");
        for i in 0..self.accepted.len() {
            println!(">> Particle {}", i);

            let acc = self.accepted[i];
            let rej = self.rejected[i];
            let total = acc + rej;

            println!(">>>> Total moves: {}", total);
            println!(">>>> Accepted moves: {} ({:.2} %) ", acc, (100.0 * acc as f64 / total as f64));
            println!(">>>> Rejected moves: {} ({:.2} %) ", rej, (100.0 * rej as f64 / total as f64));
        }

    }
}