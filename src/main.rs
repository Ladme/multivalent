// Released under MIT License.
// Copyright (c) 2023 Ladislav Bartos

mod statistics;
mod bond;
mod particle;
mod simulation;
mod parser;
mod diffusion;
mod wl;

use parser::parse_input;
use std::env;

pub const VERSION: &str = "0.4.0";
/// Arbitrary constant that is used to multiply the obtained diffusion values, so they are human-readable.
pub const DIFFUSION_MULTIPLIER: f64 = 1000.0;


fn main() {

    // get path to simulation input file from arguments
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        eprintln!("\nError. Incorrect number of arguments. Usage: {} SIMULATION_INPUT\n", args[0]);
        return;
    }

    println!("\n@@@@@@@@@@@@@@@@@@");
    println!("MULTIVALENT v{}", VERSION);
    println!("@@@@@@@@@@@@@@@@@@\n");

    println!("Reading simulation input from `{}`...", &args[1]);
    let mut system = match parse_input(&args[1]) {
        Ok(x) => x,
        Err(_) => return,
    };

    if !system.sanity_check() {
        return;
    }

    system.display();
    if !system.run() { return; }
    system.print_statistics();
}