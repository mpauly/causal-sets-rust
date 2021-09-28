use std::env;
use std::path::Path;
use std::process;

use causal_sets::Config;
use causal_sets::run_simulation;


fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() != 2 {
        eprintln!("Please provide a working directory as the single argument.");
        process::exit(1);
    }
    let working_dir = &args[1];
    let working_dir = Path::new(working_dir);
    
    let config = Config::new_from_file(&working_dir.join("parameters.yml"));
    println!("Working dir: {}", working_dir.display());
    run_simulation(config);
}
