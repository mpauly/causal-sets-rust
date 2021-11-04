use rand::prelude::*;
use rayon::prelude::*;
use rand::{distributions::Standard, Rng};
use std::fs::File;
use std::sync::{Arc, Mutex};
use std::io::prelude::*;
use typed_arena::Arena;
use yaml_rust::YamlLoader;

mod causal_set;
use crate::causal_set::{Configuration, Node, NodeCoordinate};

pub struct Config {
    grid_size: i64,
    epsilon: Vec<f64>,
    beta: Vec<f64>,
    verbose: bool,
    nr_of_sweeps: i64,
    output_path: std::path::PathBuf,
    // final_config_path: std::path::PathBuf,
}

pub struct RunConfig {
    grid_size: i64,
    epsilon: f64,
    beta: f64,
    verbose: bool,
    nr_of_sweeps: i64,
    output_file: Arc<Mutex<File>>,
}

impl Config {
    pub fn new_from_file(file_name: &std::path::Path) -> Config {
        let mut parameter_file = match File::open(file_name) {
            Ok(f) => f,
            Err(e) => panic!(
                "Error opening file parameters.yml in {}: {}",
                file_name.display(),
                e
            ),
        };
        let mut parameter_file_string = String::new();
        parameter_file
            .read_to_string(&mut parameter_file_string)
            .unwrap();

        let parameter_files = YamlLoader::load_from_str(&parameter_file_string).unwrap();
        let parameters = &parameter_files[0];

        let output_path = match parameters["output_file"].as_str() {
            Some(val) => file_name.with_file_name(val),
            None => file_name.with_file_name("output.dat"),
        };
        //let final_config = match parameters["final_config_file"].as_str() {
        //    Some(val) => file_name.with_file_name(val),
        //    None => file_name.with_file_name("final_config.dat"),
        //};
        let verbosity = match parameters["verbose"].as_bool() {
            Some(val) => val,
            None => false,
        };
        let nr_of_sweeps = match parameters["nr_of_sweeps"].as_i64() {
            Some(val) => val,
            None => 5000,
        };

        let betas = parameters["beta"].as_vec();
        let betas = match betas {
            Some(vec) => vec.iter().map(|x| x.as_f64().unwrap()).collect(),
            None => vec![parameters["beta"].as_f64().unwrap()],
        };
        let epsilons = parameters["epsilon"].as_vec();
        let epsilons = match epsilons {
            Some(vec) => vec.iter().map(|x| x.as_f64().unwrap()).collect(),
            None => vec![parameters["epsilon"].as_f64().unwrap()],
        };

        Config {
            grid_size: parameters["nr_of_points"].as_i64().unwrap(),
            epsilon: epsilons,
            beta: betas,
            verbose: verbosity,
            nr_of_sweeps: nr_of_sweeps,
            output_path: output_path,
            // final_config_path: final_config,
        }
    }
}

pub fn run_simulation(config: Config) {
    let output_file = Arc::new(Mutex::new(File::create(config.output_path).unwrap()));

    let n_string: Vec<String> = (0..=config.grid_size - 2)
        .map(|x| format!("N{}", x))
        .collect();
    let n_string = n_string.join("\t");
    {
        let mut output_file = output_file.lock().unwrap();
        match writeln!(output_file, "beta\tepsilon\tsweep\taction\t{}", n_string) {
            Ok(_) => (),
            Err(e) => panic!("Error writing to output file: {}", e),
        };
    }

    let verbosity = config.verbose;
    let grid_size = config.grid_size;
    let nr_of_sweeps = config.nr_of_sweeps;

    let runconfs: Vec<RunConfig> = itertools::iproduct!(config.beta.iter(), config.epsilon.iter()).map(|(beta, epsilon)| RunConfig {
        epsilon: *epsilon,
        beta: *beta,
        verbose: verbosity,
        grid_size: grid_size,
        nr_of_sweeps: nr_of_sweeps,
        output_file: Arc::clone(&output_file),
    }).collect();

    let _results : Result<Vec<(i64, i64)>, _> = runconfs.par_iter().map(run_simulation_for_params).collect();
}


fn run_simulation_for_params(config: &RunConfig) -> Result<(i64, i64), &str> {
    println!(
        "Starting run with N={} at eps={}, beta={}",
        config.grid_size, config.epsilon, config.beta
    );

    // following Surya 2012
    let steps_per_sweep = config.grid_size * (config.grid_size - 1) / 2;
    // setup nodes
    // we are using an arena here so that all nodes have the same lifetime and can reference each other
    let node_arena = Arena::new();
    let mut nodes: Vec<&Node> = Vec::new();
    let mut rng = rand::thread_rng();
    let u_values: Vec<i64> = (0..config.grid_size).collect();
    let mut v_values: Vec<i64> = (0..config.grid_size).collect();
    v_values.shuffle(&mut rng);
    for (i, (u, v)) in u_values.iter().zip(v_values.iter()).enumerate() {
        let node = node_arena.alloc(Node::new(i, *u, *v));
        nodes.push(node);
    }

    // make causal set
    let causal_set = Configuration::new(&mut nodes, config.grid_size);
    causal_set.position_nodes();
    let mut rng = rand::thread_rng();

    let mut old_action = causal_set.action_bd(config.epsilon);
    let mut rejected = 0;
    let mut accepted = 0;
    for sweep in 0..config.nr_of_sweeps {
        if config.verbose {
            println!(" Starting sweep {} of {}", sweep + 1, config.nr_of_sweeps);
        }
        for _step in 0..steps_per_sweep {
            let coord: NodeCoordinate = rand::random();
            let (node1, node2) = causal_set.random_nodes2();
            // here we already alter the set
            causal_set.exchange_node_coordinate(node1, node2, coord);
            let new_action = causal_set.action_bd(config.epsilon);
            let delta_action = new_action - old_action;
            // if delta_action is negative accept, otherwise sample
            if delta_action > 0.0 {
                let exponent = -config.beta * delta_action;
                let random_nr: f64 = rng.sample(Standard);
                if exponent.exp() < random_nr {
                    // reject the move and undo it
                    causal_set.exchange_node_coordinate(node1, node2, coord);
                    rejected += 1;
                } else {
                    accepted += 1;
                    old_action = new_action;
                }
            } else {
                accepted += 1;
                old_action = new_action;
            }
            //causal_set.print_set();
        }

        let path_count_string = causal_set
            .interval_count()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join("\t");
        let action_value = causal_set.action_bd(config.epsilon);
        if config.verbose {
            println!(
                "  finished - (accepted, rejected, action) = ({},{},{})",
                accepted, rejected, action_value
            );
        }
        // println!("  {}", path_count_string);
        // causal_set.print_set();
        {
            let mut output_file = config.output_file.lock().unwrap();
            match writeln!(
                output_file,
                "{}\t{}\t{}\t{}\t{}",
                config.beta, config.epsilon, sweep, action_value, path_count_string
            ) {
                Ok(_) => (),
                Err(e) => panic!("Error writing to output file: {}", e),
            };
        }
    }
    Ok((accepted, rejected))
}