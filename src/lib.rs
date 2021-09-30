use rand::prelude::*;
use rand::{distributions::Standard, Rng};
use std::fs::File;
use std::io::prelude::*;
use typed_arena::Arena;
use yaml_rust::YamlLoader;

mod causal_set;
use crate::causal_set::{Configuration, Node, NodeCoordinate};

pub struct Config {
    grid_size: i64,
    epsilon: f64,
    beta: f64,
    nr_of_sweeps: i64,
    output_path: std::path::PathBuf,
    final_config_path: std::path::PathBuf,
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
        let final_config = match parameters["final_config_file"].as_str() {
            Some(val) => file_name.with_file_name(val),
            None => file_name.with_file_name("final_config.dat"),
        };
        let nr_of_sweeps = match parameters["nr_of_sweeps"].as_i64() {
            Some(val) => val,
            None => 5000,
        };

        Config {
            grid_size: parameters["nr_of_points"].as_i64().unwrap(),
            epsilon: parameters["epsilon"].as_f64().unwrap(),
            beta: parameters["beta"].as_f64().unwrap(),
            nr_of_sweeps: nr_of_sweeps,
            output_path: output_path,
            final_config_path: final_config,
        }
    }
}

pub fn run_simulation(config: Config) {
    println!(
        "Starting run with N={} at eps={}, beta={}",
        config.grid_size, config.epsilon, config.beta
    );
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

    let mut output_file = match File::create(config.output_path) {
        Ok(f) => f,
        Err(e) => panic!("Error creating output file: {}", e),
    };

    // make causal set
    let causal_set = Configuration::new(&mut nodes, config.grid_size);
    causal_set.position_nodes();
    let mut rng = rand::thread_rng();

    // following Surya 2012
    let steps_per_sweep = config.grid_size * (config.grid_size - 1) / 2;

    let n_string: Vec<String> = (0..=config.grid_size - 2)
        .map(|x| format!("N{}", x))
        .collect();
    let n_string = n_string.join("\t");
    match writeln!(output_file, "sweep\taction\t{}", n_string) {
        Ok(_) => (),
        Err(e) => panic!("Error writing to output file: {}", e),
    };

    let mut old_action = causal_set.action_bd(config.epsilon);
    for sweep in 0..config.nr_of_sweeps {
        let mut rejected = 0;
        let mut accepted = 0;
        println!(" Starting sweep {} of {}", sweep, config.nr_of_sweeps);
        for _step in 0..steps_per_sweep {
            let coord: NodeCoordinate = rand::random();
            let (node1, node2) = causal_set.random_nodes2();
            // here we already alter the set
            causal_set.exchange_node_coordinate(node1, node2, coord);
            let new_action = causal_set.action_bd(config.epsilon);
            let delta_action = new_action - old_action;
            old_action = new_action;
            // if delta_action is negative accept, otherwise sample
            if delta_action > 0.0 {
                let exponent = -config.beta * delta_action;
                let random_nr: f64 = rng.sample(Standard);
                if exponent.exp() < random_nr {
                    // reject the move and undo it
                    // we do not construct the relations here as we only need those to compute the action
                    causal_set.exchange_node_coordinate(node1, node2, coord);
                    rejected += 1;
                }
            } else {
                accepted += 1
            }
        }

        let path_count_string = causal_set
            .interval_count()
            .iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join("\t");
        let action_value = causal_set.action_bd(config.epsilon);
        println!(
            "  finished - (accepted, rejected, action) = ({},{},{})",
            accepted, rejected, action_value
        );
        // println!("  {}", path_count_string);
        // causal_set.print_set();
        match writeln!(
            output_file,
            "{}\t{}\t{}",
            sweep, action_value, path_count_string
        ) {
            Ok(_) => (),
            Err(e) => panic!("Error writing to output file: {}", e),
        };
    }
    write_configuration_to_file(&causal_set, &config.final_config_path);
    causal_set.print_set();
    causal_set.print_nodelist();
}

fn write_configuration_to_file(config: &Configuration, output_path: &std::path::Path) {
    let mut output_file = match File::create(output_path) {
        Ok(f) => f,
        Err(e) => panic!("Error creating configuration output file: {}", e),
    };
    match writeln!(output_file, "id,u,v") {
        Ok(_) => (),
        Err(e) => panic!("Error writing to configuration output file: {}", e),
    };
    match writeln!(output_file, "{}", config.text_representation()) {
        Ok(_) => (),
        Err(e) => panic!("Error writing to configuration output file: {}", e),
    };
}
