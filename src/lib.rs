use std::fs::File;
use std::io::prelude::*;
use rand::prelude::*;
use rand::{
    distributions::{Standard},
    Rng
};
use typed_arena::Arena;
use yaml_rust::YamlLoader;

mod causal_set;
use crate::causal_set::{Node, NodeCoordinate, Configuration};

pub struct Config {
    grid_size: i64,
    epsilon: f64,
    beta: f64,
    nr_of_sweeps: i64,
    output_path: std::path::PathBuf,
}

impl Config {
    pub fn new_from_file(file_name: &std::path::Path) -> Config {
        let mut parameter_file = match File::open(file_name) {
            Ok(f)=> f,
            Err(e) => panic!("Error opening file parameters.yml in {}: {}", file_name.display(), e)
        };
        let mut parameter_file_string = String::new();
        parameter_file.read_to_string(&mut parameter_file_string).unwrap();
    
        let parameter_files = YamlLoader::load_from_str(&parameter_file_string).unwrap();
        let parameters = &parameter_files[0];
    
        let output_path = file_name.with_file_name("output.dat");
        Config{
            grid_size: parameters["nr_of_points"].as_i64().unwrap(),
            epsilon: parameters["epsilon"].as_f64().unwrap(),
            beta: parameters["beta"].as_f64().unwrap(),
            nr_of_sweeps: 5000,
            output_path: output_path,
        }
    }
}

pub fn run_simulation(config: Config) {
    println!("Starting run with N={} at eps={}, beta={}", config.grid_size, config.epsilon, config.beta);
    // setup nodes
    // we are using an arena here so that all nodes have the same lifetime and can reference each other
    let node_arena = Arena::new();
    let mut nodes: Vec<&Node> = Vec::new();
    let mut rng = rand::thread_rng();
    let u_values: Vec<i64> = (0..config.grid_size).collect();
    let mut v_values: Vec<i64> = (0..config.grid_size).collect();
    v_values.shuffle(&mut rng);
    for (i,(u, v)) in u_values.iter().zip(v_values.iter()).enumerate() {
        let node = node_arena.alloc(Node::new(i, *u, *v));
        nodes.push(node);
    }

    let mut output_file = match File::create(config.output_path) {
        Ok(f)=> f,
        Err(e) => panic!("Error creating output file: {}", e)
    };
    
    // make causal set
    let causal_set = Configuration::new(& mut nodes, config.grid_size);
    causal_set.position_nodes();
    causal_set.construct_relations();
    let mut rng = rand::thread_rng();
    
    // following Surya 2012
    let steps_per_sweep = config.grid_size * (config.grid_size - 1) / 2;

    match writeln!(output_file, "sweep\taction\tnr_edges") {
        Ok(_) => (),
        Err(e) => panic!("Error writing to output file: {}", e)
    };
    
    for sweep in 0..config.nr_of_sweeps {
        println!(" Starting sweep {} of {}", sweep, config.nr_of_sweeps);
        for step in 0..steps_per_sweep {
            let coord: NodeCoordinate = rand::random();
            let (node1, node2) = causal_set.random_nodes2();
            // here we already alter the set
            causal_set.exchange_node_coordinate(node1, node2, coord);
            causal_set.construct_relations();
            let exponent = - config.beta * causal_set.action_bd(config.epsilon);
            let random_nr: f64 = rng.sample(Standard);
            if exponent < random_nr.ln() {
                // reject the move and undo it
                // we do not construct the correction relations here as we only need those to compute the action
                causal_set.exchange_node_coordinate(node1, node2, coord);
            }
            if step % 100 == 0 {
                print!(".");
            }
        }
        println!("finished");
        causal_set.construct_relations();
        let path_count_string = causal_set.path_count().iter().map(|x| x.to_string()).collect::<Vec<String>>().join("\t");
        match writeln!(output_file, "{}\t{}\t{}", sweep, causal_set.action_bd(config.epsilon), path_count_string) {
            Ok(_) => (),
            Err(e) => panic!("Error writing to output file: {}", e)
        };
    }
    causal_set.print_set();
    causal_set.print_nodelist();

}