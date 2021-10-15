use ndarray::Array2;
use rand::prelude::*;
use rand::{
    distributions::{Distribution, Standard},
    Rng,
};
use std::cell::{Cell, RefCell};

#[derive(Debug, Clone, Copy)]
pub enum NodeCoordinate {
    U,
    V,
}

impl Distribution<NodeCoordinate> for Standard {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> NodeCoordinate {
        match rng.gen_range(0..=1) {
            0 => NodeCoordinate::U,
            _ => NodeCoordinate::V,
        }
    }
}

#[derive(Debug)]
pub struct Node {
    id: usize,
    u: Cell<i64>,
    v: Cell<i64>,
}

impl<'a> PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Node {
    pub fn new(id: usize, u: i64, v: i64) -> Node {
        Node {
            id: id,
            u: Cell::new(u),
            v: Cell::new(v),
        }
    }

    fn before(&self, second: &Node) -> bool {
        self.u > second.u && self.v > second.v
    }
    #[allow(dead_code)]
    fn after(&self, second: &Node) -> bool {
        self.u < second.u && self.v < second.v
    }
    fn related(&self, second: &Node) -> bool {
        (self.u > second.u && self.v > second.v) || (self.u < second.u && self.v < second.v)
    }
}

#[derive(Debug)]
pub struct Configuration<'a> {
    nodes: &'a mut Vec<&'a Node>,
    data: RefCell<Array2<Option<&'a Node>>>,
    grid_size: usize,
    grid_size_i: i64,
    seed: u64,
    steps: u64,
}

impl<'a> Configuration<'a> {
    pub fn new(nodes: &'a mut Vec<&'a Node>, grid_size: i64) -> Configuration<'a> {
        let config = Configuration {
            grid_size: grid_size as usize,
            grid_size_i: grid_size,
            seed: 42,
            steps: 0,
            nodes: nodes,
            data: RefCell::new(Array2::<Option<&Node>>::from_elem(
                (grid_size as usize, grid_size as usize),
                None,
            )),
        };
        config
    }

    pub fn text_representation(&self) -> String {
        self.nodes
            .iter()
            .map(|&n| format!("{},{},{}", n.id, n.u.get(), n.v.get()))
            .collect::<Vec<String>>()
            .join("\n")
    }

    pub fn position_nodes(&self) {
        for node1 in self.nodes.iter() {
            let mut data = self.data.borrow_mut();
            data[[node1.u.get() as usize, node1.v.get() as usize]] = Some(node1);
        }
    }

    pub fn print_set(&self) {
        let data = self.data.borrow();
        let nr_of_slices = 2 * self.grid_size - 1;
        for t in 0..nr_of_slices {
            let mut repr = String::new();
            let lower_half = t >= self.grid_size;

            let lower_bound = if lower_half {
                t - self.grid_size + 1
            } else {
                0
            };
            let upper_bound = if lower_half { self.grid_size } else { t + 1 };
            let white_spaces = if lower_half {
                t - self.grid_size + 1
            } else {
                self.grid_size - t
            };
            for _ in 0..white_spaces {
                repr += " "
            }
            repr += if lower_half { "\\  " } else { "/ " };
            for x in lower_bound..upper_bound {
                let u = x;
                let v = t - x;
                match data[[u, v]] {
                    Some(x) => repr += &x.id.to_string(),
                    None => repr += ".",
                }
                if t % 1 == 0 {
                    repr += " ";
                }
            }
            repr += if lower_half { " /" } else { "\\" };
            println!("{}", repr);
        }
    }

    #[allow(dead_code)]
    pub fn print_nodelist(&self) {
        for node in self.nodes.iter() {
            println!(
                "Node {} at (u,v)=({}, {}): ",
                node.id,
                node.u.get(),
                node.v.get()
            );
        }
    }

    pub fn random_nodes2(&self) -> (&'a Node, &'a Node) {
        let sample: Vec<_> = self
            .nodes
            .choose_multiple(&mut rand::thread_rng(), 2)
            .collect();
        (&sample[0], &sample[1])
    }

    pub fn exchange_node_coordinate(
        &self,
        first_node: &'a Node,
        second_node: &'a Node,
        coordinate: NodeCoordinate,
    ) {
        let mut data = self.data.borrow_mut();
        data[[first_node.u.get() as usize, first_node.v.get() as usize]] = None;
        data[[second_node.u.get() as usize, second_node.v.get() as usize]] = None;
        match coordinate {
            NodeCoordinate::U => {
                let temp = first_node.u.get();
                first_node.u.set(second_node.u.get());
                second_node.u.set(temp);
            }
            NodeCoordinate::V => {
                let temp = first_node.v.get();
                first_node.v.set(second_node.v.get());
                second_node.v.set(temp);
            }
        }
        data[[first_node.u.get() as usize, first_node.v.get() as usize]] = Some(first_node);
        data[[second_node.u.get() as usize, second_node.v.get() as usize]] = Some(second_node);
    }

    pub fn interval_cardinality(&self, first_node: &Node, second_node: &Node) -> i64 {
        if !first_node.related(second_node) {
            panic!("Only pass related nodes");
        }
        let mut nodes = 0;
        let u_min = std::cmp::min(first_node.u.get(), second_node.u.get());
        let u_max = std::cmp::max(first_node.u.get(), second_node.u.get());
        let v_min = std::cmp::min(first_node.v.get(), second_node.v.get());
        let v_max = std::cmp::max(first_node.v.get(), second_node.v.get());
        let data = self.data.borrow();
        for u in u_min..=u_max {
            for v in v_min..=v_max {
                match data[[u as usize, v as usize]] {
                    Some(_) => nodes += 1,
                    None => (),
                }
            }
        }
        // don't count the nodes themselves
        nodes - 2
    }

    pub fn interval_count(&self) -> Vec<i64> {
        let mut counts: Vec<i64> = Vec::new();
        // we are iterating over all nodes twice but are not double-counting as we are only counting towards the future
        for node1 in self.nodes.iter() {
            for node2 in self.nodes.iter() {
                if node1.before(node2) {
                    counts.push(self.interval_cardinality(node1, node2));
                }
            }
        }
        // grid_size - 1 is correct here: the vector goes from 0 to N-2
        counts.iter().fold(vec![0; self.grid_size - 1], |mut v, x| {
            v[*x as usize] += 1;
            v
        })
    }

    pub fn action_bd(&self, eps: f64) -> f64 {
        // see Eqn. (1) and (2) in 1706.06432
        let f = |n: f64, eps: f64| -> f64 {
            let temp: f64 = 1.0 - eps;
            temp.powf(n)
                * (1.0 - 2.0 * eps * n / temp + eps * eps * n * (n - 1.0) / 2.0 / temp / temp)
        };

        let counts = self.interval_count();

        let mut action = counts
            .iter()
            .enumerate()
            .map(|(n, n_n)| (*n_n as f64) * f(n as f64, eps))
            .sum();
        action = (self.grid_size_i as f64) - 2.0 * eps * action;
        action *= 4.0 * eps;
        action
    }
}
