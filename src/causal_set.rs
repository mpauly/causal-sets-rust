use ndarray::Array2;
use rand::prelude::*;
use rand::{
    distributions::{Distribution, Standard},
    Rng
};
use std::cell::{Cell,RefCell};

#[derive(Debug,Clone,Copy)]
pub enum NodeCoordinate {
    U,
    V
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
pub struct Node<'a>{
    id: usize,
    u: Cell<i64>,
    v: Cell<i64>,
    parents: RefCell<Vec<&'a Node<'a>>>,
    children: RefCell<Vec<&'a Node<'a>>>,
}

impl <'a> PartialEq for Node<'a> {
    fn eq(&self, other: &Self) -> bool { self.id == other.id }
}

impl<'a> Node<'a> {
    pub fn new(id: usize, u: i64, v:i64) -> Node<'a> {
        Node {
            id: id,
            u: Cell::new(u), 
            v: Cell::new(v),
            parents: RefCell::new(Vec::new()),
            children: RefCell::new(Vec::new())
        }
    }

    fn t(&self)-> i64 { self.u.get() + self.v.get() }
    fn x(&self)-> i64 { self.u.get() - self.v.get() }
    fn before(&self, second: &Node<'a>)->bool { 
        self.t() < second.t() 
    }
    fn after(&self, second: &Node<'a>)->bool { 
        self.t() > second.t() 
    }
    fn related(&self, second: &Node<'a>) -> bool {
        let del_x = self.x() - second.x();
        let del_t = self.t() - second.t();
        if del_x * del_x > del_t * del_t{
            return false;
        }
        true
    }

    fn add_child(&self, node: &'a Node<'a>){
        self.children.borrow_mut().push(node);
    }

    fn add_parent(&self, node: &'a Node<'a>){
        self.parents.borrow_mut().push(node);
    }

    fn find_in_future(&self, node: &'a Node<'a>, count:i64) -> Vec<i64> {
        if self == node { return vec![count]; }
        return self.children.borrow().iter().map(|c|{c.find_in_future(node, count + 1)}).flatten().collect()
    }
}


#[derive(Debug)]
pub struct Configuration<'a> {
    nodes: &'a mut Vec<&'a Node<'a>>,
    data: RefCell<Array2::<Option<&'a Node<'a>>>>,
    grid_size: usize,
    grid_size_i: i64,
    seed: u64,
    steps: u64,
}

impl<'a> Configuration<'a> {
    pub fn new(nodes: &'a mut Vec<&'a Node<'a>>, grid_size: i64) -> Configuration<'a> {
        let config = Configuration {
            grid_size: grid_size as usize,
            grid_size_i: grid_size,
            seed: 42,
            steps: 0,
            nodes: nodes,
            data: RefCell::new(Array2::<Option<&Node>>::from_elem((grid_size as usize, grid_size as usize), None)),
        };
        config
    }

    fn nodes_directly_related(&self, first_node: &'a Node<'a>, second_node: &'a Node<'a>)-> bool {
        if first_node == second_node {return false}
        if !first_node.related(second_node) {return false}
        for node3 in self.nodes.iter() {
            if first_node.related(node3) && second_node.related(node3) {
                if first_node.before(node3) && second_node.after(node3) {
                    return false
                } else if first_node.after(node3) && second_node.before(node3) {
                    return false
                }
            }
        }
        true
    }

    pub fn position_nodes(&self) {
        for node1 in self.nodes.iter() {
            let mut data = self.data.borrow_mut();
            data[[node1.u.get() as usize, node1.v.get() as usize]]= Some(node1);
        }
    }

    pub fn construct_relations(& self) {
        // the following is O(N^3) and could maybe be faster by making use of the grid data structure
        for node1 in self.nodes.iter() {
            for node2 in self.nodes.iter() {
                if self.nodes_directly_related(node1, node2) {
                    if node1.before(&node2) { 
                        node1.add_child(&node2);
                    } else {
                        node1.add_parent(&node2);
                    }   
                }             
            }
        }
    }

    pub fn print_set(& self) {  
        let data = self.data.borrow();
        let nr_of_slices = 2 * self.grid_size - 1;
        for t in 0..nr_of_slices {
            let mut repr = String::new();
            let lower_half = t>=self.grid_size;
            
            let lower_bound = if lower_half {t - self.grid_size + 1} else {0}; 
            let upper_bound = if lower_half {self.grid_size} else { t + 1 };
            let white_spaces = if lower_half {t - self.grid_size + 1} else {self.grid_size -t};
            for _ in 0..white_spaces {
                repr += " "
            }
            repr += if lower_half {"\\  "} else {"/ "};
            for x in lower_bound..upper_bound {
                let u = x;
                let v = t - x;
                match data[[u,v]] {
                    Some(x)=> repr += &x.id.to_string() ,
                    None => repr += "."
                }
                if t % 1 == 0 {repr += " ";}
            }
            repr += if lower_half {" /"} else {"\\"};
            println!("{}",repr);
        }
        
    }

    pub fn print_nodelist(& self) {  
        for node in self.nodes.iter() {
            println!("Node {} at (u,v)=({}, {}): ", node.id, node.u.get(), node.v.get());
            let parents: Vec<_> = node.parents.borrow().iter().map(|n| n.id.to_string() ).collect();
            println!("  - Parents: {}", parents.join(","));
            let children: Vec<_> = node.children.borrow().iter().map(|n| n.id.to_string() ).collect();
            println!("  - Children: {}", children.join(","));
        }
    }

    pub fn random_nodes2(& self) -> (&'a Node<'a>, &'a Node<'a>) {
        let sample: Vec<_> = self.nodes
            .choose_multiple(&mut rand::thread_rng(), 2)
            .collect();
        (& sample[0], & sample[1])
    }

    pub fn exchange_node_coordinate(&self, first_node: &'a Node<'a>, second_node: &'a Node<'a>, coordinate: NodeCoordinate) {
        let mut data = self.data.borrow_mut();
        data[[first_node.u.get() as usize, first_node.v.get() as usize]] = None;
        data[[second_node.u.get() as usize, second_node.v.get() as usize]] = None;
        match coordinate {
            NodeCoordinate::U => {
                let temp = first_node.u.get();
                first_node.u.set(second_node.u.get());
                second_node.u.set(temp);
            },
            NodeCoordinate::V => {
                let temp = first_node.v.get();
                first_node.v.set(second_node.v.get());
                second_node.v.set(temp);
            }
        }
        data[[first_node.u.get() as usize, first_node.v.get() as usize]]= Some(first_node);
        data[[second_node.u.get() as usize, second_node.v.get() as usize]]= Some(second_node);
        self.clear_children_and_parent_info();
    }

    fn clear_children_and_parent_info(&self) {
        for node in self.nodes.iter() {
            node.children.borrow_mut().clear();
            node.parents.borrow_mut().clear();
        }
    }

    pub fn path_count(&self) -> Vec<i64> {
        let mut counts : Vec<i64> = Vec::new();
        // we are iterating over all nodes twice but are not double-counting as we are only counting toward the future
        for node1 in self.nodes.iter() {
            for node2 in self.nodes.iter() {
                counts.extend(node1.find_in_future(node2, 0));
            }
        }
        // find_in_future counts intermediate edges, need to switch to intermediate nodes
        counts.iter()
            .filter(|&x| *x > 0) // get rid of the number of nodes
            .map(|x| x-1) // switch to number of intermediate nodes 
            .fold(vec![0;self.grid_size - 1], |mut v, x| {v[x as usize]+= 1; v}) // count which number occurs how often
    }

    pub fn action_bd(&self, eps:f64) -> f64 {
        // see Eqn. (1) and (2) in 1706.06432
        let f = |n: f64, eps:f64| -> f64 {
            let temp : f64 = 1.0-eps;
            temp.powf(n) * (1.0 - 2.0 * eps * n / temp + eps * eps * n * (n-1.0) / 2.0 /temp / temp) 
        };

        let counts = self.path_count();

        let mut action = counts.iter().enumerate().map(|(n, n_n)| (*n_n as f64) * f(n as f64, eps)).sum();
        action = (self.grid_size_i as f64) - 2.0 * eps * action; 
        action *= 4.0 * eps;
        action
    }
}