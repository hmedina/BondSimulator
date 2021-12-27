pub mod agent_types {
    use petgraph::prelude::*;

    #[derive(Debug)]
    pub struct AxnNode {
        pub x1_bond: Option<EdgeIndex>,
        pub x2_bond: Option<EdgeIndex>,
        pub xp_bond: Option<EdgeIndex>
    }   

    impl Default for AxnNode {
        fn default() -> Self {
            AxnNode {
               x1_bond: None,
               x2_bond: None,
               xp_bond: None
            }
        }
    }   

    #[derive(Debug)]
    pub struct ApcNode {
        pub p1_bond: Option<EdgeIndex>,
        pub p2_bond: Option<EdgeIndex>,
        pub p3_bond: Option<EdgeIndex>,
        pub pp_bond: Option<EdgeIndex>
    }   

    impl Default for ApcNode {
        fn default() -> Self {
            ApcNode {
               p1_bond: None,
               p2_bond: None,
               p3_bond: None,
               pp_bond: None
            }
        }
    }

    #[derive(Clone, Debug)]
    pub enum AgentType {
        AxnNode(Option<EdgeIndex>, Option<EdgeIndex>, Option<EdgeIndex>),
        ApcNode(Option<EdgeIndex>, Option<EdgeIndex>, Option<EdgeIndex>, Option<EdgeIndex>)
    }
}

pub mod edge_ends {
    use std::cmp::Ordering;
    use petgraph::prelude::*;

    #[derive(Clone, Copy, Debug)]
    pub struct EdgeEnds {
        // holds the node indexes for the bond termini
        // char holds the site name; this allows "parallel" edge resolution via site-name awareness
        pub a: NodeIndex,
        pub b: NodeIndex,
        pub a_s: char,
        pub b_s: char,
        pub z: EdgeIndex
    }
    
    impl PartialEq for EdgeEnds {
        fn eq(&self, other: &Self) -> bool {
            (self.a == other.a && self.b == other.b && self.a_s == other.a_s && self.b_s == other.b_s) || (self.a == other.b && self.b == other.a && self.a_s == other.b_s && self.b_s == other.a_s)
        }
    }
    
    impl Eq for EdgeEnds {}
    
    impl PartialOrd for EdgeEnds {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }
    
    impl Ord for EdgeEnds {
        fn cmp(&self, other: &Self) -> Ordering {
            if self.a == other.a { self.b.cmp(&other.b) }
            else {self.a.cmp(&other.a)}
        }
    }
}

pub mod open_ports {
    use std::collections::HashSet;
    use petgraph::prelude::NodeIndex;

    #[derive(Debug)]
    pub struct OpenPorts {
        // sets of NodeIndexes, from the universe_graph
        pub xh_free: HashSet<NodeIndex>,
        pub xt_free: HashSet<NodeIndex>,
        pub xp_free: HashSet<NodeIndex>,
        pub p1_free: HashSet<NodeIndex>,
        pub p2_free: HashSet<NodeIndex>,
        pub p3_free: HashSet<NodeIndex>,
        pub pp_free: HashSet<NodeIndex>
    }
    
    impl OpenPorts {
        pub fn update_from(&mut self, other: &Self) {
            for item in &other.xh_free {
                self.xh_free.insert(*item);
            }
            for item in &other.xt_free {
                self.xh_free.insert(*item);
            }
            for item in &other.xp_free {
                self.xh_free.insert(*item);
            }
            for item in &other.p1_free {
                self.xh_free.insert(*item);
            }
            for item in &other.p2_free {
                self.xh_free.insert(*item);
            }
            for item in &other.p3_free {
                self.xh_free.insert(*item);
            }
            for item in &other.pp_free {
                self.xh_free.insert(*item);
            }
        }

        pub fn default_axn(ix: NodeIndex) -> Self {
            let mut s: HashSet<NodeIndex> = HashSet::new();
            s.insert(ix);
            OpenPorts { 
                xh_free: s.clone(),
                xt_free: s.clone(),
                xp_free: s,
                p1_free: HashSet::new(),
                p2_free: HashSet::new(),
                p3_free: HashSet::new(),
                pp_free: HashSet::new()
            }
        }

        pub fn default_apc(ix: NodeIndex) -> Self {
            let mut s: HashSet<NodeIndex> = HashSet::new();
            s.insert(ix);
            OpenPorts {
                xh_free: HashSet::new(),
                xt_free: HashSet::new(),
                xp_free: HashSet::new(),
                p1_free: s.clone(),
                p2_free: s.clone(),
                p3_free: s.clone(),
                pp_free: s
            }
        }

        pub fn default_empty() -> Self {
            OpenPorts {
                xh_free: HashSet::new(),
                xt_free: HashSet::new(),
                xp_free: HashSet::new(),
                p1_free: HashSet::new(),
                p2_free: HashSet::new(),
                p3_free: HashSet::new(),
                pp_free: HashSet::new()
            }
        }
    }
}

pub mod edge_types {
    use std::collections::BTreeSet;
    use std::rc::Rc;
    use std::cell::RefCell;
    use crate::edge_ends::EdgeEnds;

    #[derive(Debug)]
    pub struct EdgeTypes {
        // sets of edges, defined by their ends (node id & site names)
        pub xh_xt: BTreeSet<Rc<RefCell<EdgeEnds>>>,
        pub p1_xp: BTreeSet<Rc<RefCell<EdgeEnds>>>,
        pub p2_xp: BTreeSet<Rc<RefCell<EdgeEnds>>>,
        pub p3_xp: BTreeSet<Rc<RefCell<EdgeEnds>>>,
        pub pp_pp: BTreeSet<Rc<RefCell<EdgeEnds>>>
    }
    
    impl EdgeTypes {
        pub fn update_from(&mut self, other: &Self) {
            for item in &other.xh_xt {
                self.xh_xt.insert(Rc::clone(&item));
            }
            for item in &other.p1_xp {
                self.xh_xt.insert(Rc::clone(&item));
            }
            for item in &other.p2_xp {
                self.xh_xt.insert(Rc::clone(&item));
            }
            for item in &other.p3_xp {
                self.xh_xt.insert(Rc::clone(&item));
            }
            for item in &other.pp_pp {
                self.xh_xt.insert(Rc::clone(&item));
            }
        }

        pub fn default_empty() -> Self {
            EdgeTypes {
                xh_xt: BTreeSet::new(),
                p1_xp: BTreeSet::new(),
                p2_xp: BTreeSet::new(),
                p3_xp: BTreeSet::new(),
                pp_pp: BTreeSet::new()
            }
        }
    }
}

pub mod rule_activities {
    pub struct RuleActivities {
        pub axn_axn_u_bind: MassActionTerm,
        pub axn_axn_b_bind: MassActionTerm,
        pub axn_axn_u_free: MassActionTerm,
        pub axn_axn_b_free: MassActionTerm,
        pub ap1_axn_u_bind: MassActionTerm,
        pub ap1_axn_b_bind: MassActionTerm,
        pub ap1_axn_u_free: MassActionTerm,
        pub ap1_axn_b_free: MassActionTerm,
        pub ap2_axn_u_bind: MassActionTerm,
        pub ap2_axn_b_bind: MassActionTerm,
        pub ap2_axn_u_free: MassActionTerm,
        pub ap2_axn_b_free: MassActionTerm,
        pub ap3_axn_u_bind: MassActionTerm,
        pub ap3_axn_b_bind: MassActionTerm,
        pub ap3_axn_u_free: MassActionTerm,
        pub ap3_axn_b_free: MassActionTerm,
        pub apc_apc_u_bind: MassActionTerm,
        pub apc_apc_b_bind: MassActionTerm,
        pub apc_apc_u_free: MassActionTerm,
        pub apc_apc_b_free: MassActionTerm,
    }

    pub struct MassActionTerm {
        pub mass: usize,
        pub rate: f64
    }

    pub struct RuleRates {
        pub axn_axn_u_bind: f64,
        pub axn_axn_b_bind: f64,
        pub axn_axn_u_free: f64,
        pub axn_axn_b_free: f64,
        pub ap1_axn_u_bind: f64,
        pub ap1_axn_b_bind: f64,
        pub ap1_axn_u_free: f64,
        pub ap1_axn_b_free: f64,
        pub ap2_axn_u_bind: f64,
        pub ap2_axn_b_bind: f64,
        pub ap2_axn_u_free: f64,
        pub ap2_axn_b_free: f64,
        pub ap3_axn_u_bind: f64,
        pub ap3_axn_b_bind: f64,
        pub ap3_axn_u_free: f64,
        pub ap3_axn_b_free: f64,
        pub apc_apc_u_bind: f64,
        pub apc_apc_b_bind: f64,
        pub apc_apc_u_free: f64,
        pub apc_apc_b_free: f64,
    }
}

pub mod unary_embeds {
    use petgraph::prelude::NodeIndex;
    use nalgebra::*;

    #[derive(Debug)]
    pub struct UnaryEmbeds {
        pub xh_xt: DMatrix<bool>,
        pub p1_xp: DMatrix<bool>,
        pub p2_xp: DMatrix<bool>,
        pub x3_xp: DMatrix<bool>,
        pub pp_pp: DMatrix<bool>
    }

    impl UnaryEmbeds {
        pub fn new_from_node(node: NodeIndex) -> Self {
            UnaryEmbeds {
                xh_xt: DMatrix::from_element()
            }
        }
    }
}

pub mod reaction_mixture {
    use crate::open_ports::OpenPorts;
    use crate::edge_types::EdgeTypes;
    use crate::edge_ends::EdgeEnds;
    use crate::agent_types::{AgentType, AxnNode, ApcNode};
    use crate::rule_activities::{RuleActivities, RuleRates, MassActionTerm};
    use crate::unary_embeds::UnaryEmbeds;
    use petgraph::{graph::Graph, algo::astar::astar, prelude::*};
    use rand::prelude::*;
    use std::cell::{RefCell, RefMut};
    use std::cmp::Ordering;
    use std::collections::{VecDeque, BTreeMap, BTreeSet, HashSet};
    use std::rc::Rc;

    #[derive(Debug)]
    struct MixtureSpecies {
        ports: OpenPorts,
        edges: EdgeTypes,
        size: usize,
        agent_set: BTreeSet<NodeIndex>
    }

    impl Eq for MixtureSpecies {}

    impl PartialEq for MixtureSpecies {
        fn eq(&self, other: &Self) -> bool {
            (self.agent_set == other.agent_set) &&
            (self.edges.xh_xt == other.edges.xh_xt) &&
            (self.edges.p1_xp == other.edges.p1_xp) &&
            (self.edges.p2_xp == other.edges.p2_xp) &&
            (self.edges.p3_xp == other.edges.p3_xp) &&
            (self.edges.pp_pp == other.edges.pp_pp)
        }
    }

    impl Ord for MixtureSpecies {
        fn cmp(&self, other: &Self) -> Ordering {
            if self.agent_set == other.agent_set {
                if self.edges.xh_xt == other.edges.xh_xt {
                    if self.edges.p1_xp == other.edges.p1_xp {
                        if self.edges.p2_xp == other.edges.p2_xp {
                            if self.edges.p3_xp == other.edges.p3_xp {
                                self.edges.pp_pp.cmp(&other.edges.pp_pp)
                            }
                            else {self.edges.p3_xp.cmp(&other.edges.p3_xp)}
                        }
                        else {self.edges.p2_xp.cmp(&other.edges.p2_xp)}
                    }
                    else {self.edges.p1_xp.cmp(&other.edges.p2_xp)}
                }
                else {self.edges.xh_xt.cmp(&other.edges.xh_xt)}
            }
            else {self.agent_set.cmp(&other.agent_set)}
        }
    }

    impl PartialOrd for MixtureSpecies {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }

    pub struct Mixture {
        universe_graph: Graph<AgentType, bool, Undirected>,
        species_annots: BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies>>>,   // key is agent index, value is reference-counted pointer to a reference-cell, that holds the mixture species data
        species_set: VecDeque<Rc<RefCell<MixtureSpecies>>>,                 // cache, iteration, & "final persistent owner"
        ports: OpenPorts,                                                   // struct with vector, each a list of agent indexes
        edges: EdgeTypes,
        //unary_binding_pairs: UnaryEmbeds,  //how to store the unary embeds?
        rule_activities: RuleActivities
    }
    
    impl Mixture {
        pub fn to_kappa(self) -> String {
            // printer
            let mut mixture_string: String = String::new();
            mixture_string.push_str("# synthetic mixture\n\n");
            let debug_str = format!("# mixture has {} nodes and {} edges\n", self.universe_graph.node_count(), self.universe_graph.edge_count());
            mixture_string.push_str(&debug_str);
            for this_species_ref in &self.species_set {
                let this_species = this_species_ref.borrow();
                let size_annot: String = format!("init: 1 /*{} agents*/ ", this_species.size);
                let mut species_string: String = String::new();
                species_string.push_str(&size_annot);
                for agent_id_ref in &this_species.agent_set {
                    let agent_id: NodeIndex = agent_id_ref.clone();
                    let mut s = format!("x{}:", agent_id.index());
                    let agent_str: String = match self.universe_graph.node_weight(agent_id).unwrap() {
                        AgentType::AxnNode(x1, x2, xp) => {
                            s.push_str("Axn(x1[");
                            match x1 {
                                Some(i) => s.push_str(i.index().to_string().as_str()),
                                None => s.push_str(".")
                            };
                            s.push_str("], x2[");
                            match x2 {
                                Some(i) => s.push_str(i.index().to_string().as_str()),
                                None => s.push_str(".")
                            };
                            s.push_str("], p[");
                            match xp {
                                Some(i) => s.push_str(i.index().to_string().as_str()),
                                None => s.push_str(".")
                            };
                            s.push_str("]), ");
                            s
                        }
                        AgentType::ApcNode(x1, x2, x3, pp) => {
                            s.push_str("APC(ax1[");
                            match x1 {
                                Some(i) => s.push_str(i.index().to_string().as_str()),
                                None => s.push_str(".")
                            };
                            s.push_str("], x2[");
                            match x2 {
                                Some(i) => s.push_str(i.index().to_string().as_str()),
                                None => s.push_str(".")
                            };
                            s.push_str("], x3[");
                            match x3 {
                                Some(i) => s.push_str(i.index().to_string().as_str()),
                                None => s.push_str(".")
                            };
                            s.push_str("], p[");
                            match pp {
                                Some(i) => s.push_str(i.index().to_string().as_str()),
                                None => s.push_str(".")
                            };
                            s.push_str("]), ");
                            s
                        }
                    };
                    species_string.push_str(&agent_str);
                }
                mixture_string.push_str(species_string.strip_suffix(", ").unwrap());
                mixture_string.push_str("\n");
            }
            mixture_string            
        }

        pub fn new_from_monomers(x_mass: usize, p_mass: usize, rule_rates: RuleRates) -> Mixture {
            // from initial list of species, build network
            // build list of species and each's cache
            // define initial abundances for tracking of mass action
            
            // the universe graph
            let mut net: Graph<AgentType, bool, Undirected> = Graph::with_capacity(x_mass + p_mass, x_mass * 3 + p_mass * 4);
            // the species map & global open ports caches
            let mut spc: BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies>>> = BTreeMap::new();
            let mut sps: VecDeque<Rc<RefCell<MixtureSpecies>>> = VecDeque::with_capacity(x_mass + p_mass);
            let mut opr = OpenPorts::default_empty();
            for _i in 0..x_mass {
                let node_ix = net.add_node(AgentType::AxnNode(None, None, None));
                let mut node_set = BTreeSet::new();
                assert!(node_set.insert(node_ix));
                let species_data = Rc::new(RefCell::new(MixtureSpecies{
                    ports: OpenPorts::default_axn(node_ix),
                    edges: EdgeTypes::default_empty(),
                    size: 1,
                    agent_set: node_set}));
                assert_eq!(None, spc.insert(node_ix, Rc::clone(&species_data)));
                sps.push_back(Rc::clone(&species_data));
                opr.xh_free.insert(node_ix);
                opr.xt_free.insert(node_ix);
                opr.xp_free.insert(node_ix);
            }
            for _i in 0..p_mass {
                let node_ix = net.add_node(AgentType::ApcNode(None, None, None, None));
                let mut node_set = BTreeSet::new();
                assert!(node_set.insert(node_ix));
                let species_data = Rc::new(RefCell::new(MixtureSpecies{
                    ports: OpenPorts::default_apc(node_ix),
                    edges: EdgeTypes::default_empty(),
                    size: 1,
                    agent_set: node_set}));
                assert_eq!(None, spc.insert(node_ix, Rc::clone(&species_data)));
                sps.push_back(Rc::clone(&species_data));
                opr.p1_free.insert(node_ix);
                opr.p2_free.insert(node_ix);
                opr.p3_free.insert(node_ix);
                opr.pp_free.insert(node_ix);
            }
            // edge cache
            let edg = EdgeTypes::default_empty();
            // rule activities
            let rac = RuleActivities {
                // all monomeric: there are no unary binding opportunities
                axn_axn_u_bind: MassActionTerm{mass: 0, rate: rule_rates.axn_axn_u_bind},
                ap1_axn_u_bind: MassActionTerm{mass: 0, rate: rule_rates.ap1_axn_u_bind},
                ap2_axn_u_bind: MassActionTerm{mass: 0, rate: rule_rates.ap2_axn_u_bind},
                ap3_axn_u_bind: MassActionTerm{mass: 0, rate: rule_rates.ap3_axn_u_bind},
                apc_apc_u_bind: MassActionTerm{mass: 0, rate: rule_rates.apc_apc_u_bind},
                // all monomeric: no bonds to break, cycles or no cycles
                axn_axn_u_free: MassActionTerm{mass: 0, rate: rule_rates.axn_axn_u_free},
                axn_axn_b_free: MassActionTerm{mass: 0, rate: rule_rates.axn_axn_b_free},
                ap1_axn_u_free: MassActionTerm{mass: 0, rate: rule_rates.ap1_axn_u_free},
                ap1_axn_b_free: MassActionTerm{mass: 0, rate: rule_rates.ap1_axn_b_free},
                ap2_axn_u_free: MassActionTerm{mass: 0, rate: rule_rates.ap2_axn_u_free},
                ap2_axn_b_free: MassActionTerm{mass: 0, rate: rule_rates.ap2_axn_b_free},
                ap3_axn_u_free: MassActionTerm{mass: 0, rate: rule_rates.ap3_axn_u_free},
                ap3_axn_b_free: MassActionTerm{mass: 0, rate: rule_rates.ap3_axn_b_free},
                apc_apc_u_free: MassActionTerm{mass: 0, rate: rule_rates.apc_apc_u_free},
                apc_apc_b_free: MassActionTerm{mass: 0, rate: rule_rates.apc_apc_b_free},
                // all monomeric: everything that can happen is of binary binding type
                axn_axn_b_bind: MassActionTerm{mass: x_mass * x_mass, rate: rule_rates.axn_axn_b_bind},
                ap1_axn_b_bind: MassActionTerm{mass: p_mass * x_mass, rate: rule_rates.ap1_axn_b_bind},
                ap2_axn_b_bind: MassActionTerm{mass: p_mass * x_mass, rate: rule_rates.ap2_axn_b_bind},
                ap3_axn_b_bind: MassActionTerm{mass: p_mass * x_mass, rate: rule_rates.ap3_axn_b_bind},
                apc_apc_b_bind: MassActionTerm{mass: p_mass * p_mass, rate: rule_rates.apc_apc_b_bind}
            };
            Mixture {universe_graph: net, species_annots: spc, ports: opr, edges: edg, rule_activities: rac, species_set: sps}
        }
    
        pub fn axn_axn_unary_bind(&mut self) {
            // identify pair to bind
            let a: NodeIndex = self.ports.xh_free.iter().choose(&mut rand::thread_rng()).unwrap().to_owned();
            let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&a).unwrap().borrow_mut();
            let b: NodeIndex = target_species.ports.xt_free.iter().choose(&mut rand::thread_rng()).unwrap().to_owned();
            // create edge & update caches
            let new_edge_index = self.universe_graph.add_edge(a, b, true);
            let new_edge = Rc::new(RefCell::new(EdgeEnds{a: a, b: b, a_s: 'h', b_s: 't', z: new_edge_index}));
            let node_a: &mut AgentType = self.universe_graph.node_weight_mut(a).unwrap();
            if let AgentType::AxnNode(x1, _, _) = node_a {*x1 = Some(new_edge_index)}
            let node_b = self.universe_graph.node_weight_mut(b).unwrap();
            if let AgentType::AxnNode(_, x2, _) = node_b {*x2 = Some(new_edge_index)}
            target_species.edges.xh_xt.insert(Rc::clone(&new_edge));
            self.edges.xh_xt.insert(Rc::clone(&new_edge));
            self.ports.xh_free.remove(&a);
            self.ports.xt_free.remove(&b);
            target_species.ports.xh_free.remove(&a);
            target_species.ports.xt_free.remove(&b);
            // helper closure for bond resilience updates
            let mut resilience_iterate_typed = |edge_iter: &BTreeSet<std::rc::Rc<std::cell::RefCell<EdgeEnds>>>| {
                for boxed_edge in edge_iter {
                    let edge = boxed_edge.borrow();
                    let edge_id = self.universe_graph.find_edge(edge.a, edge.b).unwrap();
                    // check bond's weight, which is a boolean, marking if the bond is a known cycle-member
                    if ! self.universe_graph.edge_weight(edge_id).unwrap() {
                        let mut counterfactual_graph = self.universe_graph.clone();
                        counterfactual_graph.remove_edge(edge_id);
                        let path = astar(&counterfactual_graph, edge.a, |finish| finish == edge.b, |_| 1, |_| 0);
                        match path {
                            Some(_) => {self.universe_graph.update_edge(edge.a, edge.b, true);},
                            None => ()
                        }
                    }
                }
            };
            // update resilience of that species' edges
            resilience_iterate_typed(&target_species.edges.xh_xt);
            resilience_iterate_typed(&target_species.edges.p1_xp);
            resilience_iterate_typed(&target_species.edges.p2_xp);
            resilience_iterate_typed(&target_species.edges.p3_xp);
            resilience_iterate_typed(&target_species.edges.pp_pp);
            // update rule activities
            self.rule_activities.axn_axn_u_bind.mass -= 1;
            self.rule_activities.axn_axn_u_free.mass -= 1;
        }
    
        pub fn axn_axn_binary_bind(&mut self) {
            // identify pair
            let a: NodeIndex = self.ports.xh_free.iter().choose(&mut rand::thread_rng()).unwrap().to_owned();
            let invalid_targets: &HashSet<NodeIndex> = &self.species_annots.get(&a).unwrap().borrow().ports.xt_free.to_owned();
            let b: NodeIndex = self.ports.xt_free.difference(invalid_targets).into_iter().choose(&mut rand::thread_rng()).unwrap().to_owned();
            // define which will be the "host" species, receiving the cache data from "eaten" species
            let (host_index, eaten_index) = if self.species_annots.get(&a).unwrap().borrow().size >= self.species_annots.get(&b).unwrap().borrow().size {(a, b)} else {(b, a)};
            // update rule activities
            self.rule_activities.axn_axn_b_free.mass += 1;  //the new bond
            //self.rule_activities.axn_axn_u_free.mass      // bi -> uni binding does not change the number of cycles in the system
            self.rule_activities.axn_axn_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.xh_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xt_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.xh_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xt_free.len());    // potential bindings that were binary, are now unary
            self.rule_activities.axn_axn_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.xh_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xt_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.xh_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xt_free.len());    // and so one counter increases by the other's loss
            self.rule_activities.ap1_axn_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap1_axn_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap2_axn_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap2_axn_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap3_axn_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap3_axn_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.apc_apc_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free.len());
            self.rule_activities.apc_apc_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free.len());
            // merge caches held by the species; uses interior-mutability pattern
            self.species_annots.get(&host_index).unwrap().borrow_mut().ports.update_from(&self.species_annots.get(&eaten_index).unwrap().borrow().ports);
            self.species_annots.get(&host_index).unwrap().borrow_mut().edges.update_from(&self.species_annots.get(&eaten_index).unwrap().borrow().edges);
            self.species_annots.get(&host_index).unwrap().borrow_mut().size += self.species_annots.get(&eaten_index).unwrap().borrow().size;
            self.species_annots.get(&host_index).unwrap().borrow_mut().agent_set.append(&mut self.species_annots.get(&eaten_index).unwrap().borrow().agent_set.clone());
            // create edge & update terminus trackers
            let new_edge_index = self.universe_graph.add_edge(a, b, false);
            let new_edge = Rc::new(RefCell::new(EdgeEnds{a: a, b:b, a_s: 'h', b_s: 't', z: new_edge_index}));
            self.species_annots.get_mut(&host_index).unwrap().borrow_mut().edges.xh_xt.insert(Rc::clone(&new_edge));
            self.species_annots.get_mut(&host_index).unwrap().borrow_mut().ports.xh_free.remove(&a);
            self.species_annots.get_mut(&host_index).unwrap().borrow_mut().ports.xt_free.remove(&b);
            self.edges.xh_xt.insert(Rc::clone(&new_edge));
            self.ports.xh_free.remove(&a);
            self.ports.xt_free.remove(&b);
            // update species annotation tracker
            let spec_ix: usize = self.species_set.binary_search(self.species_annots.get(&eaten_index).unwrap()).unwrap();
            self.species_set.remove(spec_ix);
            //assert!(self.species_set.remove());
            let indexes_of_eaten: BTreeSet<NodeIndex> = self.species_annots.get(&eaten_index).unwrap().borrow().agent_set.clone();
            for agent_index in indexes_of_eaten {
                let new_ref = Rc::clone(&self.species_annots.get(&host_index).unwrap());
                self.species_annots.entry(agent_index).and_modify(|e| {*e = new_ref});
            }
            // update bond identifier trackers
            let node_a_weight = self.universe_graph.node_weight_mut(a).unwrap();
            match node_a_weight {
                AgentType::AxnNode(x1, _, _) => *x1 = Some(new_edge_index),
                AgentType::ApcNode(_, _, _, _) => panic!("Node weight for an Axin head matched and APC!")
            };
            let node_b_weight = self.universe_graph.node_weight_mut(b).unwrap();
            match node_b_weight {
                AgentType::AxnNode(_, x2, _) => *x2 = Some(new_edge_index),
                AgentType::ApcNode(_, _, _, _) => panic!("Node weight for an Axin tail matched an APC!")
            };
        }
    }
}