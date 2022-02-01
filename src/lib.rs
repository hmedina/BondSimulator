pub mod agent_types {
    use petgraph::prelude::*;

    #[derive(Debug)]
    /**
     * Structure for an Axn type-node. The fields hold either a None when unbound, or a 
     * `Some(EdgeIndex)` when bound. One field per binding site.
    */
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
    /**
     * Structure for an APC type-node. The fields hold either a None when unbound, or a
     * `Some(EdgeIndex)` when bound. One field per binding site.
    */
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
    /** 
     * Structure whose fields store the NodeIndexes of the binding pair in `a` & `b`, the site type
     * (to resolve parallel edges), and the EdgeIndex in `z` from the UniverseGraph.
    */
    pub struct EdgeEnds {
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
    use std::collections::{HashSet, BTreeSet};
    use petgraph::prelude::NodeIndex;

    #[derive(Debug)]
    /** 
     * Structure of HashSets, one per site type. The HashSets hold NodeIndexes, from the UniverseGraph.
     */
    pub struct OpenPorts {
        pub xh_free: HashSet<NodeIndex>,
        pub xt_free: HashSet<NodeIndex>,
        pub xp_free: HashSet<NodeIndex>,
        pub p1_free: HashSet<NodeIndex>,
        pub p2_free: HashSet<NodeIndex>,
        pub p3_free: HashSet<NodeIndex>,
        pub pp_free: HashSet<NodeIndex>
    }
    
    impl OpenPorts {
        /**
         * Update every site-specific HashSet with the elements from the other.
        */
        pub fn update_from(&mut self, other: &Self) {
            for item in &other.xh_free {
                self.xh_free.insert(*item);
            }
            for item in &other.xt_free {
                self.xt_free.insert(*item);
            }
            for item in &other.xp_free {
                self.xp_free.insert(*item);
            }
            for item in &other.p1_free {
                self.p1_free.insert(*item);
            }
            for item in &other.p2_free {
                self.p2_free.insert(*item);
            }
            for item in &other.p3_free {
                self.p3_free.insert(*item);
            }
            for item in &other.pp_free {
                self.pp_free.insert(*item);
            }
        }

        /**
         * Partition the site-specific HashSets using a set of NodeIndexes, the "ejected" nodes.
         * A free site belonging to one of the agents in the ejected set, gets ejected and returned
         * in a new OpenPort object. Sites not belonging to agents in the ejected set are retained.
        */
        pub fn eject_where(&mut self, ejected_nodes: &BTreeSet<NodeIndex>) -> Self {
            // this would be simplified by "drain_filter" moving out of the Nightly-only
            let mut xh_ejected: HashSet<NodeIndex> = HashSet::new();
            let mut xt_ejected: HashSet<NodeIndex> = HashSet::new();
            let mut xp_ejected: HashSet<NodeIndex> = HashSet::new();
            let mut p1_ejected: HashSet<NodeIndex> = HashSet::new();
            let mut p2_ejected: HashSet<NodeIndex> = HashSet::new();
            let mut p3_ejected: HashSet<NodeIndex> = HashSet::new();
            let mut pp_ejected: HashSet<NodeIndex> = HashSet::new();
            let mut xh_kept: HashSet<NodeIndex> = HashSet::new();
            let mut xt_kept: HashSet<NodeIndex> = HashSet::new();
            let mut xp_kept: HashSet<NodeIndex> = HashSet::new();
            let mut p1_kept: HashSet<NodeIndex> = HashSet::new();
            let mut p2_kept: HashSet<NodeIndex> = HashSet::new();
            let mut p3_kept: HashSet<NodeIndex> = HashSet::new();
            let mut pp_kept: HashSet<NodeIndex> = HashSet::new();
            for item in self.xh_free.drain() {
                if ejected_nodes.contains(&item) {assert!(xh_ejected.insert(item))}
                else {assert!(xh_kept.insert(item))}
            }
            for item in self.xt_free.drain() {
                if ejected_nodes.contains(&item) {assert!(xt_ejected.insert(item))}
                else {assert!(xt_kept.insert(item))}
            }
            for item in self.xp_free.drain() {
                if ejected_nodes.contains(&item) {assert!(xp_ejected.insert(item))}
                else {assert!(xp_kept.insert(item))}
            }
            for item in self.p1_free.drain() {
                if ejected_nodes.contains(&item) {assert!(p1_ejected.insert(item))}
                else {assert!(p1_kept.insert(item))}
            }
            for item in self.p2_free.drain() {
                if ejected_nodes.contains(&item) {assert!(p2_ejected.insert(item))}
                else {assert!(p2_kept.insert(item))}
            }
            for item in self.p3_free.drain() {
                if ejected_nodes.contains(&item) {assert!(p3_ejected.insert(item))}
                else {assert!(p3_kept.insert(item))}
            }
            for item in self.pp_free.drain() {
                if ejected_nodes.contains(&item) {assert!(pp_ejected.insert(item))}
                else {assert!(pp_kept.insert(item))}
            }
            // update self
            self.xh_free = xh_kept;
            self.xt_free = xt_kept;
            self.xp_free = xp_kept;
            self.p1_free = p1_kept;
            self.p2_free = p2_kept;
            self.p3_free = p3_kept;
            self.pp_free = pp_kept;
            // return new
            OpenPorts {
                xh_free: xh_ejected,
                xt_free: xt_ejected,
                xp_free: xp_ejected,
                p1_free: p1_ejected,
                p2_free: p2_ejected,
                p3_free: p3_ejected,
                pp_free: pp_ejected
            }
        }

        /**
         * When introduced as a monomer, an Axn-type species contains only itself in the Axn-type
         * free-site sets. It contains nothing in the APC-type free-site sets, as there are no
         * APCs in this species.
        */
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

        /**
         * When introduced as a monomer, an APC-type species contains only itself in the APC-type
         * free-site sets. It contains nothing in the Axn-type free-site sets, as there are no
         * Axns in this species.
        */
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

        /**
         * When constructing the universal tracker, pre-allocate based on the mass sizes, leaving
         * the values empty.
        */
        pub fn default_empty(axn_mass: usize, apc_mass: usize) -> Self {
            OpenPorts {
                xh_free: HashSet::with_capacity(axn_mass),
                xt_free: HashSet::with_capacity(axn_mass),
                xp_free: HashSet::with_capacity(axn_mass),
                p1_free: HashSet::with_capacity(apc_mass),
                p2_free: HashSet::with_capacity(apc_mass),
                p3_free: HashSet::with_capacity(apc_mass),
                pp_free: HashSet::with_capacity(apc_mass)
            }
        }
    }
}

pub mod edge_types {
    use std::collections::BTreeSet;
    use std::rc::Rc;
    use std::cell::RefCell;
    use petgraph::prelude::*;
    use crate::edge_ends::EdgeEnds;

    #[derive(Debug)]
    /** 
     * Structure of BTreeSets, one per bond type, holding <Rc<RefCell<EdgeEnds>>>.
     */
    pub struct EdgeTypes {
        pub xh_xt: BTreeSet<Rc<RefCell<EdgeEnds>>>,
        pub p1_xp: BTreeSet<Rc<RefCell<EdgeEnds>>>,
        pub p2_xp: BTreeSet<Rc<RefCell<EdgeEnds>>>,
        pub p3_xp: BTreeSet<Rc<RefCell<EdgeEnds>>>,
        pub pp_pp: BTreeSet<Rc<RefCell<EdgeEnds>>>
    }
    
    impl EdgeTypes {
        /**
         * Update this instance with the elements from the other EdgeType structure.
        */
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

        /**
         * Partition the site-specific BTreeSets using a set of NodeIndexes, the "ejected" nodes.
         * A bond where both agents belong to the ejected set, gets ejected and returned
         * in a new EdgeType object. Bonds where both agents do not belong to agents in the
         * ejected set are retained. If one agent is ejected but another retained, this function
         * panics, as it should not be breaking cycles as a side-effect.
        */
        pub fn eject_where(&mut self, ejected_nodes: &BTreeSet<NodeIndex>) -> Self {
            // this would be simplified by "drain_filter" moving out of the Nightly-only
            let mut xh_xt_ejected: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut p1_xp_ejected: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut p2_xp_ejected: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut p3_xp_ejected: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut pp_pp_ejected: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut xh_xt_kept: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut p1_xp_kept: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut p2_xp_kept: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut p3_xp_kept: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            let mut pp_pp_kept: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
            for item in &self.xh_xt {
                if ejected_nodes.contains(&item.borrow().a) && ejected_nodes.contains(&item.borrow().b) {assert!(xh_xt_ejected.insert(Rc::clone(item)))}
                else if !ejected_nodes.contains(&item.borrow().a) && !ejected_nodes.contains(&item.borrow().b) {assert!(xh_xt_kept.insert(Rc::clone(item)))}
                else {panic!{"This binary break would sever more than one bond!"}}
            }
            for item in &self.p1_xp {
                if ejected_nodes.contains(&item.borrow().a) && ejected_nodes.contains(&item.borrow().b) {assert!(p1_xp_ejected.insert(Rc::clone(item)))}
                else if !ejected_nodes.contains(&item.borrow().a) && !ejected_nodes.contains(&item.borrow().b) {assert!(p1_xp_kept.insert(Rc::clone(item)))}
                else {panic!{"This binary break would sever more than one bond!"}}
            }
            for item in &self.p2_xp {
                if ejected_nodes.contains(&item.borrow().a) && ejected_nodes.contains(&item.borrow().b) {assert!(p2_xp_ejected.insert(Rc::clone(item)))}
                else if !ejected_nodes.contains(&item.borrow().a) && !ejected_nodes.contains(&item.borrow().b) {assert!(p2_xp_kept.insert(Rc::clone(item)))}
                else {panic!{"This binary break would sever more than one bond!"}}
            }
            for item in &self.p3_xp {
                if ejected_nodes.contains(&item.borrow().a) && ejected_nodes.contains(&item.borrow().b) {assert!(p3_xp_ejected.insert(Rc::clone(item)))}
                else if !ejected_nodes.contains(&item.borrow().a) && !ejected_nodes.contains(&item.borrow().b) {assert!(p3_xp_kept.insert(Rc::clone(item)))}
                else {panic!{"This binary break would sever more than one bond!"}}
            }
            for item in &self.pp_pp {
                if ejected_nodes.contains(&item.borrow().a) && ejected_nodes.contains(&item.borrow().b) {assert!(pp_pp_ejected.insert(Rc::clone(item)))}
                else if !ejected_nodes.contains(&item.borrow().a) && !ejected_nodes.contains(&item.borrow().b) {assert!(pp_pp_kept.insert(Rc::clone(item)))}
                else {panic!{"This binary break would sever more than one bond!"}}
            }
            // update self
            self.xh_xt = xh_xt_kept;
            self.p1_xp = p1_xp_kept;
            self.p2_xp = p2_xp_kept;
            self.p3_xp = p3_xp_kept;
            self.pp_pp = pp_pp_kept;
            // return new
            EdgeTypes {
                xh_xt: xh_xt_ejected,
                p1_xp: p1_xp_ejected,
                p2_xp: p2_xp_ejected,
                p3_xp: p3_xp_ejected,
                pp_pp: pp_pp_ejected
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
    use std::fmt;

    #[derive(Debug)]
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

    impl fmt::Display for RuleActivities {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "Axn-Axn\n\tu bind: {}\n\tb bind: {}\n\tu free: {}\n\tb free: {}\nAp1_Axn\n\tu bind: {}\n\tb bind: {}\n\tu free: {}\n\tb free: {}\nAp2_Axn\n\tu bind: {}\n\tb bind: {}\n\tu free: {}\n\tb free: {}\nAp3_Axn\n\tu bind: {}\n\tb bind: {}\n\tu free: {}\n\tb free: {}\nApc_Apc\n\tu bind: {}\n\tb bind: {}\n\tu free: {}\n\tb free: {}\n", self.axn_axn_u_bind, self.axn_axn_b_bind, self.axn_axn_u_free, self.axn_axn_b_free, self.ap1_axn_u_bind, self.ap1_axn_b_bind, self.ap1_axn_u_free, self.ap1_axn_b_free, self.ap2_axn_u_bind, self.ap2_axn_b_bind, self.ap2_axn_u_free, self.ap2_axn_b_free, self.ap3_axn_u_bind, self.ap3_axn_b_bind, self.ap3_axn_u_free, self.ap3_axn_b_free, self.apc_apc_u_bind, self.apc_apc_b_bind, self.apc_apc_u_free, self.apc_apc_b_free)
        }
    }

    #[derive(Debug)]
    pub struct MassActionTerm {
        pub mass: usize,
        pub rate: f64
    }

    impl fmt::Display for MassActionTerm {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "|{}| * {}", self.mass, self.rate)
        }
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
    use nalgebra::*;

    #[derive(Debug)]
    /** 
     * Structure of boolean matrices, one per bond type. The value flags an interaction as allowed
     * via unary interactions (i.e. the agents belong to the same connected component).
    */
    pub struct UnaryEmbeds {
        pub xh_xt: DMatrix<bool>,
        pub p1_xp: DMatrix<bool>,
        pub p2_xp: DMatrix<bool>,
        pub p3_xp: DMatrix<bool>,
        pub pp_pp: DMatrix<bool>
    }

    impl UnaryEmbeds {
        pub fn new_from_masses(axn_mass: usize, apc_mass: usize) -> Self {
            UnaryEmbeds {
                xh_xt: DMatrix::from_element(axn_mass, axn_mass, false),
                p1_xp: DMatrix::from_element(apc_mass, axn_mass, false),
                p2_xp: DMatrix::from_element(apc_mass, axn_mass, false),
                p3_xp: DMatrix::from_element(apc_mass, axn_mass, false),
                pp_pp: DMatrix::from_element(apc_mass, apc_mass, false),
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
    use petgraph::{graph::Graph, algo::astar::astar, prelude::*, visit::Dfs};
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
        edge_index_map: BTreeMap<EdgeIndex, Rc<RefCell<EdgeEnds>>>,         // used to update Z when removing bonds from the universe graph
        unary_binding_pairs: UnaryEmbeds,
        rule_activities: RuleActivities
    }
    
    impl Mixture {
        pub fn to_kappa(&self) -> String {
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
            let mut net: Graph<AgentType, bool, Undirected> = Graph::with_capacity(x_mass + p_mass, x_mass * 3 + p_mass * 4);   // the universe graph
            let mut spc: BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies>>> = BTreeMap::new();                                    // the global agent -> species map
            let mut sps: VecDeque<Rc<RefCell<MixtureSpecies>>> = VecDeque::with_capacity(x_mass + p_mass);                      // the species set
            let edg = EdgeTypes::default_empty();                                                                               // the edge-type tracker
            let uem = UnaryEmbeds::new_from_masses(x_mass, p_mass);                                                             // unary embedding cache
            let eim = BTreeMap::new();                                                                                          // the global edge-index -> EdgeEnds map
            let mut opr = OpenPorts::default_empty(x_mass, p_mass);                                                             // the open ports tracker
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
            let rac = RuleActivities {
                axn_axn_u_bind: MassActionTerm{mass: 0, rate: rule_rates.axn_axn_u_bind}, // all monomeric: there are no unary binding opportunities
                ap1_axn_u_bind: MassActionTerm{mass: 0, rate: rule_rates.ap1_axn_u_bind},
                ap2_axn_u_bind: MassActionTerm{mass: 0, rate: rule_rates.ap2_axn_u_bind},
                ap3_axn_u_bind: MassActionTerm{mass: 0, rate: rule_rates.ap3_axn_u_bind},
                apc_apc_u_bind: MassActionTerm{mass: 0, rate: rule_rates.apc_apc_u_bind},
                axn_axn_u_free: MassActionTerm{mass: 0, rate: rule_rates.axn_axn_u_free}, // all monomeric: no cycles to open
                ap1_axn_u_free: MassActionTerm{mass: 0, rate: rule_rates.ap1_axn_u_free},
                ap2_axn_u_free: MassActionTerm{mass: 0, rate: rule_rates.ap2_axn_u_free},
                ap3_axn_u_free: MassActionTerm{mass: 0, rate: rule_rates.ap3_axn_u_free},
                apc_apc_u_free: MassActionTerm{mass: 0, rate: rule_rates.apc_apc_u_free},
                axn_axn_b_free: MassActionTerm{mass: 0, rate: rule_rates.axn_axn_b_free}, // all monomeric: no complexes to break apart
                ap1_axn_b_free: MassActionTerm{mass: 0, rate: rule_rates.ap1_axn_b_free},
                ap2_axn_b_free: MassActionTerm{mass: 0, rate: rule_rates.ap2_axn_b_free},
                ap3_axn_b_free: MassActionTerm{mass: 0, rate: rule_rates.ap3_axn_b_free},
                apc_apc_b_free: MassActionTerm{mass: 0, rate: rule_rates.apc_apc_b_free},
                axn_axn_b_bind: MassActionTerm{mass: x_mass * (x_mass-1), rate: rule_rates.axn_axn_b_bind}, // all monomeric: everything that can happen is of binary binding type
                ap1_axn_b_bind: MassActionTerm{mass: p_mass * x_mass, rate: rule_rates.ap1_axn_b_bind},
                ap2_axn_b_bind: MassActionTerm{mass: p_mass * x_mass, rate: rule_rates.ap2_axn_b_bind},
                ap3_axn_b_bind: MassActionTerm{mass: p_mass * x_mass, rate: rule_rates.ap3_axn_b_bind},
                apc_apc_b_bind: MassActionTerm{mass: p_mass * (p_mass-1), rate: rule_rates.apc_apc_b_bind}
            };
            // bringing it all together
            Mixture {universe_graph: net, species_annots: spc, ports: opr, edges: edg, rule_activities: rac, species_set: sps, unary_binding_pairs: uem, edge_index_map: eim}
        }
    
        pub fn axn_axn_unary_bind(&mut self, head_node: NodeIndex, tail_node: NodeIndex) {
            // sanity check for unary-ness
            assert_eq!(self.species_annots.get(&head_node), self.species_annots.get(&tail_node), "This head and tail nodes do not already belong to the same species; can't unary bind them!");
            let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_node).unwrap().borrow_mut();
            // create edge & update caches
            let new_edge_index = self.universe_graph.add_edge(head_node, tail_node, true);
            let new_edge = Rc::new(RefCell::new(EdgeEnds{a: head_node, b: tail_node, a_s: 'h', b_s: 't', z: new_edge_index}));
            let node_a: &mut AgentType = self.universe_graph.node_weight_mut(head_node).unwrap();
            match node_a {
                AgentType::AxnNode(x1, _, _) => match x1 {
                    Some(_) => panic!("This Axn head was already bound!"),
                    None => *x1 = Some(new_edge_index)
                }
                AgentType::ApcNode(_, _, _, _) => panic!("This matched an APC node instead of a free-headed Axin node!")
            }
            let node_b = self.universe_graph.node_weight_mut(tail_node).unwrap();
            match node_b {
                AgentType::AxnNode(_, x2, _) => match x2 {
                    Some(_) => panic!("This Axn tail was already bound!"),
                    None => *x2 = Some(new_edge_index)
                }
                AgentType::ApcNode(_, _, _, _) => panic!("This matched an APC node instead of a free-tailed Axin node!")
            }
            target_species.edges.xh_xt.insert(Rc::clone(&new_edge));
            self.edge_index_map.insert(new_edge_index, Rc::clone(&new_edge));
            self.edges.xh_xt.insert(Rc::clone(&new_edge));
            self.ports.xh_free.remove(&head_node);
            self.ports.xt_free.remove(&tail_node);
            target_species.ports.xh_free.remove(&head_node);
            target_species.ports.xt_free.remove(&tail_node);
            // helper closure for bond resilience updates
            let mut resilience_iterate_typed = |edge_iter: &BTreeSet<std::rc::Rc<std::cell::RefCell<EdgeEnds>>>| -> usize {
                let mut bonds_cyclized: usize = 0;
                for boxed_edge in edge_iter {
                    let edge = boxed_edge.borrow();
                    let edge_id = self.universe_graph.find_edge(edge.a, edge.b).unwrap();
                    // check bond's weight, which is a boolean, marking if the bond is a known cycle-member
                    if ! self.universe_graph.edge_weight(edge_id).unwrap() {
                        let mut counterfactual_graph = self.universe_graph.clone();
                        counterfactual_graph.remove_edge(edge_id);
                        let path = astar(&counterfactual_graph, edge.a, |finish| finish == edge.b, |_| 1, |_| 0);
                        match path {
                            Some(_) => {self.universe_graph.update_edge(edge.a, edge.b, true); bonds_cyclized += 1;},
                            None => ()
                        }
                    }
                }
                bonds_cyclized
            };
            // update resilience of that species' edges
            let xh_xt_cyclized: usize = resilience_iterate_typed(&target_species.edges.xh_xt);
            let p1_xp_cyclized: usize = resilience_iterate_typed(&target_species.edges.p1_xp);
            let p2_xp_cyclized: usize = resilience_iterate_typed(&target_species.edges.p2_xp);
            let p3_xp_cyclized: usize = resilience_iterate_typed(&target_species.edges.p3_xp);
            let pp_pp_cyclized: usize = resilience_iterate_typed(&target_species.edges.pp_pp);
            // update rule activities
            let binary_combinatorics_lost_heads: usize = self.ports.xt_free.len() - target_species.ports.xt_free.len();
            let binary_combinatorics_lost_tails: usize = self.ports.xh_free.len() - target_species.ports.xh_free.len();
            let binary_combinatorics_lost: usize = binary_combinatorics_lost_tails + binary_combinatorics_lost_heads;
            let unary_embeds_made_here: usize = 1;
            let unary_combinatorics_lost_head: usize = target_species.ports.xt_free.len();
            let unary_combinatorics_lost_tail: usize = target_species.ports.xh_free.len();
            let unary_combinatorics_lost: usize = unary_combinatorics_lost_head + unary_combinatorics_lost_tail;
            self.rule_activities.axn_axn_b_free.mass -= xh_xt_cyclized;
            self.rule_activities.axn_axn_u_free.mass += unary_embeds_made_here + xh_xt_cyclized;     //the new bond
            self.rule_activities.axn_axn_b_bind.mass -= binary_combinatorics_lost;
            self.rule_activities.axn_axn_u_bind.mass -= unary_embeds_made_here + unary_combinatorics_lost;
            //  the rest are binary -> unary conversions
            self.rule_activities.ap1_axn_u_free.mass += p1_xp_cyclized;
            self.rule_activities.ap1_axn_b_free.mass -= p1_xp_cyclized;
            self.rule_activities.ap2_axn_u_free.mass += p2_xp_cyclized;
            self.rule_activities.ap2_axn_b_free.mass -= p2_xp_cyclized;
            self.rule_activities.ap3_axn_u_free.mass += p3_xp_cyclized;
            self.rule_activities.ap3_axn_b_free.mass -= p3_xp_cyclized;
            self.rule_activities.apc_apc_u_free.mass += pp_pp_cyclized;
            self.rule_activities.apc_apc_b_free.mass -= pp_pp_cyclized;
            // update unaries
            // this could be replaced with a for loop, iterating over the respective heads & tails for
            //  this specific species, which may be faster than the nalgebra's fill_row | fill_column
            self.unary_binding_pairs.xh_xt.fill_row(head_node.index(), false);
            self.unary_binding_pairs.xh_xt.fill_column(tail_node.index(), false);
        }
    
        pub fn axn_axn_binary_bind(&mut self, head_node: NodeIndex, tail_node: NodeIndex) {
            // sanity check for binary-ness
            assert_ne!(self.species_annots.get(&head_node), self.species_annots.get(&tail_node), "This head and tail nodes already belong to the same species; can't binary bind them!");
            // define which will be the "host" species, receiving the cache data from "eaten" species
            let (host_index, eaten_index) = if self.species_annots.get(&head_node).unwrap().borrow().size >= self.species_annots.get(&tail_node).unwrap().borrow().size {(head_node, tail_node)} else {(tail_node, head_node)};
            // merge caches held by the species; uses interior-mutability pattern
            assert!(self.species_annots.get(&head_node).unwrap().borrow_mut().ports.xh_free.remove(&head_node), "This head node was not listed as free in-species already, can't bind to it!");
            assert!(self.species_annots.get(&tail_node).unwrap().borrow_mut().ports.xt_free.remove(&tail_node), "This tail node was not listed as free in-species already, can't bind to it!");
            assert!(self.ports.xh_free.remove(&head_node), "This head node was not globally listed as free already, can't bind to it!");
            assert!(self.ports.xt_free.remove(&tail_node), "This tail node was not globally listed as free already, can't bind to it!");
            // update rule activities
            //  for the Axn-Axn case, we can ignore the _free'ing activities
            //  for the other bond types, as we are not modifying those bond
            //  counts, ergo their unbinding counts
            let binary_combinatorics_lost_heads: usize = self.ports.xt_free.len() - self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free.len();
            let binary_combinatorics_lost_tails: usize = self.ports.xh_free.len() - self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free.len();
            let binary_combinatorics_lost: usize = binary_combinatorics_lost_heads + binary_combinatorics_lost_tails;
            let transformed_to_unary: usize = (self.species_annots.get(&host_index).unwrap().borrow().ports.xh_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xt_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.xh_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xt_free.len());
            let binary_embeds_made_here: usize = 1;     // the new bond
            // since self-bonds are not allowed in this model,
            //  this special-cases the Axn-Axn unary treatment,
            //  discounting the monomer self-binding
            let unary_combinatorics_lost_head: usize = 
                if self.species_annots.get(&head_node).unwrap().borrow().size > 1 {
                    self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free.len()
                }
                else {0};
            let unary_combinatorics_lost_tail: usize =
                if self.species_annots.get(&tail_node).unwrap().borrow().size > 1 {
                    self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free.len()
                }
                else {0};
            let unary_combinatorics_lost = unary_combinatorics_lost_head + unary_combinatorics_lost_tail;
            self.rule_activities.axn_axn_b_free.mass += binary_embeds_made_here;                    // the new bond
            //self.rule_activities.axn_axn_u_free.mass      // does not apply
            self.rule_activities.axn_axn_b_bind.mass -= transformed_to_unary + binary_combinatorics_lost + binary_embeds_made_here;
            self.rule_activities.axn_axn_u_bind.mass += transformed_to_unary;                       // broken into two operations to avoid underflowing
            self.rule_activities.axn_axn_u_bind.mass -= unary_combinatorics_lost;                   // when the former is 0, but the latter is non-zero
            //  the rest are unary -> binary conversions
            self.rule_activities.ap1_axn_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap1_axn_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap2_axn_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap2_axn_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap3_axn_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap3_axn_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.apc_apc_b_bind.mass -= (self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free.len());
            self.rule_activities.apc_apc_u_bind.mass += (self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free.len()) + (self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free.len());
            self.species_annots.get(&host_index).unwrap().borrow_mut().ports.update_from(&self.species_annots.get(&eaten_index).unwrap().borrow().ports);
            self.species_annots.get(&host_index).unwrap().borrow_mut().edges.update_from(&self.species_annots.get(&eaten_index).unwrap().borrow().edges);
            self.species_annots.get(&host_index).unwrap().borrow_mut().size += self.species_annots.get(&eaten_index).unwrap().borrow().size;
            self.species_annots.get(&host_index).unwrap().borrow_mut().agent_set.append(&mut self.species_annots.get(&eaten_index).unwrap().borrow().agent_set.clone());
            // create & update edge trackers
            let new_edge_index = self.universe_graph.add_edge(head_node, tail_node, false);
            let new_edge = Rc::new(RefCell::new(EdgeEnds{a: head_node, b: tail_node, a_s: 'h', b_s: 't', z: new_edge_index}));
            self.edge_index_map.insert(new_edge_index, Rc::clone(&new_edge));
            self.species_annots.get(&host_index).unwrap().borrow_mut().edges.xh_xt.insert(Rc::clone(&new_edge));
            self.edges.xh_xt.insert(Rc::clone(&new_edge));
            // update species annotation tracker
            let spec_ix: usize = self.species_set.binary_search(self.species_annots.get(&eaten_index).unwrap()).unwrap();
            self.species_set.remove(spec_ix);
            let indexes_of_eaten: BTreeSet<NodeIndex> = self.species_annots.get(&eaten_index).unwrap().borrow().agent_set.clone();
            for agent_index in indexes_of_eaten {
                let new_ref = Rc::clone(&self.species_annots.get(&host_index).unwrap());
                self.species_annots.entry(agent_index).and_modify(|e| {*e = new_ref});
            }
            // update bond identifier trackers
            let node_a_weight = self.universe_graph.node_weight_mut(head_node).unwrap();
            match node_a_weight {
                AgentType::AxnNode(x1, _, _) => match x1 {
                    Some(_) => panic!("This Axn head was already bound!"),
                    None => *x1 = Some(new_edge_index),
                }
                AgentType::ApcNode(_, _, _, _) => panic!("Node weight for an Axn head matched and APC!")
            };
            let node_b_weight = self.universe_graph.node_weight_mut(tail_node).unwrap();
            match node_b_weight {
                AgentType::AxnNode(_, x2, _) => match x2 {
                    Some(_) => panic!("This Axn tail was already bound!"),
                    None => *x2 = Some(new_edge_index),
                }
                AgentType::ApcNode(_, _, _, _) => panic!("Node weight for an Axn tail matched an APC!")
            };
            // update unaries
            for xh_port in &self.species_annots.get(&host_index).unwrap().borrow().ports.xh_free{
                for xt_port in &self.species_annots.get(&eaten_index).unwrap().borrow().ports.xt_free {
                    self.unary_binding_pairs.xh_xt[(xh_port.index(), xt_port.index())] = true;
                }
            }
            for xt_port in &self.species_annots.get(&host_index).unwrap().borrow().ports.xt_free{
                for xh_port in &self.species_annots.get(&eaten_index).unwrap().borrow().ports.xh_free {
                    self.unary_binding_pairs.xh_xt[(xh_port.index(), xt_port.index())] = true;
                }
            }
            for xp_port in &self.species_annots.get(&eaten_index).unwrap().borrow().ports.xp_free {
                for p1_port in &self.species_annots.get(&host_index).unwrap().borrow().ports.p1_free{
                    self.unary_binding_pairs.p1_xp[(p1_port.index(), xp_port.index())] = true;
                }
                for p2_port in &self.species_annots.get(&host_index).unwrap().borrow().ports.p2_free{
                    self.unary_binding_pairs.p2_xp[(p2_port.index(), xp_port.index())] = true;
                }
                for p3_port in &self.species_annots.get(&host_index).unwrap().borrow().ports.p3_free{
                    self.unary_binding_pairs.p3_xp[(p3_port.index(), xp_port.index())] = true;
                }
            }
            for xp_port in &self.species_annots.get(&host_index).unwrap().borrow().ports.xp_free {
                for p1_port in &self.species_annots.get(&eaten_index).unwrap().borrow().ports.p1_free{
                    self.unary_binding_pairs.p1_xp[(p1_port.index(), xp_port.index())] = true;
                }
                for p2_port in &self.species_annots.get(&eaten_index).unwrap().borrow().ports.p2_free{
                    self.unary_binding_pairs.p2_xp[(p2_port.index(), xp_port.index())] = true;
                }
                for p3_port in &self.species_annots.get(&eaten_index).unwrap().borrow().ports.p3_free{
                    self.unary_binding_pairs.p3_xp[(p3_port.index(), xp_port.index())] = true;
                }
            }
            for p1_port in &self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free{
                for p2_port in &self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free {
                    self.unary_binding_pairs.xh_xt[(p1_port.index(), p2_port.index())] = true;
                }
            }
            for p1_port in &self.species_annots.get(&eaten_index).unwrap().borrow().ports.pp_free{
                for p2_port in &self.species_annots.get(&host_index).unwrap().borrow().ports.pp_free {
                    self.unary_binding_pairs.xh_xt[(p1_port.index(), p2_port.index())] = true;
                }
            }
        }

        pub fn axn_axn_unary_unbind(&mut self, target_edge: Rc<RefCell<EdgeEnds>>) {
            let EdgeEnds{a: head_node, b: tail_node, a_s: head_type, b_s: tail_type, z: edge_index} = *target_edge.borrow();
            // sanity check for reslient bond
            assert!(self.universe_graph.edge_weight(edge_index).unwrap(), "This bond is not flagged as resilient! Breaking it would not yield a still-connected component.");
            match head_type {
                'h' => (),
                _ => panic!("Did not find an Axn head-type annotation!")
            }
            match tail_type {
                't' => (),
                _ => panic!("Did not find an Axn tail-type annotation!")
            }
            let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_node).unwrap().borrow_mut();
            let old_last_edge_index = self.universe_graph.edge_indices().last().unwrap();
            self.universe_graph.remove_edge(edge_index).unwrap();
            let swapped_edge = self.edge_index_map.remove(&old_last_edge_index).unwrap();
            if old_last_edge_index != edge_index {
                // remove edge & update caches; keep edge index tracker up to date
                swapped_edge.borrow_mut().z = edge_index;
                let swapped_colateral_a = self.universe_graph.node_weight_mut(swapped_edge.borrow().a).unwrap();
                match swapped_colateral_a {
                    AgentType::AxnNode(x1, x2, x3) => match (&x1, &x2, &x3) {
                        (Some(this_index), _, _) if this_index == &old_last_edge_index => *x1 = Some(edge_index),
                        (_, Some(this_index), _) if this_index == &old_last_edge_index => *x2 = Some(edge_index),
                        (_, _, Some(this_index)) if this_index == &old_last_edge_index => *x3 = Some(edge_index),
                        _ => panic!("This node, an Axn, did not have the expected old edge index for updating!")
                    }
                    AgentType::ApcNode(p1, p2, p3, p4) => match (&p1, &p2, &p3, &p4) {
                        (Some(this_index), _, _, _) if this_index == &old_last_edge_index => *p1 = Some(edge_index),
                        (_, Some(this_index), _, _) if this_index == &old_last_edge_index => *p2 = Some(edge_index),
                        (_, _, Some(this_index), _) if this_index == &old_last_edge_index => *p3 = Some(edge_index),
                        (_, _, _, Some(this_index)) if this_index == &old_last_edge_index => *p4 = Some(edge_index),
                        _ => panic!("This node, an APC, did not have the expected old edge index for updating!")
                    }
                }
                let swapped_colateral_b = self.universe_graph.node_weight_mut(swapped_edge.borrow().b).unwrap();
                match swapped_colateral_b {
                    AgentType::AxnNode(x1, x2, x3) => match (&x1, &x2, &x3) {
                        (Some(this_index), _, _) if this_index == &old_last_edge_index => *x1 = Some(edge_index),
                        (_, Some(this_index), _) if this_index == &old_last_edge_index => *x2 = Some(edge_index),
                        (_, _, Some(this_index)) if this_index == &old_last_edge_index => *x3 = Some(edge_index),
                        _ => panic!("This node, an Axn, did not have the expected old edge index for updating!")
                    }
                    AgentType::ApcNode(p1, p2, p3, p4) => match (&p1, &p2, &p3, &p4) {
                        (Some(this_index), _, _, _) if this_index == &old_last_edge_index => *p1 = Some(edge_index),
                        (_, Some(this_index), _, _) if this_index == &old_last_edge_index => *p2 = Some(edge_index),
                        (_, _, Some(this_index), _) if this_index == &old_last_edge_index => *p3 = Some(edge_index),
                        (_, _, _, Some(this_index)) if this_index == &old_last_edge_index => *p4 = Some(edge_index),
                        _ => panic!("This node, an APC, did not have the expected old edge index for updating!")
                    }
                }
                *self.edge_index_map.get_mut(&edge_index).unwrap() = Rc::clone(&swapped_edge);
            }
            // sanity checks for bond breaking; update bond indexes on the agent-type
            let node_a: &mut AgentType = self.universe_graph.node_weight_mut(head_node).unwrap();
            match node_a {
                AgentType::AxnNode(x1, _, _) => match x1 {
                    Some(_) => *x1 = None,
                    None => panic!("This Axn head was not already bound!")
                }
                AgentType::ApcNode(_, _, _, _) => panic!("This matched an APC node instead of a bond-headed Axin node!")
            }
            let node_b = self.universe_graph.node_weight_mut(tail_node).unwrap();
            match node_b {
                AgentType::AxnNode(_, x2, _) => match x2 {
                    Some(_) => *x2 = None,
                    None => panic!("This Axn tail was not already bound!")
                }
                AgentType::ApcNode(_, _, _, _) => panic!("This matched an APC node instead of a bond-tailed Axin node!")
            }
            target_species.edges.xh_xt.remove(&target_edge);
            self.edges.xh_xt.remove(&target_edge);
            assert!(self.ports.xh_free.insert(head_node), "This Axn head was already listed as free prior to this unbinding!");
            assert!(self.ports.xt_free.insert(tail_node), "This Axn tail was already listed as free prior to this unbinding!");
            assert!(target_species.ports.xh_free.insert(head_node), "This Axn head was already listed as free prior to this unbinding!");
            assert!(target_species.ports.xt_free.insert(tail_node), "This Axn tail was already listed as free prior to this unbinding!");
            // helper closure for bond resilience updates
            let mut resilience_iterate_typed = |edge_iter: &BTreeSet<std::rc::Rc<std::cell::RefCell<EdgeEnds>>>| -> usize {
                let mut bonds_decyclized: usize = 0;
                for boxed_edge in edge_iter {
                    let edge = boxed_edge.borrow();
                    // check bond's weight, which is a boolean, marking if the bond is a known cycle-member
                    if *self.universe_graph.edge_weight(edge.z).unwrap() {
                        let mut counterfactual_graph = self.universe_graph.clone();
                        counterfactual_graph.remove_edge(edge.z);
                        let path = astar(&counterfactual_graph, edge.a, |finish| finish == edge.b, |_| 1, |_| 0);
                        match path {
                            Some(_) => (),
                            None => {self.universe_graph.update_edge(edge.a, edge.b, false); bonds_decyclized += 1;}
                        }
                    }
                }
                bonds_decyclized
            };
            // update resilience of that species' edges
            let xh_xt_decyclized: usize = resilience_iterate_typed(&target_species.edges.xh_xt);
            let p1_xp_decyclized: usize = resilience_iterate_typed(&target_species.edges.p1_xp);
            let p2_xp_decyclized: usize = resilience_iterate_typed(&target_species.edges.p2_xp);
            let p3_xp_decyclized: usize = resilience_iterate_typed(&target_species.edges.p3_xp);
            let pp_pp_decyclized: usize = resilience_iterate_typed(&target_species.edges.pp_pp);
            // update rule activities
            let binary_combinatorics_gained_heads: usize = self.ports.xt_free.len() - target_species.ports.xt_free.len();
            let binary_combinatorics_gained_tails: usize = self.ports.xh_free.len() - target_species.ports.xh_free.len();
            let binary_combinatorics_gained: usize = binary_combinatorics_gained_tails + binary_combinatorics_gained_heads;
            let unary_embeds_made_here: usize = 1;
            let unary_combinatorics_gained_head: usize = target_species.ports.xt_free.len();
            let unary_combinatorics_gained_tail: usize = target_species.ports.xh_free.len();
            let unary_combinatorics_gained: usize = unary_combinatorics_gained_head * unary_combinatorics_gained_tail;
            self.rule_activities.axn_axn_b_free.mass += xh_xt_decyclized;
            self.rule_activities.axn_axn_u_free.mass -= unary_embeds_made_here + xh_xt_decyclized;
            self.rule_activities.axn_axn_b_bind.mass += binary_combinatorics_gained;
            self.rule_activities.axn_axn_u_bind.mass += unary_combinatorics_gained;
            //  the rest are binary -> unary conversions
            self.rule_activities.ap1_axn_u_free.mass -= p1_xp_decyclized;
            self.rule_activities.ap1_axn_b_free.mass += p1_xp_decyclized;
            self.rule_activities.ap2_axn_u_free.mass -= p2_xp_decyclized;
            self.rule_activities.ap2_axn_b_free.mass += p2_xp_decyclized;
            self.rule_activities.ap3_axn_u_free.mass -= p3_xp_decyclized;
            self.rule_activities.ap3_axn_b_free.mass += p3_xp_decyclized;
            self.rule_activities.apc_apc_u_free.mass -= pp_pp_decyclized;
            self.rule_activities.apc_apc_b_free.mass += pp_pp_decyclized;
            // update unary binding pair: one bond removed, combinatorics now possible
            for xh_port in &target_species.ports.xh_free {
                for xt_port in &target_species.ports.xt_free {
                    self.unary_binding_pairs.xh_xt[(xh_port.index(), xt_port.index())] = true;
                }
            }
        }

        pub fn axn_axn_binary_unbind(&mut self, target_edge: Rc<RefCell<EdgeEnds>>) {
            let EdgeEnds{a: head_node, b: tail_node, a_s: head_type, b_s: tail_type, z: edge_index} = *target_edge.borrow();
            // sanity check for non-reslient bond
            assert!(! self.universe_graph.edge_weight(edge_index).unwrap(), "This bond is flagged as resilient! Breaking it would yield a still-connected component.");
            match head_type {
                'h' => (),
                _ => panic!("Did not find an Axn head-type annotation!")
            }
            match tail_type {
                't' => (),
                _ => panic!("Did not find an Axn tail-type annotation!")
            }
            let old_last_edge_index = self.universe_graph.edge_indices().last().unwrap();
            self.universe_graph.remove_edge(edge_index).unwrap();
            let swapped_edge = self.edge_index_map.remove(&old_last_edge_index).unwrap();
            if old_last_edge_index != edge_index {
                // remove edge & update caches; keep edge index tracker up to date
                swapped_edge.borrow_mut().z = edge_index;
                let swapped_colateral_a = self.universe_graph.node_weight_mut(swapped_edge.borrow().a).unwrap();
                match swapped_colateral_a {
                    AgentType::AxnNode(x1, x2, x3) => match (&x1, &x2, &x3) {
                        (Some(this_index), _, _) if this_index == &old_last_edge_index => *x1 = Some(edge_index),
                        (_, Some(this_index), _) if this_index == &old_last_edge_index => *x2 = Some(edge_index),
                        (_, _, Some(this_index)) if this_index == &old_last_edge_index => *x3 = Some(edge_index),
                        _ => panic!("This node, an Axn, did not have the expected old edge index for updating!")
                    }
                    AgentType::ApcNode(p1, p2, p3, p4) => match (&p1, &p2, &p3, &p4) {
                        (Some(this_index), _, _, _) if this_index == &old_last_edge_index => *p1 = Some(edge_index),
                        (_, Some(this_index), _, _) if this_index == &old_last_edge_index => *p2 = Some(edge_index),
                        (_, _, Some(this_index), _) if this_index == &old_last_edge_index => *p3 = Some(edge_index),
                        (_, _, _, Some(this_index)) if this_index == &old_last_edge_index => *p4 = Some(edge_index),
                        _ => panic!("This node, an APC, did not have the expected old edge index for updating!")
                    }
                }
                let swapped_colateral_b = self.universe_graph.node_weight_mut(swapped_edge.borrow().b).unwrap();
                match swapped_colateral_b {
                    AgentType::AxnNode(x1, x2, x3) => match (&x1, &x2, &x3) {
                        (Some(this_index), _, _) if this_index == &old_last_edge_index => *x1 = Some(edge_index),
                        (_, Some(this_index), _) if this_index == &old_last_edge_index => *x2 = Some(edge_index),
                        (_, _, Some(this_index)) if this_index == &old_last_edge_index => *x3 = Some(edge_index),
                        _ => panic!("This node, an Axn, did not have the expected old edge index for updating!")
                    }
                    AgentType::ApcNode(p1, p2, p3, p4) => match (&p1, &p2, &p3, &p4) {
                        (Some(this_index), _, _, _) if this_index == &old_last_edge_index => *p1 = Some(edge_index),
                        (_, Some(this_index), _, _) if this_index == &old_last_edge_index => *p2 = Some(edge_index),
                        (_, _, Some(this_index), _) if this_index == &old_last_edge_index => *p3 = Some(edge_index),
                        (_, _, _, Some(this_index)) if this_index == &old_last_edge_index => *p4 = Some(edge_index),
                        _ => panic!("This node, an APC, did not have the expected old edge index for updating!")
                    }
                }
                *self.edge_index_map.get_mut(&edge_index).unwrap() = swapped_edge;
            }
            self.edges.xh_xt.remove(&target_edge);
            self.species_annots.get(&head_node).unwrap().borrow_mut().edges.xh_xt.remove(&target_edge);
            // discover what node indexes ended where; the loop terminates once the smallest graph is mapped-out
            let mut dfs_head = Dfs::new(&self.universe_graph, head_node);
            let mut dfs_tail = Dfs::new(&self.universe_graph, tail_node);
            let mut head_graph_indexes: BTreeSet<NodeIndex> = BTreeSet::new();
            let mut tail_graph_indexes: BTreeSet<NodeIndex> = BTreeSet::new();
            let mut nx_head = dfs_head.next(&self.universe_graph);
            let mut nx_tail = dfs_tail.next(&self.universe_graph);
            loop {
                match (nx_head, nx_tail) {
                    (Some(in_head_found), Some(in_tail_found)) => {
                        head_graph_indexes.insert(in_head_found);
                        tail_graph_indexes.insert(in_tail_found);
                        nx_head = dfs_head.next(&self.universe_graph);
                        nx_tail = dfs_tail.next(&self.universe_graph);
                    },
                    (None, Some(_)) => {
                        // we've mapped-out the head node's graph
                        tail_graph_indexes = self.species_annots.get(&tail_node).unwrap().borrow().agent_set.difference(&head_graph_indexes).cloned().collect();
                        break
                    },
                    (Some(_), None) => {
                        // we've mapped-out the tail node's graph
                        head_graph_indexes = self.species_annots.get(&head_node).unwrap().borrow().agent_set.difference(&tail_graph_indexes).cloned().collect();
                        break
                    },
                    (None, None) => break
                }
            }
            let (ejected_mark, ejected_indexes, retained_mark, retained_indexes) =  // ToDo should this be reversed? Iterate over small & lookup over big instead?
                if head_graph_indexes.len() <= tail_graph_indexes.len()
                    {(head_node, head_graph_indexes, tail_node, tail_graph_indexes)}
                else 
                    {(tail_node, tail_graph_indexes, head_node, head_graph_indexes)};
            // create new species cache, in-place modify old species cache
            let ejected_ports: OpenPorts = self.species_annots.get(&retained_mark).unwrap().borrow_mut().ports.eject_where(&ejected_indexes);
            let ejected_edges: EdgeTypes = self.species_annots.get(&retained_mark).unwrap().borrow_mut().edges.eject_where(&ejected_indexes);
            self.species_annots.get(&retained_mark).unwrap().borrow_mut().size -= &ejected_indexes.len();
            self.species_annots.get(&retained_mark).unwrap().borrow_mut().agent_set = retained_indexes;
            let new_species = Rc::new(RefCell::new(MixtureSpecies{
                ports: ejected_ports,
                edges: ejected_edges,
                size: ejected_indexes.len(),
                agent_set: ejected_indexes.clone()
            }));
            // update the species annotation tracker
            for node_index in &ejected_indexes {
                self.species_annots.entry(*node_index).and_modify(|e| {*e = Rc::clone(&new_species)});
            }
            // sanity checks for bond breaking; update bond indexes on the agent-type
            let node_a: &mut AgentType = self.universe_graph.node_weight_mut(head_node).unwrap();
            match node_a {
                AgentType::AxnNode(x1, _, _) => match x1 {
                    Some(_) => *x1 = None,
                    None => panic!("This Axn head was not already bound!")
                }
                AgentType::ApcNode(_, _, _, _) => panic!("This matched an APC node instead of a bond-headed Axin node!")
            }
            let node_b = self.universe_graph.node_weight_mut(tail_node).unwrap();
            match node_b {
                AgentType::AxnNode(_, x2, _) => match x2 {
                    Some(_) => *x2 = None,
                    None => panic!("This Axn tail was not already bound!")
                }
                AgentType::ApcNode(_, _, _, _) => panic!("This matched an APC node instead of a bond-tailed Axin node!")
            }
            // update the species set, open port trackers, global edge tracker
            self.species_set.push_back(Rc::clone(&new_species));
            assert!(self.ports.xh_free.insert(head_node), "This Axn head was already listed as free prior to this unbinding!");
            assert!(self.ports.xt_free.insert(tail_node), "This Axn tail was already listed as free prior to this unbinding!");
            assert!(self.species_annots.get(&head_node).unwrap().borrow_mut().ports.xh_free.insert(head_node), "This Axn head was already listed as free prior to this unbinding!");
            assert!(self.species_annots.get(&tail_node).unwrap().borrow_mut().ports.xt_free.insert(tail_node), "This Axn tail was already listed as free prior to this unbinding!");
            // update unary trackers
            for xh_port in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.xh_free {
                for xt_port in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xt_free {
                    self.unary_binding_pairs.xh_xt[(xh_port.index(), xt_port.index())] = false;
                }
            }
            for xt_port in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.xt_free {
                for xh_port in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xh_free {
                    self.unary_binding_pairs.xh_xt[(xh_port.index(), xt_port.index())] = false;
                }
            }
            for xp_port in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free {
                for p1_port in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p1_free{
                    self.unary_binding_pairs.p1_xp[(p1_port.index(), xp_port.index())] = false;
                }
                for p2_port in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p2_free{
                    self.unary_binding_pairs.p2_xp[(p2_port.index(), xp_port.index())] = false;
                }
                for p3_port in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p3_free{
                    self.unary_binding_pairs.p3_xp[(p3_port.index(), xp_port.index())] = false;
                }
            }
            for xp_port in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free {
                for p1_port in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.p1_free{
                    self.unary_binding_pairs.p1_xp[(p1_port.index(), xp_port.index())] = false;
                }
                for p2_port in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.p2_free{
                    self.unary_binding_pairs.p2_xp[(p2_port.index(), xp_port.index())] = false;
                }
                for p3_port in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.p3_free{
                    self.unary_binding_pairs.p3_xp[(p3_port.index(), xp_port.index())] = false;
                }
            }
            for p1_port in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.pp_free{
                for p2_port in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.pp_free {
                    self.unary_binding_pairs.xh_xt[(p1_port.index(), p2_port.index())] = false;
                }
            }
            for p1_port in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.pp_free{
                for p2_port in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.pp_free {
                    self.unary_binding_pairs.xh_xt[(p1_port.index(), p2_port.index())] = false;
                }
            }
            // update self.rule_activities
            let binary_embeds_made_here: usize = 1;
            // since self-bonds are not allowed in this model,
            //  this special-cases the Axn-Axn unary treatment,
            //  discounting the monomer self-binding
            let unary_combinatorics_gained_head: usize = 
                if self.species_annots.get(&head_node).unwrap().borrow().size > 1 {
                    self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free.len()
                }
                else {0};
            let unary_combinatorics_gained_tail: usize = 
                if self.species_annots.get(&tail_node).unwrap().borrow().size > 1 { 
                    self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free.len()
                } else {0};
            let unary_combinatorics_gained: usize = unary_combinatorics_gained_head + unary_combinatorics_gained_tail;
            let transformed_to_binary: usize = 
                (self.species_annots.get(&tail_node).unwrap().borrow().ports.xt_free.len() * self.species_annots.get(&head_node).unwrap().borrow().ports.xh_free.len()) + 
                (self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free.len() * self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free.len()) -
                binary_embeds_made_here;                                                                                                    // correct for the bond we just broke
            // for the Axn-Axn case, we can ignore the _free'ing activities
            //  for the other bond types, as we are not modifying those bond
            //  counts, ergo their unbinding counts
            let binary_combinatorics_gained_head: usize = self.ports.xt_free.len() - self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free.len() - 1;    // the bond that would bind the tail & head nodes will
            let binary_combinatorics_gained_tail: usize = self.ports.xh_free.len() - self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free.len();        // be present in both; substract 1 to correct count
            let binary_combinatorics_gained: usize = binary_combinatorics_gained_head + binary_combinatorics_gained_tail;
            println!("bis gained {}, unis gained {}, transformed to binary {}", binary_combinatorics_gained, unary_combinatorics_gained, transformed_to_binary);
            self.rule_activities.axn_axn_b_free.mass -= binary_embeds_made_here;                                                            // the bond we just broke
            //self.rule_activities.axn_axn_u_free.mass      // does not apply
            self.rule_activities.axn_axn_b_bind.mass += binary_combinatorics_gained + transformed_to_binary;
            self.rule_activities.axn_axn_u_bind.mass += unary_combinatorics_gained;             // broken into two operations to avoid underflowing
            self.rule_activities.axn_axn_u_bind.mass -= transformed_to_binary;                  // when the former is 0, but the latter is non-zero
            //  the rest are unary -> binary conversions
            self.rule_activities.ap1_axn_b_bind.mass += (self.species_annots.get(&retained_mark).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap1_axn_u_bind.mass -= (self.species_annots.get(&retained_mark).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p1_free.len() * self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap2_axn_b_bind.mass += (self.species_annots.get(&retained_mark).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap2_axn_u_bind.mass -= (self.species_annots.get(&retained_mark).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p2_free.len() * self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap3_axn_b_bind.mass += (self.species_annots.get(&retained_mark).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.ap3_axn_u_bind.mass -= (self.species_annots.get(&retained_mark).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.len()) + (self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p3_free.len() * self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free.len());
            self.rule_activities.apc_apc_b_bind.mass += (self.species_annots.get(&retained_mark).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&ejected_mark).unwrap().borrow().ports.pp_free.len()) + (self.species_annots.get(&ejected_mark).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&retained_mark).unwrap().borrow().ports.pp_free.len());
            self.rule_activities.apc_apc_u_bind.mass -= (self.species_annots.get(&retained_mark).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&ejected_mark).unwrap().borrow().ports.pp_free.len()) + (self.species_annots.get(&ejected_mark).unwrap().borrow().ports.pp_free.len() * self.species_annots.get(&retained_mark).unwrap().borrow().ports.pp_free.len());
        }

        pub fn print_rule_activities(&self) {
            println!("{}", self.rule_activities)
        }
    }
}