pub mod agent_types {
    use petgraph::prelude::*;

    #[derive(Debug)]
    /**
     * Structure for an Axn type-node. The fields hold either a `None` when unbound, or a 
     * `Some(EdgeIndex)` when bound. One field per binding site.
    */
    pub struct AxnNode {
        pub x1_bond: Option<EdgeIndex>,
        pub x2_bond: Option<EdgeIndex>,
        pub xp_bond: Option<EdgeIndex>
    }   

    /**
     * By default, Axn nodes are created fully unbound.
    */
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
     * Structure for an APC type-node. The fields hold either a `None` when unbound, or a
     * `Some(EdgeIndex)` when bound. One field per binding site.
    */
    pub struct ApcNode {
        pub p1_bond: Option<EdgeIndex>,
        pub p2_bond: Option<EdgeIndex>,
        pub p3_bond: Option<EdgeIndex>,
        pub pp_bond: Option<EdgeIndex>
    }   

    /**
     * By default, APC nodes are created fully unbound.
    */
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
    use std::fmt;
    use std::cmp::Ordering;
    use petgraph::prelude::*;

    #[derive(Clone, Copy, Debug)]
    /** 
     * Structure whose fields store the NodeIndexes of the binding pair in `a` & `b`, and the
     * EdgeIndex in `z` from the UniverseGraph, if any. The bond is read in an ordered fashion,
     * and this avoids issues with parallel edges, as the NodeIndexes are globally unique, and
     * the contact map does not have parallel edges when read in agent-name-alphabetic order.
    */
    pub struct EdgeEnds {
        pub a: NodeIndex,
        pub b: NodeIndex,
        pub z: Option<EdgeIndex>
    }
    
    impl PartialEq for EdgeEnds {
        fn eq(&self, other: &Self) -> bool {
            (self.z == other.z) && ( (self.a == other.a && self.b == other.b) || (self.a == other.b && self.b == other.a) )
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

    impl fmt::Display for EdgeEnds {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            let str_rep = match self.z {
                Some(i) => format!("#{} @ {} -> {}", i.index(), self.a.index(), self.b.index()),
                None => format!("{} -> {}", self.a.index(), self.b.index())
            };
            write!(f, "{}", str_rep)
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
    use std::fmt;
    use std::rc::Rc;
    use std::cell::RefCell;
    use petgraph::prelude::*;
    use crate::edge_ends::EdgeEnds;

    #[derive(Clone, Debug, PartialEq)]
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

    impl fmt::Display for EdgeTypes {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            let mut str_rep: String = String::new();
            let helper_closure = |input_set: &BTreeSet<Rc<RefCell<EdgeEnds>>>| -> String {
                let mut temp_str = String::new();
                for boxed_edge in input_set {
                    temp_str.push_str(&format!("\n\t{}", boxed_edge.borrow()))
                }
                temp_str
            };
            str_rep.push_str("xh_xt:");
            str_rep.push_str(&helper_closure(&self.xh_xt));
            str_rep.push_str("\np1_xp:");
            str_rep.push_str(&helper_closure(&self.p1_xp));
            str_rep.push_str("\np2_xp:");
            str_rep.push_str(&helper_closure(&self.p2_xp));
            str_rep.push_str("\np3_xp:");
            str_rep.push_str(&helper_closure(&self.p3_xp));
            str_rep.push_str("\npp_pp:");
            str_rep.push_str(&helper_closure(&self.pp_pp));
            write!(f, "{}\n", str_rep)
        }
    }
}

pub mod rule_activities {
    use std::fmt;

    #[derive(Clone, Debug, PartialEq)]
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

    impl RuleActivities {
        pub fn calculate_rule_activities(&self) -> Vec<f64>{
            vec![
                self.axn_axn_u_bind.calculate_activity(),
                self.axn_axn_b_bind.calculate_activity(),
                self.axn_axn_u_free.calculate_activity(),
                self.axn_axn_b_free.calculate_activity(),
                self.ap1_axn_u_bind.calculate_activity(),
                self.ap1_axn_b_bind.calculate_activity(),
                self.ap1_axn_u_free.calculate_activity(),
                self.ap1_axn_b_free.calculate_activity(),
                self.ap2_axn_u_bind.calculate_activity(),
                self.ap2_axn_b_bind.calculate_activity(),
                self.ap2_axn_u_free.calculate_activity(),
                self.ap2_axn_b_free.calculate_activity(),
                self.ap3_axn_u_bind.calculate_activity(),
                self.ap3_axn_b_bind.calculate_activity(),
                self.ap3_axn_u_free.calculate_activity(),
                self.ap3_axn_b_free.calculate_activity(),
                self.apc_apc_u_bind.calculate_activity(),
                self.apc_apc_b_bind.calculate_activity(),
                self.apc_apc_u_free.calculate_activity(),
                self.apc_apc_b_free.calculate_activity()
            ]
        }

        pub fn calculate_system_activity(&self) -> f64 {
            self.calculate_rule_activities().iter().sum()
        }
    }

    #[derive(Clone, Debug, PartialEq)]
    pub struct MassActionTerm {
        pub mass: usize,
        pub rate: f64,
    }

    impl fmt::Display for MassActionTerm {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "|{}| * {} => {}", self.mass, self.rate, self.calculate_activity())
        }
    }

    impl MassActionTerm {
        pub fn calculate_activity(&self) -> f64 {
            self.mass as f64 * self.rate
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

pub mod possible_bond_embeds {
    use itertools::Itertools;
    use std::{fmt, collections::HashSet};
    use petgraph::prelude::NodeIndex;

    #[derive(Clone, Debug, PartialEq)]
    /** 
     * Structure of hashsets, one per bond type. The value contained in them is a tuple of 
     * NodeIndexes, specifying the oriented bond's origin & destination. These are pairs
     * that could bind.
    */
    pub struct PossibleBondEmbeds {
        pub xh_xt: HashSet<(NodeIndex, NodeIndex)>,
        pub p1_xp: HashSet<(NodeIndex, NodeIndex)>,
        pub p2_xp: HashSet<(NodeIndex, NodeIndex)>,
        pub p3_xp: HashSet<(NodeIndex, NodeIndex)>,
        pub pp_pp: HashSet<(NodeIndex, NodeIndex)>
    }

    impl PossibleBondEmbeds {
        /**
         * Used to track unary instances; constructor allocates maximum required space but leaves
         * the sets empty.
        */
        pub fn new_from_masses_unary(axn_mass: usize, apc_mass: usize) -> Self {
            PossibleBondEmbeds {
                xh_xt: HashSet::with_capacity(axn_mass * axn_mass),
                p1_xp: HashSet::with_capacity(apc_mass * axn_mass),
                p2_xp: HashSet::with_capacity(apc_mass * axn_mass),
                p3_xp: HashSet::with_capacity(apc_mass * axn_mass),
                pp_pp: HashSet::with_capacity(apc_mass * apc_mass),
            }
        }

        /**
         * Used to track binary instances; constructor allocates maximum required space, and
         * populates the sets with the appropriate NodeIndex combinatorics, cleaning up 
         * self-bonds afterwards.
        */
        pub fn new_from_masses_binary(axn_mass: usize, apc_mass: usize) -> Self {
            let x_mass: u32 = axn_mass as u32;
            let p_mass: u32 = apc_mass as u32;
            let mut xh_xt = (0..x_mass).map(|x| NodeIndex::from(x)).cartesian_product((0..x_mass).map(|x| NodeIndex::from(x))).collect::<HashSet<(NodeIndex, NodeIndex)>>();
            let p1_xp = (x_mass..(x_mass + p_mass)).map(|x| NodeIndex::from(x)).cartesian_product((0..x_mass).map(|x| NodeIndex::from(x))).collect::<HashSet<(NodeIndex, NodeIndex)>>();
            let p2_xp = (x_mass..(x_mass + p_mass)).map(|x| NodeIndex::from(x)).cartesian_product((0..x_mass).map(|x| NodeIndex::from(x))).collect::<HashSet<(NodeIndex, NodeIndex)>>();
            let p3_xp = (x_mass..(x_mass + p_mass)).map(|x| NodeIndex::from(x)).cartesian_product((0..x_mass).map(|x| NodeIndex::from(x))).collect::<HashSet<(NodeIndex, NodeIndex)>>();
            let pp_pp_vect_combs = (x_mass..(x_mass + p_mass)).map(|x| NodeIndex::from(x)).combinations(2);
            let mut pp_pp: HashSet<(NodeIndex, NodeIndex)> = HashSet::with_capacity(apc_mass * apc_mass / 2);
            for mut elem in pp_pp_vect_combs {
                let node_b: NodeIndex = elem.pop().unwrap();
                let node_a: NodeIndex = elem.pop().unwrap();
                pp_pp.insert((node_a, node_b));
            }
            // cleanup of self-bonds
            for this_axn in 0..x_mass {xh_xt.remove(&(NodeIndex::from(this_axn), NodeIndex::from(this_axn)));}
            PossibleBondEmbeds {xh_xt, p1_xp, p2_xp, p3_xp, pp_pp}
        }
    }

    impl fmt::Display for PossibleBondEmbeds {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            let mut str_rep: String = String::new();
            let helper_closure = |input_set: &HashSet<(NodeIndex, NodeIndex)>| -> String {
                let mut temp_str = String::new();
                let mut tuple_vec = input_set.into_iter().cloned().collect::<Vec<(NodeIndex, NodeIndex)>>();
                tuple_vec.sort();
                for (side_a, side_b) in tuple_vec {
                    temp_str.push_str(&format!("\n\t({}, {})", side_a.index(), side_b.index()));
                }
                temp_str
            };
            str_rep.push_str("xh_xt:");
            str_rep.push_str(&helper_closure(&self.xh_xt));
            str_rep.push_str("\np1_xp:");
            str_rep.push_str(&helper_closure(&self.p1_xp));
            str_rep.push_str("\np2_xp:");
            str_rep.push_str(&helper_closure(&self.p2_xp));
            str_rep.push_str("\np3_xp:");
            str_rep.push_str(&helper_closure(&self.p3_xp));
            str_rep.push_str("\npp_pp:");
            str_rep.push_str(&helper_closure(&self.pp_pp));
            write!(f, "{}\n", str_rep)
        }
    }
}

pub mod reaction_mixture {
    use crate::open_ports::OpenPorts;
    use crate::edge_types::EdgeTypes;
    use crate::edge_ends::EdgeEnds;
    use crate::agent_types::AgentType;
    use crate::rule_activities::{RuleActivities, RuleRates, MassActionTerm};
    use crate::possible_bond_embeds::PossibleBondEmbeds;
    use petgraph::{graph::Graph, algo::astar::astar, prelude::*, visit::Dfs};
    use rand::{distributions::WeightedIndex, prelude::*};
    use std::cell::{RefCell, RefMut};
    use std::cmp::Ordering;
    use std::collections::{VecDeque, BTreeMap, BTreeSet, HashSet};
    use std::rc::Rc;

    #[derive(Debug)]
    struct MixtureSpecies {
        ports: OpenPorts,
        edges: EdgeTypes,
        size: usize,
        agent_set: BTreeSet<NodeIndex>,
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
            if self.size == other.size {
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
            else {self.size.cmp(&other.size)}
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
        species_set: VecDeque<Rc<RefCell<MixtureSpecies>>>,                 // used for iteration, pretty printing, & "final persistent owner"
        ports: OpenPorts,                                                   // struct with vector, each a list of agent indexes
        cycle_edges: EdgeTypes,                                             // to store edges whose breaking does not yield fission of its complex
        tree_edges: EdgeTypes,                                              // to store edges whose breaking yields fission of its complex
        edge_index_map: BTreeMap<EdgeIndex, Rc<RefCell<EdgeEnds>>>,         // used to update Z when removing bonds from the universe graph
        unary_binding_pairs: PossibleBondEmbeds,
        binary_binding_pairs: PossibleBondEmbeds,
        rule_activities: RuleActivities,
        simulated_time: f64,
        simulated_events: usize,
        simulator_rng: ThreadRng
    }
    
    impl Mixture {
        fn advance_simulation_metrics(&mut self) {
            self.simulated_events += 1;
            let tau: f64 = 1.0 / self.rule_activities.calculate_system_activity();
            let rdy: f64 = (1.0 / self.simulator_rng.gen::<f64>()).ln();
            self.simulated_time += tau * rdy;
        }

        /**
         * Update the mass trackers, based on the potential embedding trackers.
        */
        fn update_mass_actions(&mut self) {
            self.rule_activities.axn_axn_b_free.mass = self.tree_edges.xh_xt.len();
            self.rule_activities.axn_axn_u_free.mass = self.cycle_edges.xh_xt.len();
            self.rule_activities.axn_axn_b_bind.mass = self.binary_binding_pairs.xh_xt.len();
            self.rule_activities.axn_axn_u_bind.mass = self.unary_binding_pairs.xh_xt.len();
            self.rule_activities.ap1_axn_b_free.mass = self.tree_edges.p1_xp.len();
            self.rule_activities.ap1_axn_u_free.mass = self.cycle_edges.p1_xp.len();
            self.rule_activities.ap1_axn_b_bind.mass = self.binary_binding_pairs.p1_xp.len();
            self.rule_activities.ap1_axn_u_bind.mass = self.unary_binding_pairs.p1_xp.len();
            self.rule_activities.ap2_axn_b_free.mass = self.tree_edges.p2_xp.len();
            self.rule_activities.ap2_axn_u_free.mass = self.cycle_edges.p2_xp.len();
            self.rule_activities.ap2_axn_b_bind.mass = self.binary_binding_pairs.p2_xp.len();
            self.rule_activities.ap2_axn_u_bind.mass = self.unary_binding_pairs.p2_xp.len();
            self.rule_activities.ap3_axn_b_free.mass = self.tree_edges.p3_xp.len();
            self.rule_activities.ap3_axn_u_free.mass = self.cycle_edges.p3_xp.len();
            self.rule_activities.ap3_axn_b_bind.mass = self.binary_binding_pairs.p3_xp.len();
            self.rule_activities.ap3_axn_u_bind.mass = self.unary_binding_pairs.p3_xp.len();
            self.rule_activities.apc_apc_b_free.mass = self.tree_edges.pp_pp.len();
            self.rule_activities.apc_apc_u_free.mass = self.cycle_edges.pp_pp.len();
            self.rule_activities.apc_apc_b_bind.mass = self.binary_binding_pairs.pp_pp.len();
            self.rule_activities.apc_apc_u_bind.mass = self.unary_binding_pairs.pp_pp.len();
        }

        pub fn choose_and_apply_next_rule(&mut self) {
            let dist = WeightedIndex::new(self.rule_activities.calculate_rule_activities()).unwrap();
            let chosen_rule: usize = dist.sample(&mut self.simulator_rng);
            match chosen_rule {
                0 => self.pick_targets_axn_axn_unary_bind(),
                1 => self.pick_targets_axn_axn_binary_bind(),
                2 => self.pick_targets_axn_axn_unary_unbind(),
                3 => self.pick_targets_axn_axn_binary_unbind(),
                _ => panic!("We are not here yet!")
            }
            self.advance_simulation_metrics();
        }

        fn pick_targets_axn_axn_unary_bind(&mut self) {
            let (target_a, target_b): (NodeIndex, NodeIndex) = self.unary_binding_pairs.xh_xt.iter().choose(&mut self.simulator_rng).unwrap().clone();
            //println!("Applying Axn-Axn Unary Bind to {} and {}", target_a.index(), target_b.index());
            self.axn_axn_unary_bind(target_a, target_b);
        }

        fn pick_targets_axn_axn_binary_bind(&mut self) {
            let (target_a, target_b): (NodeIndex, NodeIndex) = self.binary_binding_pairs.xh_xt.iter().choose(&mut self.simulator_rng).unwrap().clone();
            //println!("Applying Axn-Axn Binary Bind to {} and {}", target_a.index(), target_b.index());
            self.axn_axn_binary_bind(target_a, target_b);
        }

        fn pick_targets_axn_axn_unary_unbind(&mut self) {
            let target_edge: Rc<RefCell<EdgeEnds>> = self.cycle_edges.xh_xt.iter().choose(&mut self.simulator_rng).unwrap().clone();
            //println!("Applying Axn-Axn Unary Unbind to {}", target_edge.borrow());
            self.axn_axn_unary_unbind(target_edge);
        }

        fn pick_targets_axn_axn_binary_unbind(&mut self) {
            let target_edge: Rc<RefCell<EdgeEnds>> = self.tree_edges.xh_xt.iter().choose(&mut self.simulator_rng).unwrap().clone();
            //println!("Applying Axn-Axn Binary Unbind to {}", target_edge.borrow());
            self.axn_axn_binary_unbind(target_edge);
        }

        /**
         * Yield a snapshot representation that can be printed or compared against others in a
         * kappa-aware way.
        */
        pub fn to_kappa(&self) -> String {
            let mut mixture_string: String = String::new();
            mixture_string.push_str(&format!("// Snapshot [Event:{}]\n", self.simulated_events));
            mixture_string.push_str(&format!("// Mixture has {} nodes and {} edges\n", self.universe_graph.node_count(), self.universe_graph.edge_count()));
            mixture_string.push_str(&format!("%def: \"T0\" \"{}\"\n\n", self.simulated_time));
            let mut cloned_species = self.species_set.clone();
            cloned_species.make_contiguous().sort_by(|a, b| b.cmp(a));
            for this_species_ref in cloned_species {
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

        /**
         * Create a new mixture from a desired number of monomers, and desired rule rates.
        */
        pub fn new_from_monomers(x_mass: usize, p_mass: usize, rule_rates: RuleRates) -> Mixture {
            let mut net: Graph<AgentType, bool, Undirected> = Graph::with_capacity(x_mass + p_mass, x_mass * 3 + p_mass * 4);   // the universe graph
            let mut spc: BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies>>> = BTreeMap::new();                                    // the global agent -> species map
            let mut sps: VecDeque<Rc<RefCell<MixtureSpecies>>> = VecDeque::with_capacity(x_mass + p_mass);                      // the species set
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
            let bbp = PossibleBondEmbeds::new_from_masses_binary(x_mass, p_mass);
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
                axn_axn_b_bind: MassActionTerm{mass: bbp.xh_xt.len(), rate: rule_rates.axn_axn_b_bind},
                ap1_axn_b_bind: MassActionTerm{mass: bbp.p1_xp.len(), rate: rule_rates.ap1_axn_b_bind},
                ap2_axn_b_bind: MassActionTerm{mass: bbp.p2_xp.len(), rate: rule_rates.ap2_axn_b_bind},
                ap3_axn_b_bind: MassActionTerm{mass: bbp.p3_xp.len(), rate: rule_rates.ap3_axn_b_bind},
                apc_apc_b_bind: MassActionTerm{mass: bbp.pp_pp.len(), rate: rule_rates.apc_apc_b_bind}
            };
            // bringing it all together
            Mixture {universe_graph: net, species_annots: spc, ports: opr, rule_activities: rac, species_set: sps,
                cycle_edges: EdgeTypes::default_empty(), tree_edges: EdgeTypes::default_empty(),  
                unary_binding_pairs: PossibleBondEmbeds::new_from_masses_unary(x_mass, p_mass),
                binary_binding_pairs: bbp,
                edge_index_map: BTreeMap::new(), simulator_rng: rand::thread_rng(),
                simulated_events: 0, simulated_time: 0.0}
        }
    
        /**
         * Create a bond that creates a cycle, connecting two agents that were already connected
         * via some other components.
        */
        pub fn axn_axn_unary_bind(&mut self, head_node: NodeIndex, tail_node: NodeIndex) {
            // sanity check for unary-ness
            assert_eq!(self.species_annots.get(&head_node), self.species_annots.get(&tail_node), "This head and tail nodes ({}, {}) do not already belong to the same species; can't unary bind them!", head_node.index(), tail_node.index());
            {
                let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_node).unwrap().borrow_mut();
                // create edge & update caches
                let new_edge_index = self.universe_graph.add_edge(head_node, tail_node, true);
                let new_edge = Rc::new(RefCell::new(EdgeEnds{a: head_node, b: tail_node, z: Some(new_edge_index)}));
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
                assert!(self.cycle_edges.xh_xt.insert(Rc::clone(&new_edge)));
                // bindings no longer possible due to occupied sites
                let mut binary_combinatorics_lost: HashSet<(NodeIndex, NodeIndex)> = HashSet::new();
                for tail in &self.ports.xt_free.difference(&target_species.ports.xt_free).cloned().collect::<Vec<NodeIndex>>() {
                    match self.binary_binding_pairs.xh_xt.take(&(head_node, *tail)) {
                        Some((a_head, a_tail)) => {assert!(binary_combinatorics_lost.insert((a_head, a_tail)), "This lost pair already cached into the binary lost tracker!")}
                        None => ()
                    }
                }
                for head in &self.ports.xh_free.difference(&target_species.ports.xh_free).cloned().collect::<Vec<NodeIndex>>() {
                    match self.binary_binding_pairs.xh_xt.take(&(*head, tail_node)) {
                        Some((a_head, a_tail)) => {assert!(binary_combinatorics_lost.insert((a_head, a_tail)), "This lost pair already cached into the binary lost racker!")},
                        None => ()
                    }
                }
                let mut unary_combinatorics_lost: HashSet<(NodeIndex, NodeIndex)> = HashSet::new();
                for tail in &target_species.ports.xt_free {
                    match self.unary_binding_pairs.xh_xt.take(&(head_node, *tail)) {
                        Some((a_head, a_tail)) => {assert!(unary_combinatorics_lost.insert((a_head, a_tail)), "This lost pair already cached into the unary lost tracker!")},
                        None => ()
                    }
                }
                for head in &target_species.ports.xh_free {
                    match self.unary_binding_pairs.xh_xt.take(&(*head, tail_node)) {
                        Some((a_head, a_tail)) => {assert!(unary_combinatorics_lost.insert((a_head, a_tail)), "This lost pair already cached into the unary los tracker!")},
                        None => ()
                    }
                }
                assert!(self.ports.xh_free.remove(&head_node));
                assert!(self.ports.xt_free.remove(&tail_node));
                assert!(target_species.ports.xh_free.remove(&head_node));
                assert!(target_species.ports.xt_free.remove(&tail_node));
                // helper closure for bond resilience updates
                let mut resilience_iterate_typed = |edge_iter: &BTreeSet<std::rc::Rc<std::cell::RefCell<EdgeEnds>>>| -> BTreeSet<Rc<RefCell<EdgeEnds>>> {
                    let mut bonds_cyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
                    for boxed_edge in edge_iter {
                        let edge = boxed_edge.borrow();
                        // check bond's weight, which is a boolean, marking if the bond is a known cycle-member
                        if ! self.universe_graph.edge_weight(edge.z.unwrap()).unwrap() {
                            let mut counterfactual_graph = self.universe_graph.clone();
                            counterfactual_graph.remove_edge(edge.z.unwrap());
                            let path = astar(&counterfactual_graph, edge.a, |finish| finish == edge.b, |_| 1, |_| 0);
                            match path {
                                Some(_) => {self.universe_graph.update_edge(edge.a, edge.b, true); assert!(bonds_cyclized.insert(Rc::clone(&boxed_edge)), "This edge already seen; can't insert it into the set!")},
                                None => ()
                            }
                        }
                    }
                    bonds_cyclized
                };
                // update resilience of that species' edges; update edge trackers
                let mut xh_xt_cyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.xh_xt);
                let mut p1_xp_cyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.p1_xp);
                let mut p2_xp_cyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.p2_xp);
                let mut p3_xp_cyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.p3_xp);
                let mut pp_pp_cyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.pp_pp);
                let tree_edges_retained_xh_xt = &self.tree_edges.xh_xt - &xh_xt_cyclized;
                let tree_edges_retained_p1_xp = &self.tree_edges.p1_xp - &p1_xp_cyclized;
                let tree_edges_retained_p2_xp = &self.tree_edges.p2_xp - &p2_xp_cyclized;
                let tree_edges_retained_p3_xp = &self.tree_edges.p3_xp - &p3_xp_cyclized;
                let tree_edges_retained_pp_pp = &self.tree_edges.pp_pp - &pp_pp_cyclized;
                self.cycle_edges.xh_xt.append(&mut xh_xt_cyclized);
                self.cycle_edges.p1_xp.append(&mut p1_xp_cyclized);
                self.cycle_edges.p2_xp.append(&mut p2_xp_cyclized);
                self.cycle_edges.p3_xp.append(&mut p3_xp_cyclized);
                self.cycle_edges.pp_pp.append(&mut pp_pp_cyclized);
                self.tree_edges.xh_xt = tree_edges_retained_xh_xt;
                self.tree_edges.p1_xp = tree_edges_retained_p1_xp;
                self.tree_edges.p2_xp = tree_edges_retained_p2_xp;
                self.tree_edges.p3_xp = tree_edges_retained_p3_xp;
                self.tree_edges.pp_pp = tree_edges_retained_pp_pp;
            }
            self.update_mass_actions();
        }
    
        /**
         * Create a bond that joins two species into one.
        */
        pub fn axn_axn_binary_bind(&mut self, head_node: NodeIndex, tail_node: NodeIndex) {
            // sanity check for binary-ness
            assert_ne!(self.species_annots.get(&head_node), self.species_annots.get(&tail_node), "This head and tail nodes ({}, {}) already belong to the same species; can't binary bind them!", head_node.index(), tail_node.index());
            let (host_index, eaten_index) = if self.species_annots.get(&head_node).unwrap().borrow().size >= self.species_annots.get(&tail_node).unwrap().borrow().size {(head_node, tail_node)} else {(tail_node, head_node)};
            // bindings no longer possible due to occupied sites
            let mut unary_combinatorics_lost: HashSet<(NodeIndex, NodeIndex)> = HashSet::with_capacity(self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free.len() + self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free.len());
            for tail in &self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free {
                match self.unary_binding_pairs.xh_xt.take(&(head_node, *tail)) {
                    Some((a_head, a_tail)) => {assert!(unary_combinatorics_lost.insert((a_head, a_tail)), "This lost pair already cached into the unary lost tracker!")},
                    None => ()
                }
            }
            for head in &self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free {
                match self.unary_binding_pairs.xh_xt.take(&(*head, tail_node)) {
                    Some((a_head, a_tail)) => {assert!(unary_combinatorics_lost.insert((a_head, a_tail)), "This lost pair already ached into the unary lost tracker!")},
                    None => ()
                }
            }
            let mut binary_combinatorics_lost: HashSet<(NodeIndex, NodeIndex)> = HashSet::new();
            for tail in &self.ports.xt_free.difference(&self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free).cloned().collect::<Vec<NodeIndex>>() {
                match self.binary_binding_pairs.xh_xt.take(&(head_node, *tail)) {
                    Some((a_head, a_tail)) => {assert!(binary_combinatorics_lost.insert((a_head, a_tail)), "This lost pair already cached into the binary lost tracker!")},
                    None => ()
                }
            }
            for head in &self.ports.xh_free.difference(&self.species_annots.get(&head_node).unwrap().borrow().ports.xh_free).cloned().collect::<Vec<NodeIndex>>() {
                match self.binary_binding_pairs.xh_xt.take(&(*head, tail_node)) {
                    Some((a_head, a_tail)) => {assert!(binary_combinatorics_lost.insert((a_head, a_tail)), "This lost pair already cached into the binary lost tracker!")},
                    None => ()
                }
            }
            assert!(self.species_annots.get(&head_node).unwrap().borrow_mut().ports.xh_free.remove(&head_node), "This head node was not listed as free in-species already, can't bind to it!");
            assert!(self.species_annots.get(&tail_node).unwrap().borrow_mut().ports.xt_free.remove(&tail_node), "This tail node was not listed as free in-species already, can't bind to it!");
            assert!(self.ports.xh_free.remove(&head_node), "This head node was not globally listed as free already, can't bind to it!");
            assert!(self.ports.xt_free.remove(&tail_node), "This tail node was not globally listed as free already, can't bind to it!");
            // binding opportunities that were binary, but will now be unary
            for free_head in self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free.union(&self.species_annots.get(&head_node).unwrap().borrow().ports.xh_free).into_iter() {
                for free_tail in self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free.union(&self.species_annots.get(&tail_node).unwrap().borrow().ports.xt_free).into_iter() {
                    match self.binary_binding_pairs.xh_xt.take(&(*free_head, *free_tail)) {
                        Some((a_head, a_tail)) => {assert!(self.unary_binding_pairs.xh_xt.insert((a_head, a_tail)), "This pair already cached as unary!")},
                        None => ()
                    }
                }
            }
            for free_head in self.species_annots.get(&tail_node).unwrap().borrow().ports.p1_free.union(&self.species_annots.get(&head_node).unwrap().borrow().ports.p1_free).into_iter() {
                for free_tail in self.species_annots.get(&head_node).unwrap().borrow().ports.xp_free.union(&self.species_annots.get(&tail_node).unwrap().borrow().ports.xp_free).into_iter() {
                    match self.binary_binding_pairs.p1_xp.take(&(*free_head, *free_tail)) {
                        Some((a_head, a_tail)) => {assert!(self.unary_binding_pairs.p1_xp.insert((a_head, a_tail)), "This pair already cached as unary!")},
                        None => ()
                    }
                }
            }
            for free_head in self.species_annots.get(&tail_node).unwrap().borrow().ports.p2_free.union(&self.species_annots.get(&head_node).unwrap().borrow().ports.p2_free).into_iter() {
                for free_tail in self.species_annots.get(&head_node).unwrap().borrow().ports.xp_free.union(&self.species_annots.get(&tail_node).unwrap().borrow().ports.xp_free).into_iter() {
                    match self.binary_binding_pairs.p2_xp.take(&(*free_head, *free_tail)) {
                        Some((a_head, a_tail)) => {assert!(self.unary_binding_pairs.p2_xp.insert((a_head, a_tail)), "This pair already cached as unary!")},
                        None => ()
                    }
                }
            }
            for free_head in self.species_annots.get(&tail_node).unwrap().borrow().ports.p3_free.union(&self.species_annots.get(&head_node).unwrap().borrow().ports.p3_free).into_iter() {
                for free_tail in self.species_annots.get(&head_node).unwrap().borrow().ports.xp_free.union(&self.species_annots.get(&tail_node).unwrap().borrow().ports.xp_free).into_iter() {
                    match self.binary_binding_pairs.p3_xp.take(&(*free_head, *free_tail)) {
                        Some((a_head, a_tail)) => {assert!(self.unary_binding_pairs.p3_xp.insert((a_head, a_tail)), "This pair already cached as unary!")},
                        None => ()
                    }
                }
            }
            for free_head in &self.species_annots.get(&tail_node).unwrap().borrow().ports.pp_free {
                for free_tail in &self.species_annots.get(&head_node).unwrap().borrow().ports.pp_free {
                    match self.binary_binding_pairs.pp_pp.take(&(*free_head, *free_tail)) {
                        Some((a_head, a_tail)) => {assert!(self.unary_binding_pairs.pp_pp.insert((a_head, a_tail)), "This pair already cached as unary!")},
                        None => ()
                    }
                }
            }
            self.species_annots.get(&host_index).unwrap().borrow_mut().ports.update_from(&self.species_annots.get(&eaten_index).unwrap().borrow().ports);
            self.species_annots.get(&host_index).unwrap().borrow_mut().edges.update_from(&self.species_annots.get(&eaten_index).unwrap().borrow().edges);
            self.species_annots.get(&host_index).unwrap().borrow_mut().size += self.species_annots.get(&eaten_index).unwrap().borrow().size;
            self.species_annots.get(&host_index).unwrap().borrow_mut().agent_set.append(&mut self.species_annots.get(&eaten_index).unwrap().borrow().agent_set.clone());
            // create & update edge trackers
            let new_edge_index = self.universe_graph.add_edge(head_node, tail_node, false);
            let new_edge = Rc::new(RefCell::new(EdgeEnds{a: head_node, b: tail_node, z: Some(new_edge_index)}));
            self.edge_index_map.insert(new_edge_index, Rc::clone(&new_edge));
            self.species_annots.get(&host_index).unwrap().borrow_mut().edges.xh_xt.insert(Rc::clone(&new_edge));
            self.tree_edges.xh_xt.insert(Rc::clone(&new_edge));
            // update species annotation tracker
            self.species_set.make_contiguous().sort();
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
            self.update_mass_actions();
        }

        /**
         * Break a bond that is a cycle member; this does not partition one species into two, it
         * simply reduces the internal connectivity of said species.
        */
        pub fn axn_axn_unary_unbind(&mut self, target_edge: Rc<RefCell<EdgeEnds>>) {
            if let EdgeEnds{a: head_node, b: tail_node, z: Some(edge_index)} = *target_edge.borrow() { 
                // sanity check for reslient bond
                assert!(self.universe_graph.edge_weight(edge_index).unwrap(), "This bond ({}) is not flagged as resilient; breaking it would not yield a still-connected component!", target_edge.borrow());
                let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_node).unwrap().borrow_mut();
                let old_last_edge_index = self.universe_graph.edge_indices().last().unwrap();
                self.universe_graph.remove_edge(edge_index).unwrap();
                let swapped_edge = self.edge_index_map.remove(&old_last_edge_index).unwrap();
                if old_last_edge_index != edge_index {
                    // remove edge & update caches; keep edge index tracker up to date
                    swapped_edge.borrow_mut().z = Some(edge_index);
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
                assert!(self.cycle_edges.xh_xt.remove(&target_edge));
                assert!(self.ports.xh_free.insert(head_node), "This Axn head was already listed as free prior to this unbinding!");
                assert!(self.ports.xt_free.insert(tail_node), "This Axn tail was already listed as free prior to this unbinding!");
                assert!(target_species.ports.xh_free.insert(head_node), "This Axn head was already listed as free prior to this unbinding!");
                assert!(target_species.ports.xt_free.insert(tail_node), "This Axn tail was already listed as free prior to this unbinding!");
                // bindings now possible due to free'd sites
                for tail in &self.ports.xt_free.difference(&target_species.ports.xt_free).cloned().collect::<Vec<NodeIndex>>() {
                    assert!(self.binary_binding_pairs.xh_xt.insert((head_node, *tail)), "This gained pair already cached into the binary tracker!")}
                for head in &self.ports.xh_free.difference(&target_species.ports.xh_free).cloned().collect::<Vec<NodeIndex>>() {
                    assert!(self.binary_binding_pairs.xh_xt.insert((*head, tail_node)), "This gained pair already cached into the binary tracker!")}
                for tail in &target_species.ports.xt_free {
                    assert!(self.unary_binding_pairs.xh_xt.insert((head_node, *tail)), "This gained pair already cached into the unary tracker!")}
                self.unary_binding_pairs.xh_xt.remove(&(head_node, tail_node));     // hack to remove pair that will be visited from the other side in the next loop
                for head in &target_species.ports.xh_free {
                    assert!(self.unary_binding_pairs.xh_xt.insert((*head, tail_node)), "This gained pair already cached into the unary tracker!")}
                // helper closure for bond resilience updates
                let mut resilience_iterate_typed = |edge_iter: &BTreeSet<std::rc::Rc<std::cell::RefCell<EdgeEnds>>>| -> BTreeSet<Rc<RefCell<EdgeEnds>>> {
                    let mut bonds_decyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = BTreeSet::new();
                    for boxed_edge in edge_iter {
                        let edge = boxed_edge.borrow();
                        // check bond's weight, which is a boolean, marking if the bond is a known cycle-member
                        if *self.universe_graph.edge_weight(edge.z.unwrap()).unwrap() {
                            let mut counterfactual_graph = self.universe_graph.clone();
                            counterfactual_graph.remove_edge(edge.z.unwrap());
                            let path = astar(&counterfactual_graph, edge.a, |finish| finish == edge.b, |_| 1, |_| 0);
                            match path {
                                Some(_) => (),
                                None => {self.universe_graph.update_edge(edge.a, edge.b, false); assert!(bonds_decyclized.insert(Rc::clone(&boxed_edge)), "This edge already seen; can't insert it into the set!")}
                            }
                        }
                    }
                    bonds_decyclized
                };
                // update resilience of that species' edges
                let mut xh_xt_decyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.xh_xt);
                let mut p1_xp_decyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.p1_xp);
                let mut p2_xp_decyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.p2_xp);
                let mut p3_xp_decyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.p3_xp);
                let mut pp_pp_decyclized: BTreeSet<Rc<RefCell<EdgeEnds>>> = resilience_iterate_typed(&target_species.edges.pp_pp);
                let cycle_edges_retained_xh_xt = &self.cycle_edges.xh_xt - &xh_xt_decyclized;
                let cycle_edges_retained_p1_xp = &self.cycle_edges.p1_xp - &p1_xp_decyclized;
                let cycle_edges_retained_p2_xp = &self.cycle_edges.p2_xp - &p2_xp_decyclized;
                let cycle_edges_retained_p3_xp = &self.cycle_edges.p3_xp - &p3_xp_decyclized;
                let cycle_edges_retained_pp_pp = &self.cycle_edges.pp_pp - &pp_pp_decyclized;
                self.tree_edges.xh_xt.append(&mut xh_xt_decyclized);
                self.tree_edges.p1_xp.append(&mut p1_xp_decyclized);
                self.tree_edges.p2_xp.append(&mut p2_xp_decyclized);
                self.tree_edges.p3_xp.append(&mut p3_xp_decyclized);
                self.tree_edges.pp_pp.append(&mut pp_pp_decyclized);
                self.cycle_edges.xh_xt = cycle_edges_retained_xh_xt;
                self.cycle_edges.p1_xp = cycle_edges_retained_p1_xp;
                self.cycle_edges.p2_xp = cycle_edges_retained_p2_xp;
                self.cycle_edges.p3_xp = cycle_edges_retained_p3_xp;
                self.cycle_edges.pp_pp = cycle_edges_retained_pp_pp;
            } else {panic!("This edge did not have a bond index!")}
            self.update_mass_actions();
        }

        /**
         * Break a bond that is not a cycle member; this partitions one species into two.
        */
        pub fn axn_axn_binary_unbind(&mut self, target_edge: Rc<RefCell<EdgeEnds>>) {
            if let EdgeEnds{a: head_node, b: tail_node, z: Some(edge_index)} = *target_edge.borrow() {
                // sanity check for non-reslient bond
                assert!(! self.universe_graph.edge_weight(edge_index).unwrap(), "This bond ({}) is flagged as resilient; breaking it would yield a still-connected component!", target_edge.borrow());
                let old_last_edge_index = self.universe_graph.edge_indices().last().unwrap();
                self.universe_graph.remove_edge(edge_index).unwrap();
                let swapped_edge = self.edge_index_map.remove(&old_last_edge_index).unwrap();
                if old_last_edge_index != edge_index {
                    // remove edge & update caches; keep edge index tracker up to date
                    swapped_edge.borrow_mut().z = Some(edge_index);
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
                self.tree_edges.xh_xt.remove(&target_edge);
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
                // binding opportunities that were unary, but will now be binary
                for free_head in self.species_annots.get(&retained_mark).unwrap().borrow().ports.xh_free.union(&self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xh_free).into_iter() {
                    for free_tail in self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xt_free.union(&self.species_annots.get(&retained_mark).unwrap().borrow().ports.xt_free).into_iter() {
                        match self.unary_binding_pairs.xh_xt.take(&(*free_head, *free_tail)) {
                            Some((a_head, a_tail)) => {assert!(self.binary_binding_pairs.xh_xt.insert((a_head, a_tail)), "This pair already cached as binary!")},
                            None => ()
                        }
                    }
                }
                for free_head in self.species_annots.get(&retained_mark).unwrap().borrow().ports.p1_free.union(&self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p1_free).into_iter() {
                    for free_tail in self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.union(&self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free).into_iter() {
                        match self.unary_binding_pairs.p1_xp.take(&(*free_head, *free_tail)) {
                            Some((a_head, a_tail)) => {assert!(self.binary_binding_pairs.p1_xp.insert((a_head, a_tail)), "This pair already cached as binary!")},
                            None => ()
                        }
                    }
                }
                for free_head in self.species_annots.get(&retained_mark).unwrap().borrow().ports.p2_free.union(&self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p2_free).into_iter() {
                    for free_tail in self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.union(&self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free).into_iter() {
                        match self.unary_binding_pairs.p2_xp.take(&(*free_head, *free_tail)) {
                            Some((a_head, a_tail)) => {assert!(self.binary_binding_pairs.p2_xp.insert((a_head, a_tail)), "This pair already cached as binary!")},
                            None => ()
                        }
                    }
                }
                for free_head in self.species_annots.get(&retained_mark).unwrap().borrow().ports.p3_free.union(&self.species_annots.get(&ejected_mark).unwrap().borrow().ports.p3_free).into_iter() {
                    for free_tail in self.species_annots.get(&ejected_mark).unwrap().borrow().ports.xp_free.union(&self.species_annots.get(&retained_mark).unwrap().borrow().ports.xp_free).into_iter() {
                        match self.unary_binding_pairs.p3_xp.take(&(*free_head, *free_tail)) {
                            Some((a_head, a_tail)) => {assert!(self.binary_binding_pairs.p3_xp.insert((a_head, a_tail)), "This pair already cached as binary!")},
                            None => ()
                        }
                    }
                }
                for free_head in &self.species_annots.get(&retained_mark).unwrap().borrow().ports.pp_free {
                    for free_tail in &self.species_annots.get(&ejected_mark).unwrap().borrow().ports.pp_free {
                        match self.unary_binding_pairs.pp_pp.take(&(*free_head, *free_tail)) {
                            Some((a_head, a_tail)) => {assert!(self.binary_binding_pairs.pp_pp.insert((a_head, a_tail)), "This pair already cached as binary!")},
                            None => ()
                        }
                    }
                }
                assert!(self.ports.xh_free.insert(head_node), "This Axn head was already listed as free prior to this unbinding!");
                assert!(self.ports.xt_free.insert(tail_node), "This Axn tail was already listed as free prior to this unbinding!");
                assert!(self.species_annots.get(&head_node).unwrap().borrow_mut().ports.xh_free.insert(head_node), "This Axn head was already listed as free prior to this unbinding!");
                assert!(self.species_annots.get(&tail_node).unwrap().borrow_mut().ports.xt_free.insert(tail_node), "This Axn tail was already listed as free prior to this unbinding!");
                // bindings now possible due to the free'd sites
                for tail in &self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free {
                    assert!(self.unary_binding_pairs.xh_xt.insert((head_node, *tail)), "This gained pair already cached into the unary tracker")}
                for head in &self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free {
                    assert!(self.unary_binding_pairs.xh_xt.insert((*head, tail_node)), "This gained pair already cached into the unary tracker")}
                for tail in &self.ports.xt_free.difference(&self.species_annots.get(&head_node).unwrap().borrow().ports.xt_free).cloned().collect::<Vec<NodeIndex>>() {
                    assert!(self.binary_binding_pairs.xh_xt.insert((head_node, *tail)), "This gained pair already cached into the binary tracker")}
                self.binary_binding_pairs.xh_xt.remove(&(head_node, tail_node));    // hack to remove pair that will be visited from the other side in the next loop
                for head in &self.ports.xh_free.difference(&self.species_annots.get(&tail_node).unwrap().borrow().ports.xh_free).cloned().collect::<Vec<NodeIndex>>() {
                    assert!(self.binary_binding_pairs.xh_xt.insert((*head, tail_node)), "This gained pair already cached into the binary tracker")}
                // self-edge cleanup
                self.binary_binding_pairs.xh_xt.remove(&(head_node, head_node));
                self.binary_binding_pairs.xh_xt.remove(&(tail_node, tail_node));
                self.unary_binding_pairs.xh_xt.remove(&(head_node, head_node));
                self.unary_binding_pairs.xh_xt.remove(&(tail_node, tail_node));
            } else {panic!("This edge did not have a bond index!")}
            self.update_mass_actions();
        }

        /** Get a new copy of the mixture's rule activities. */
        pub fn current_rule_activities(&self) -> RuleActivities {
            self.rule_activities.clone()
        }

        /** Pretty print the mixture's rule activities. */
        pub fn print_rule_activities(&self) {
            println!("Rule activities:\n{}", self.rule_activities)
        }

        /** 
         * Get a new copy of the mixture's unary binding pairs. These are bonds that could occur,
         * but would not change the number of species in the mixture, occuring between agents
         * already attached to some species.
        */
        pub fn current_unary_binding_pairs(&self) -> PossibleBondEmbeds {
            self.unary_binding_pairs.clone()
        }

        /** 
         * Pretty print the mixture's unary binding pairs. These are bonds that could occur, but
         * would not change the number of species in the mixture, occuring between agents already
         * attached to some species.
        */
        pub fn print_unary_binding_pairs(&self) {
            println!("Unary binding pairs:\n{}", self.unary_binding_pairs)
        }

        /** 
         * Get a new copy of the mixture's binary binding pairs. These are bonds that could occur,
         * and would reduce by one the number of species in the mixture, fusing two species into
         * one.
        */
        pub fn current_binary_binding_pairs(&self) -> PossibleBondEmbeds {
            self.binary_binding_pairs.clone()
        }

        /** 
         * Pretty print the mixture's binary binding pairs. These are bonds that could occur, and
         * would reduce by one the number of species in the mixture, fusing two species into one.
        */
        pub fn print_binary_binding_pairs(&self) {
            println!("Binary binding pairs:\n{}", self.binary_binding_pairs)
        }

        /** 
         * Get a new copy of the mixture's bonds, who are cycle members. These are bonds that
         * exist, and whose breaking would not lead to a complex fragmenting into two complexes,
         * but just to a less inter-connected complex.
         */
        pub fn current_cycle_edges(&self) -> EdgeTypes {
            self.cycle_edges.clone()
        }

        /** 
         * Pretty print the mixture's bonds, who are cycle members. These are bonds that exist,
         * and whose breaking would not lead to a complex fragmenting into two complexes, but
         * just to a less inter-connected complex.
         */
        pub fn print_cycle_edges(&self) {
            println!("Edges in cycles:\n{}", self.cycle_edges)
        }

        /** 
         * Get a new copy of the mixture's bonds, who are not cycle members. These are bonds that
         * exist, and whose breaking would lead to a complex fragmenting into two complexes.
        */
        pub fn current_tree_edges(&self) -> EdgeTypes {
            self.tree_edges.clone()
        }

        /** 
         * Pretty print the mixture's bonds, who are not cycle members. These are bonds that
         * exist, and whose breaking would lead to a complex fragmenting into two complexes.
        */
        pub fn print_tree_edges(&self) {
            println!("Edges in trees:\n{}", self.tree_edges)
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{rule_activities::RuleRates, reaction_mixture::Mixture, edge_ends::EdgeEnds};
    use petgraph::prelude::{NodeIndex, EdgeIndex};
    use std::{rc::Rc, cell::RefCell};

    #[test]
    fn simple_reversability() {
        let my_rates = RuleRates {
            axn_axn_u_bind: 1.0,
            axn_axn_b_bind: 1.0,
            axn_axn_u_free: 1.0,
            axn_axn_b_free: 1.0,
            ap1_axn_u_bind: 1.0,
            ap1_axn_b_bind: 1.0,
            ap1_axn_u_free: 1.0,
            ap1_axn_b_free: 1.0,
            ap2_axn_u_bind: 1.0,
            ap2_axn_b_bind: 1.0,
            ap2_axn_u_free: 1.0,
            ap2_axn_b_free: 1.0,
            ap3_axn_u_bind: 1.0,
            ap3_axn_b_bind: 1.0,
            ap3_axn_u_free: 1.0,
            ap3_axn_b_free: 1.0,
            apc_apc_u_bind: 1.0,
            apc_apc_b_bind: 1.0,
            apc_apc_u_free: 1.0,
            apc_apc_b_free: 1.0,
        };
        let mut my_mix = Mixture::new_from_monomers(5, 5, my_rates);
        let u_b_0 = my_mix.current_unary_binding_pairs();
        let b_b_0 = my_mix.current_binary_binding_pairs();
        let u_f_0 = my_mix.current_cycle_edges();
        let b_f_0 = my_mix.current_tree_edges();
        
        my_mix.axn_axn_binary_bind(NodeIndex::from(0), NodeIndex::from(3));
        let u_b_1 = my_mix.current_unary_binding_pairs();
        let b_b_1 = my_mix.current_binary_binding_pairs();
        let u_f_1 = my_mix.current_cycle_edges();
        let b_f_1 = my_mix.current_tree_edges();

        my_mix.axn_axn_binary_bind(NodeIndex::from(2), NodeIndex::from(0));
        let u_b_2 = my_mix.current_unary_binding_pairs();
        let b_b_2 = my_mix.current_binary_binding_pairs();
        let u_f_2 = my_mix.current_cycle_edges();
        let b_f_2 = my_mix.current_tree_edges();
        
        my_mix.axn_axn_binary_bind(NodeIndex::from(1), NodeIndex::from(2));
        
        my_mix.axn_axn_binary_unbind(Rc::new(RefCell::new(EdgeEnds{
            a: NodeIndex::new(1), b: NodeIndex::new(2), z: Some(EdgeIndex::new(2))})));
        let u_b_4 = my_mix.current_unary_binding_pairs();
        let b_b_4 = my_mix.current_binary_binding_pairs();
        let u_f_4 = my_mix.current_cycle_edges();
        let b_f_4 = my_mix.current_tree_edges();

        assert_eq!(u_b_2, u_b_4);
        assert_eq!(b_b_2, b_b_4);
        assert_eq!(u_f_2, u_f_4);
        assert_eq!(b_f_2, b_f_4);
        
        my_mix.axn_axn_binary_unbind(Rc::new(RefCell::new(EdgeEnds{
            a: NodeIndex::new(2), b: NodeIndex::new(0), z: Some(EdgeIndex::new(1))})));
        let u_b_5 = my_mix.current_unary_binding_pairs();
        let b_b_5 = my_mix.current_binary_binding_pairs();
        let u_f_5 = my_mix.current_cycle_edges();
        let b_f_5 = my_mix.current_tree_edges();
        
        assert_eq!(u_b_1, u_b_5);
        assert_eq!(b_b_1, b_b_5);
        assert_eq!(u_f_1, u_f_5);
        assert_eq!(b_f_1, b_f_5);
    
        my_mix.axn_axn_binary_unbind(Rc::new(RefCell::new(EdgeEnds{
            a: NodeIndex::new(0), b: NodeIndex::new(3), z: Some(EdgeIndex::new(0))})));
        let u_b_6 = my_mix.current_unary_binding_pairs();
        let b_b_6 = my_mix.current_binary_binding_pairs();
        let u_f_6 = my_mix.current_cycle_edges();
        let b_f_6 = my_mix.current_tree_edges();

        assert_eq!(u_b_0, u_b_6);
        assert_eq!(b_b_0, b_b_6);
        assert_eq!(u_f_0, u_f_6);
        assert_eq!(b_f_0, b_f_6);
    }

    #[test]
    #[should_panic(expected = "This head and tail nodes (3, 0) already belong to the same species; can't binary bind them!")]
    fn unary_binding_embed_panics_binary_binding() {
        let my_rates = RuleRates {
            axn_axn_u_bind: 1.0,
            axn_axn_b_bind: 1.0,
            axn_axn_u_free: 1.0,
            axn_axn_b_free: 1.0,
            ap1_axn_u_bind: 0.0,
            ap1_axn_b_bind: 0.0,
            ap1_axn_u_free: 0.0,
            ap1_axn_b_free: 0.0,
            ap2_axn_u_bind: 0.0,
            ap2_axn_b_bind: 0.0,
            ap2_axn_u_free: 0.0,
            ap2_axn_b_free: 0.0,
            ap3_axn_u_bind: 0.0,
            ap3_axn_b_bind: 0.0,
            ap3_axn_u_free: 0.0,
            ap3_axn_b_free: 0.0,
            apc_apc_u_bind: 0.0,
            apc_apc_b_bind: 0.0,
            apc_apc_u_free: 0.0,
            apc_apc_b_free: 0.0,
        };
        let mut my_mix = Mixture::new_from_monomers(5, 5, my_rates);
        my_mix.axn_axn_binary_bind(NodeIndex::from(0), NodeIndex::from(1));
        my_mix.axn_axn_binary_bind(NodeIndex::from(1), NodeIndex::from(2));
        my_mix.axn_axn_binary_bind(NodeIndex::from(2), NodeIndex::from(3));
        my_mix.axn_axn_binary_bind(NodeIndex::from(3), NodeIndex::from(0));
    }

    #[test]
    #[should_panic(expected = "This head and tail nodes (3, 4) do not already belong to the same species; can't unary bind them!")]
    fn binary_binding_embed_panics_unary_binding() {
        let my_rates = RuleRates {
            axn_axn_u_bind: 1.0,
            axn_axn_b_bind: 1.0,
            axn_axn_u_free: 1.0,
            axn_axn_b_free: 1.0,
            ap1_axn_u_bind: 0.0,
            ap1_axn_b_bind: 0.0,
            ap1_axn_u_free: 0.0,
            ap1_axn_b_free: 0.0,
            ap2_axn_u_bind: 0.0,
            ap2_axn_b_bind: 0.0,
            ap2_axn_u_free: 0.0,
            ap2_axn_b_free: 0.0,
            ap3_axn_u_bind: 0.0,
            ap3_axn_b_bind: 0.0,
            ap3_axn_u_free: 0.0,
            ap3_axn_b_free: 0.0,
            apc_apc_u_bind: 0.0,
            apc_apc_b_bind: 0.0,
            apc_apc_u_free: 0.0,
            apc_apc_b_free: 0.0,
        };
        let mut my_mix = Mixture::new_from_monomers(6, 1, my_rates);
        my_mix.axn_axn_binary_bind(NodeIndex::from(0), NodeIndex::from(1));
        my_mix.axn_axn_binary_bind(NodeIndex::from(1), NodeIndex::from(2));
        my_mix.axn_axn_binary_bind(NodeIndex::from(2), NodeIndex::from(3));
        my_mix.axn_axn_unary_bind(NodeIndex::from(3), NodeIndex::from(4));
    }

    #[test]
    #[should_panic(expected = "This bond (#1 @ 1 -> 2) is flagged as resilient; breaking it would yield a still-connected component!")]
    fn unary_unbinding_embed_panics_binary_unbinding() {
        let my_rates = RuleRates {
            axn_axn_u_bind: 1.0,
            axn_axn_b_bind: 1.0,
            axn_axn_u_free: 1.0,
            axn_axn_b_free: 1.0,
            ap1_axn_u_bind: 0.0,
            ap1_axn_b_bind: 0.0,
            ap1_axn_u_free: 0.0,
            ap1_axn_b_free: 0.0,
            ap2_axn_u_bind: 0.0,
            ap2_axn_b_bind: 0.0,
            ap2_axn_u_free: 0.0,
            ap2_axn_b_free: 0.0,
            ap3_axn_u_bind: 0.0,
            ap3_axn_b_bind: 0.0,
            ap3_axn_u_free: 0.0,
            ap3_axn_b_free: 0.0,
            apc_apc_u_bind: 0.0,
            apc_apc_b_bind: 0.0,
            apc_apc_u_free: 0.0,
            apc_apc_b_free: 0.0,
        };
        let mut my_mix = Mixture::new_from_monomers(5, 5, my_rates);
        my_mix.axn_axn_binary_bind(NodeIndex::from(0), NodeIndex::from(1));
        my_mix.axn_axn_binary_bind(NodeIndex::from(1), NodeIndex::from(2));
        my_mix.axn_axn_binary_bind(NodeIndex::from(2), NodeIndex::from(3));
        my_mix.axn_axn_unary_bind(NodeIndex::from(3), NodeIndex::from(0));
        my_mix.axn_axn_binary_unbind(Rc::new(RefCell::new(EdgeEnds{a: NodeIndex::new(1), b: NodeIndex::new(2), z: Some(EdgeIndex::new(1))})))
    }

    #[test]
    #[should_panic(expected = "This bond (#2 @ 2 -> 3) is not flagged as resilient; breaking it would not yield a still-connected component!")]
    fn binary_unbinding_embed_panics_unary_unbinding() {
        let my_rates = RuleRates {
            axn_axn_u_bind: 1.0,
            axn_axn_b_bind: 1.0,
            axn_axn_u_free: 1.0,
            axn_axn_b_free: 1.0,
            ap1_axn_u_bind: 0.0,
            ap1_axn_b_bind: 0.0,
            ap1_axn_u_free: 0.0,
            ap1_axn_b_free: 0.0,
            ap2_axn_u_bind: 0.0,
            ap2_axn_b_bind: 0.0,
            ap2_axn_u_free: 0.0,
            ap2_axn_b_free: 0.0,
            ap3_axn_u_bind: 0.0,
            ap3_axn_b_bind: 0.0,
            ap3_axn_u_free: 0.0,
            ap3_axn_b_free: 0.0,
            apc_apc_u_bind: 0.0,
            apc_apc_b_bind: 0.0,
            apc_apc_u_free: 0.0,
            apc_apc_b_free: 0.0,
        };
        let mut my_mix = Mixture::new_from_monomers(6, 5, my_rates);
        my_mix.axn_axn_binary_bind(NodeIndex::from(0), NodeIndex::from(1));
        my_mix.axn_axn_binary_bind(NodeIndex::from(1), NodeIndex::from(2));
        my_mix.axn_axn_binary_bind(NodeIndex::from(2), NodeIndex::from(3));
        my_mix.axn_axn_binary_bind(NodeIndex::from(3), NodeIndex::from(4));
        my_mix.axn_axn_binary_bind(NodeIndex::from(4), NodeIndex::from(5));
        my_mix.axn_axn_unary_unbind(Rc::new(RefCell::new(EdgeEnds{a: NodeIndex::new(2), b: NodeIndex::new(3), z: Some(EdgeIndex::new(2))})))
    }
}