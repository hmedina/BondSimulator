#![feature(hash_drain_filter, btree_drain_filter)]
#![warn(clippy::all)]

/**
 * Module for storing stuff to build the simulator; rules & abundances & rates & bond-types...
*/
pub mod building_blocks {
    use std::collections::BTreeMap;

    /**
     * A resource vector, using a map, where the key is an agent name, and the value is a protomer
     * abundance.
    */
    pub struct ProtomerResources <'a> {
        pub map: BTreeMap<&'a str, usize>
    }

    /**
     * Structure to group & name the four rate constants of a binding interaction, the agents &
     * sites involved in that interaction, and a user-memorable name for the binding interaction.
    */
    pub struct InteractionData <'a> {
        pub name: Option<&'a str>,
        pub bind_unary: RateFudge,
        pub bind_binary: RateFudge,
        pub free_unary: RateFudge,
        pub free_binary: RateFudge,
        pub agent_a: &'a str,
        pub agent_b: &'a str,
        pub site_a: &'a str,
        pub site_b: &'a str
    }

    /**
     * A base rate constant, plus some fudge factors.
    */
    pub struct RateFudge {
        pub base: f64,
        pub multiplicative: Option<f64>,
        pub additive: Option<f64>
    }
}


/**
 * Module for primitive types & objects.
*/
pub mod primitives {
    use petgraph::prelude::{EdgeIndex, NodeIndex};
    use std::{collections::BTreeMap, fmt};


    /**
     * Represents an agent in the reaction mixture.
     * ```
     * use axin_apc_simulator::primitives::Agent;
     * let foo = Agent::new_from_signature("Bob", vec!["s1", "site_3", "s2"]);
     * assert_eq!("Bob(s1[.], s2[.], site_3[.])", format!("{}", foo));
     * ```
     * Supports agent identifiers, and bond state.
     * ```
     * use axin_apc_simulator::primitives::Agent;
     * use petgraph::prelude::NodeIndex;
     * use std::collections::BTreeMap;
     * let bar = Agent{name: "José", id: Some(NodeIndex::new(5)), sites: BTreeMap::new(("x", None), ("ÿ", EdgeIndex::new(3)), ("z", None))};
     * assert_eq!("x5:José(x[.], ÿ[3], z[.])", format!("{}", bar));
     * ```
    */
    #[derive(Clone, Debug)]
    pub struct Agent <'a> {
        pub name: &'a str,
        pub id: Option<NodeIndex>,
        pub sites: BTreeMap<&'a str, Option<EdgeIndex>>
    }

    impl <'a> Agent <'a> {
        pub fn new_from_signature(agent_name: &'a str, site_names: Vec<&'a str>) -> Agent<'a> {
            Agent {
                name: agent_name,
                id: None,
                sites: BTreeMap::from_iter(site_names.iter().map(|n| (*n, None)))
            }
        }
    }

    impl fmt::Display for Agent <'_> {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            let mut str_store: Vec<String> = Vec::new();
            for (s_name, b_state) in self.sites.iter() {
                str_store.push( match b_state {
                    Some(b) => format!("{}[{}]", s_name, b.index()),
                    None => format!("{}[.]", s_name)
                } )
            }
            let prefix: String = match self.id {
                Some(id) => format!("x{}:", id.index()),
                None => String::new()
            };
            write!(f, "{}{}({})", prefix, self.name, str_store.join(", "))
        }
    }


    /**
     * Represents the landing site of a bond. Since site names are unique within the agent 
     * namespace, and agent names are unique within the global namespace, this pair is unique
     * globally.
    */
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct AgentSite <'a> {
        pub agent: &'a str,
        pub site: &'a str
    }

    impl fmt::Display for AgentSite <'_> {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "{}.{}", self.agent, self.site)
        }
    }


    /**
     * Represents a bond type: a pair of agents and their sites. Constructor orders alphabetically
     * by the agent names, then by site names. This order is observed for comparison and sorting.
     * ```
     * use axin_apc_simulator::primitives::BondType;
     * let foo = BondType::new("Mango", "stone", "Earth", "ground");
     * assert_eq!("Earth(ground[1]), Mango(stone[1])", format!("{}", foo));
     * ```
    */
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct BondType <'a> {
        pub pair_1: AgentSite <'a>,
        pub pair_2: AgentSite <'a>,
        pub arity: InteractionArity,
        pub direction: InteractionDirection
    }

    impl fmt::Display for BondType <'_> {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "{}({}[1]), {}({}[1])", self.pair_1.agent, self.pair_1.site, self.pair_2.agent, self.pair_2.site)
        }
    }

    impl <'a> BondType <'a> {
        pub fn new(agent_a: &'a str, site_a: &'a str, agent_b: &'a str, site_b: &'a str, arity: InteractionArity, direction: InteractionDirection) -> Self {
            if agent_a < agent_b {
                BondType{pair_1: AgentSite{agent: agent_a, site: site_a}, pair_2: AgentSite{agent: agent_b, site: site_b}, arity, direction}
            } else if agent_b < agent_a {
                BondType{pair_2: AgentSite{agent: agent_a, site: site_a}, pair_1: AgentSite{agent: agent_b, site: site_b}, arity, direction}
            } else if site_a < site_b {
                BondType{pair_1: AgentSite{agent: agent_a, site: site_a}, pair_2: AgentSite{agent: agent_b, site: site_b}, arity, direction}
            } else {
                BondType{pair_2: AgentSite{agent: agent_a, site: site_a}, pair_1: AgentSite{agent: agent_b, site: site_b}, arity, direction}
            }
        }
    }

    
    /** 
     * Structure for representing a bond embed, either current or possible, along with the 
     * corresponding edge index, if any. The agent indexes correspond to the agent identifiers.
    */
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct BondEmbed <'a> {
        pub a_index: NodeIndex,
        pub b_index: NodeIndex,
        pub z: Option<EdgeIndex>,
        pub bond_type: BondType <'a>
    }
    
    impl <'a> fmt::Display for BondEmbed <'a> {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            let str_rep = match self.z {
                Some(i) => format!("x{}:{}({}[{}]) <-> x{}:{}({}[{}])", self.a_index.index(), self.bond_type.pair_1.agent, self.bond_type.pair_1.site, i.index(), self.b_index.index(), self.bond_type.pair_2.agent, self.bond_type.pair_2.site, i.index()),
                None => format!("x{}:{}({}[.]) <-> x{}:{}({}[.])", self.a_index.index(), self.bond_type.pair_1.agent, self.bond_type.pair_1.site, self.b_index.index(), self.bond_type.pair_2.agent, self.bond_type.pair_2.site)
            };
            write!(f, "{}", str_rep)
        }
    }


    /**
     * Marks if an interaction occurs within a complex, or between complexes. Used for type'ing
     * transformations, embed trackers and their constructors.
    */
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub enum InteractionArity {
        Unary,
        Binary,
    }

    impl fmt::Display for InteractionArity {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                InteractionArity::Unary => write!(f, "Unary"),
                InteractionArity::Binary => write!(f, "Binary"),
            }
        }
    }


    /**
     * Marks if an interaction creates a bond or frees a pair of sites for binding. Used for 
     * type'ing transformations, embed trackers and their constructors.
    */
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub enum InteractionDirection {
        Bind,
        Free
    }

    impl fmt::Display for InteractionDirection {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            match self {
                InteractionDirection::Bind => write!(f, "Bind"),
                InteractionDirection::Free => write!(f, "Free"),
            }
        }
    }
}

/**
 * Module for various collectors of primitives.
*/
mod collectors {
    use itertools::Itertools;
    use petgraph::prelude::{NodeIndex};
    use rand::{distributions::WeightedIndex, prelude::*};
    use std::{cell::RefCell, collections::{HashSet, BTreeSet, BTreeMap}, cmp::Ordering, fmt, rc::Rc};
    use crate::building_blocks::{InteractionData, ProtomerResources, RateFudge};
    use crate::primitives::{AgentSite, BondEmbed, BondType, InteractionArity, InteractionDirection};
    

    /**
     * Annotation for a specific species in the mixture.
    */
    #[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct MixtureSpecies <'a> {
        pub size: usize,
        pub agent_set: BTreeSet<NodeIndex>,
        pub edges: BTreeSet<Rc<RefCell<BondEmbed<'a>>>>,
        pub ports: OpenPorts <'a>,
    }


    /** 
     * Map of HashSets, one per site type. The HashSets hold NodeIndexes, from the UniverseGraph.
     * The key of the map is the site name.
     */
    #[derive(Debug, PartialEq, Eq)]
    pub struct OpenPorts <'a> {
        pub map: BTreeMap<AgentSite<'a>, HashSet<NodeIndex>>
    }

    impl <'a> fmt::Display for OpenPorts <'a> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            let mut str_store: Vec<String> = Vec::new();
            for (agent_site, agent_indexes) in self.map.iter() {
                str_store.push(format!("\t{}\n{}", agent_site, agent_indexes.iter().map(|x| x.index().to_string()).collect::<Vec<String>>().join(", ")));
            }
            write!(f, "{}", str_store.join("\n"))
        }
    }

    impl <'a> Ord for OpenPorts <'a> {
        fn cmp(&self, other: &Self) -> Ordering {
            let mut zipper = self.map.keys().zip(other.map.keys());
            let eval: Ordering = loop {
                match zipper.next() {
                    Some((self_key, other_key)) => {
                        if self_key != other_key {break self_key.cmp(other_key)}
                        else {
                            let self_set: Vec<&NodeIndex> = self.map.get(self_key).unwrap().iter().collect();
                            let other_set: Vec<&NodeIndex> = self.map.get(other_key).unwrap().iter().collect();
                            if self_set != other_set {break self_set.cmp(&other_set)}
                        }
                    },
                    None => {
                        if self.map.len() != other.map.len() {break self.map.len().cmp(&other.map.len())}   // if they were equal up to the point one returned empty, compare length of key-arrays
                        else {break Ordering::Equal}                                                        // if they ran out at the same time, with equal keys so far, and equal arrays, they are equal
                    }
                };
            };
            eval
        }
    }

    impl <'a> PartialOrd for OpenPorts <'a> {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }
    
    impl <'a> OpenPorts <'a> {
        /**
         * Update every site-specific HashSet with the elements from the other, consuming it.
        */
        pub fn update_from(&mut self, other: &Self) {
            for (agent_site, site_data) in other.map.iter() {
                self.map.entry(*agent_site).or_insert_with(HashSet::new).extend(site_data)
            }
        }

        /**
         * Partition the site-specific HashSets using a set of NodeIndexes, the "ejected" nodes.
         * A free site belonging to one of the agents in the ejected set, gets ejected and returned
         * in a new OpenPort object. Sites not belonging to agents in the ejected set are retained.
        */
        pub fn eject_where(&mut self, ejected_nodes: &BTreeSet<NodeIndex>) -> Self {
            let mut ejected_ports: BTreeMap<AgentSite<'a>, HashSet<NodeIndex>> = BTreeMap::new();
            for (agent_site, site_data) in self.map.iter_mut() {
                let sink: HashSet<NodeIndex> = site_data.drain_filter(|n| ejected_nodes.contains(n)).collect();
                ejected_ports.insert(*agent_site, sink);
            }
            OpenPorts{map: ejected_ports}
        }

        /**
         * Constructor for the smallest species.
        */
        pub fn new_monomer(agent_sites: Vec<AgentSite<'a>>) -> Self {
            OpenPorts {
                map: BTreeMap::from_iter(agent_sites.iter().map(|n| (*n, HashSet::new())))
            }
        }

        /**
         * Constructor for the universe tracker.
        */
        pub fn new_universe(agent_sites: Vec<AgentSite<'a>>, site_abundances: Vec<usize>) -> Self {
            let mut map: BTreeMap<AgentSite<'a>, HashSet<NodeIndex>> = BTreeMap::new();
            for (agent_site, site_abundance) in agent_sites.iter().zip(site_abundances) {
                map.insert(*agent_site, HashSet::with_capacity(site_abundance));
            }
            OpenPorts {map}
        }
    }

    


    /**
     * Convenience structure for calculating the activity of interactions, while keeping the bias
     * metrics & calculation separate.
    */
    #[derive(Debug, PartialEq, PartialOrd)]
    pub struct MassActionParameters {
        /// The abundance of some reactant
        mass: usize,

        /// The rate constant of the interaction
        rate: f64,

        /// If given, includes a factor that multiplies by the complex size(s)
        bias_m: Option<f64>,

        /// If given, includes a summand that adds a term with a product of the complex size(s)
        bias_a: Option<f64>
    }

    impl MassActionParameters {
        /**
         * Activity is calculated as:
         * 
         * `α = rate * mass [* bias_m * a_size * b_size] [+ bias_a * a_size * b_size]`
         * - `bias_m` is the optional multiplicative bias for the complex size,
         * - `bias_a` is the optional additive bias for the complex size,
         * - `a_size` is the size of the left agent's complex
         * - `b_size` is the size of the right agent's complex
         * 
         * With both biases omitted:
         * `α = rate * mass * 1 + 0`
        */
        pub fn calculate_activity(&self, a_size: usize, b_size: usize) -> f64 {
            let m: f64 = match self.bias_m {
                Some(f) => f * a_size as f64 * b_size as f64,
                None => 1.0
            };
            let a: f64 = match self.bias_a {
                Some(f) => f * a_size as f64 * b_size as f64,
                None => 0.0
            };
            self.mass as f64 * self.rate * m + a
        }
    }


    /** 
     * Collects the interacting pairs, either bound already, or eligible for a future binding
     * event.
     */
    #[derive(Debug, PartialEq, PartialOrd)]
    pub struct InteractingTracker <'a> {
        /// The set of interacting pairs.
        pub set: BTreeSet<Rc<RefCell<BondEmbed<'a>>>>,

        /// The parameters governing the activity of this interaction.
        activity_parameters: MassActionParameters,

        /// The activity of this interaction.
        pub activity_calculated: f64,
    }
    
    impl <'a> InteractingTracker <'a> {
        /**
         * Update set with the elements from the other.
        */
        pub fn update_from(&mut self, other: &mut Self, species_annots: &BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies<'a>>>> ) {
            self.set.append(&mut other.set);
            self.update_activity(species_annots)
        }

        /**
         * Partition the site-specific BTreeSets using a set of NodeIndexes, the "ejected" nodes.
         * A bond where both agents belong to the ejected set, gets ejected and returned
         * in a new object. Bonds where both agents do not belong to agents in the
         * ejected set are retained. If one agent is ejected but another retained, this function
         * panics, as it should not be breaking cycles as a side-effect.
        */
        pub fn eject_where(&mut self, ejected_nodes: &BTreeSet<NodeIndex>, species_annots: &BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies<'a>>>>) -> Self {
            let sink: BTreeSet<Rc<RefCell<BondEmbed>>> = self.set.drain_filter( |b| ejected_nodes.contains(&b.borrow().a_index) && ejected_nodes.contains(&b.borrow().b_index)).collect();
            self.update_activity(species_annots);
            let mut ejected_pairs = InteractingTracker{
                activity_calculated: 0.0,
                activity_parameters: MassActionParameters{
                    mass: sink.len(), 
                    rate: self.activity_parameters.rate, 
                    bias_m: self.activity_parameters.bias_m, 
                    bias_a: self.activity_parameters.bias_a},
                set: sink
            };
            ejected_pairs.update_activity(species_annots);
            ejected_pairs
        }

        /**
         * Constructor used when starting from a monomeric initial state. For interactions of 
         * `InteractionArity::Binary` and `InteractionDirection::Bind`, it will generate the 
         * appropriate combinatorics for the possible bond types. For the other classes, it will be
         * empty, but pre-allocated.
        */
        pub fn new_tracker(interaction: &InteractionData<'a>, resources: &ProtomerResources, arity: InteractionArity, direction: InteractionDirection) -> Self {
            let mut index_offset = 0;
            let bond_type = BondType::new(interaction.agent_a, interaction.site_a, interaction.agent_b, interaction.site_b, arity, direction);  // this BondType construction also sorts their names
            let a_agent = bond_type.pair_1.agent;
            let b_agent = bond_type.pair_2.agent;
            let a_site = bond_type.pair_1.site;
            let b_site = bond_type.pair_2.site;
            let mut set: BTreeSet<Rc<RefCell<BondEmbed>>> = BTreeSet::new();
            let (activity_parameters, activity_calculated) = match (arity, direction) {
                (InteractionArity::Binary, InteractionDirection::Bind) => {
                    let a_indexes: Vec<NodeIndex> = (index_offset..(index_offset + resources.map.get(a_agent).unwrap())).map(|x| NodeIndex::from(x as u32)).collect();
                    index_offset += resources.map.get(a_agent).unwrap();
                    if a_agent != b_agent {
                        // heterogenous dimerization, using cartesian product
                        let b_indexes: Vec<NodeIndex> = (index_offset..(index_offset + resources.map.get(b_agent).unwrap())).map(|x| NodeIndex::from(x as u32)).collect();
                        index_offset += resources.map.get(b_agent).unwrap();
                        for (a_ix, b_ix) in a_indexes.iter().cartesian_product(b_indexes.iter()) {
                            let comb = BondEmbed{a_index: *a_ix, b_index: *b_ix, z: None, bond_type};
                            assert!(set.insert(Rc::new(RefCell::new(comb))));
                        }
                    }
                    else if a_site != b_site {
                        // homogenous head-to-tail concatenation, using 2-permutations
                        for v_perms in a_indexes.iter().permutations(2) {
                            let comb = BondEmbed{a_index: *v_perms[0], b_index: *v_perms[1], z: None, bond_type};
                            assert!(set.insert(Rc::new(RefCell::new(comb))));
                        }
                    }
                    else {
                        // symmetric dimerization, using 2-combinations
                        for v_combs in a_indexes.iter().combinations(2) {
                            let comb = BondEmbed{a_index: *v_combs[0], b_index: *v_combs[1], z: None, bond_type};
                            assert!(set.insert(Rc::new(RefCell::new(comb))));
                        }
                    }
                    // calculate initial activity, accounting for fudge factors
                    let m: f64 = match interaction.bind_binary.multiplicative {
                        Some(f) => f * (set.len()^2) as f64,
                        None => 1.0
                    };
                    let a: f64 = match interaction.bind_binary.additive {
                        Some(f) => f * (set.len()^2) as f64,
                        None => 0.0
                    };
                    let act: f64 = set.len() as f64 * interaction.bind_binary.base * m + a;
                    (MassActionParameters{mass: set.len(), rate: interaction.bind_binary.base, bias_m: interaction.bind_binary.multiplicative, bias_a: interaction.bind_binary.additive}, act)
                },
                (InteractionArity::Unary, InteractionDirection::Bind) => {
                    (MassActionParameters{mass: 0, rate: interaction.bind_unary.base, bias_m: interaction.bind_unary.multiplicative, bias_a: interaction.bind_unary.additive}, 0.0)
                },
                (InteractionArity::Binary, InteractionDirection::Free) => {
                    (MassActionParameters{mass: 0, rate: interaction.free_binary.base, bias_m: interaction.free_binary.multiplicative, bias_a: interaction.free_binary.additive}, 0.0)
                },
                (InteractionArity::Unary, InteractionDirection::Free) => {
                    (MassActionParameters{mass: 0, rate: interaction.free_unary.base, bias_m: interaction.free_unary.multiplicative, bias_a: interaction.free_unary.additive}, 0.0)
                }
            };
            InteractingTracker {set, activity_parameters, activity_calculated}
        }

        /**
         * Helper method to build a vector of float values with the activity of each interacting
         * pair using their sizes and the activity parameters for the interaction type.
        */
        fn build_embed_weights(&self, species_annots: &BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies<'a>>>>) -> Vec<f64> {
            let helper_closure = |pair: &Rc<RefCell<BondEmbed>>| -> f64 {
                let borrowed_pair = pair.borrow();
                let index_a: NodeIndex = borrowed_pair.a_index;
                let index_b: NodeIndex = borrowed_pair.b_index;
                self.activity_parameters.calculate_activity(species_annots.get(&index_a).unwrap().borrow().size, species_annots.get(&index_b).unwrap().borrow().size)
            };
            let embed_weights: Vec<f64> = self.set.iter().map(helper_closure).collect::<Vec<f64>>();
            embed_weights
        }

        /**
         * Update self's activity, using optionally a function to bias the weight from each
         * binding pair.
        */
        pub fn update_activity(&mut self, species_annots: &BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies<'a>>>>) {
            self.activity_calculated = self.build_embed_weights(species_annots).iter().sum();
        }

        /**
         * Using complex sizes in the weighting, pick a pair of targets.
        */
        pub fn pick_targets(&self, rng: ThreadRng, species_annots: &BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies<'a>>>>) -> Rc<RefCell<BondEmbed<'a>>> {
            let dist = WeightedIndex::new(self.build_embed_weights(species_annots)).unwrap();
            let chosen_ix: usize = dist.sample(&mut rng);
            *self.set.iter().nth(chosen_ix).unwrap()
        }
    }

    impl <'a> fmt::Display for InteractingTracker <'a> {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            let bond_instances_closure = |input_set: &BTreeSet<Rc<RefCell<BondEmbed>>>| -> String {
                input_set.iter().map(|b| format!("{}", b.borrow())).collect::<Vec<String>>().join("\n\t")
            };
            write!(f, "{}", bond_instances_closure(&self.set))
        }
    }
}


pub mod reaction_mixture {
    use crate::primitives::{Agent, AgentSite, BondType, BondEmbed, InteractionArity, InteractionDirection};
    use crate::collectors::{OpenPorts, MixtureSpecies, InteractingTracker};
    use petgraph::{graph::Graph, algo::astar::astar, prelude::*, visit::Dfs};
    use rand::{distributions::WeightedIndex, prelude::*};
    use std::cell::{RefCell, RefMut, Ref};
    use std::collections::{VecDeque, BTreeMap, BTreeSet};
    use std::rc::Rc;




    pub struct Mixture <'a> {
        /// The graph that tracks connectivity; used for cycle detection & species partitioning.
        universe_graph: Graph<Agent<'a>, bool, Undirected>,
        
        /// Maps the species cache to each node in the graph.
        species_annots: BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies<'a>>>>,
        
        /// Used for iteration for pretty-printing as a Kappa snapshot.
        species_set: VecDeque<Rc<RefCell<MixtureSpecies<'a>>>>,

        /// Stores a global cache of open ports.
        ports: OpenPorts<'a>,
        
        /// Stores the interactions in the system.
        interactions: BTreeMap<&'a BondType<'a>, InteractingTracker<'a>>,

        /// Maps an EdgeIndex to a BondEmbed, for quick updating when deleting a bond invalidates
        /// the last EdgeIndex in the graph.
        edge_index_map: BTreeMap<EdgeIndex, Rc<RefCell<BondEmbed<'a>>>>,

        simulated_time: f64,
        simulated_events: usize,

        /// Cache for the random number generator.
        simulator_rng: ThreadRng
    }
    
    impl <'a> Mixture <'a> {

        /**
         * Advance simulation time and event number tracker.
        */
        fn advance_simulation_metrics(&mut self) {
            self.simulated_events += 1;
            let mut tau: f64 = 0.0;
            for (i_type, i_data) in self.interactions {
                tau += i_data.activity_calculated;
            }
            let rdy: f64 = (1.0 / self.simulator_rng.gen::<f64>()).ln();
            self.simulated_time += tau * rdy;
        }

        /**
         * Based on the activities of the interactions, choose one to realize.
        */
        pub fn choose_and_apply_next_rule(&mut self) {
            let dist = WeightedIndex::new(self.interactions.iter().map(|i| i.1.activity_calculated)).unwrap();
            let chosen_index: usize = dist.sample(&mut self.simulator_rng);
            let chosen_bond_type: &BondType = self.interactions.iter().nth(chosen_index).unwrap().0;
            let i = self.interactions.get_mut(chosen_bond_type).unwrap();
            let bond: Rc<RefCell<BondEmbed>> = i.pick_targets(self.simulator_rng, &self.species_annots);
            match (chosen_bond_type.arity, chosen_bond_type.direction) {
                (InteractionArity::Unary, InteractionDirection::Bind) => self.unary_bind(bond),
                (InteractionArity::Unary, InteractionDirection::Free) => self.unary_free(bond),
                (InteractionArity::Binary, InteractionDirection::Bind) => self.binary_bind(bond),
                (InteractionArity::Binary, InteractionDirection::Free) => self.binary_free(bond)
            };
            self.advance_simulation_metrics();
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
        pub fn unary_bind(&mut self, target_bond: Rc<RefCell<BondEmbed<'a>>>) {
            let head_index: NodeIndex = target_bond.borrow().a_index;
            let tail_index: NodeIndex = target_bond.borrow().b_index;
            let bond_type: BondType<'a> = target_bond.borrow_mut().bond_type;
            // sanity check for unary-ness
            assert_eq!(self.species_annots.get(&head_index), self.species_annots.get(&tail_index), "This head and tail nodes ({}, {}) do not already belong to the same species; can't unary bind them!", head_index.index(), tail_index.index());
            let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_index).unwrap().borrow_mut();
            // create edge & update caches
            assert!(self.interactions.get(&bond_type).unwrap().set.remove(&target_bond), "The target bond {} was not already in the expected tracker for {}!", target_bond.borrow(), bond_type);
            let new_edge_index = self.universe_graph.add_edge(head_index, tail_index, true);
            bond_type.direction = InteractionDirection::Free;
            target_bond.borrow_mut().z = Some(new_edge_index);
            let agent_a: &mut Agent = self.universe_graph.node_weight_mut(head_index).unwrap();
            let agent_b: &mut Agent = self.universe_graph.node_weight_mut(tail_index).unwrap();
            let site_a = agent_a.sites.get_mut(bond_type.pair_1.site).unwrap();
            match site_a {
                Some(_) => panic!("Node {} was already bound at {}!", agent_a, bond_type.pair_1.site),
                None => *site_a = Some(new_edge_index)
            }
            let site_b = agent_b.sites.get_mut(bond_type.pair_2.site).unwrap();
            match site_b {
                Some(_) => panic!("Node {} was already bound at {}!", agent_b, bond_type.pair_2.site),
                None => *site_b = Some(new_edge_index)
            }
            target_species.edges.insert(Rc::clone(&target_bond));
            self.edge_index_map.insert(new_edge_index, Rc::clone(&target_bond));
            self.interactions.get_mut(&bond_type).unwrap().set.insert(Rc::clone(&target_bond));
            // bindings no longer possible due to occupied sites
            let mut unary_combinatorics_lost: BTreeSet<BondEmbed> = BTreeSet::new();
            for tail in target_species.ports.map.get(&bond_type.pair_2).unwrap().iter().collect::<Vec<&NodeIndex>>() {
                let b = self.interactions.get_mut(&bond_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: head_index, b_index: *tail, z: None, bond_type}))).unwrap();
                assert!(unary_combinatorics_lost.insert(*b.borrow()), "This lost pair ({}, {}) already cached into the unary lost tracker!", head_index.index(), tail.index())
            }
            for head in target_species.ports.map.get(&bond_type.pair_1).unwrap().iter().collect::<Vec<&NodeIndex>>() {
                let b = self.interactions.get_mut(&bond_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: *head, b_index: tail_index, z: None, bond_type}))).unwrap();
                assert!(unary_combinatorics_lost.insert(*b.borrow()), "This lost pair ({}, {}) already cached into the unary lost tracker!", head.index(), tail_index.index())
            }
            let mut binary_combinatorics_lost: BTreeSet<BondEmbed> = BTreeSet::new();
            let purged_type = BondType{arity: InteractionArity::Binary, ..bond_type};
            for tail in &self.ports.map.get(&bond_type.pair_2).unwrap().difference(&target_species.ports.map.get(&bond_type.pair_2).unwrap()).cloned().collect::<Vec<NodeIndex>>() {
                let b = self.interactions.get_mut(&purged_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: head_index, b_index: *tail, z: None, bond_type: purged_type}))).unwrap();
                assert!(binary_combinatorics_lost.insert(*b.borrow()), "This lost pair ({}, {}) already cached into the binary lost tracker!", head_index.index(), tail.index())
            }
            for head in &self.ports.map.get(&bond_type.pair_1).unwrap().difference(&target_species.ports.map.get(&bond_type.pair_1).unwrap()).cloned().collect::<Vec<NodeIndex>>() {
                let b = self.interactions.get_mut(&purged_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: *head, b_index: tail_index, z: None, bond_type: purged_type}))).unwrap();
                assert!(binary_combinatorics_lost.insert(*b.borrow()), "This lost pair ({}, {}) already cached into the binary lost racker!", head.index(), tail_index.index())
            }
            assert!(self.ports.map.get_mut(&bond_type.pair_1).unwrap().remove(&head_index));
            assert!(self.ports.map.get_mut(&bond_type.pair_2).unwrap().remove(&tail_index));
            assert!(target_species.ports.map.get_mut(&bond_type.pair_1).unwrap().remove(&head_index));
            assert!(target_species.ports.map.get_mut(&bond_type.pair_2).unwrap().remove(&tail_index));
            // update resiliency
            // move things from self.interactions that match Free & Binary (aka tree edges) into
            // the corresponding
            // things in self.interactions that match Free & Unary (aka cycle edges)
            let mut newly_cyclized_bond = |boxed_bond: &Rc<RefCell<BondEmbed>>| -> bool {
                let bond_embed = boxed_bond.borrow();
                if ! self.universe_graph.edge_weight(bond_embed.z.unwrap()).unwrap() {  // check bond's weight, which is a boolean, marking if the bond is an already known cycle-member
                    let mut counterfactual_graph = self.universe_graph.clone();
                    counterfactual_graph.remove_edge(bond_embed.z.unwrap());
                    let path = astar(&counterfactual_graph, bond_embed.a_index, |finish| finish == bond_embed.b_index, |_| 1, |_| 0);
                    match path {
                        Some(_) => {self.universe_graph.update_edge(bond_embed.a_index, bond_embed.b_index, true); true}
                        None => (false)
                    }
                } else {false}      // if the bond's weight was already true, then it is not newly cyclized
            };
            for (i_type, i_data) in self.interactions {
                if i_type.arity == InteractionArity::Binary && i_type.direction == InteractionDirection::Free {
                    let mut newly_cyclized_bonds: BTreeSet<Rc<RefCell<BondEmbed>>> = i_data.set.drain_filter(newly_cyclized_bond).collect();
                    let target_container = BondType{arity: InteractionArity::Unary, ..*i_type};
                    self.interactions.get_mut(&target_container).unwrap().set.append(&mut newly_cyclized_bonds);
                }
            }
        }
    
        /**
         * Create a bond that joins two species into one.
        */
        pub fn binary_bind(&mut self, target_bond: Rc<RefCell<BondEmbed<'a>>>) {
            let head_index: NodeIndex = target_bond.borrow().a_index;
            let tail_index: NodeIndex = target_bond.borrow().b_index;
            let bond_type: BondType<'a> = target_bond.borrow_mut().bond_type;
            // sanity check for binary-ness
            assert_ne!(self.species_annots.get(&head_index), self.species_annots.get(&tail_index), "Nodes ({}, {}) already belong to the same species; can't binary bind them!", head_index.index(), tail_index.index());
            let head_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_index).unwrap().borrow_mut();
            let tail_species: RefMut<MixtureSpecies> = self.species_annots.get(&tail_index).unwrap().borrow_mut();
            // bindings no longer possible due to occupied sites
            let mut binary_combinatorics_lost: BTreeSet<BondEmbed> = BTreeSet::new();
            for open_tail in &self.ports.map.get(&bond_type.pair_2).unwrap().difference(&head_species.ports.map.get(&bond_type.pair_2).unwrap()).cloned().collect::<Vec<NodeIndex>>() {
                let b = self.interactions.get_mut(&bond_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: head_index, b_index: *open_tail, z: None, bond_type}))).unwrap();
                assert!(binary_combinatorics_lost.insert(*b.borrow()), "This lost binary pair ({}, {}) already cached into the binary lost tracker!", head_index.index(), open_tail.index())
            }
            for open_head in &self.ports.map.get(&bond_type.pair_1).unwrap().difference(&head_species.ports.map.get(&bond_type.pair_1).unwrap()).cloned().collect::<Vec<NodeIndex>>() {
                let b = self.interactions.get_mut(&bond_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: *open_head, b_index: tail_index, z: None, bond_type}))).unwrap();
                assert!(binary_combinatorics_lost.insert(*b.borrow()), "This lost binary pair ({}, {}) already cached into the binary lost tracker!", open_head.index(), tail_index.index())
            }
            let mut unary_combinatorics_lost: BTreeSet<BondEmbed> = BTreeSet::new();
            let purged_type = BondType{arity: InteractionArity::Unary, ..bond_type};
            for open_tail in head_species.ports.map.get(&bond_type.pair_2).unwrap().iter().collect::<Vec<&NodeIndex>>() {
                let b = self.interactions.get_mut(&purged_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: head_index, b_index: *open_tail, z: None, bond_type: purged_type}))).unwrap();
                assert!(unary_combinatorics_lost.insert(*b.borrow()), "This lost unary pair ({}, {}) already cached into the unary lost tracker!", head_index.index(), open_tail.index())
            }
            for open_head in tail_species.ports.map.get(&bond_type.pair_1).unwrap().iter().collect::<Vec<&NodeIndex>>() {
                let b = self.interactions.get_mut(&purged_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: *open_head, b_index: tail_index, z: None, bond_type: purged_type}))).unwrap();
                assert!(unary_combinatorics_lost.insert(*b.borrow()), "This lost unary pair ({}, {}) already cached into the unary lost tracker!", open_head.index(), tail_index.index())
            }
            assert!(self.ports.map.get_mut(&bond_type.pair_1).unwrap().remove(&head_index), "This head node was not globally listed as free already, can't bind to it!");
            assert!(self.ports.map.get_mut(&bond_type.pair_2).unwrap().remove(&tail_index), "This tail node was not globally listed as free already, can't bind to it!");
            assert!(head_species.ports.map.get_mut(&bond_type.pair_1).unwrap().remove(&head_index), "Port {} on index {} was not listed as free in-species already, can't bind to it!", bond_type.pair_1, head_index.index());
            assert!(tail_species.ports.map.get_mut(&bond_type.pair_2).unwrap().remove(&tail_index), "Port {} on index {} was not listed as free in-species already, can't bind to it!", bond_type.pair_2, tail_index.index());
            // binding opportunities that were binary, but will now be unary
            for some_binary_bond_type in self.interactions.keys().filter(|b| b.arity == InteractionArity::Binary && b.direction == InteractionDirection::Bind) {
                let some_unary_bond_type: BondType = BondType{arity: InteractionArity::Unary, ..**some_binary_bond_type};
                for open_side_a in tail_species.ports.map.get(&some_binary_bond_type.pair_1).unwrap().union(&head_species.ports.map.get(&some_binary_bond_type.pair_1).unwrap()).into_iter() {
                    for open_side_b in tail_species.ports.map.get(&some_binary_bond_type.pair_2).unwrap().union(&head_species.ports.map.get(&some_binary_bond_type.pair_2).unwrap()).into_iter() {
                        let transformed_embed = self.interactions.get_mut(some_binary_bond_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: *open_side_a, b_index: *open_side_b, z: None, bond_type: **some_binary_bond_type}))).unwrap();
                        assert!(self.interactions.get_mut(&some_unary_bond_type).unwrap().set.insert(transformed_embed), "This embed {} already cached as {}!", *transformed_embed.borrow(), some_unary_bond_type);
                    }
                }
            }
            // recast from "head & tail" to "host & eaten"; iterate over smaller sets while merging caches
            let (host_index, host_port_type, eaten_index, eaten_port_type) = if self.species_annots.get(&head_index).unwrap().borrow().size >= self.species_annots.get(&tail_index).unwrap().borrow().size {
                (head_index, bond_type.pair_1, tail_index, bond_type.pair_2)
            } else {
                (tail_index, bond_type.pair_2, head_index, bond_type.pair_1)
            };
            let mut host_species: RefMut<MixtureSpecies> = self.species_annots.get(&host_index).unwrap().borrow_mut();
            let mut eaten_species: Rc<RefCell<MixtureSpecies>> = *self.species_annots.get(&eaten_index).unwrap();
            let mut node_host: &mut Agent = self.universe_graph.node_weight_mut(host_index).unwrap();
            let mut node_eaten: &mut Agent = self.universe_graph.node_weight_mut(eaten_index).unwrap();
            // create & update edge trackers
            assert!(self.interactions.get(&bond_type).unwrap().set.remove(&target_bond), "The target bond {} was not already in the expected tracker for {}!", target_bond.borrow(), bond_type);
            let new_edge_index = self.universe_graph.add_edge(head_index, tail_index, false);
            bond_type.direction = InteractionDirection::Free;
            target_bond.borrow_mut().z = Some(new_edge_index);
            let agent_a: &mut Agent = self.universe_graph.node_weight_mut(head_index).unwrap();
            let agent_b: &mut Agent = self.universe_graph.node_weight_mut(tail_index).unwrap();
            let site_a = agent_a.sites.get_mut(bond_type.pair_1.site).unwrap();
            match site_a {
                Some(_) => panic!("Node {} was already bound at {}!", agent_a, bond_type.pair_1.site),
                None => *site_a = Some(new_edge_index)
            }
            let site_b = agent_b.sites.get_mut(bond_type.pair_2.site).unwrap();
            match site_b {
                Some(_) => panic!("Node {} was already bound at {}!", agent_b, bond_type.pair_2.site),
                None => *site_b = Some(new_edge_index)
            }
            host_species.edges.insert(Rc::clone(&target_bond));
            self.edge_index_map.insert(new_edge_index, Rc::clone(&target_bond));
            self.interactions.get_mut(&bond_type).unwrap().set.insert(Rc::clone(&target_bond));
            // update species annotation tracker
            self.species_set.make_contiguous().sort();
            let spec_ix: usize = self.species_set.binary_search(&eaten_species).unwrap();
            self.species_set.remove(spec_ix);
            let indexes_of_eaten: BTreeSet<NodeIndex> = eaten_species.borrow().agent_set.clone();
            for agent_index in indexes_of_eaten {
                let new_ref = Rc::clone(&self.species_annots.get(&host_index).unwrap());
                self.species_annots.entry(agent_index).and_modify(|e| {*e = new_ref});
            }
            host_species.ports.update_from(&eaten_species.borrow_mut().ports);
            host_species.edges.append(&mut eaten_species.borrow_mut().edges);
            host_species.size += eaten_species.borrow().size;
            host_species.agent_set.append(&mut eaten_species.borrow_mut().agent_set);
        }

        /**
         * Break a bond that is a cycle member; this does not partition one species into two, it
         * simply reduces the internal connectivity of said species.
        */
        pub fn unary_free(&mut self, target_bond: Rc<RefCell<BondEmbed<'a>>>) {
            let head_node: NodeIndex = target_bond.borrow().a_index;
            let tail_node: NodeIndex = target_bond.borrow().b_index;
            let bond_type: BondType<'a> = target_bond.borrow().bond_type;
            let edge_index: EdgeIndex = target_bond.borrow().z.unwrap();
            // sanity check for reslient bond
            assert!(self.universe_graph.edge_weight(edge_index).unwrap(), "This bond ({}) is not flagged as resilient; breaking it would not yield a still-connected component!", target_bond.borrow());
            let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_node).unwrap().borrow_mut();
            let old_last_edge_index = self.universe_graph.edge_indices().last().unwrap();
            self.universe_graph.remove_edge(edge_index).unwrap();
            let swapped_edge = self.edge_index_map.remove(&old_last_edge_index).unwrap();
            // deal with the last-edge-index invalidation on edge deletion
            // aka, the collaterals
            // we need to update the weight of Agents in the universe_graph,
            // and use the edge_index_map to update the BondEmbeds
            if old_last_edge_index != edge_index {
                swapped_edge.borrow_mut().z = Some(edge_index);
                let collateral_bond_type: BondType = swapped_edge.borrow().bond_type;
                let collateral_agent_a: &mut Agent = self.universe_graph.node_weight_mut(swapped_edge.borrow().a_index).unwrap();
                let collateral_agent_b: &mut Agent = self.universe_graph.node_weight_mut(swapped_edge.borrow().b_index).unwrap();
                let collateral_site_a = collateral_agent_a.sites.get(collateral_bond_type.pair_1.site).unwrap();
                let collateral_site_b = collateral_agent_b.sites.get(collateral_bond_type.pair_2.site).unwrap();
                match collateral_site_a {
                    Some(this_index) if this_index == &old_last_edge_index => *this_index = edge_index,
                    _ => panic!("Collateral node {} did not have the expected edge index {} for updating to {}!", collateral_agent_a, old_last_edge_index.index(), edge_index.index())
                }
                match collateral_site_b {
                    Some(this_index) if this_index == &old_last_edge_index => *this_index = edge_index,
                    _ => panic!("Collateral node {} did not have the expected edge index {} for updating to {}!", collateral_agent_b, old_last_edge_index.index(), edge_index.index())
                }
                *self.edge_index_map.get_mut(&edge_index).unwrap() = swapped_edge;
            }
            // sanity checks for bond breaking; then update bond indexes on the agent
            let node_a: &mut Agent = self.universe_graph.node_weight_mut(head_node).unwrap();
            let site_a = node_a.sites.get_mut(bond_type.pair_1.site).unwrap();
            site_a = match site_a {
                Some(some_edge) if *some_edge == edge_index => &mut None,
                _ => panic!("This node {} was not already bound at expected bond {}!", node_a, edge_index.index())
            };
            let node_b: &mut Agent = self.universe_graph.node_weight_mut(tail_node).unwrap();
            let site_b = node_b.sites.get_mut(bond_type.pair_2.site).unwrap();
            site_b = match site_b {
                Some(some_edge) if *some_edge == edge_index => &mut None,
                _ => panic!("This node {} was not already bound at expected bond {}!", node_b, edge_index.index())
            };
            assert!(target_species.edges.remove(&target_bond), "Bond {} was not already present as an available option in the species tracker!", target_bond.borrow());
            assert!(self.interactions.get_mut(&bond_type).unwrap().set.remove(&target_bond), "Bond {} was not already present as an available option!", target_bond.borrow());
            assert!(self.ports.map.get_mut(&bond_type.pair_1).unwrap().insert(head_node), "Site {} on node {} was already listed as free prior to unbinding!", bond_type.pair_1, node_a);
            assert!(self.ports.map.get_mut(&bond_type.pair_2).unwrap().insert(tail_node), "Site {} on node {} was already listed as free prior to unbinding!", bond_type.pair_2, node_b);
            assert!(target_species.ports.map.get_mut(&bond_type.pair_1).unwrap().insert(head_node), "Site {} on node {} was already listed as free prior to this unbinding on the species cache!", bond_type.pair_1, node_a);
            assert!(target_species.ports.map.get_mut(&bond_type.pair_2).unwrap().insert(tail_node), "Site {} on node {} was already listed as free prior to this unbinding on the species cache!", bond_type.pair_2, node_b);
            // bindings now possible due to free'd sites
            let binary_binding_subtype = BondType{arity: InteractionArity::Binary, direction: InteractionDirection::Bind, ..bond_type};
            for tail in &self.ports.map.get(&bond_type.pair_2).unwrap().difference(&target_species.ports.map.get(&bond_type.pair_2).unwrap()).cloned().collect::<Vec<NodeIndex>>() {
                assert!(self.interactions.get_mut(&binary_binding_subtype).unwrap().set.insert(Rc::new(RefCell::new(BondEmbed{a_index: head_node, b_index: *tail, z: None, bond_type: binary_binding_subtype}))), "This bond embed was already cached into the binary tracker!")
            }
            for head in &self.ports.map.get(&bond_type.pair_1).unwrap().difference(&target_species.ports.map.get(&bond_type.pair_1).unwrap()).cloned().collect::<Vec<NodeIndex>>() {
                assert!(self.interactions.get_mut(&binary_binding_subtype).unwrap().set.insert(Rc::new(RefCell::new(BondEmbed{a_index: *head, b_index: tail_node, z: None, bond_type: binary_binding_subtype}))), "This bond embed was already cached into the binary tracker!")
            }
            let unary_binding_subtype = BondType{arity: InteractionArity::Unary, direction: InteractionDirection::Bind, ..bond_type};
            for tail in &target_species.ports.map.get(&bond_type.pair_2).unwrap().iter().cloned().collect::<Vec<NodeIndex>>() {
                assert!(self.interactions.get_mut(&unary_binding_subtype).unwrap().set.insert(Rc::new(RefCell::new(BondEmbed{a_index: head_node, b_index: *tail, z: None, bond_type: unary_binding_subtype}))), "This bond embed was already cached into the unary tracker!")
            }
            self.interactions.get_mut(&unary_binding_subtype).unwrap().set.remove(&Rc::new(RefCell::new(BondEmbed{a_index: head_node, b_index: tail_node, z: None, bond_type: unary_binding_subtype}))); // hack to remove the pair that will be visited from the other side in the next loop
            for head in &target_species.ports.map.get(&bond_type.pair_1).unwrap().iter().cloned().collect::<Vec<NodeIndex>>() {
                assert!(self.interactions.get_mut(&unary_binding_subtype).unwrap().set.insert(Rc::new(RefCell::new(BondEmbed{a_index: *head, b_index: tail_node, z: None, bond_type: unary_binding_subtype}))), "This bond embed was already cached into the unary tracker!")
            }
            // update resiliency
            // move things from self.interactions that match Free & Unary (aka cycle edges) into
            // the corresponding
            // things in self.interactions that match Free & Binary (aka tree edges)
            let mut newly_treed_bond = |boxed_bond: &Rc<RefCell<BondEmbed>>| -> bool {
                let bond_embed = boxed_bond.borrow();
                if *self.universe_graph.edge_weight(bond_embed.z.unwrap()).unwrap() {
                    let mut counterfactual_graph = self.universe_graph.clone();
                    counterfactual_graph.remove_edge(bond_embed.z.unwrap());
                    let path = astar(&counterfactual_graph, bond_embed.a_index, |finish| finish == bond_embed.b_index, |_| 1, |_| 0);
                    match path {
                        Some(_) => (false),
                        None => {self.universe_graph.update_edge(bond_embed.a_index, bond_embed.b_index, true); true}
                    }
                } else {false}      // if the bond's weight was already false, then it was not tree'd
            };
            for (i_type, i_data) in self.interactions {
                if i_type.arity == InteractionArity::Unary && i_type.direction == InteractionDirection::Free {
                    let mut newly_treed_bonds: BTreeSet<Rc<RefCell<BondEmbed>>> = i_data.set.drain_filter(newly_treed_bond).collect();
                    let target_container = BondType{arity: InteractionArity::Binary, ..*i_type};
                    self.interactions.get_mut(&target_container).unwrap().set.append(&mut newly_treed_bonds);
                }
            }
        }

        /**
         * Break a bond that is not a cycle member; this partitions one species into two.
        */
        pub fn binary_free(&mut self, target_bond: Rc<RefCell<BondEmbed<'a>>>) {
            let head_node: NodeIndex = target_bond.borrow().a_index;
            let tail_node: NodeIndex = target_bond.borrow().b_index;
            let bond_type: BondType<'a> = target_bond.borrow().bond_type;
            let edge_index: EdgeIndex = target_bond.borrow().z.unwrap();
            // sanity check for non-reslient bond
            assert!(! self.universe_graph.edge_weight(edge_index).unwrap(), "This bond ({}) is flagged as resilient; breaking it would yield a still-connected component!", target_bond.borrow());
            let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_node).unwrap().borrow_mut();
            let old_last_edge_index = self.universe_graph.edge_indices().last().unwrap();
            self.universe_graph.remove_edge(edge_index).unwrap();
            let swapped_edge = self.edge_index_map.remove(&old_last_edge_index).unwrap();
            // deal with the last-edge-index invalidation on edge deletion
            // aka, the collaterals
            // we need to update the weight of Agents in the universe_graph,
            // and use the edge_index_map to update the BondEmbeds
            if old_last_edge_index != edge_index {
                swapped_edge.borrow_mut().z = Some(edge_index);
                let collateral_bond_type: BondType = swapped_edge.borrow().bond_type;
                let collateral_agent_a: &mut Agent = self.universe_graph.node_weight_mut(swapped_edge.borrow().a_index).unwrap();
                let collateral_agent_b: &mut Agent = self.universe_graph.node_weight_mut(swapped_edge.borrow().b_index).unwrap();
                let collateral_site_a = collateral_agent_a.sites.get(collateral_bond_type.pair_1.site).unwrap();
                let collateral_site_b = collateral_agent_b.sites.get(collateral_bond_type.pair_2.site).unwrap();
                match collateral_site_a {
                    Some(this_index) if this_index == &old_last_edge_index => *this_index = edge_index,
                    _ => panic!("Collateral node {} did not have the expected edge index {} for updating to {}!", collateral_agent_a, old_last_edge_index.index(), edge_index.index())
                }
                match collateral_site_b {
                    Some(this_index) if this_index == &old_last_edge_index => *this_index = edge_index,
                    _ => panic!("Collateral node {} did not have the expected edge index {} for updating to {}!", collateral_agent_b, old_last_edge_index.index(), edge_index.index())
                }
                *self.edge_index_map.get_mut(&edge_index).unwrap() = swapped_edge;
            }
            // sanity checks for bond breaking; update bond indexes on the agent-type
            let node_a: &mut Agent = self.universe_graph.node_weight_mut(head_node).unwrap();
            let site_a = node_a.sites.get_mut(bond_type.pair_1.site).unwrap();
            site_a = match site_a {
                Some(some_edge) if *some_edge == edge_index => &mut None,
                _ => panic!("This node {} was not already bound at expected bond {}!", node_a, edge_index.index())
            };
            let node_b: &mut Agent = self.universe_graph.node_weight_mut(tail_node).unwrap();
            let site_b = node_b.sites.get_mut(bond_type.pair_2.site).unwrap();
            site_b = match site_b {
                Some(some_edge) if *some_edge == edge_index => &mut None,
                _ => panic!("This node {} was not already bound at expected bond {}!", node_b, edge_index.index())
            };
            assert!(target_species.edges.remove(&target_bond), "Bond {} was not already present as an available option in the species tracker!", target_bond.borrow());
            assert!(self.interactions.get_mut(&bond_type).unwrap().set.remove(&target_bond), "Bond {} was not already present as an available option!", target_bond.borrow());
            assert!(self.ports.map.get_mut(&bond_type.pair_1).unwrap().insert(head_node), "Site {} on node {} was already listed as free prior to unbinding!", bond_type.pair_1, node_a);
            assert!(self.ports.map.get_mut(&bond_type.pair_2).unwrap().insert(tail_node), "Site {} on node {} was already listed as free prior to unbinding!", bond_type.pair_2, node_b);
            assert!(target_species.ports.map.get_mut(&bond_type.pair_1).unwrap().insert(head_node), "Site {} on node {} was already listed as free prior to this unbinding on the species cache!", bond_type.pair_1, node_a);
            assert!(target_species.ports.map.get_mut(&bond_type.pair_2).unwrap().insert(tail_node), "Site {} on node {} was already listed as free prior to this unbinding on the species cache!", bond_type.pair_2, node_b);
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
                        // we've mapped-out the head node's graph; tail node's graph is the difference of original set minus this mapped-out set
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
            let (ejected_mark, ejected_indexes, retained_mark, retained_indexes) =  // Iterate over small & lookup over big
                if head_graph_indexes.len() >= tail_graph_indexes.len()
                    {(head_node, head_graph_indexes, tail_node, tail_graph_indexes)}
                else 
                    {(tail_node, tail_graph_indexes, head_node, head_graph_indexes)};
            // create new species cache, in-place modify old species cache
            let retained_species: RefMut<MixtureSpecies> = self.species_annots.get_mut(&retained_mark).unwrap().borrow_mut();
            let ejected_ports: OpenPorts = retained_species.ports.eject_where(&ejected_indexes);
            let ejected_edges: BTreeSet<Rc<RefCell<BondEmbed>>> = retained_species.edges.drain_filter(|b| ejected_indexes.contains(&b.borrow().a_index)).collect();
            retained_species.size -= &ejected_indexes.len();
            retained_species.agent_set = retained_indexes;
            let new_species = Rc::new(RefCell::new(MixtureSpecies{
                ports: ejected_ports,
                edges: ejected_edges,
                size: ejected_indexes.len(),
                agent_set: ejected_indexes.clone()
            }));
            // update the species trackers, per node index & global
            for node_index in &ejected_indexes {
                self.species_annots.entry(*node_index).and_modify(|e| {*e = Rc::clone(&new_species)});
            }
            self.species_set.push_back(Rc::clone(&new_species));
            // binding opportunities that were unary, but will now be binary
            for some_unary_bond_type in self.interactions.keys().filter(|b| b.arity == InteractionArity::Unary && b.direction == InteractionDirection::Bind) {
                let some_binary_bond_type: BondType = BondType{arity: InteractionArity::Binary, ..**some_unary_bond_type};
                for open_side_a in retained_species.ports.map.get(&some_binary_bond_type.pair_1).unwrap().union(&new_species.borrow().ports.map.get(&some_binary_bond_type.pair_1).unwrap()).into_iter() {
                    for open_side_b in retained_species.ports.map.get(&some_binary_bond_type.pair_2).unwrap().union(&new_species.borrow().ports.map.get(&some_binary_bond_type.pair_2).unwrap()).into_iter() {
                        let transformed_embed = self.interactions.get_mut(some_unary_bond_type).unwrap().set.take(&Rc::new(RefCell::new(BondEmbed{a_index: *open_side_a, b_index: *open_side_b, z: None, bond_type: **some_unary_bond_type}))).unwrap();
                        assert!(self.interactions.get_mut(&some_binary_bond_type).unwrap().set.insert(transformed_embed), "This embed {} already cached as {}!", *transformed_embed.borrow(), some_binary_bond_type);
                    }
                }
            }
            // binding opportunities now possible due to free'd sites
            let head_species: Ref<MixtureSpecies> = self.species_annots.get(&head_node).unwrap().borrow();
            let tail_species: Ref<MixtureSpecies> = self.species_annots.get(&head_node).unwrap().borrow();
            let binary_binding_subtype = BondType{arity: InteractionArity::Binary, direction: InteractionDirection::Bind, ..bond_type};
            for tail in &self.ports.map.get(&bond_type.pair_2).unwrap().difference(&tail_species.ports.map.get(&bond_type.pair_2).unwrap()).cloned().collect::<Vec<NodeIndex>>() {
                assert!(self.interactions.get_mut(&binary_binding_subtype).unwrap().set.insert(Rc::new(RefCell::new(BondEmbed{a_index: head_node, b_index: *tail, z: None, bond_type: binary_binding_subtype}))), "This bond embed was already cached into the binary tracker!")
            }
            self.interactions.get_mut(&binary_binding_subtype).unwrap().set.remove(&Rc::new(RefCell::new(BondEmbed{a_index: head_node, b_index: tail_node, z: None, bond_type: binary_binding_subtype}))); // hack to remove the pair that will be visited from the other side in the next loop
            for head in &self.ports.map.get(&bond_type.pair_1).unwrap().difference(&head_species.ports.map.get(&bond_type.pair_1).unwrap()).cloned().collect::<Vec<NodeIndex>>() {
                assert!(self.interactions.get_mut(&binary_binding_subtype).unwrap().set.insert(Rc::new(RefCell::new(BondEmbed{a_index: *head, b_index: tail_node, z: None, bond_type: binary_binding_subtype}))), "This bond embed was already cached into the binary tracker!")
            }
            let unary_binding_subtype = BondType{arity: InteractionArity::Unary, direction: InteractionDirection::Bind, ..bond_type};
            for tail in &head_species.ports.map.get(&bond_type.pair_2).unwrap().iter().cloned().collect::<Vec<NodeIndex>>() {
                assert!(self.interactions.get_mut(&unary_binding_subtype).unwrap().set.insert(Rc::new(RefCell::new(BondEmbed{a_index: head_node, b_index: *tail, z: None, bond_type: unary_binding_subtype}))), "This bond embed was already cached into the unary tracker!")
            }
            for head in &tail_species.ports.map.get(&bond_type.pair_1).unwrap().iter().cloned().collect::<Vec<NodeIndex>>() {
                assert!(self.interactions.get_mut(&unary_binding_subtype).unwrap().set.insert(Rc::new(RefCell::new(BondEmbed{a_index: *head, b_index: tail_node, z: None, bond_type: unary_binding_subtype}))), "This bond embed was already cached into the unary tracker!")
            }
            // self-edge cleanup; these may be present, or not, depending on the mixture state
            self.interactions.get_mut(&binary_binding_subtype).unwrap().set.remove(&Rc::new(RefCell::new(BondEmbed{a_index: head_node, b_index: tail_node, z: None, bond_type: binary_binding_subtype})));
            self.interactions.get_mut(&binary_binding_subtype).unwrap().set.remove(&Rc::new(RefCell::new(BondEmbed{a_index: tail_node, b_index: head_node, z: None, bond_type: binary_binding_subtype})));
            self.interactions.get_mut(&unary_binding_subtype).unwrap().set.remove(&Rc::new(RefCell::new(BondEmbed{a_index: head_node, b_index: tail_node, z: None, bond_type: unary_binding_subtype})));
            self.interactions.get_mut(&unary_binding_subtype).unwrap().set.remove(&Rc::new(RefCell::new(BondEmbed{a_index: tail_node, b_index: head_node, z: None, bond_type: unary_binding_subtype})));
        }
    /*
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
    */
    }
}
/*
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
*/