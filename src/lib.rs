#![feature(hash_drain_filter, btree_drain_filter)]
#![warn(clippy::all)]

/**
 * Module for storing stuff to build the simulator; rules & abundances & rates & bond-types...
*/
pub mod building_blocks {
    use crate::primitives::AgentSite;
    use std::error::Error;
    use std::collections::{BTreeMap, BTreeSet};
    use std::fs::File;
    use std::io::BufReader;
    use std::path::Path;
    use serde::Deserialize;
    

    /**
     * A resource vector, using a map, where the key is an agent name, and the value is a protomer
     * abundance.
    */
    #[derive(Debug, Deserialize)]
    pub struct ProtomerResources {
        pub map: BTreeMap<String, usize>
    }

    impl ProtomerResources {
        pub fn from_json<P: AsRef<Path>>(path: P) -> Result<ProtomerResources, Box<dyn Error>> {
            let file = File::open(path)?;
            let reader = BufReader::new(file);
            let pr = serde_json::from_reader(reader)?;
            Ok(pr)
        }

        /// Total amount of protomers in the system.
        pub fn total_mass(&self) -> usize {
            self.map.values().sum()
        }
    }


    /**
     * Structure to group & name the four rate constants of a binding interaction, the agents &
     * sites involved in that interaction, and a user-memorable name for the binding interaction.
    */
    #[derive(Debug, Deserialize)]
    pub struct InteractionData {
        pub name: Option<String>,
        pub bind_unary: RateFudge,
        pub bind_binary: RateFudge,
        pub free_unary: RateFudge,
        pub free_binary: RateFudge,
        pub agent_a: String,
        pub agent_b: String,
        pub site_a: String,
        pub site_b: String
    }

    impl InteractionData {
        pub fn generate_agent_sites(&self) -> (AgentSite, AgentSite) {
            let pair_a = AgentSite{agent: &self.agent_a, site: &self.site_a};
            let pair_b = AgentSite{agent: &self.agent_b, site: &self.site_b};
            if pair_a <= pair_b { (pair_a, pair_b) }
            else { (pair_b, pair_a) }
        }
    }


    /**
     * Container, used for persistent storage and easier deserialization of the interaction set.
    */
    #[derive(Debug, Deserialize)]
    pub struct AllInteractionData {
        pub interactions: Vec<InteractionData>
    }

    impl AllInteractionData {
        /**
         * Read a JSON-formatted text file with interaction data. See files in `/examples`.
        */
        pub fn from_json<P: AsRef<Path>>(path: P) -> Result<AllInteractionData, Box<dyn Error>> {
            let file = File::open(path)?;
            let reader = BufReader::new(file);
            let id = serde_json::from_reader(reader)?;
            Ok(id)
        }

        /**
         * Yield a collection of the AgentSite types from the interaction data supplied.
        */
        pub fn generate_agent_sites(&self) -> BTreeSet<AgentSite> {
            let mut collector: BTreeSet<AgentSite> = BTreeSet::new();
            for r in self.interactions.iter() {
                let pair_a: AgentSite = AgentSite{agent: &r.agent_a, site: &r.site_a};
                let pair_b: AgentSite = AgentSite{agent: &r.agent_b, site: &r.site_b};
                collector.insert(pair_a);
                collector.insert(pair_b);
            }
            collector
        }
    }


    /**
     * A base rate constant, plus some fudge factors.
    */
    #[derive(Debug, Deserialize, PartialEq)]
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
    use std::{collections::{BTreeMap, BTreeSet}, fmt};
    use petgraph::prelude::{EdgeIndex, NodeIndex};
    

    /**
     * Represents an agent in the reaction mixture.
    */
    #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct Agent <'a> {
        pub name: &'a str,
        pub id: Option<NodeIndex>,
        pub sites: BTreeMap<&'a str, Option<EdgeIndex>>
    }

    impl <'a> Agent <'a> {
        pub fn new_from_agent_sites(agent_name: &'a str, sites: &BTreeSet<AgentSite<'a>>) -> Agent<'a> {
            Agent {
                name: agent_name,
                id: None,
                sites: BTreeMap::from_iter(sites.iter().map(|n| (n.site, None)))
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
     * Represents an interaction type: a pair of agents, their sites, the operation's direction,
     * and the operation's arity.
     * Constructor orders alphabetically by the agent names, then by site names. This order is
     * observed for comparison and sorting.
     * ```
     * use axin_apc_simulator::primitives::{AgentSite, BondType, InteractionArity, InteractionDirection};
     * let foo = BondType::new(AgentSite{agent: "Mango", site: "stone"}, AgentSite{agent: "Earth", site: "ground"}, InteractionArity::Binary, InteractionDirection::Bind);
     * assert_eq!("Earth(ground[.]) + Mango(stone[.]) -> Earth(ground[0]) , Mango(stone[0])", format!("{}", foo));
     *
     * /// The interactions are displayed as:
     * 
     * let un_bi = BondType::new(AgentSite{agent: "A", site: "b"}, AgentSite{agent: "C", site: "d"}, InteractionArity::Unary, InteractionDirection::Bind);
     * assert_eq!("A(b[.]) , C(d[.]) -> A(b[0]) , C(d[0])", format!("{}", un_bi));
     * 
     * let un_fr = BondType::new(AgentSite{agent: "A", site: "b"}, AgentSite{agent: "C", site: "d"}, InteractionArity::Unary, InteractionDirection::Free);
     * assert_eq!("A(b[0]) , C(d[0]) -> A(b[.]) , C(d[.])", format!("{}", un_fr));
     *
     * let bi_bi = BondType::new(AgentSite{agent: "A", site: "b"}, AgentSite{agent: "C", site: "d"}, InteractionArity::Binary, InteractionDirection::Bind);
     * assert_eq!("A(b[.]) + C(d[.]) -> A(b[0]) , C(d[0])", format!("{}", bi_bi));
     * 
     * let bi_fr = BondType::new(AgentSite{agent: "A", site: "b"}, AgentSite{agent: "C", site: "d"}, InteractionArity::Binary, InteractionDirection::Free);
     * assert_eq!("A(b[0]) , C(d[0]) -> A(b[.]) + C(d[.])", format!("{}", bi_fr));
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
            write!(f, "{}", self.as_string_decorated(None, None, None))
        }
    }

    impl <'a> BondType <'a> {
        pub fn new(pair_a: AgentSite<'a>, pair_b: AgentSite<'a>, arity: InteractionArity, direction: InteractionDirection) -> Self {
            if pair_a <= pair_b {
                BondType{pair_1: pair_a, pair_2: pair_b, arity, direction}
            } else {
                BondType{pair_1: pair_b, pair_2: pair_a, arity, direction}
            }
        }

        fn as_string_decorated(&self, agent_1_index: Option<NodeIndex>, agent_2_index: Option<NodeIndex>, bond_index: Option<EdgeIndex>) -> String {
            let b: String = match bond_index {
                Some(i) => i.index().to_string(),
                None => String::from("0")
            };
            let (pre_bond, post_bond, pre_join, post_join) = match self.direction {
                InteractionDirection::Bind => match self.arity {
                    InteractionArity::Binary => (".", &b[..], "+", ","),
                    InteractionArity::Unary => (".", &b[..], ",", ",")
                }
                InteractionDirection::Free => match self.arity {
                    InteractionArity::Binary => (&b[..], ".", ",", "+"),
                    InteractionArity::Unary => (&b[..], ".", ",", ",")
                }
            };
            let prefix_1: String = match agent_1_index {
                Some(i) => format!("x{}:", i.index()),
                None => String::from("")
            };
            let prefix_2: String = match agent_2_index {
                Some(i) => format!("x{}:", i.index()),
                None => String::from("")
            };
            let lhs = format!("{}{}({}[{}]) {} {}{}({}[{}])", prefix_1, self.pair_1.agent, self.pair_1.site, pre_bond, pre_join, prefix_2, self.pair_2.agent, self.pair_2.site, pre_bond);
            let rhs = format!("{}{}({}[{}]) {} {}{}({}[{}])", prefix_1, self.pair_1.agent, self.pair_1.site, post_bond, post_join, prefix_2, self.pair_2.agent, self.pair_2.site, post_bond);
            format!("{} -> {}", &lhs, &rhs)
        }
    }

    
    /** 
     * Structure for representing a bond embed, either current or possible, along with the 
     * corresponding edge index, if any. The agent indexes correspond to the agent identifiers.
    */
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub struct BondEmbed <'a> {
        /// Index of the left agent, from the universe_graph.
        pub a_index: NodeIndex,

        /// Index of the right agent, from the universe_graph.
        pub b_index: NodeIndex,

        /// The type of the interaction.
        pub bond_type: BondType <'a>,

        /// For interactions currently bound, the index from the universe_graph.
        /// 
        /// This index is invalidated when it is the last index in the graph and a bond deletion
        /// operation occurs; hence why BondEmbeds are used behind a smart pointer 
        /// Rc(RefCell(BondEmbed)), so the indexes can be updated as "collateral" in the event
        /// cycle that invalidates one. This also means the sorting evaluation should not consider
        /// this item. Using the derived versions of PartialOrd will use `z` in sorting; however in
        /// practice, no embed tracker should have more than one entry with the same `a_index`,
        /// `b_index`, and `bond_type`, and so the derived sorting should never end up utilizing
        /// the information on `z`.
        pub z: Option<EdgeIndex>,
    }
    
    impl <'a> fmt::Display for BondEmbed <'a> {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "{}", self.bond_type.as_string_decorated(Some(self.a_index), Some(self.b_index), self.z))
        }
    }


    /**
     * Marks if an interaction occurs within a complex, or between complexes. Used for type'ing
     * transformations, embed trackers and their constructors.
    */
    #[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
    pub enum InteractionArity {
        Binary,
        Unary,
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
    use crate::building_blocks::{InteractionData, ProtomerResources};
    use crate::primitives::{Agent, AgentSite, BondEmbed, BondType, InteractionArity, InteractionDirection};
    

    /**
     * Annotation for a specific species in the mixture.
    */
    #[derive(Debug, PartialEq, Eq)]
    pub struct MixtureSpecies <'a> {
        pub agent_set: BTreeSet<Rc<RefCell<Agent<'a>>>>,
        pub edges: BTreeSet<Rc<RefCell<BondEmbed<'a>>>>,
        pub ports: OpenPorts <'a>,
    }

    impl <'a> fmt::Display for MixtureSpecies <'a> {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(f, "{}", self.agent_set.iter().map(|a| format!("{}", a.borrow())).join(", "))
        }
    }

    // Manually derived to use agent_set length as most significant sorting attribute;
    // idem for PartialOrd. This allows snapshots to print sorted species by size. Since the
    // VecDeque will accumulate the biggest species on one end, fetching it will be efficient.
    // Sorting a close-to-sorted VecDeque will also be efficient.
    impl <'a> Ord for MixtureSpecies <'a> {
        fn cmp(&self, other: &Self) -> Ordering {
            if self.agent_set.len() == other.agent_set.len() {
                if self.agent_set == other.agent_set {
                    if self.edges == other.edges {
                        self.ports.cmp(&other.ports)
                    } else {self.edges.cmp(&other.edges) }
                } else { self.agent_set.cmp(&other.agent_set) }
            } else { self.agent_set.len().cmp(&other.agent_set.len()) }
        }
    }

    impl <'a> PartialOrd for MixtureSpecies <'a> {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }


    /** 
     * Map of HashSets, one per site type. The HashSets hold NodeIndexes, from the UniverseGraph.
     * The key of the map is the site name.
     */
    #[derive(Clone, Debug, PartialEq, Eq)]
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
        pub fn new_monomer(agent_sites: &BTreeSet<AgentSite<'a>>, index: &NodeIndex) -> Self {
            OpenPorts {
                map: BTreeMap::from_iter(agent_sites.iter().map(|n| (*n, HashSet::from_iter([*index]))))
            }
        }

        /**
         * Constructor for the universe tracker.
        */
        pub fn new_universe(agent_sites: &BTreeSet<AgentSite<'a>>, site_abundances: &ProtomerResources) -> Self {
            let mut map: BTreeMap<AgentSite<'a>, HashSet<NodeIndex>> = BTreeMap::new();
            for agent_site in agent_sites {
                map.insert(*agent_site, HashSet::with_capacity(*site_abundances.map.get(agent_site.agent).unwrap()));
            }
            OpenPorts {map}
        }
    }

    
    /**
     * Convenience structure for calculating the activity of interactions, while keeping the bias
     * metrics & calculation separate.
    */
    #[derive(Clone, Debug, PartialEq, PartialOrd)]
    pub struct MassActionParameters {
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
        pub fn calculate_activity(&self, mass: usize, a_size: usize, b_size: usize) -> f64 {
            let m: f64 = match self.bias_m {
                Some(f) => f * a_size as f64 * b_size as f64,
                None => 1.0
            };
            let a: f64 = match self.bias_a {
                Some(f) => f * a_size as f64 * b_size as f64,
                None => 0.0
            };
            mass as f64 * self.rate * m + a
        }
    }


    /** 
     * Collects the interacting pairs, either bound already, or eligible for a future binding
     * event.
     */
    #[derive(Clone, Debug, PartialEq, PartialOrd)]
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
         * Constructor used when starting from a monomeric initial state. For interactions of 
         * `InteractionArity::Binary` and `InteractionDirection::Bind`, it will generate the 
         * appropriate combinatorics for the possible bond types. For the other classes, it will be
         * empty, but pre-allocated.
        */
        pub fn new_tracker(interaction: &'a InteractionData, resources: &'a ProtomerResources, arity: InteractionArity, direction: InteractionDirection) -> Self {
            let (pair_a, pair_b) = interaction.generate_agent_sites();          // this provides sorted names, which is key for the tandem-build of this tracker with the universe_graph keeping indexing in sync
            let mut index_offset: usize = resources.map.range(..String::from(pair_a.agent)).map(|i| *i.1).sum();   // cumulative sum of previous agents
            let bond_type = BondType::new(pair_a, pair_b, arity, direction);
            let mut set: BTreeSet<Rc<RefCell<BondEmbed>>> = BTreeSet::new();
            let (activity_parameters, activity_calculated) = match (arity, direction) {
                (InteractionArity::Binary, InteractionDirection::Bind) => {
                    let a_indexes: Vec<NodeIndex> = (index_offset..(index_offset + resources.map.get(pair_a.agent).unwrap())).map(|x| NodeIndex::from(x as u32)).collect();
                    index_offset += resources.map.get(pair_a.agent).unwrap();
                    if pair_a.agent != pair_b.agent {
                        // heterogenous dimerization, using cartesian product
                        let b_indexes: Vec<NodeIndex> = (index_offset..(index_offset + resources.map.get(pair_b.agent).unwrap())).map(|x| NodeIndex::from(x as u32)).collect();
                        index_offset += resources.map.get(pair_b.agent).unwrap();
                        for (a_ix, b_ix) in a_indexes.iter().cartesian_product(b_indexes.iter()) {
                            let comb = BondEmbed{a_index: *a_ix, b_index: *b_ix, z: None, bond_type};
                            assert!(set.insert(Rc::new(RefCell::new(comb))));
                        }
                    }
                    else if pair_a.site != pair_b.site {
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
                    (MassActionParameters{rate: interaction.bind_binary.base, bias_m: interaction.bind_binary.multiplicative, bias_a: interaction.bind_binary.additive}, act)
                },
                (InteractionArity::Unary, InteractionDirection::Bind) => {
                    (MassActionParameters{rate: interaction.bind_unary.base, bias_m: interaction.bind_unary.multiplicative, bias_a: interaction.bind_unary.additive}, 0.0)
                },
                (InteractionArity::Binary, InteractionDirection::Free) => {
                    (MassActionParameters{rate: interaction.free_binary.base, bias_m: interaction.free_binary.multiplicative, bias_a: interaction.free_binary.additive}, 0.0)
                },
                (InteractionArity::Unary, InteractionDirection::Free) => {
                    (MassActionParameters{rate: interaction.free_unary.base, bias_m: interaction.free_unary.multiplicative, bias_a: interaction.free_unary.additive}, 0.0)
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
                self.activity_parameters.calculate_activity(self.set.len(), species_annots.get(&index_a).unwrap().borrow().agent_set.len(), species_annots.get(&index_b).unwrap().borrow().agent_set.len())
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
         * Using complex sizes in the weighting, pick a target embed.
        */
        pub fn pick_targets(&self, rng: &mut ThreadRng, species_annots: &BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies<'a>>>>) -> Rc<RefCell<BondEmbed<'a>>> {
            let dist = WeightedIndex::new(self.build_embed_weights(species_annots)).unwrap();
            let chosen_ix: usize = dist.sample(rng);
            let chosen_embed: Rc<RefCell<BondEmbed>> = Rc::clone(self.set.iter().nth(chosen_ix).unwrap());
            // To avoid mutating self, this returns a new reference to the Rc<RefCell<BondEmbed>>,
            // and so the mixture transformation methods get a reference of their own
            chosen_embed
        }
    }

    impl <'a> fmt::Display for InteractingTracker <'a> {
        fn fmt(&self, f:&mut fmt::Formatter<'_>) -> fmt::Result {
            if self.set.is_empty() { write!(f, "[empty interaction set]") }
            else { write!(f, "{}", self.set.iter().map(|b| format!("{}", b.borrow())).collect::<Vec<String>>().join("\n")) }
        }
    }    
}


/**
 * Module for the reaction mixture & its transformations.
*/
pub mod reaction_mixture {
    use crate::building_blocks::{ProtomerResources, AllInteractionData};
    use crate::primitives::{Agent, AgentSite, BondType, BondEmbed, InteractionArity, InteractionDirection};
    use crate::collectors::{OpenPorts, MixtureSpecies, InteractingTracker};
    use std::cell::{RefCell, RefMut, Ref};
    use std::collections::{VecDeque, BTreeMap, BTreeSet, HashSet};
    use std::fs::File;
    use std::io::prelude::*;
    use std::path::Path;
    use std::rc::Rc;
    use petgraph::{graph::Graph, algo::astar::astar, prelude::*, visit::Dfs};
    use rand::{distributions::WeightedIndex, prelude::*};
    use uuid::Uuid;



    #[derive(Clone)]
    pub struct Mixture <'a> {
        /// The graph that tracks connectivity; used for cycle detection & species partitioning.
        universe_graph: Graph<Rc<RefCell<Agent<'a>>>, bool, Undirected>,
        
        /// Maps the species cache to each node in the graph.
        species_annots: BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies<'a>>>>,
        
        /// Used for iteration for pretty-printing as a Kappa snapshot.
        species_set: VecDeque<Rc<RefCell<MixtureSpecies<'a>>>>,

        /// Stores a global cache of open ports.
        ports: OpenPorts<'a>,
        
        /// Stores the interactions in the system.
        interactions: BTreeMap<BondType<'a>, InteractingTracker<'a>>,

        /// Maps an EdgeIndex to a BondEmbed, for quick updating when deleting a bond invalidates
        /// the last EdgeIndex in the graph.
        edge_index_map: BTreeMap<EdgeIndex, Rc<RefCell<BondEmbed<'a>>>>,

        simulated_time: f64,
        simulated_events: usize,
        uuid: Uuid,

        /// Cache for the random number generator.
        simulator_rng: ThreadRng,
    }
    
    impl <'a> Mixture <'a> {

        /**
         * Advance simulation time and event number tracker.
        */
        fn advance_simulation_metrics(&mut self) {
            self.simulated_events += 1;
            let mut tau: f64 = 0.0;
            for i_data in self.interactions.values() {
                tau += i_data.activity_calculated;
            }
            let rdy: f64 = (1.0 / self.simulator_rng.gen::<f64>()).ln();
            self.simulated_time += tau * rdy;
        }

        /**
         * Based on the activities of the interactions, choose one to realize.
        */
        pub fn choose_and_apply_next_rule(&mut self, time_p: Option<f64>, event_p: Option<usize>, dir: &str) {
            // if next event would be past the observation period, dump snapshot
            let time_pre: f64 = self.simulated_time();
            self.advance_simulation_metrics();
            let time_pst: f64 = self.simulated_time();
            if let Some(period) = time_p {
                let cycle_pre = time_pre.div_euclid(period);
                let cycle_pst = time_pst.div_euclid(period);
                if cycle_pre != cycle_pst {
                    self.snapshot_to_file(dir, Some(cycle_pst * period), Some(self.simulated_events - 1))
                }
            };
            // perform event
            let dist = WeightedIndex::new(self.interactions.iter().map(|i| i.1.activity_calculated)).unwrap();
            let chosen_index: usize = dist.sample(&mut self.simulator_rng);
            let chosen_bond_type: BondType = *self.interactions.iter().nth(chosen_index).unwrap().0;
            let i = self.interactions.get_mut(&chosen_bond_type).unwrap();
            let bond: Rc<RefCell<BondEmbed>> = i.pick_targets(&mut self.simulator_rng, &self.species_annots);
            println!("Event: {:>5},\tTime: {:.5},\tchosen: {}", self.simulated_events, self.simulated_time, bond.borrow());
            match (chosen_bond_type.arity, chosen_bond_type.direction) {
                (InteractionArity::Unary, InteractionDirection::Bind) => self.unary_bind(bond),
                (InteractionArity::Unary, InteractionDirection::Free) => self.unary_free(bond),
                (InteractionArity::Binary, InteractionDirection::Bind) => self.binary_bind(bond),
                (InteractionArity::Binary, InteractionDirection::Free) => self.binary_free(bond)
            }
            // check if event period is satisfied; if so dump snapshot
            if event_p.is_some() && (self.event_number() % event_p.unwrap() == 0) {
                self.snapshot_to_file(dir, None, None)
            }
        }

        /**
         * Simulate a given amount of events, optionally dumping snapshots with some event
         * periodicity.
        */
        pub fn simulate_up_to_event(&mut self, max_event: usize, snap_p: Option<usize>, dir: &str) {
            while self.event_number() <= max_event {
                self.choose_and_apply_next_rule(None, snap_p, dir);
            }
        }

        /**
         * Simulate a given amount of events, optionally dumping snapshots with some event
         * periodicity.
        */
        pub fn simulate_up_to_time(&mut self, max_time: f64, snap_p: Option<f64>, dir: &str) {
            while self.simulated_time() <= max_time {
                self.choose_and_apply_next_rule(snap_p, None, dir);
            }
        }

        /**
         * Yield a snapshot representation that can be printed or compared against others in a
         * kappa-aware way.
        */
        pub fn to_kappa(&self, time_override: Option<f64>, event_override: Option<usize>) -> String {
            let mut mixture_vec: Vec<String> = Vec::with_capacity(self.species_set.len() + 3);
            let event_to_print: usize = match event_override {
                Some(e) => e,
                None => self.simulated_events
            };
            mixture_vec.push(format!("// Snapshot [Event:{}]\n// \"uuid\" : \"{}\"", event_to_print, self.uuid));
            mixture_vec.push(format!("// Mixture has {} protomers and {} bonds in {} species", self.universe_graph.node_count(), self.universe_graph.edge_count(), self.species_set.len()));
            let time_to_print: f64 = match time_override {
                Some(t) => t,
                None => self.simulated_time
            };
            mixture_vec.push(format!("%def: \"T0\" \"{}\"\n", time_to_print));
            let mut cloned_species = self.species_set.clone();
            cloned_species.make_contiguous().sort_by(|a, b| b.cmp(a));
            for this_species in cloned_species {
                mixture_vec.push(format!("init: 1 /*{} agents*/ {}", this_species.borrow().agent_set.len(), this_species.borrow()));
            }
            let mut ka_str: String = mixture_vec.join("\n");
            ka_str.push('\n');
            ka_str
        }

        /**
         * Save to file in kappa-snapshot format.
        */
        pub fn snapshot_to_file(&self, dir: &str, time_override: Option<f64>, event_override: Option<usize>) {
            let file_name: String = format!("{}/snapshot_{}.ka", dir, self.event_number());
            let path = Path::new(&file_name[..]);
            let display = path.display();
            let mut file = match File::create(&path) {
                Err(why) => panic!("Could not create {}: {}", display, why),
                Ok(file) => file
            };
            if let Err(why) = file.write_all(self.to_kappa(time_override, event_override).as_bytes()) { panic!("Could not write to {}: {}", display, why) }
        }

        /**
         * Create a new mixture from a desired number of monomers, and desired rule rates.
        */
        pub fn new_monomeric_from(raw_abundances: &'a ProtomerResources, raw_interactions: &'a AllInteractionData) -> Mixture<'a> {
            let mut universe_graph: Graph<Rc<RefCell<Agent>>, bool, Undirected> = Graph::with_capacity(raw_abundances.total_mass(), raw_abundances.total_mass() ^ 2);  // ToDo: better estimation for max. edge count
            let mut species_annots: BTreeMap<NodeIndex, Rc<RefCell<MixtureSpecies>>> = BTreeMap::new();
            let mut species_set: VecDeque<Rc<RefCell<MixtureSpecies>>> = VecDeque::with_capacity(raw_abundances.total_mass());
            let mut ports: OpenPorts = OpenPorts::new_universe(&raw_interactions.generate_agent_sites(), raw_abundances);
            let mut interactions: BTreeMap<BondType<'a>, InteractingTracker<'a>> = BTreeMap::new();
            let edge_index_map: BTreeMap<EdgeIndex, Rc<RefCell<BondEmbed<'a>>>> = BTreeMap::new();
            let simulated_time: f64 = 0.0;
            let simulated_events: usize = 0;
            let uuid: Uuid = Uuid::new_v4();
            let simulator_rng = rand::thread_rng();
            // iterate over raw_abundances to populate universe_graph, species_annots, species_set, and ports
            for (this_protomer, this_abundance) in &raw_abundances.map {
                let site_set: BTreeSet<AgentSite> = raw_interactions.generate_agent_sites().drain_filter(|a| a.agent == this_protomer).collect();
                for _ix in 0..*this_abundance {
                    let this_agent: Rc<RefCell<Agent>> = Rc::new(RefCell::new(Agent::new_from_agent_sites(this_protomer, &site_set)));
                    let node_ix: NodeIndex = universe_graph.add_node(Rc::clone(&this_agent));
                    this_agent.borrow_mut().id = Some(node_ix);
                    let new_species = Rc::new(RefCell::new(MixtureSpecies {
                        agent_set: BTreeSet::from_iter([Rc::clone(&this_agent)]),
                        edges: BTreeSet::new(),
                        ports: OpenPorts::new_monomer(&site_set, &node_ix)
                    }));
                    species_annots.insert(node_ix, Rc::clone(&new_species));
                    species_set.push_back(Rc::clone(&new_species));
                    for some_site in &site_set {
                        ports.map.entry(*some_site).or_insert_with(HashSet::new).insert(node_ix);
                    }
                }
            }
            // iterate over raw_interactions to populate interactions
            for some_interaction in &raw_interactions.interactions {
                let (pair_a, pair_b) = some_interaction.generate_agent_sites();
                for ari in [InteractionArity::Binary, InteractionArity::Unary] {
                    for dir in [InteractionDirection::Bind, InteractionDirection::Free] {
                        let bond_type = BondType::new(pair_a, pair_b, ari, dir);
                        let mut tracker = InteractingTracker::new_tracker(some_interaction, raw_abundances, ari, dir);
                        tracker.update_activity(&species_annots);
                        interactions.insert(bond_type, tracker);
                    }
                }
            }
            Mixture{universe_graph, species_annots, species_set, ports, interactions, edge_index_map, simulated_time, simulated_events, uuid, simulator_rng}
        }
    
        /**
         * Create a bond that creates a cycle, connecting two agents that were already connected
         * via some other components.
        */
        pub fn unary_bind(&mut self, target_bond: Rc<RefCell<BondEmbed<'a>>>) {
            let head_index: NodeIndex = target_bond.borrow().a_index;
            let tail_index: NodeIndex = target_bond.borrow().b_index;
            let mut bond_type: BondType<'a> = target_bond.borrow_mut().bond_type;
            assert!(bond_type.arity == InteractionArity::Unary, "Incorrect arity in interaction type! Expected {}, got {}", InteractionArity::Unary, bond_type.arity);
            assert!(bond_type.direction == InteractionDirection::Bind, "Incorrect direction in interaction type! Expected {}, got {}", InteractionDirection::Bind, bond_type.direction);
            // sanity check for unary-ness
            assert_eq!(self.species_annots.get(&head_index), self.species_annots.get(&tail_index), "This head and tail nodes ({}, {}) do not already belong to the same species; can't unary bind them!", head_index.index(), tail_index.index());
            {
                let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_index).unwrap().borrow_mut();
                // create edge & update caches
                assert!(self.interactions.get_mut(&bond_type).unwrap().set.remove(&target_bond), "The target bond {} was not already in the expected tracker for {}!", target_bond.borrow(), bond_type);
                let new_edge_index = self.universe_graph.add_edge(head_index, tail_index, true);
                bond_type.direction = InteractionDirection::Free;
                target_bond.borrow_mut().z = Some(new_edge_index);
                target_bond.borrow_mut().bond_type = bond_type;
                {                
                    let mut agent_a: RefMut<Agent> = self.universe_graph.node_weight_mut(head_index).unwrap().borrow_mut();
                    let site_a = agent_a.sites.get_mut(bond_type.pair_1.site).unwrap();
                    match site_a {
                        Some(_) => panic!("Node {} was already bound at {}!", agent_a, bond_type.pair_1.site),
                        None => *site_a = Some(new_edge_index)
                    }
                }
                {
                    let mut agent_b: RefMut<Agent> = self.universe_graph.node_weight_mut(tail_index).unwrap().borrow_mut();
                    let site_b = agent_b.sites.get_mut(bond_type.pair_2.site).unwrap();
                    match site_b {
                        Some(_) => panic!("Node {} was already bound at {}!", agent_b, bond_type.pair_2.site),
                        None => *site_b = Some(new_edge_index)
                    }
                }
                target_species.edges.insert(Rc::clone(&target_bond));
                self.edge_index_map.insert(new_edge_index, Rc::clone(&target_bond));
                self.interactions.get_mut(&bond_type).unwrap().set.insert(Rc::clone(&target_bond));
                // bindings no longer possible due to occupied sites
                let binary_binding_types: Vec<BondType> = self.interactions.keys().filter(|b| b.arity == InteractionArity::Binary && b.direction == InteractionDirection::Bind).cloned().collect();
                for some_bin_bond_type in binary_binding_types {
                    let some_uni_bond_type: BondType = BondType{arity: InteractionArity::Unary, ..some_bin_bond_type};
                    let bi_combinatorics_lost: BTreeSet<Rc<RefCell<BondEmbed>>>;
                    let un_combinatorics_lost: BTreeSet<Rc<RefCell<BondEmbed>>>;
                    (bi_combinatorics_lost, un_combinatorics_lost) = if some_bin_bond_type.pair_1 == some_bin_bond_type.pair_2 {
                        (
                            self.interactions.get_mut(&some_bin_bond_type).unwrap().set.drain_filter(|b| (b.borrow().a_index == head_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_1) || (b.borrow().b_index == tail_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_2) || (b.borrow().a_index == tail_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_2) || (b.borrow().b_index == head_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_1)).collect(),
                            self.interactions.get_mut(&some_uni_bond_type).unwrap().set.drain_filter(|b| (b.borrow().a_index == head_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_1) || (b.borrow().b_index == tail_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_2) || (b.borrow().a_index == tail_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_2) || (b.borrow().b_index == head_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_1)).collect()
                        )
                    } else {
                        (
                            self.interactions.get_mut(&some_bin_bond_type).unwrap().set.drain_filter(|b| (b.borrow().a_index == head_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_1) || (b.borrow().b_index == tail_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_2)).collect(),
                            self.interactions.get_mut(&some_uni_bond_type).unwrap().set.drain_filter(|b| (b.borrow().a_index == head_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_1) || (b.borrow().b_index == tail_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_2)).collect()
                        )
                    };
                }
                assert!(self.ports.map.get_mut(&bond_type.pair_1).unwrap().remove(&head_index));
                assert!(self.ports.map.get_mut(&bond_type.pair_2).unwrap().remove(&tail_index));
                assert!(target_species.ports.map.get_mut(&bond_type.pair_1).unwrap().remove(&head_index));
                assert!(target_species.ports.map.get_mut(&bond_type.pair_2).unwrap().remove(&tail_index));
            }
            // update resiliency
            // move things from self.interactions that match Free & Binary (aka tree edges) into
            // the corresponding
            // things in self.interactions that match Free & Unary (aka cycle edges)
            let is_newly_cyclized_bond = |boxed_bond: &Rc<RefCell<BondEmbed>>| -> bool {
                let bond_embed = boxed_bond.borrow();
                if ! self.universe_graph.edge_weight(bond_embed.z.unwrap()).unwrap() {  // check bond's weight, which is a boolean, marking if the bond is an already known cycle-member
                    let mut counterfactual_graph = self.universe_graph.clone();
                    counterfactual_graph.remove_edge(bond_embed.z.unwrap());
                    let path = astar(&counterfactual_graph, bond_embed.a_index, |finish| finish == bond_embed.b_index, |_| 1, |_| 0);
                    path.is_some()
                } else {false}      // if the bond's weight was already true, then it is not newly cyclized
            };
            // double for-loop to maintain single-mutable-borrow of self.interactions
            let mut temp_container: BTreeMap<BondType, BTreeSet<Rc<RefCell<BondEmbed>>>> = BTreeMap::new();
            for (i_type, i_data) in self.interactions.iter_mut() {
                if i_type.arity == InteractionArity::Binary && i_type.direction == InteractionDirection::Free {
                    let newly_cyclized_bonds: BTreeSet<Rc<RefCell<BondEmbed>>> = i_data.set.drain_filter(is_newly_cyclized_bond).collect();
                    let target_type = BondType{arity: InteractionArity::Unary, ..*i_type};
                    for this_embed in &newly_cyclized_bonds {
                        this_embed.borrow_mut().bond_type = target_type;
                    }
                    temp_container.insert(target_type, newly_cyclized_bonds);
                }
            }
            for (target_type, mut conv_data) in temp_container {
                for conv_bond in &conv_data {
                    let edge_weight = self.universe_graph.edge_weight_mut(conv_bond.borrow().z.unwrap()).unwrap();
                    *edge_weight = true;
                }
                self.interactions.get_mut(&target_type).unwrap().set.append(&mut conv_data);
            }
            // update rule activities
            for i_data in self.interactions.values_mut() {
                i_data.update_activity(&self.species_annots);
            }
        }
    
        /**
         * Create a bond that joins two species into one.
        */
        pub fn binary_bind(&mut self, target_bond: Rc<RefCell<BondEmbed<'a>>>) {
            let head_index: NodeIndex = target_bond.borrow().a_index;
            let tail_index: NodeIndex = target_bond.borrow().b_index;
            let mut bond_type: BondType<'a> = target_bond.borrow_mut().bond_type;
            assert!(bond_type.arity == InteractionArity::Binary, "Incorrect arity in interaction type! Expected {}, got {}", InteractionArity::Binary, bond_type.arity);
            assert!(bond_type.direction == InteractionDirection::Bind, "Incorrect direction in interaction type! Expected {}, got {}", InteractionDirection::Bind, bond_type.direction);
            // sanity check for binary-ness
            assert_ne!(self.species_annots.get(&head_index).unwrap(), self.species_annots.get(&tail_index).unwrap(), "Nodes ({}, {}) already belong to the same species; can't binary bind them!", head_index.index(), tail_index.index());
            {
                let mut head_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_index).unwrap().borrow_mut();
                let mut tail_species: RefMut<MixtureSpecies> = self.species_annots.get(&tail_index).unwrap().borrow_mut();
                assert!(self.ports.map.get_mut(&bond_type.pair_1).unwrap().remove(&head_index), "Port {} on index {} was not globally listed as free already, can't bind to it!", bond_type.pair_1, head_index.index());
                assert!(self.ports.map.get_mut(&bond_type.pair_2).unwrap().remove(&tail_index), "Port {} on index {} was not globally listed as free already, can't bind to it!", bond_type.pair_2, tail_index.index());
                assert!(head_species.ports.map.get_mut(&bond_type.pair_1).unwrap().remove(&head_index), "Port {} on index {} was not listed as free in-species already, can't bind to it!", bond_type.pair_1, head_index.index());
                assert!(tail_species.ports.map.get_mut(&bond_type.pair_2).unwrap().remove(&tail_index), "Port {} on index {} was not listed as free in-species already, can't bind to it!", bond_type.pair_2, tail_index.index());
                self.interactions.get_mut(&bond_type).unwrap().set.remove(&Rc::clone(&target_bond));
                let binary_binding_types: Vec<BondType> = self.interactions.keys().filter(|b| b.arity == InteractionArity::Binary && b.direction == InteractionDirection::Bind).cloned().collect();
                for some_bin_bond_type in binary_binding_types {
                    let some_uni_bond_type: BondType = BondType{arity: InteractionArity::Unary, ..some_bin_bond_type};
                    // binding opportunities that were binary, but will now be unary
                    let mut transformed_embeds: BTreeSet<Rc<RefCell<BondEmbed>>>;
                    transformed_embeds = self.interactions.get_mut(&some_bin_bond_type).unwrap().set.drain_filter(
                        |e| ( head_species.agent_set.contains(self.universe_graph.node_weight(e.borrow().a_index).unwrap()) && tail_species.agent_set.contains(self.universe_graph.node_weight(e.borrow().b_index).unwrap()) ) ||
                            ( head_species.agent_set.contains(self.universe_graph.node_weight(e.borrow().b_index).unwrap()) && tail_species.agent_set.contains(self.universe_graph.node_weight(e.borrow().a_index).unwrap()) )
                    ).collect();
                    for embed in transformed_embeds.iter() {
                        embed.borrow_mut().bond_type = some_uni_bond_type;
                    }
                    self.interactions.get_mut(&some_uni_bond_type).unwrap().set.append(&mut transformed_embeds);
                    // bindings no longer possible due to occupied sites
                    let bi_combinatorics_lost: BTreeSet<Rc<RefCell<BondEmbed>>>;
                    let un_combinatorics_lost: BTreeSet<Rc<RefCell<BondEmbed>>>;
                    (bi_combinatorics_lost, un_combinatorics_lost) = if some_bin_bond_type.pair_1 == some_bin_bond_type.pair_2 {
                        (
                            self.interactions.get_mut(&some_bin_bond_type).unwrap().set.drain_filter(|b| (b.borrow().a_index == head_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_1) || (b.borrow().b_index == tail_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_2) || (b.borrow().a_index == tail_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_2) || (b.borrow().b_index == head_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_1)).collect(),
                            self.interactions.get_mut(&some_uni_bond_type).unwrap().set.drain_filter(|b| (b.borrow().a_index == head_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_1) || (b.borrow().b_index == tail_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_2) || (b.borrow().a_index == tail_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_2) || (b.borrow().b_index == head_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_1)).collect()
                        )
                    } else {
                        (
                            self.interactions.get_mut(&some_bin_bond_type).unwrap().set.drain_filter(|b| (b.borrow().a_index == head_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_1) || (b.borrow().b_index == tail_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_2)).collect(),
                            self.interactions.get_mut(&some_uni_bond_type).unwrap().set.drain_filter(|b| (b.borrow().a_index == head_index && b.borrow().bond_type.pair_1 == target_bond.borrow().bond_type.pair_1) || (b.borrow().b_index == tail_index && b.borrow().bond_type.pair_2 == target_bond.borrow().bond_type.pair_2)).collect()
                        )
                    };
                }
            }
            // recast from "head & tail" to "host & eaten"; iterate over smaller sets while merging caches
            let (host_index, eaten_index) = if self.species_annots.get(&head_index).unwrap().borrow().agent_set.len() >= self.species_annots.get(&tail_index).unwrap().borrow().agent_set.len() {
                (head_index, tail_index)
            } else {
                (tail_index, head_index)
            };
            let host_species: Rc<RefCell<MixtureSpecies>> = Rc::clone(self.species_annots.get(&host_index).unwrap());
            let eaten_species: Rc<RefCell<MixtureSpecies>> = Rc::clone(self.species_annots.get(&eaten_index).unwrap());
            // create & update edge trackers
            let new_edge_index = self.universe_graph.add_edge(head_index, tail_index, false);
            bond_type.direction = InteractionDirection::Free;
            target_bond.borrow_mut().z = Some(new_edge_index);
            target_bond.borrow_mut().bond_type = bond_type;
            {
                let mut agent_a: RefMut<Agent> = self.universe_graph.node_weight(head_index).unwrap().borrow_mut();
                let site_a = agent_a.sites.get_mut(bond_type.pair_1.site).unwrap();
                match site_a {
                    Some(_) => panic!("Node {} was already bound at {}!", agent_a, bond_type.pair_1.site),
                    None => *site_a = Some(new_edge_index)
                }
                let mut agent_b: RefMut<Agent> = self.universe_graph.node_weight(tail_index).unwrap().borrow_mut();
                let site_b = agent_b.sites.get_mut(bond_type.pair_2.site).unwrap();
                match site_b {
                    Some(_) => panic!("Node {} was already bound at {}!", agent_b, bond_type.pair_2.site),
                    None => *site_b = Some(new_edge_index)
                }
            }
            host_species.borrow_mut().edges.insert(Rc::clone(&target_bond));
            self.edge_index_map.insert(new_edge_index, Rc::clone(&target_bond));
            self.interactions.get_mut(&bond_type).unwrap().set.insert(Rc::clone(&target_bond));
            // update species annotation tracker
            self.species_set.make_contiguous().sort();
            let spec_ix: usize = self.species_set.binary_search(&eaten_species).unwrap();
            self.species_set.remove(spec_ix);
            let indexes_of_eaten: BTreeSet<NodeIndex> = eaten_species.borrow().agent_set.iter().map(|a| a.borrow().id.unwrap()).collect();
            for agent_index in indexes_of_eaten {
                self.species_annots.entry(agent_index).and_modify(|e| {*e = Rc::clone(&host_species)});
            }
            host_species.borrow_mut().ports.update_from(&eaten_species.borrow_mut().ports);
            host_species.borrow_mut().edges.append(&mut eaten_species.borrow_mut().edges);
            host_species.borrow_mut().agent_set.append(&mut eaten_species.borrow_mut().agent_set);
            // update rule activities
            for i_data in self.interactions.values_mut() {
                i_data.update_activity(&self.species_annots);
            }
        }

        /**
         * Break a bond that is a cycle member; this does not partition one species into two, it
         * simply reduces the internal connectivity of said species.
        */
        pub fn unary_free(&mut self, target_bond: Rc<RefCell<BondEmbed<'a>>>) {
            let head_index: NodeIndex = target_bond.borrow().a_index;
            let tail_index: NodeIndex = target_bond.borrow().b_index;
            let bond_type: BondType<'a> = target_bond.borrow().bond_type;
            let edge_index: EdgeIndex = target_bond.borrow().z.unwrap();
            assert!(bond_type.arity == InteractionArity::Unary, "Incorrect arity in interaction type! Expected {}, got {}", InteractionArity::Unary, bond_type.arity);
            assert!(bond_type.direction == InteractionDirection::Free, "Incorrect direction in interaction type! Expected {}, got {}", InteractionDirection::Free, bond_type.direction);
            // sanity check for reslient bond
            assert!(self.universe_graph.edge_weight(edge_index).unwrap(), "This bond ({}) is not flagged as resilient; breaking it would not yield a still-connected component!", target_bond.borrow());
            {
                let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_index).unwrap().borrow_mut();
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
                    let mut collateral_agent_a: RefMut<Agent> = self.universe_graph.node_weight(swapped_edge.borrow().a_index).unwrap().borrow_mut();
                    let mut collateral_agent_b: RefMut<Agent> = self.universe_graph.node_weight(swapped_edge.borrow().b_index).unwrap().borrow_mut();
                    let collateral_site_a = collateral_agent_a.sites.get_mut(collateral_bond_type.pair_1.site).unwrap();
                    let collateral_site_b = collateral_agent_b.sites.get_mut(collateral_bond_type.pair_2.site).unwrap();
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

                {
                    // sanity checks for bond breaking; then update bond indexes on the agent
                    let mut node_a: RefMut<Agent> = self.universe_graph.node_weight(head_index).unwrap().borrow_mut();
                    let site_a = node_a.sites.get_mut(bond_type.pair_1.site).unwrap();
                    match site_a {
                        Some(some_edge) if *some_edge == edge_index => *site_a = None,
                        _ => panic!("This node {} was not already bound at expected bond {}!", node_a, edge_index.index())
                    };
                    let mut node_b: RefMut<Agent> = self.universe_graph.node_weight(tail_index).unwrap().borrow_mut();
                    let site_b = node_b.sites.get_mut(bond_type.pair_2.site).unwrap();
                    match site_b {
                        Some(some_edge) if *some_edge == edge_index => *site_b = None,
                        _ => panic!("This node {} was not already bound at expected bond {}!", node_b, edge_index.index())
                    };
                    assert!(target_species.edges.remove(&target_bond), "Bond {} was not already present as an available option in the species tracker!", target_bond.borrow());
                    assert!(self.interactions.get_mut(&bond_type).unwrap().set.remove(&target_bond), "Bond {} was not already present as an available option!", target_bond.borrow());
                    assert!(self.ports.map.get_mut(&bond_type.pair_1).unwrap().insert(head_index), "Site {} on node {} was already listed as free prior to unbinding!", bond_type.pair_1, node_a);
                    assert!(self.ports.map.get_mut(&bond_type.pair_2).unwrap().insert(tail_index), "Site {} on node {} was already listed as free prior to unbinding!", bond_type.pair_2, node_b);
                    assert!(target_species.ports.map.get_mut(&bond_type.pair_1).unwrap().insert(head_index), "Site {} on node {} was already listed as free prior to this unbinding on the species cache!", bond_type.pair_1, node_a);
                    assert!(target_species.ports.map.get_mut(&bond_type.pair_2).unwrap().insert(tail_index), "Site {} on node {} was already listed as free prior to this unbinding on the species cache!", bond_type.pair_2, node_b);
                }

                let unary_binding_types: Vec<BondType> = self.interactions.keys().filter(|b| b.arity == InteractionArity::Unary && b.direction == InteractionDirection::Bind).cloned().collect();
                for some_uni_bond_type in unary_binding_types {
                    let some_bin_bond_type: BondType = BondType{arity: InteractionArity::Binary, ..some_uni_bond_type};
                    // bindings now possible due to free'd sites
                    let mut heads_in_target: HashSet<NodeIndex> = 
                        if target_species.ports.map.contains_key(&some_uni_bond_type.pair_1) {
                            target_species.ports.map.get(&some_uni_bond_type.pair_1).unwrap().iter().copied().collect()
                        } else { HashSet::new() };
                    let mut tails_in_target: HashSet<NodeIndex> = 
                        if target_species.ports.map.contains_key(&some_uni_bond_type.pair_2) {
                            target_species.ports.map.get(&some_uni_bond_type.pair_2).unwrap().iter().copied().collect()
                        } else { HashSet::new() };
                    let heads_outside_target: HashSet<NodeIndex> = self.ports.map.get(&some_bin_bond_type.pair_1).unwrap().difference(&heads_in_target).copied().collect();
                    let tails_outside_target: HashSet<NodeIndex> = self.ports.map.get(&some_bin_bond_type.pair_2).unwrap().difference(&tails_in_target).copied().collect();
                    // binary combinatorics
                    let mut bin_combinatorics_gained: BTreeSet<Rc<RefCell<BondEmbed>>> = BTreeSet::new();
                    if bond_type.pair_1 == some_bin_bond_type.pair_1 {
                        for tail in tails_outside_target {
                            let bond_embed = BondEmbed{a_index: head_index, b_index: tail, bond_type: some_bin_bond_type, z: None};
                            bin_combinatorics_gained.insert(Rc::new(RefCell::new(bond_embed)));
                        }
                    }
                    if bond_type.pair_2 == some_bin_bond_type.pair_2 {
                        for head in heads_outside_target {
                            let bond_embed = BondEmbed{a_index: head, b_index: tail_index, bond_type: some_bin_bond_type, z: None};
                            bin_combinatorics_gained.insert(Rc::new(RefCell::new(bond_embed)));
                        }
                    }
                    // unary combinatorics
                    let mut uni_combinatorics_gained: BTreeSet<Rc<RefCell<BondEmbed>>> = BTreeSet::new();
                    if bond_type.pair_1 == some_uni_bond_type.pair_1 {
                        if bond_type.pair_1 == bond_type.pair_2 {tails_in_target.remove(&head_index);}
                        if bond_type.pair_1.agent == bond_type.pair_2.agent {tails_in_target.remove(&head_index);}
                        for tail in tails_in_target {
                            let bond_embed = BondEmbed{a_index: head_index, b_index: tail, bond_type: some_uni_bond_type, z: None};
                            uni_combinatorics_gained.insert(Rc::new(RefCell::new(bond_embed)));
                        }
                    }
                    if bond_type.pair_2 == some_uni_bond_type.pair_2 {
                        if bond_type.pair_1 == bond_type.pair_2 {heads_in_target.remove(&tail_index);}
                        if bond_type.pair_1.agent == bond_type.pair_2.agent {heads_in_target.remove(&tail_index);}
                        for head in heads_in_target {
                            let bond_embed = BondEmbed{a_index: head, b_index: tail_index, bond_type: some_uni_bond_type, z: None};
                            uni_combinatorics_gained.insert(Rc::new(RefCell::new(bond_embed)));
                        }
                    }
                    self.interactions.get_mut(&some_bin_bond_type).unwrap().set.append(&mut bin_combinatorics_gained);
                    self.interactions.get_mut(&some_uni_bond_type).unwrap().set.append(&mut uni_combinatorics_gained);
                }
            }

            // update resiliency
            // move things from self.interactions that match Free & Unary (aka cycle edges) into
            // the corresponding
            // things in self.interactions that match Free & Binary (aka tree edges)
            let is_newly_treed_bond = |boxed_bond: &Rc<RefCell<BondEmbed>>| -> bool {
                let bond_embed = boxed_bond.borrow();
                if *self.universe_graph.edge_weight(bond_embed.z.unwrap()).unwrap() {
                    let mut counterfactual_graph = self.universe_graph.clone();
                    counterfactual_graph.remove_edge(bond_embed.z.unwrap());
                    let path = astar(&counterfactual_graph, bond_embed.a_index, |finish| finish == bond_embed.b_index, |_| 1, |_| 0);
                    path.is_none()
                } else {false}      // if the bond's weight was already false, then it was not tree'd
            };
            // double for-loop to maintain single-mutable-borrow of self.interactions
            let mut temp_container: BTreeMap<BondType, BTreeSet<Rc<RefCell<BondEmbed>>>> = BTreeMap::new();
            for (i_type, i_data) in self.interactions.iter_mut() {
                if i_type.arity == InteractionArity::Unary && i_type.direction == InteractionDirection::Free {
                    let newly_treed_bonds: BTreeSet<Rc<RefCell<BondEmbed>>> = i_data.set.drain_filter(is_newly_treed_bond).collect();
                    let target_type = BondType{arity: InteractionArity::Binary, ..*i_type};
                    for this_embed in &newly_treed_bonds {
                        this_embed.borrow_mut().bond_type = target_type;
                    }
                    temp_container.insert(target_type, newly_treed_bonds);
                }
            }
            for (target_type, mut conv_data) in temp_container {
                for conv_bond in &conv_data {
                    let edge_weight = self.universe_graph.edge_weight_mut(conv_bond.borrow().z.unwrap()).unwrap();
                    *edge_weight = false;
                }
                self.interactions.get_mut(&target_type).unwrap().set.append(&mut conv_data);
            }
            // update rule activities
            for i_data in self.interactions.values_mut() {
                i_data.update_activity(&self.species_annots);
            }
        }

        /**
         * Break a bond that is not a cycle member; this partitions one species into two.
        */
        pub fn binary_free(&mut self, target_bond: Rc<RefCell<BondEmbed<'a>>>) {
            let head_index: NodeIndex = target_bond.borrow().a_index;
            let tail_index: NodeIndex = target_bond.borrow().b_index;
            let bond_type: BondType<'a> = target_bond.borrow().bond_type;
            let edge_index: EdgeIndex = target_bond.borrow().z.unwrap();
            assert!(bond_type.arity == InteractionArity::Binary, "Incorrect arity in interaction type! Expected {}, got {}", InteractionArity::Binary, bond_type.arity);
            assert!(bond_type.direction == InteractionDirection::Free, "Incorrect direction in interaction type! Expected {}, got {}", InteractionDirection::Free, bond_type.direction);
            // sanity check for non-reslient bond
            assert!(! self.universe_graph.edge_weight(edge_index).unwrap(), "This bond ({}) is flagged as resilient; breaking it would yield a still-connected component!", target_bond.borrow());
            {
                // operations on target_species
                let mut target_species: RefMut<MixtureSpecies> = self.species_annots.get(&head_index).unwrap().borrow_mut();
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
                    let mut collateral_agent_a: RefMut<Agent> = self.universe_graph.node_weight(swapped_edge.borrow().a_index).unwrap().borrow_mut();
                    let mut collateral_agent_b: RefMut<Agent> = self.universe_graph.node_weight(swapped_edge.borrow().b_index).unwrap().borrow_mut();
                    let collateral_site_a = collateral_agent_a.sites.get_mut(collateral_bond_type.pair_1.site).unwrap();
                    let collateral_site_b = collateral_agent_b.sites.get_mut(collateral_bond_type.pair_2.site).unwrap();
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
                let mut node_a: RefMut<Agent> = self.universe_graph.node_weight(head_index).unwrap().borrow_mut();
                let site_a = node_a.sites.get_mut(bond_type.pair_1.site).unwrap();
                match site_a {
                    Some(some_edge) if *some_edge == edge_index => *site_a = None,
                    _ => panic!("This node {} was not already bound at expected bond {}!", node_a, edge_index.index())
                };
                let mut node_b: RefMut<Agent> = self.universe_graph.node_weight(tail_index).unwrap().borrow_mut();
                let site_b = node_b.sites.get_mut(bond_type.pair_2.site).unwrap();
                match site_b {
                    Some(some_edge) if *some_edge == edge_index => *site_b = None,
                    _ => panic!("This node {} was not already bound at expected bond {}!", node_b, edge_index.index())
                };
                assert!(target_species.edges.remove(&target_bond), "Bond {} was not already present as an available option in the species tracker!", target_bond.borrow());
                assert!(self.interactions.get_mut(&bond_type).unwrap().set.remove(&target_bond), "Bond {} was not already present as an available option!", target_bond.borrow());
                assert!(self.ports.map.get_mut(&bond_type.pair_1).unwrap().insert(head_index), "Site {} on node {} was already listed as free prior to unbinding!", bond_type.pair_1, node_a);
                assert!(self.ports.map.get_mut(&bond_type.pair_2).unwrap().insert(tail_index), "Site {} on node {} was already listed as free prior to unbinding!", bond_type.pair_2, node_b);
                assert!(target_species.ports.map.get_mut(&bond_type.pair_1).unwrap().insert(head_index), "Site {} on node {} was already listed as free prior to this unbinding on the species cache!", bond_type.pair_1, node_a);
                assert!(target_species.ports.map.get_mut(&bond_type.pair_2).unwrap().insert(tail_index), "Site {} on node {} was already listed as free prior to this unbinding on the species cache!", bond_type.pair_2, node_b);
            }
            // discover what node indexes ended where; the loop terminates once the smallest graph is mapped-out
            let mut dfs_head = Dfs::new(&self.universe_graph, head_index);
            let mut dfs_tail = Dfs::new(&self.universe_graph, tail_index);
            let mut head_graph_indexes: BTreeSet<NodeIndex> = BTreeSet::new();
            let mut tail_graph_indexes: BTreeSet<NodeIndex> = BTreeSet::new();
            let head_graph_agents: BTreeSet<Rc<RefCell<Agent<'a>>>>;    // created from iterator at below loop's termination
            let tail_graph_agents: BTreeSet<Rc<RefCell<Agent<'a>>>>;    // created from iterator at below loop's termination
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
                        head_graph_agents = head_graph_indexes.iter().map(|i| Rc::clone(self.universe_graph.node_weight(*i).unwrap())).collect::<BTreeSet<Rc<RefCell<Agent>>>>();
                        tail_graph_agents = self.species_annots.get(&tail_index).unwrap().borrow().agent_set.difference(&head_graph_agents).cloned().collect();
                        tail_graph_indexes = tail_graph_agents.iter().map(|i| i.borrow().id.unwrap()).collect();
                        break
                    },
                    (Some(_), None) => {
                        // we've mapped-out the tail node's graph
                        tail_graph_agents = tail_graph_indexes.iter().map(|i| Rc::clone(self.universe_graph.node_weight(*i).unwrap())).collect::<BTreeSet<Rc<RefCell<Agent>>>>();
                        head_graph_agents = self.species_annots.get(&head_index).unwrap().borrow().agent_set.difference(&tail_graph_agents).cloned().collect();
                        head_graph_indexes = head_graph_agents.iter().map(|i| i.borrow().id.unwrap()).collect();
                        break
                    },
                    (None, None) => {
                        head_graph_agents = head_graph_indexes.iter().map(|i| Rc::clone(self.universe_graph.node_weight(*i).unwrap())).collect::<BTreeSet<Rc<RefCell<Agent>>>>();
                        tail_graph_agents = tail_graph_indexes.iter().map(|i| Rc::clone(self.universe_graph.node_weight(*i).unwrap())).collect::<BTreeSet<Rc<RefCell<Agent>>>>();
                        break
                    }
                }
            }
            let (ejected_agents, ejected_indexes, retained_mark, retained_agents, retained_indexes) =  // Iterate over small & lookup over big
                if head_graph_agents.len() >= tail_graph_agents.len()
                    {(head_graph_agents, head_graph_indexes, tail_index, tail_graph_agents, tail_graph_indexes)}
                else 
                    {(tail_graph_agents, tail_graph_indexes, head_index, head_graph_agents, head_graph_indexes)};
            // create new species cache, in-place modify old species cache
            let ejected_ports: OpenPorts = self.species_annots.get(&retained_mark).unwrap().borrow_mut().ports.eject_where(&ejected_indexes);
            let ejected_edges: BTreeSet<Rc<RefCell<BondEmbed>>> = self.species_annots.get(&retained_mark).unwrap().borrow_mut().edges.drain_filter(|b| ejected_indexes.contains(&b.borrow().a_index)).collect();
            self.species_annots.get(&retained_mark).unwrap().borrow_mut().agent_set = retained_agents;
            let new_species = Rc::new(RefCell::new(MixtureSpecies{
                ports: ejected_ports,
                edges: ejected_edges,
                agent_set: ejected_agents
            }));
            // update the species trackers, per node index & global
            for node_index in &ejected_indexes {
                self.species_annots.entry(*node_index).and_modify(|e| {*e = Rc::clone(&new_species)});
            }
            self.species_set.push_back(Rc::clone(&new_species));
            // rephrasing operations to use types: head & tail
            let head_species: Ref<MixtureSpecies> = self.species_annots.get(&head_index).unwrap().borrow();
            let tail_species: Ref<MixtureSpecies> = self.species_annots.get(&tail_index).unwrap().borrow();
            let unary_binding_types: Vec<BondType> = self.interactions.keys().filter(|b| b.arity == InteractionArity::Unary && b.direction == InteractionDirection::Bind).cloned().collect();
            for some_uni_bond_type in unary_binding_types {
                let some_bin_bond_type: BondType = BondType{arity: InteractionArity::Binary, ..some_uni_bond_type};
                // binding opportunities that were unary, but will now be binary
                let mut transformed_embeds: BTreeSet<Rc<RefCell<BondEmbed>>>;
                transformed_embeds = self.interactions.get_mut(&some_uni_bond_type).unwrap().set.drain_filter(
                    |e| ( new_species.borrow().agent_set.contains(self.universe_graph.node_weight(e.borrow().a_index).unwrap()) && self.species_annots.get(&retained_mark).unwrap().borrow().agent_set.contains(self.universe_graph.node_weight(e.borrow().b_index).unwrap()) ) ||
                        ( new_species.borrow().agent_set.contains(self.universe_graph.node_weight(e.borrow().b_index).unwrap()) && self.species_annots.get(&retained_mark).unwrap().borrow().agent_set.contains(self.universe_graph.node_weight(e.borrow().a_index).unwrap()) )
                ).collect();
                for embed in transformed_embeds.iter() {
                    embed.borrow_mut().bond_type = some_bin_bond_type;
                }
                self.interactions.get_mut(&some_bin_bond_type).unwrap().set.append(&mut transformed_embeds);
                // bindings now possible due to free'd sites
                let mut heads_in_tail: HashSet<NodeIndex> = 
                    if tail_species.ports.map.contains_key(&some_uni_bond_type.pair_1) {
                        tail_species.ports.map.get(&some_uni_bond_type.pair_1).unwrap().iter().copied().collect()
                    } else { HashSet::new() };
                let mut tails_in_head: HashSet<NodeIndex> = 
                    if head_species.ports.map.contains_key(&some_uni_bond_type.pair_2) {
                        head_species.ports.map.get(&some_uni_bond_type.pair_2).unwrap().iter().copied().collect()
                    } else { HashSet::new() };
                let heads_outside_tail: HashSet<NodeIndex> = self.ports.map.get(&some_bin_bond_type.pair_1).unwrap().difference(&heads_in_tail).copied().collect();
                let tails_outside_head: HashSet<NodeIndex> = self.ports.map.get(&some_bin_bond_type.pair_2).unwrap().difference(&tails_in_head).copied().collect();
                // binary combinatorics
                let mut bin_combinatorics_gained: BTreeSet<Rc<RefCell<BondEmbed>>> = BTreeSet::new();
                if bond_type.pair_1 == some_bin_bond_type.pair_1 {
                    for tail in tails_outside_head {
                        let bond_embed = BondEmbed{a_index: head_index, b_index: tail, bond_type: some_bin_bond_type, z: None};
                        bin_combinatorics_gained.insert(Rc::new(RefCell::new(bond_embed)));
                    }
                }
                if bond_type.pair_2 == some_bin_bond_type.pair_2 {
                    for head in heads_outside_tail {
                        let bond_embed = BondEmbed{a_index: head, b_index: tail_index, bond_type: some_bin_bond_type, z: None};
                        bin_combinatorics_gained.insert(Rc::new(RefCell::new(bond_embed)));
                    }
                }
                // unary combinatorics
                let mut uni_combinatorics_gained: BTreeSet<Rc<RefCell<BondEmbed>>> = BTreeSet::new();
                if bond_type.pair_1 == some_uni_bond_type.pair_1 {
                    tails_in_head.remove(&head_index);
                    for tail in tails_in_head {
                        let bond_embed = BondEmbed{a_index: head_index, b_index: tail, bond_type: some_uni_bond_type, z: None};
                        assert!(uni_combinatorics_gained.insert(Rc::new(RefCell::new(bond_embed))), "Embed {} already cached into the unary-gained temp. tracker!", bond_embed);
                    }
                }
                if bond_type.pair_2 == some_uni_bond_type.pair_2 {
                    heads_in_tail.remove(&tail_index);
                    for head in heads_in_tail {
                        let bond_embed = BondEmbed{a_index: head, b_index: tail_index, bond_type: some_uni_bond_type, z: None};
                        assert!(uni_combinatorics_gained.insert(Rc::new(RefCell::new(bond_embed))), "Embed {} already cached into the unary-gained temp. tracker!", bond_embed);
                    }
                }
                self.interactions.get_mut(&some_bin_bond_type).unwrap().set.append(&mut bin_combinatorics_gained);
                self.interactions.get_mut(&some_uni_bond_type).unwrap().set.append(&mut uni_combinatorics_gained);
            }
            // update rule activities
            for i_data in self.interactions.values_mut() {
                i_data.update_activity(&self.species_annots);
            }
        }

        /** 
         * Pretty print the mixture's rule activities. 
        */
        pub fn print_activities(&self) {
            let mut str_vec: Vec<String> = Vec::new();
            for (i_name, i_data) in &self.interactions {
                str_vec.push(format!("{:.4e}\t\t{}", i_data.activity_calculated, i_name))
            }
            println!("Rule activities:\n{}", str_vec.join("\n"))
        }

        /**
         * Pretty print the mixture's possible transformations (aka embeddings).
        */
        pub fn print_all_transformation_posibilities(&self) {
            let str_vec: Vec<String> = self.interactions.iter().map(|i| format!("{} ({:.4e}):\n{}", i.0, i.1.activity_calculated, i.1)).collect();
            println!("{}", str_vec.join("\n\n"))
        }

        /**
         * Pretty print the mixture's possible unary binding transformations (aka embeddings).
        */
        pub fn print_unary_bind_transformation_posibilities(&self) {
            let str_vec: Vec<String> = self.interactions.iter().filter(|n| n.0.arity == InteractionArity::Unary && n.0.direction == InteractionDirection::Bind).map(|i| format!("{}", i.1)).collect();
            println!("{}", str_vec.join("\n\n"))
        }

        /**
         * Pretty print the mixture's possible unary free'ing transformations (aka embeddings).
        */
        pub fn print_unary_free_transformation_posibilities(&self) {
            let str_vec: Vec<String> = self.interactions.iter().filter(|n| n.0.arity == InteractionArity::Unary && n.0.direction == InteractionDirection::Free).map(|i| format!("{}", i.1)).collect();
            println!("{}", str_vec.join("\n\n"))
        }

        /**
         * Pretty print the mixture's possible binary binding transformations (aka embeddings).
        */
        pub fn print_binary_bind_transformation_posibilities(&self) {
            let str_vec: Vec<String> = self.interactions.iter().filter(|n| n.0.arity == InteractionArity::Binary && n.0.direction == InteractionDirection::Bind).map(|i| format!("{}", i.1)).collect();
            println!("{}", str_vec.join("\n\n"))
        }

        /**
         * Pretty print the mixture's possible binary free'ing transformations (aka embeddings).
        */
        pub fn print_binary_free_transformation_posibilities(&self) {
            let str_vec: Vec<String> = self.interactions.iter().filter(|n| n.0.arity == InteractionArity::Binary && n.0.direction == InteractionDirection::Free).map(|i| format!("{}", i.1)).collect();
            println!("{}", str_vec.join("\n\n"))
        }

        /**
         * Return number of simulated events so far.
        */
        pub fn event_number(&self) -> usize {
            self.simulated_events
        }

        /**
         * Return amount of simulated time so far.
        */
        pub fn simulated_time(&self) -> f64 {
            self.simulated_time
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::building_blocks::{AllInteractionData, ProtomerResources, RateFudge};
    
    use petgraph::prelude::{NodeIndex, EdgeIndex};
    use std::{rc::Rc, cell::RefCell};

    #[test]
    fn read_protomer_resources() {
        let r = ProtomerResources::from_json("examples/APC_HUMAN--AXIN1_HUMAN--intestinal_crypt/Abundances.json").unwrap();
        assert_eq!(*r.map.get("Axin").unwrap(), 135000);
        assert_eq!(*r.map.get("APC").unwrap(), 270000);
    }

    #[test]
    fn read_all_interaction_data() {
        let i = AllInteractionData::from_json("examples/APC_HUMAN--AXIN1_HUMAN--intestinal_crypt/Interactions.json").unwrap();
        assert_eq!(i.interactions.len(), 5);
        assert_eq!(i.interactions[0].agent_a, "Axin");
        assert_eq!(i.interactions[0].agent_b, "Axin");
        assert_eq!(i.interactions[0].site_a, "head");
        assert_eq!(i.interactions[0].site_b, "tail");
        assert_eq!(i.interactions[0].name, Some(String::from("head to tail polymerization")));
        assert_eq!(i.interactions[0].free_unary.base, 1.0e-2);
        assert_eq!(i.interactions[0].free_unary.multiplicative, None);
        assert_eq!(i.interactions[0].free_unary.additive, None);
        assert_eq!(i.interactions[0].bind_unary, RateFudge{base: 1.0e-2, multiplicative: None, additive: None});
    }
    
}