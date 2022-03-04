use axin_apc_simulator::building_blocks::{AllInteractionData, ProtomerResources};
use axin_apc_simulator::reaction_mixture::Mixture;
use clap::Parser;


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
struct Args {
    /// Path to file with protomer abundances
    #[clap(short='a', long, default_value = "./examples/APC_HUMAN--AXIN1_HUMAN--intestinal_crypt/Abundances_toy.json")]
    abundances_path: String,

    /// Path to file with interaction data
    #[clap(short='i', long, default_value = "./examples/APC_HUMAN--AXIN1_HUMAN--intestinal_crypt/Interactions.json")]
    interactions_path: String,

    /// Number of events to simulate
    #[clap(short='e', long, default_value_t = 1)]
    requested_events: usize
}

fn main() {
    let args = Args::parse();
    let raw_interactions: AllInteractionData = AllInteractionData::from_json(args.interactions_path).unwrap();
    let raw_abundances: ProtomerResources = ProtomerResources::from_json(args.abundances_path).unwrap();
   
    let mut my_mix = Mixture::new_monomeric_from(&raw_abundances, &raw_interactions);
    println!("{}", my_mix.to_kappa());
    my_mix.choose_and_apply_next_rule();
    println!("{}", my_mix.to_kappa());
    
    
    //for _ in 0..args.requested_events {
    //    my_mix.choose_and_apply_next_rule();
    //}
    //println!("{}", my_mix.to_kappa());
}