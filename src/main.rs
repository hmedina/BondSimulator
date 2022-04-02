use axin_apc_simulator::building_blocks::{AllInteractionData, ProtomerResources};
use axin_apc_simulator::reaction_mixture::Mixture;
use std::fs;
use std::io::ErrorKind;
use std::time::SystemTime;
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

    /// Maximum number of events to simulate
    #[clap(short='e', long)]
    requested_events: Option<usize>,

    /// Maximum time to simulate
    #[clap(short='t', long)]
    requested_time: Option<f64>,

    /// Path to file where snapshots should be saved
    #[clap(short='o', long, default_value = "./output")]
    snapshot_dir: String,

    /// Snapshot-taking period, in time units
    #[clap(long)]
    snap_time_p: Option<f64>,

    /// Snapshot-taking period, in event number
    #[clap(long)]
    snap_event_p: Option<usize>,
}

fn main() {
    let args = Args::parse();
    let raw_interactions: AllInteractionData = AllInteractionData::from_json(args.interactions_path).unwrap();
    let raw_abundances: ProtomerResources = ProtomerResources::from_json(args.abundances_path).unwrap();
    match fs::create_dir(&args.snapshot_dir) {
        Err(why) => match why.kind() {
            ErrorKind::AlreadyExists => {println!("Directory {} already exists; reusing it.", args.snapshot_dir)},
            _ => println!("Could not create directory {}: {}", args.snapshot_dir, why)
        }
        Ok(_) => {}
    }

    match (args.requested_events, args.requested_time) {
        (None, Some(m_time)) => {
            let init_time = SystemTime::now();
            println!("Starting initialization...");
            let mut my_mix = Mixture::new_monomeric_from(&raw_abundances, &raw_interactions);
            println!("Initialization complete! Duration was {} seconds; simulation starting...", init_time.elapsed().unwrap().as_secs());
            my_mix.simulate_up_to_time(m_time, args.snap_time_p, Some(&args.snapshot_dir));
            println!("Simulation complete! Duration was {} seconds.", init_time.elapsed().unwrap().as_secs());
        },
        (Some(m_event), None) => {
            let init_time = SystemTime::now();
            println!("Starting initialization...");
            let mut my_mix = Mixture::new_monomeric_from(&raw_abundances, &raw_interactions);
            println!("Initialization complete! Duration was {} seconds; simulation starting...", init_time.elapsed().unwrap().as_secs());
            my_mix.simulate_up_to_event(m_event, args.snap_event_p, Some(&args.snapshot_dir));
            println!("Simulation complete! Duration was {} seconds.", init_time.elapsed().unwrap().as_secs());
        },
        (None, None) => {println!("Error; supply a run-length of either time (-t) or events (-e).")},
        (Some(_), Some(_)) => {panic!("Can't run in both time and event mode! Supply one or the other.")},
    }
}