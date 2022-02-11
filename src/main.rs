use clap::Parser;
use axin_apc_simulator::{rule_activities::RuleRates, reaction_mixture::Mixture};


#[derive(Parser, Debug)]
#[clap(author, version, about, long_about=None)]
struct Args {
    /// Number of Axn-type protomers in the initial mix.
    #[clap(short='x', long, default_value_t = 0)]
    axn: usize,

    /// Number of APC-type protomers in the initial mix.
    #[clap(short='p', long, default_value_t = 0)]
    apc: usize,

    /// Number of events to simulate.
    #[clap(short, long, default_value_t = 10)]
    events: usize,

    /// Axn-Axn unary binding rate constant.
    #[clap(long, default_value_t = 0.0)]
    axn_axn_u_bind: f64,

    /// Axn_Axn binary binding rate constant.
    #[clap(long, default_value_t = 0.0)]
    axn_axn_b_bind: f64,

    /// Axn-Axn unary unbinding rate constant.
    #[clap(long, default_value_t = 0.0)]
    axn_axn_u_free: f64,

    /// Axn_Axn binary unbinding rate constant.
    #[clap(long, default_value_t = 0.0)]
    axn_axn_b_free: f64
}

fn main() {
    let args = Args::parse();

    let my_rates = RuleRates {
        axn_axn_u_bind: args.axn_axn_u_bind,
        axn_axn_b_bind: args.axn_axn_b_bind,
        axn_axn_u_free: args.axn_axn_u_free,
        axn_axn_b_free: args.axn_axn_b_free,
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
    let mut my_mix = Mixture::new_from_monomers(args.axn, args.apc, my_rates);
    for _ in 0..args.events {
        my_mix.choose_and_apply_next_rule();
    }
    println!("{}", my_mix.to_kappa());
}