use axin_apc_simulator::{rule_activities::RuleRates, reaction_mixture::Mixture};
use petgraph::prelude::*;


fn main() {
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
    let mut my_mix = Mixture::new_from_monomers(3, 2, my_rates);
    my_mix.axn_axn_binary_bind(NodeIndex::from(0), NodeIndex::from(1));
    println!("{}", my_mix.to_kappa());
    my_mix.axn_axn_binary_bind(NodeIndex::from(1), NodeIndex::from(2));
    println!("{}", my_mix.to_kappa());
    my_mix.axn_axn_unary_bind(NodeIndex::from(2), NodeIndex::from(0));
    println!("{}", my_mix.to_kappa());
}