use axin_apc_simulator::{rule_activities::RuleRates, reaction_mixture::Mixture};


fn main() {
    let my_rates = RuleRates {
        axn_axn_u_bind: 1.0e-2,
        axn_axn_b_bind: 1.0e-4,
        axn_axn_u_free: 1.0e-2,
        axn_axn_b_free: 1.0e-2,
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
    let mut my_mix = Mixture::new_from_monomers(1000, 0, my_rates);
    for _ in 0..1000 {
        my_mix.choose_and_apply_next_rule();
    }
    println!("{}", my_mix.to_kappa());
}