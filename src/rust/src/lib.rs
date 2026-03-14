roxido_registration!();
use roxido::*;

use slice_sampler::univariate::doubling::{
    TuningParameters as DoublingParameters, univariate_slice_sampler_doubling_and_shrinkage,
};
use slice_sampler::univariate::stepping_out::{
    TuningParameters as SteppingOutParameters, univariate_slice_sampler_stepping_out_and_shrinkage,
};

#[roxido]
fn univariate_slice_sampler_stepping_out(
    x: f64,
    target: &RFunction,
    on_log_scale: bool,
    width: f64,
) {
    let tuning_parameters = SteppingOutParameters::new().width(width);
    let result = univariate_slice_sampler_stepping_out_and_shrinkage(
        x,
        |x| target.call1(x.to_r(pc), pc).stop().as_scalar().stop().f64(),
        on_log_scale,
        &tuning_parameters,
        &mut None,
    );
    let rval = RList::with_names(&["value", "n_evaluations"], pc);
    rval.set(0, result.0.to_r(pc)).stop();
    rval.set(1, (result.1 as i32).to_r(pc)).stop();
    rval
}

#[roxido]
fn univariate_slice_sampler_doubling(x: f64, target: &RFunction, on_log_scale: bool, width: f64) {
    let tuning_parameters = DoublingParameters::new().width(width);
    let result = univariate_slice_sampler_doubling_and_shrinkage(
        x,
        |x| target.call1(x.to_r(pc), pc).stop().as_scalar().stop().f64(),
        on_log_scale,
        &tuning_parameters,
        &mut None,
    );
    let rval = RList::with_names(&["value", "n_evaluations"], pc);
    rval.set(0, result.0.to_r(pc)).stop();
    rval.set(1, (result.1 as i32).to_r(pc)).stop();
    rval
}
