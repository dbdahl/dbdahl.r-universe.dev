roxido_registration!();

use crp_distr::{log_crp_integer_partition_pmf_with_log_s, sample_crp_integer_partition};
use lorenz_ip_distr::{
    CubicSpline, Lorenz, LorenzFromPartition, LorenzISpline, LorenzIp, LorenzLinear,
    TadpoleLorenzIpBuilder, TidalLorenzIpBuilder,
    counting::{
        bell_number_log, enumerate_integer_partitions as rust_enumerate_integer_partitions,
        enumerate_set_partitions, integer_partition_count_log,
        integer_partition_count_with_parts_log, multinomial_partition_log,
        stirling_second_kind_log,
    },
    gini_coefficient as rust_gini_coefficient,
};
use rand::{RngExt, SeedableRng};
use rand_pcg::Pcg64Mcg;
use roxido::*;
use std::collections::BTreeMap;

// ============================================================================
// Lorenz curve external pointer support
// ============================================================================

/// Create a LorenzLinear from partition and wrap in external pointer.
fn new_lorenz_linear(partition: &[f64]) -> Result<LorenzLinear, lorenz_ip_distr::LorenzError> {
    // Sort partition (poorest to richest)
    let mut sorted = partition.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    LorenzLinear::from_partition(&sorted)
}

/// Create a LorenzISpline from partition with parameters and wrap in external pointer.
fn new_lorenz_ispline(
    partition: &[f64],
    n_interior_knots: Option<usize>,
    degree: usize,
    lambda: f64,
) -> Result<LorenzISpline, lorenz_ip_distr::LorenzError> {
    // Sort partition (poorest to richest)
    let mut sorted = partition.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Auto-compute n_interior_knots if not specified
    let knots = n_interior_knots.unwrap_or_else(|| {
        let k = sorted.len();
        (k.saturating_sub(2)).clamp(1, 10)
    });

    LorenzISpline::from_partition_with(&sorted, knots, degree, lambda)
}

#[roxido]
fn new_lorenz_linear_extptr(partition: &RVector) {
    let partition_f64 = partition.to_f64(pc);
    let partition_vec: Vec<f64> = partition_f64.slice().to_vec();

    if partition_vec.is_empty() {
        stop!("Partition cannot be empty.");
    }
    if partition_vec.iter().any(|&x| x <= 0.0) {
        stop!("All partition elements must be positive.");
    }

    let curve = new_lorenz_linear(&partition_vec).stop();

    // Store partition in tag for reconstruction after serialization
    let tag = RVector::from_slice(&partition_vec, pc);
    let ptr = RExternalPtr::encode_full(curve, tag, true, pc);
    ptr.set_class(["lorenz_linear", "lorenz"].to_r(pc));
    ptr
}

#[roxido]
fn new_lorenz_ispline_extptr(
    partition: &RVector,
    n_interior_knots: &RObject,
    degree: i32,
    lambda: f64,
) {
    let partition_f64 = partition.to_f64(pc);
    let partition_vec: Vec<f64> = partition_f64.slice().to_vec();

    if partition_vec.is_empty() {
        stop!("Partition cannot be empty.");
    }
    if partition_vec.iter().any(|&x| x <= 0.0) {
        stop!("All partition elements must be positive.");
    }
    if degree < 0 {
        stop!("Degree must be non-negative.");
    }
    if lambda < 0.0 {
        stop!("Lambda must be non-negative.");
    }

    let knots_opt: Option<usize> = if n_interior_knots.is_null() {
        None
    } else {
        let knots_vec = n_interior_knots
            .as_vector()
            .stop_str("n_interior_knots should be a scalar");
        let knots_i32 = knots_vec.to_i32(pc);
        let k = knots_i32.get(0).stop();
        if k < 1 {
            stop!("n_interior_knots must be at least 1.");
        }
        Some(k as usize)
    };

    let curve = new_lorenz_ispline(&partition_vec, knots_opt, degree as usize, lambda).stop();

    // Store partition and parameters in tag for reconstruction
    // Tag format: [partition..., n_interior_knots_or_nan, degree, lambda]
    let k = partition_vec.len();
    let tag = RVector::<f64>::new(k + 3, pc);
    for (i, &v) in partition_vec.iter().enumerate() {
        tag.set(i, v).stop();
    }
    tag.set(k, knots_opt.map(|x| x as f64).unwrap_or(f64::NAN))
        .stop();
    tag.set(k + 1, degree as f64).stop();
    tag.set(k + 2, lambda).stop();

    let ptr = RExternalPtr::encode_full(curve, tag, true, pc);
    ptr.set_class(["lorenz_ispline", "lorenz"].to_r(pc));
    ptr
}

/// Resurrect a LorenzLinear from an external pointer with validation.
fn resurrect_lorenz_linear(params: &mut RExternalPtr) {
    let classes = params.get_class();
    if classes.get(0) != Ok("lorenz_linear") {
        stop!("Expected a 'lorenz_linear' object.");
    }

    let f = |tag: &RObject| {
        let vec = tag.as_vector().stop();
        let float_vec = vec.as_f64().stop();
        let partition: Vec<f64> = float_vec.slice().to_vec();
        new_lorenz_linear(&partition).stop()
    };
    params.reencode(f);
}

/// Resurrect a LorenzISpline from an external pointer with validation.
fn resurrect_lorenz_ispline(params: &mut RExternalPtr) {
    let classes = params.get_class();
    if classes.get(0) != Ok("lorenz_ispline") {
        stop!("Expected a 'lorenz_ispline' object.");
    }

    let f = |tag: &RObject| {
        let vec = tag.as_vector().stop();
        let float_vec = vec.as_f64().stop();
        let slice = float_vec.slice();
        if slice.len() < 4 {
            stop!("Invalid tag structure in external pointer.");
        }
        let k = slice.len() - 3;
        let partition: Vec<f64> = slice[..k].to_vec();
        let n_interior_knots = if slice[k].is_nan() {
            None
        } else {
            Some(slice[k] as usize)
        };
        let degree = slice[k + 1] as usize;
        let lambda = slice[k + 2];
        new_lorenz_ispline(&partition, n_interior_knots, degree, lambda).stop()
    };
    params.reencode(f);
}

// ============================================================================
// Cubic spline external pointer support
// ============================================================================

/// Create a CubicSpline from x and y vectors.
fn new_cubic_spline(
    x: &[f64],
    y: &[f64],
) -> Result<CubicSpline, lorenz_ip_distr::CubicSplineError> {
    CubicSpline::new(x.to_vec(), y.to_vec())
}

#[roxido]
fn new_cubic_spline_extptr(x: &RVector, y: &RVector) {
    let x_f64 = x.to_f64(pc);
    let y_f64 = y.to_f64(pc);
    let x_vec: Vec<f64> = x_f64.slice().to_vec();
    let y_vec: Vec<f64> = y_f64.slice().to_vec();

    if x_vec.len() != y_vec.len() {
        stop!("x and y must have the same length.");
    }
    if x_vec.is_empty() {
        stop!("x and y cannot be empty.");
    }

    let spline = new_cubic_spline(&x_vec, &y_vec).stop();

    // Store x and y in tag for reconstruction after serialization
    let tag = RList::new(2, pc);
    tag.set(0, RVector::from_slice(&x_vec, pc)).stop();
    tag.set(1, RVector::from_slice(&y_vec, pc)).stop();

    let ptr = RExternalPtr::encode_full(spline, tag, true, pc);
    ptr.set_class(["cubic_spline"].to_r(pc));
    ptr
}

/// Resurrect a CubicSpline from an external pointer with validation.
fn resurrect_cubic_spline(params: &mut RExternalPtr) {
    let classes = params.get_class();
    if classes.get(0) != Ok("cubic_spline") {
        stop!("Expected a 'cubic_spline' object.");
    }

    let f = |tag: &RObject| {
        let list = tag.as_list().stop();
        let x_vec = list.get(0).stop().as_vector().stop().as_f64().stop();
        let y_vec = list.get(1).stop().as_vector().stop().as_f64().stop();
        let x: Vec<f64> = x_vec.slice().to_vec();
        let y: Vec<f64> = y_vec.slice().to_vec();
        new_cubic_spline(&x, &y).stop()
    };
    params.reencode(f);
}

#[roxido]
fn cubic_spline_evaluate(params: &mut RExternalPtr, x: &RVector) {
    resurrect_cubic_spline(params);
    let spline = params.decode_ref::<CubicSpline>();
    let x_f64 = x.to_f64(pc);
    let result = RVector::<f64>::new(x_f64.len(), pc);
    for (i, &xi) in x_f64.slice().iter().enumerate() {
        let yi = spline.evaluate(xi).unwrap_or(f64::NAN);
        result.set(i, yi).stop();
    }
    result
}

#[roxido]
fn print_cubic_spline(params: &mut RExternalPtr) {
    let classes = params.get_class();
    if classes.get(0) != Ok("cubic_spline") {
        stop!("Expected a 'cubic_spline' object.");
    }

    let tag = params.tag();
    let list = tag.as_list().stop();
    let x_vec = list.get(0).stop().as_vector().stop().as_f64().stop();
    let n_knots = x_vec.len();

    resurrect_cubic_spline(params);
    let spline = params.decode_ref::<CubicSpline>();
    let n_segments = spline.num_segments();

    rprintln!(
        "Natural cubic spline ({} control points, {} segments)",
        n_knots,
        n_segments
    );
}

// ============================================================================
// Lorenz curve operations
// ============================================================================

#[roxido]
fn lorenz_evaluate_linear(params: &mut RExternalPtr, x: &RVector) {
    resurrect_lorenz_linear(params);
    let curve = params.decode_ref::<LorenzLinear>();
    let x_f64 = x.to_f64(pc);
    let result = RVector::<f64>::new(x_f64.len(), pc);
    for (i, &xi) in x_f64.slice().iter().enumerate() {
        let yi = curve.evaluate(xi).unwrap_or(f64::NAN);
        result.set(i, yi).stop();
    }
    result
}

#[roxido]
fn lorenz_evaluate_ispline(params: &mut RExternalPtr, x: &RVector) {
    resurrect_lorenz_ispline(params);
    let curve = params.decode_ref::<LorenzISpline>();
    let x_f64 = x.to_f64(pc);
    let result = RVector::<f64>::new(x_f64.len(), pc);
    for (i, &xi) in x_f64.slice().iter().enumerate() {
        let yi = curve.evaluate(xi).unwrap_or(f64::NAN);
        result.set(i, yi).stop();
    }
    result
}

#[roxido]
fn lorenz_gini_linear(params: &mut RExternalPtr) {
    resurrect_lorenz_linear(params);
    let curve = params.decode_ref::<LorenzLinear>();
    curve.gini_coefficient()
}

#[roxido]
fn lorenz_gini_ispline(params: &mut RExternalPtr) {
    resurrect_lorenz_ispline(params);
    let curve = params.decode_ref::<LorenzISpline>();
    curve.gini_coefficient()
}

#[roxido]
fn lorenz_weights_linear(params: &mut RExternalPtr, n_clusters: i32) {
    if n_clusters <= 0 {
        stop!("n_clusters must be positive.");
    }
    resurrect_lorenz_linear(params);
    let curve = params.decode_ref::<LorenzLinear>();
    let k = usize::try_from(n_clusters).stop();
    let weights = curve.weights(k).stop();
    let result = RVector::<f64>::new(weights.len(), pc);
    for (i, &w) in weights.iter().enumerate() {
        result.set(i, w).stop();
    }
    result
}

#[roxido]
fn lorenz_weights_ispline(params: &mut RExternalPtr, n_clusters: i32) {
    if n_clusters <= 0 {
        stop!("n_clusters must be positive.");
    }
    resurrect_lorenz_ispline(params);
    let curve = params.decode_ref::<LorenzISpline>();
    let k = usize::try_from(n_clusters).stop();
    let weights = curve.weights(k).stop();
    let result = RVector::<f64>::new(weights.len(), pc);
    for (i, &w) in weights.iter().enumerate() {
        result.set(i, w).stop();
    }
    result
}

#[roxido]
fn print_lorenz_linear(params: &mut RExternalPtr) {
    let classes = params.get_class();
    if classes.get(0) != Ok("lorenz_linear") {
        stop!("Expected a 'lorenz_linear' object.");
    }

    let tag = params.tag();
    let vec = tag.as_vector().stop();
    let float_vec = vec.as_f64().stop();
    let k = float_vec.len();

    resurrect_lorenz_linear(params);
    let curve = params.decode_ref::<LorenzLinear>();
    let gini = curve.gini_coefficient();

    rprintln!(
        "Piecewise linear Lorenz curve (k = {}, Gini = {:.4})",
        k,
        gini
    );
}

#[roxido]
fn print_lorenz_ispline(params: &mut RExternalPtr) {
    let classes = params.get_class();
    if classes.get(0) != Ok("lorenz_ispline") {
        stop!("Expected a 'lorenz_ispline' object.");
    }

    let tag = params.tag();
    let vec = tag.as_vector().stop();
    let float_vec = vec.as_f64().stop();
    let slice = float_vec.slice();
    let k = slice.len() - 3;
    let degree = slice[k + 1] as i32;

    resurrect_lorenz_ispline(params);
    let curve = params.decode_ref::<LorenzISpline>();
    let gini = curve.gini_coefficient();

    rprintln!(
        "I-spline Lorenz curve (k = {}, degree = {}, Gini = {:.4})",
        k,
        degree,
        gini
    );
}

// ============================================================================
// Standalone Gini coefficient function
// ============================================================================

#[roxido]
fn gini_coefficient(partition: &RVector) {
    let partition_f64 = partition.to_f64(pc);
    let partition_vec: Vec<f64> = partition_f64.slice().to_vec();

    if partition_vec.is_empty() {
        stop!("Partition cannot be empty.");
    }

    match rust_gini_coefficient(&partition_vec) {
        Ok(g) => g,
        Err(e) => stop!("{}", e),
    }
}

// ============================================================================
// Modified rlorenzip and dlorenzip with unified target parameter
// ============================================================================

/// Enum for target parameter (weights or Lorenz curve)
enum TargetParam {
    Values(Vec<f64>),
    LorenzLinear(LorenzLinear),
    LorenzISpline(LorenzISpline),
}

/// Enum for concentration parameter (values or cubic spline)
enum ConcentrationParam {
    Values(Vec<f64>),
    Curve(CubicSpline),
}

/// Enum for kernel parameters (values or cubic spline)
enum KernelParam {
    Values(Vec<f64>),
    Curve(CubicSpline),
}

/// Determine the number of clusters from all sources and validate consistency.
///
/// Sources of k in priority order:
/// 1. Length of target weights vector (if target is a vector)
/// 2. Length of concentration values vector + 1 (if concentration is a vector with len > 1)
/// 3. Length of log_skew values vector + 1 (if log_skew is a vector with len > 1)
/// 4. Length of tail_shape values vector + 1 (if tail_shape is a vector with len > 1)
/// 5. Explicit n_clusters parameter
///
/// All provided sources must be consistent.
fn determine_k(
    target: &RObject,
    concentration: &RObject,
    log_skew: &RObject,
    tail_shape: &RObject,
    n_clusters: &RObject,
    pc: &Pc,
) -> usize {
    let mut k_sources: Vec<(&str, usize)> = Vec::new();

    // Check target if it's a vector
    if target.is_vector()
        && let Ok(target_vec) = target.as_vector()
    {
        let target_f64 = target_vec.to_f64(pc);
        let len = target_f64.len();
        if len > 0 {
            k_sources.push(("target", len));
        }
    }

    // Check concentration if it's a vector with length > 1 (implies k = len + 1)
    if concentration.is_vector()
        && let Ok(conc_vec) = concentration.as_vector()
    {
        let conc_f64 = conc_vec.to_f64(pc);
        let len = conc_f64.len();
        if len > 1 {
            k_sources.push(("concentration", len + 1));
        }
    }

    // Check log_skew if it's a vector with length > 1 (implies k = len + 1)
    if !log_skew.is_null()
        && log_skew.is_vector()
        && let Ok(log_skew_vec) = log_skew.as_vector()
    {
        let log_skew_f64 = log_skew_vec.to_f64(pc);
        let len = log_skew_f64.len();
        if len > 1 {
            k_sources.push(("log_skew", len + 1));
        }
    }

    // Check tail_shape if it's a vector with length > 1 (implies k = len + 1)
    if !tail_shape.is_null()
        && tail_shape.is_vector()
        && let Ok(tail_shape_vec) = tail_shape.as_vector()
    {
        let tail_shape_f64 = tail_shape_vec.to_f64(pc);
        let len = tail_shape_f64.len();
        if len > 1 {
            k_sources.push(("tail_shape", len + 1));
        }
    }

    // Check n_clusters if provided
    if !n_clusters.is_null() {
        let n_clusters_vec = n_clusters
            .as_vector()
            .stop_str("n_clusters should be a scalar");
        let n_clusters_i32 = n_clusters_vec.to_i32(pc);
        let k = n_clusters_i32.get(0).stop();
        if k <= 0 {
            stop!("n_clusters must be positive.");
        }
        k_sources.push(("n_clusters", k as usize));
    }

    // Validate all sources are consistent
    if k_sources.is_empty() {
        stop!(
            "Cannot determine number of clusters: provide n_clusters or a vector for target/concentration/log_skew/tail_shape"
        );
    }

    let k = k_sources[0].1;
    for (source, k_val) in &k_sources[1..] {
        if *k_val != k {
            stop!(
                "Inconsistent number of clusters: {} implies k={}, but {} implies k={}",
                k_sources[0].0,
                k,
                source,
                k_val
            );
        }
    }

    k
}

/// Parse the target parameter into values or curve
fn parse_target(target: &RObject, pc: &Pc) -> TargetParam {
    if target.is_null() {
        stop!("'target' must be specified.");
    }

    // Check if target is a numeric vector (weights)
    if target.is_vector()
        && let Ok(target_vec) = target.as_vector()
    {
        let target_f64 = target_vec.to_f64(pc);
        let mut w: Vec<f64> = target_f64.slice().to_vec();
        // Auto-sort weights to non-decreasing order (the paper's convention).
        w.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return TargetParam::Values(w);
    }

    // Otherwise, treat as Lorenz curve object
    let lorenz_ptr = target
        .as_external_ptr()
        .stop_str("target should be a numeric vector or a lorenz_linear/lorenz_ispline object");

    let classes = lorenz_ptr.get_class();
    let class0 = classes.get(0).stop();

    if class0 == "lorenz_linear" {
        let tag = lorenz_ptr.tag();
        let vec = tag.as_vector().stop();
        let float_vec = vec.as_f64().stop();
        let partition: Vec<f64> = float_vec.slice().to_vec();
        let curve = new_lorenz_linear(&partition).stop();
        TargetParam::LorenzLinear(curve)
    } else if class0 == "lorenz_ispline" {
        let tag = lorenz_ptr.tag();
        let vec = tag.as_vector().stop();
        let float_vec = vec.as_f64().stop();
        let slice = float_vec.slice();
        if slice.len() < 4 {
            stop!("Invalid tag structure in external pointer.");
        }
        let part_len = slice.len() - 3;
        let partition: Vec<f64> = slice[..part_len].to_vec();
        let n_interior_knots = if slice[part_len].is_nan() {
            None
        } else {
            Some(slice[part_len] as usize)
        };
        let degree = slice[part_len + 1] as usize;
        let lambda = slice[part_len + 2];
        let curve = new_lorenz_ispline(&partition, n_interior_knots, degree, lambda).stop();
        TargetParam::LorenzISpline(curve)
    } else {
        stop!("target must be a numeric vector or a 'lorenz_linear'/'lorenz_ispline' object.");
    }
}

/// Parse the concentration parameter into values or curve
fn parse_concentration(concentration: &RObject, k: usize, pc: &Pc) -> ConcentrationParam {
    // Check if it's a cubic spline
    if !concentration.is_vector()
        && let Ok(spline_ptr) = concentration.as_external_ptr()
    {
        let classes = spline_ptr.get_class();
        if classes.get(0) == Ok("cubic_spline") {
            let tag = spline_ptr.tag();
            let list = tag.as_list().stop();
            let x_vec = list.get(0).stop().as_vector().stop().as_f64().stop();
            let y_vec = list.get(1).stop().as_vector().stop().as_f64().stop();
            let x: Vec<f64> = x_vec.slice().to_vec();
            let y: Vec<f64> = y_vec.slice().to_vec();
            let spline = new_cubic_spline(&x, &y).stop();
            return ConcentrationParam::Curve(spline);
        }
    }

    // Otherwise, treat as numeric vector
    let concentration_vec = concentration
        .as_vector()
        .stop_str("concentration should be a vector or cubic_spline");
    let concentration = concentration_vec.to_f64(pc);
    let mut concentration = concentration.slice().to_vec();
    if concentration.len() == 1 {
        let expected = k.saturating_sub(1);
        concentration.resize(expected, *concentration.first().unwrap());
    }
    if concentration.len() != k.saturating_sub(1) {
        stop!(
            "concentration vector length ({}) does not match k-1={} for k={}",
            concentration.len(),
            k.saturating_sub(1),
            k
        );
    }
    ConcentrationParam::Values(concentration)
}

/// Parse a kernel parameter into values or curve
fn parse_kernel_param(param: &RObject, k: usize, label: &str, pc: &Pc) -> Option<KernelParam> {
    if param.is_null() {
        return None;
    }

    // Check if it's a cubic spline
    if !param.is_vector()
        && let Ok(spline_ptr) = param.as_external_ptr()
    {
        let classes = spline_ptr.get_class();
        if classes.get(0) == Ok("cubic_spline") {
            let tag = spline_ptr.tag();
            let list = tag.as_list().stop();
            let x_vec = list.get(0).stop().as_vector().stop().as_f64().stop();
            let y_vec = list.get(1).stop().as_vector().stop().as_f64().stop();
            let x: Vec<f64> = x_vec.slice().to_vec();
            let y: Vec<f64> = y_vec.slice().to_vec();
            let spline = new_cubic_spline(&x, &y).stop();
            return Some(KernelParam::Curve(spline));
        }
    }

    // Otherwise, treat as numeric vector
    let param_vec = param
        .as_vector()
        .stop_str(&format!("{} should be a vector or cubic_spline", label));
    let param_vals = param_vec.to_f64(pc);
    let mut param_vals = param_vals.slice().to_vec();
    if param_vals.len() == 1 {
        let expected = k.saturating_sub(1);
        param_vals.resize(expected, *param_vals.first().unwrap());
    }
    if param_vals.len() != k.saturating_sub(1) {
        stop!(
            "{} vector length ({}) does not match k-1={} for k={}",
            label,
            param_vals.len(),
            k.saturating_sub(1),
            k
        );
    }
    Some(KernelParam::Values(param_vals))
}

fn parse_buffer_c(buffer_c: &RObject, pc: &Pc) -> f64 {
    if buffer_c.is_null() {
        return 1.0;
    }
    let vec = buffer_c
        .as_vector()
        .stop_str("buffer_c should be a numeric scalar");
    let values = vec.to_f64(pc);
    if values.len() != 1 {
        stop!("buffer_c must be a single numeric value.");
    }
    let c = values.get(0).stop();
    if !c.is_finite() || c <= 0.0 {
        stop!("buffer_c must be positive and finite.");
    }
    c
}

fn parse_partition_i64(partition: &RObject, pc: &Pc) -> Vec<i64> {
    if partition.is_matrix() {
        stop!("partition must be a vector.");
    }
    let vec = partition
        .as_vector()
        .stop_str("partition should be a numeric vector");
    let vals = vec.to_f64(pc);
    let slice = vals.slice();
    if slice.is_empty() {
        stop!("partition cannot be empty.");
    }
    let mut out = Vec::with_capacity(slice.len());
    for (i, &x) in slice.iter().enumerate() {
        if !x.is_finite() {
            stop!("partition values must be finite.");
        }
        let rounded = x.round();
        if (x - rounded).abs() > 1e-8 {
            stop!(
                "partition values must be integers; partition[{}]={}",
                i + 1,
                x
            );
        }
        out.push(rounded as i64);
    }
    out
}

fn finalize_data_frame(list: &mut RList, n_rows: usize, pc: &Pc) {
    let class = ["data.frame"].to_r(pc);
    list.set_class(class);
    if n_rows == 0 {
        let row_names = RVector::<i32>::new(0, pc);
        list.set_attribute(RSymbol::rownames(), row_names);
        return;
    }
    let row_names = RVector::<i32>::new(2, pc);
    row_names.set(0, R::na_i32()).stop();
    row_names
        .set(1, -(n_rows as i32))
        .stop_str("failed to set row names");
    list.set_attribute(RSymbol::rownames(), row_names);
}

fn assumptions_output<'a, K: lorenz_ip_distr::KernelFactory>(
    dist: &LorenzIp<K>,
    pc: &'a Pc,
) -> &'a mut RList {
    let report = dist.check_assumptions();
    let reasons: Vec<&str> = report.reasons.iter().map(|s| s.as_str()).collect();
    let reasons_r = reasons.as_slice().to_r(pc);
    let out = RList::with_names(
        &[
            "interior_probabilities",
            "first_mean_within_support",
            "satisfied",
            "reasons",
        ],
        pc,
    );
    out.set(0, report.interior_probabilities.to_r(pc)).stop();
    out.set(1, report.first_mean_within_support.to_r(pc)).stop();
    out.set(2, report.satisfied().to_r(pc)).stop();
    out.set(3, reasons_r).stop();
    out
}

fn diagnostics_path_output<'a, K: lorenz_ip_distr::KernelFactory>(
    dist: &LorenzIp<K>,
    partition: &[i64],
    kernel_is_tadpole: bool,
    pc: &'a Pc,
) -> &'a mut RList {
    let diag = dist.diagnostics(partition).stop();

    let n_steps = diag.steps.len();
    let step_vec = RVector::<i32>::new(n_steps, pc);
    let lower_vec = RVector::<i32>::new(n_steps, pc);
    let upper_vec = RVector::<i32>::new(n_steps, pc);
    let width_vec = RVector::<i32>::new(n_steps, pc);
    let mu_star_vec = RVector::<f64>::new(n_steps, pc);
    let mu_min_vec = RVector::<f64>::new(n_steps, pc);
    let mu_max_vec = RVector::<f64>::new(n_steps, pc);
    let mu_used_vec = RVector::<f64>::new(n_steps, pc);
    let delta_vec = RVector::<f64>::new(n_steps, pc);
    let buffer_vec = RVector::<f64>::new(n_steps, pc);
    let band_lower_vec = RVector::<f64>::new(n_steps, pc);
    let band_upper_vec = RVector::<f64>::new(n_steps, pc);
    let clamped_vec = RVector::<i32>::new(n_steps, pc);
    let degenerate_vec = RVector::<i32>::new(n_steps, pc);
    let mu_star_feasible_vec = RVector::<i32>::new(n_steps, pc);
    let kappa_vec = RVector::<f64>::new(n_steps, pc);
    let concentration_vec = RVector::<f64>::new(n_steps, pc);
    let log_skew_vec = RVector::<f64>::new(n_steps, pc);
    let log_skew_eff_vec = RVector::<f64>::new(n_steps, pc);
    let log_skew_adjusted_vec = RVector::<i32>::new(n_steps, pc);
    let mut tail_shape_vec = if kernel_is_tadpole {
        Some(RVector::<f64>::new(n_steps, pc))
    } else {
        None
    };

    for (i, step) in diag.steps.iter().enumerate() {
        step_vec.set(i, (step.index + 1) as i32).stop();
        lower_vec.set(i, step.lower as i32).stop();
        upper_vec.set(i, step.upper as i32).stop();
        width_vec.set(i, step.width as i32).stop();
        mu_star_vec.set(i, step.mu_star).stop();
        mu_min_vec.set(i, step.mu_min).stop();
        mu_max_vec.set(i, step.mu_max).stop();
        mu_used_vec.set(i, step.mu_used).stop();
        delta_vec.set(i, step.delta).stop();
        buffer_vec.set(i, step.buffer).stop();
        band_lower_vec.set(i, step.band_lower).stop();
        band_upper_vec.set(i, step.band_upper).stop();
        clamped_vec
            .set(i, if step.clamped { R::TRUE() } else { R::FALSE() })
            .stop();
        degenerate_vec
            .set(
                i,
                if step.degenerate {
                    R::TRUE()
                } else {
                    R::FALSE()
                },
            )
            .stop();
        mu_star_feasible_vec
            .set(
                i,
                if step.mu_star_feasible {
                    R::TRUE()
                } else {
                    R::FALSE()
                },
            )
            .stop();
        kappa_vec.set(i, step.kappa).stop();
        concentration_vec.set(i, step.concentration).stop();
        log_skew_vec.set(i, step.log_skew).stop();
        log_skew_eff_vec.set(i, step.log_skew_eff).stop();
        log_skew_adjusted_vec
            .set(
                i,
                if step.log_skew_adjusted {
                    R::TRUE()
                } else {
                    R::FALSE()
                },
            )
            .stop();
        if let Some(tail_shape_vec) = tail_shape_vec.as_mut() {
            tail_shape_vec.set(i, step.tail_shape).stop();
        }
    }

    let mut names = vec![
        "step",
        "lower",
        "upper",
        "width",
        "mu_star",
        "mu_min",
        "mu_max",
        "mu_used",
        "delta",
        "buffer",
        "band_lower",
        "band_upper",
        "clamped",
        "degenerate",
        "mu_star_feasible",
        "kappa",
        "concentration",
        "log_skew",
        "log_skew_eff",
        "log_skew_adjusted",
    ];
    if kernel_is_tadpole {
        names.push("tail_shape");
    }

    let steps = RList::with_names(&names, pc);
    let mut idx = 0;
    steps.set(idx, step_vec).stop();
    idx += 1;
    steps.set(idx, lower_vec).stop();
    idx += 1;
    steps.set(idx, upper_vec).stop();
    idx += 1;
    steps.set(idx, width_vec).stop();
    idx += 1;
    steps.set(idx, mu_star_vec).stop();
    idx += 1;
    steps.set(idx, mu_min_vec).stop();
    idx += 1;
    steps.set(idx, mu_max_vec).stop();
    idx += 1;
    steps.set(idx, mu_used_vec).stop();
    idx += 1;
    steps.set(idx, delta_vec).stop();
    idx += 1;
    steps.set(idx, buffer_vec).stop();
    idx += 1;
    steps.set(idx, band_lower_vec).stop();
    idx += 1;
    steps.set(idx, band_upper_vec).stop();
    idx += 1;
    steps.set(idx, clamped_vec).stop();
    idx += 1;
    steps.set(idx, degenerate_vec).stop();
    idx += 1;
    steps.set(idx, mu_star_feasible_vec).stop();
    idx += 1;
    steps.set(idx, kappa_vec).stop();
    idx += 1;
    steps.set(idx, concentration_vec).stop();
    idx += 1;
    steps.set(idx, log_skew_vec).stop();
    idx += 1;
    steps.set(idx, log_skew_eff_vec).stop();
    idx += 1;
    steps.set(idx, log_skew_adjusted_vec).stop();
    idx += 1;
    if let Some(tail_shape_vec) = tail_shape_vec {
        steps.set(idx, tail_shape_vec).stop();
    }
    finalize_data_frame(steps, n_steps, pc);

    let summary = RList::with_names(&["any_clamped", "max_abs_delta"], pc);
    summary.set(0, diag.any_clamped.to_r(pc)).stop();
    summary.set(1, diag.max_abs_delta.to_r(pc)).stop();
    let assumptions = assumptions_output(dist, pc);

    let out = RList::with_names(&["steps", "summary", "assumptions"], pc);
    out.set(0, steps).stop();
    out.set(1, summary).stop();
    out.set(2, assumptions).stop();
    out
}

fn diagnostics_mc_output<'a, K: lorenz_ip_distr::KernelFactory>(
    dist: &LorenzIp<K>,
    kernel_is_tadpole: bool,
    n_samples: usize,
    rng: &mut Pcg64Mcg,
    pc: &'a Pc,
) -> &'a mut RList {
    let diag = dist.diagnostics_mc(rng, n_samples).stop();

    let n_steps = diag.steps.len();
    let step_vec = RVector::<i32>::new(n_steps, pc);
    let clamped_rate_vec = RVector::<f64>::new(n_steps, pc);
    let degenerate_rate_vec = RVector::<f64>::new(n_steps, pc);
    let mean_delta_vec = RVector::<f64>::new(n_steps, pc);
    let mean_abs_delta_vec = RVector::<f64>::new(n_steps, pc);
    let max_abs_delta_vec = RVector::<f64>::new(n_steps, pc);
    let kappa_vec = RVector::<f64>::new(n_steps, pc);
    let concentration_vec = RVector::<f64>::new(n_steps, pc);
    let log_skew_vec = RVector::<f64>::new(n_steps, pc);
    let mut tail_shape_vec = if kernel_is_tadpole {
        Some(RVector::<f64>::new(n_steps, pc))
    } else {
        None
    };
    let log_skew_values = dist.log_skew();
    let tail_shape_values = dist.tail_shape();

    for (i, step) in diag.steps.iter().enumerate() {
        step_vec.set(i, (step.index + 1) as i32).stop();
        clamped_rate_vec.set(i, step.clamped_rate).stop();
        degenerate_rate_vec.set(i, step.degenerate_rate).stop();
        mean_delta_vec.set(i, step.mean_delta).stop();
        mean_abs_delta_vec.set(i, step.mean_abs_delta).stop();
        max_abs_delta_vec.set(i, step.max_abs_delta).stop();
        kappa_vec.set(i, dist.kappas()[step.index]).stop();
        concentration_vec
            .set(i, dist.concentration()[step.index])
            .stop();
        let log_skew_val = log_skew_values
            .get(step.index)
            .copied()
            .unwrap_or(R::na_f64());
        log_skew_vec.set(i, log_skew_val).stop();
        if let Some(tail_shape_vec) = tail_shape_vec.as_mut() {
            let tail_shape_val = tail_shape_values
                .get(step.index)
                .copied()
                .unwrap_or(R::na_f64());
            tail_shape_vec.set(i, tail_shape_val).stop();
        }
    }

    let mut names = vec![
        "step",
        "clamped_rate",
        "degenerate_rate",
        "mean_delta",
        "mean_abs_delta",
        "max_abs_delta",
        "kappa",
        "concentration",
        "log_skew",
    ];
    if kernel_is_tadpole {
        names.push("tail_shape");
    }

    let steps = RList::with_names(&names, pc);
    let mut idx = 0;
    steps.set(idx, step_vec).stop();
    idx += 1;
    steps.set(idx, clamped_rate_vec).stop();
    idx += 1;
    steps.set(idx, degenerate_rate_vec).stop();
    idx += 1;
    steps.set(idx, mean_delta_vec).stop();
    idx += 1;
    steps.set(idx, mean_abs_delta_vec).stop();
    idx += 1;
    steps.set(idx, max_abs_delta_vec).stop();
    idx += 1;
    steps.set(idx, kappa_vec).stop();
    idx += 1;
    steps.set(idx, concentration_vec).stop();
    idx += 1;
    steps.set(idx, log_skew_vec).stop();
    idx += 1;
    if let Some(tail_shape_vec) = tail_shape_vec {
        steps.set(idx, tail_shape_vec).stop();
    }
    finalize_data_frame(steps, n_steps, pc);

    let summary = RList::with_names(
        &[
            "any_clamped_rate",
            "mean_max_abs_delta",
            "max_abs_delta",
            "n_samples",
        ],
        pc,
    );
    summary.set(0, diag.any_clamped_rate.to_r(pc)).stop();
    summary.set(1, diag.mean_max_abs_delta.to_r(pc)).stop();
    summary.set(2, diag.max_abs_delta.to_r(pc)).stop();
    summary.set(3, (diag.n_samples as i32).to_r(pc)).stop();

    let sample_max = RVector::from_slice(&diag.sample_max_abs_delta, pc);
    let assumptions = assumptions_output(dist, pc);

    let out = RList::with_names(
        &["steps", "summary", "sample_max_abs_delta", "assumptions"],
        pc,
    );
    out.set(0, steps).stop();
    out.set(1, summary).stop();
    out.set(2, sample_max).stop();
    out.set(3, assumptions).stop();
    out
}

fn update_fixed_k(
    fixed_k: &mut Option<usize>,
    fixed_source: &mut Option<&'static str>,
    source: &'static str,
    k: usize,
) {
    if let Some(existing) = *fixed_k {
        if existing != k {
            stop!(
                "Inconsistent number of clusters: {} implies k={}, but {} implies k={}",
                fixed_source.unwrap_or("previous source"),
                existing,
                source,
                k
            );
        }
    } else {
        *fixed_k = Some(k);
        *fixed_source = Some(source);
    }
}

fn parse_log_k_weights(k_weights: &RObject, log_k_weights: bool, pc: &Pc) -> Vec<f64> {
    if k_weights.is_null() {
        stop!("k_weights must be provided.");
    }
    let vec = k_weights
        .as_vector()
        .stop_str("k_weights should be a numeric vector");
    let vec_f64 = vec.to_f64(pc);
    let slice = vec_f64.slice();
    if slice.is_empty() {
        stop!("k_weights cannot be empty.");
    }

    let mut log_weights = Vec::with_capacity(slice.len());
    for &w in slice.iter() {
        if w.is_nan() {
            stop!("k_weights must be finite.");
        }
        if log_k_weights {
            if w.is_infinite() && w.is_sign_positive() {
                stop!("log_k_weights cannot contain +Inf.");
            }
            log_weights.push(w);
        } else {
            if !w.is_finite() {
                stop!("k_weights must be finite.");
            }
            if w < 0.0 {
                stop!("k_weights must be non-negative.");
            }
            if w == 0.0 {
                log_weights.push(f64::NEG_INFINITY);
            } else {
                log_weights.push(w.ln());
            }
        }
    }

    log_weights
}

fn log_sum_exp(values: &[f64]) -> f64 {
    if values.is_empty() {
        return f64::NEG_INFINITY;
    }
    let mut max_val = f64::NEG_INFINITY;
    for &v in values {
        if v > max_val {
            max_val = v;
        }
    }
    if !max_val.is_finite() {
        return f64::NEG_INFINITY;
    }
    let mut sum = 0.0;
    for &v in values {
        if v.is_finite() {
            sum += (v - max_val).exp();
        }
    }
    max_val + sum.ln()
}

fn normalize_log_weights(log_weights: &[f64], k_max: usize) -> Vec<f64> {
    if k_max == 0 {
        stop!("k_weights has no support for k <= 0.");
    }
    let slice = &log_weights[..k_max];
    let log_norm = log_sum_exp(slice);
    if !log_norm.is_finite() {
        stop!("k_weights has no positive mass for k <= {}.", k_max);
    }
    slice.iter().map(|&w| w - log_norm).collect()
}

fn normalize_log_weights_allow_empty(log_weights: &[f64], k_max: usize) -> Vec<f64> {
    if k_max == 0 {
        return Vec::new();
    }
    let slice = &log_weights[..k_max];
    let log_norm = log_sum_exp(slice);
    if !log_norm.is_finite() {
        return vec![f64::NEG_INFINITY; k_max];
    }
    slice.iter().map(|&w| w - log_norm).collect()
}

fn log_probs_to_cdf(log_probs: &[f64]) -> Vec<f64> {
    let mut cdf = Vec::with_capacity(log_probs.len());
    let mut sum = 0.0;
    for &lp in log_probs {
        if lp.is_finite() {
            sum += lp.exp();
        }
        cdf.push(sum);
    }
    if sum <= 0.0 {
        stop!("Internal error: invalid k prior normalization.");
    }
    for v in &mut cdf {
        *v /= sum;
    }
    cdf
}

fn sample_k_from_cdf(cdf: &[f64], rng: &mut Pcg64Mcg) -> usize {
    let u: f64 = rng.random();
    for (i, &p) in cdf.iter().enumerate() {
        if u < p {
            return i + 1;
        }
    }
    cdf.len()
}

fn log_q_n<K: lorenz_ip_distr::KernelFactory>(
    partition: &[i64],
    log_p_k_n1: &[f64],
    dists_n1: &[Option<LorenzIp<K>>],
    n_items: i32,
) -> f64 {
    if partition.is_empty() {
        return f64::NEG_INFINITY;
    }
    let k = partition.len();
    let n_plus_one = n_items as f64 + 1.0;

    let mut counts: BTreeMap<i64, usize> = BTreeMap::new();
    for &val in partition {
        *counts.entry(val).or_insert(0) += 1;
    }

    let mut log_terms = Vec::with_capacity(counts.len() + 1);
    let c1 = *counts.get(&1).unwrap_or(&0) as f64;

    let k0 = k + 1;
    if k0 <= log_p_k_n1.len()
        && log_p_k_n1[k0 - 1].is_finite()
        && let Some(dist) = dists_n1.get(k0).and_then(|d| d.as_ref())
    {
        let mut parent = Vec::with_capacity(k0);
        parent.extend_from_slice(partition);
        parent.push(1);
        parent.sort_unstable();
        let log_pmf = compute_log_pmf(dist, &parent).stop();
        let log_weight = ((c1 + 1.0) / n_plus_one).ln();
        log_terms.push(log_p_k_n1[k0 - 1] + log_pmf + log_weight);
    }

    if k <= log_p_k_n1.len() && log_p_k_n1[k - 1].is_finite() {
        let dist = dists_n1
            .get(k)
            .and_then(|d| d.as_ref())
            .stop_str("internal error: missing parent distribution");
        for (&s, &c_s) in counts.iter() {
            if c_s == 0 {
                continue;
            }
            let mut parent = partition.to_vec();
            let pos = parent
                .iter()
                .position(|&v| v == s)
                .stop_str("internal error: missing part size");
            parent[pos] += 1;
            parent.sort_unstable();
            let log_pmf = compute_log_pmf(dist, &parent).stop();
            let c_s1 = *counts.get(&(s + 1)).unwrap_or(&0) as f64;
            let weight = (c_s1 + 1.0) * (s as f64 + 1.0) / n_plus_one;
            log_terms.push(log_p_k_n1[k - 1] + log_pmf + weight.ln());
        }
    }

    log_sum_exp(&log_terms)
}

fn log_crp_size_profile(
    partition: &[usize],
    n_items: usize,
    concentration: f64,
    discount: f64,
) -> f64 {
    if partition.is_empty() {
        return f64::NEG_INFINITY;
    }

    let mut log_terms = 0.0;
    for &n_i in partition {
        if n_i == 0 {
            return f64::NEG_INFINITY;
        }
        let m = n_i - 1;
        let log_factor = log_rising_factorial_general(1.0 - discount, m);
        if !log_factor.is_finite() {
            return f64::NEG_INFINITY;
        }
        log_terms += log_factor;
    }

    let log_c = log_concentration_factor(concentration, discount, partition.len());
    if !log_c.is_finite() {
        return f64::NEG_INFINITY;
    }
    let log_denom = log_rising_factorial(concentration, n_items.saturating_sub(1));
    let log_eppf = log_c + log_terms - log_denom;
    let log_count = multinomial_partition_log(n_items, partition);
    log_eppf + log_count
}

fn get_lorenz_log_pmf<K, F>(
    k: usize,
    partition: &[i64],
    dists: &mut [Option<LorenzIp<K>>],
    attempted: &mut [bool],
    build_dist: &mut F,
) -> f64
where
    K: lorenz_ip_distr::KernelFactory,
    F: FnMut(usize) -> Option<LorenzIp<K>>,
{
    if k >= dists.len() {
        return f64::NEG_INFINITY;
    }
    if !attempted[k] {
        attempted[k] = true;
        dists[k] = build_dist(k);
    }
    if let Some(dist) = dists[k].as_ref() {
        compute_log_pmf(dist, partition).unwrap_or(f64::NEG_INFINITY)
    } else {
        f64::NEG_INFINITY
    }
}

fn downsample_mc<K, F>(
    n_items: i32,
    n_samples: usize,
    log_p_k_n: &[f64],
    cdf_k_n: &[f64],
    log_p_k_n1: &[f64],
    mut build_dist: F,
) -> (f64, f64)
where
    K: lorenz_ip_distr::KernelFactory,
    F: FnMut(i32, usize) -> LorenzIp<K>,
{
    let k_max_n = log_p_k_n.len();
    let k_max_n1 = log_p_k_n1.len();

    let mut dists_n: Vec<Option<LorenzIp<K>>> = (0..=k_max_n).map(|_| None).collect();
    for k in 1..=k_max_n {
        if log_p_k_n[k - 1].is_finite() {
            dists_n[k] = Some(build_dist(n_items, k));
        }
    }

    let mut dists_n1: Vec<Option<LorenzIp<K>>> = (0..=k_max_n1).map(|_| None).collect();
    for k in 1..=k_max_n1 {
        if log_p_k_n1[k - 1].is_finite() {
            dists_n1[k] = Some(build_dist(n_items + 1, k));
        }
    }

    let seed = R::random_bytes::<16>();
    let mut rng = Pcg64Mcg::from_seed(seed);

    let mut tv_sum = 0.0;
    let mut kl_sum = 0.0;

    for _ in 0..n_samples {
        let k = sample_k_from_cdf(cdf_k_n, &mut rng);
        let dist_n = dists_n
            .get(k)
            .and_then(|d| d.as_ref())
            .stop_str("internal error: missing distribution for sampled k");
        let partition = dist_n.sample(&mut rng).stop();
        let log_pmf = compute_log_pmf(dist_n, &partition).stop();
        let log_p = log_p_k_n[k - 1] + log_pmf;
        let log_q = log_q_n(&partition, log_p_k_n1, &dists_n1, n_items);
        let ratio = (log_q - log_p).exp();
        tv_sum += (1.0 - ratio).abs();
        kl_sum += log_p - log_q;
    }

    let n_samples_f = n_samples as f64;
    let tv = 0.5 * tv_sum / n_samples_f;
    let kl = kl_sum / n_samples_f;
    (tv, kl)
}

#[roxido]
fn rlorenzip(
    n: usize,
    n_items: i32,
    target: &RObject,
    concentration: &RObject,
    log_skew: &RObject,
    tail_shape: &RObject,
    buffer_c: &RObject,
    n_clusters: &RObject,
) {
    // Determine k from all sources
    let kernel_is_tadpole = !tail_shape.is_null();
    let k = determine_k(target, concentration, log_skew, tail_shape, n_clusters, pc);

    // Parse parameters
    let target_param = parse_target(target, pc);
    let concentration_param = parse_concentration(concentration, k, pc);
    let buffer_c = parse_buffer_c(buffer_c, pc);
    // Output matrix (n rows x k cols)
    let out = RMatrix::<i32>::from_value(0, n, k, pc);

    if k == 0 {
        return out;
    }

    // RNG seeded from R's RNG
    let seed = R::random_bytes::<16>();
    let mut rng = Pcg64Mcg::from_seed(seed);

    // Build distribution based on kernel family
    if !kernel_is_tadpole {
        // TiDaL variant
        let mut builder = TidalLorenzIpBuilder::new(n_items.into());

        // Set target (weights or curve)
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };

        // Set concentration (values or curve)
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };

        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }

        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();

        for i in 0..n {
            let sample = dist.sample(&mut rng).stop();
            if sample.len() != k {
                stop!(
                    "Sampled partition length {} does not match expected {}",
                    sample.len(),
                    k
                );
            }
            for (j, &val) in sample.iter().enumerate() {
                out.set(i, j, i32::try_from(val).stop()).stop();
            }
        }
    } else {
        // TaDPoLe variant
        let mut builder = TadpoleLorenzIpBuilder::new(n_items.into());

        // Set target (weights or curve)
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };

        // Set concentration (values or curve)
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };

        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }
        if let Some(param) = parse_kernel_param(tail_shape, k, "tail_shape", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.tail_shape_values(values).stop(),
                KernelParam::Curve(spline) => builder.tail_shape_curve(spline).stop(),
            };
        }

        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();

        for i in 0..n {
            let sample = dist.sample(&mut rng).stop();
            if sample.len() != k {
                stop!(
                    "Sampled partition length {} does not match expected {}",
                    sample.len(),
                    k
                );
            }
            for (j, &val) in sample.iter().enumerate() {
                out.set(i, j, i32::try_from(val).stop()).stop();
            }
        }
    }
    out
}

/// Helper function to compute log probabilities for a partition using a LorenzIp distribution.
fn compute_log_pmf<K: lorenz_ip_distr::KernelFactory>(
    dist: &LorenzIp<K>,
    partition: &[i64],
) -> Result<f64, lorenz_ip_distr::LorenzIpError> {
    dist.log_pmf(partition)
}

#[roxido]
fn dlorenzip(
    x: &RObject,
    n_items: i32,
    target: &RObject,
    concentration: &RObject,
    log_skew: &RObject,
    tail_shape: &RObject,
    buffer_c: &RObject,
    log: bool,
    n_clusters: &RObject,
) {
    // Determine k from all sources
    let kernel_is_tadpole = !tail_shape.is_null();
    let k = determine_k(target, concentration, log_skew, tail_shape, n_clusters, pc);

    if k == 0 {
        // Return empty vector for empty weights
        return RVector::<f64>::from_value(f64::NAN, 0, pc);
    }

    // Parse parameters
    let target_param = parse_target(target, pc);
    let concentration_param = parse_concentration(concentration, k, pc);
    let buffer_c = parse_buffer_c(buffer_c, pc);
    // Determine the number of partitions and get the data
    // x can be either a vector (single partition) or a matrix (multiple partitions)
    let (n_partitions, mut x_data): (usize, Vec<i64>) = if x.is_matrix() {
        let x_mat = x.as_matrix().stop_str("x should be a matrix");
        let x_mat = x_mat.to_i32(pc);
        let nrow = x_mat.nrow();
        let ncol = x_mat.ncol();
        if ncol != k {
            stop!(
                "Number of columns ({}) does not match number of clusters ({})",
                ncol,
                k
            );
        }
        // Convert to i64 and transpose (R is column-major, we need row-major)
        let slice = x_mat.slice();
        let mut data = Vec::with_capacity(nrow * ncol);
        for i in 0..nrow {
            for j in 0..ncol {
                data.push(slice[j * nrow + i] as i64);
            }
        }
        (nrow, data)
    } else {
        // Treat as a vector (single partition)
        let x_vec = x.as_vector().stop_str("x should be a vector or matrix");
        let x_vec = x_vec.to_i32(pc);
        if x_vec.len() != k {
            stop!(
                "Length of x ({}) does not match number of clusters ({})",
                x_vec.len(),
                k
            );
        }
        let data: Vec<i64> = x_vec.slice().iter().map(|&v| v as i64).collect();
        (1, data)
    };

    if k > 1 {
        for i in 0..n_partitions {
            x_data[i * k..(i + 1) * k].sort_unstable();
        }
    }

    // Output vector
    let out = RVector::from_value(f64::NAN, n_partitions, pc);

    // Build Lorenz-IP distribution and compute probabilities
    if !kernel_is_tadpole {
        // TiDaL variant
        let mut builder = TidalLorenzIpBuilder::new(n_items.into());

        // Set target (weights or curve)
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };

        // Set concentration (values or curve)
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };

        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }

        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();

        for i in 0..n_partitions {
            let partition = &x_data[i * k..(i + 1) * k];
            let log_prob = compute_log_pmf(&dist, partition).unwrap_or(f64::NEG_INFINITY);
            let result = if log { log_prob } else { log_prob.exp() };
            out.set(i, result).stop();
        }
    } else {
        // TaDPoLe variant
        let mut builder = TadpoleLorenzIpBuilder::new(n_items.into());

        // Set target (weights or curve)
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };

        // Set concentration (values or curve)
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };

        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }
        if let Some(param) = parse_kernel_param(tail_shape, k, "tail_shape", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.tail_shape_values(values).stop(),
                KernelParam::Curve(spline) => builder.tail_shape_curve(spline).stop(),
            };
        }

        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();

        for i in 0..n_partitions {
            let partition = &x_data[i * k..(i + 1) * k];
            let log_prob = compute_log_pmf(&dist, partition).unwrap_or(f64::NEG_INFINITY);
            let result = if log { log_prob } else { log_prob.exp() };
            out.set(i, result).stop();
        }
    }
    out
}

#[roxido]
fn lorenzip_diagnostics_path(
    partition: &RObject,
    n_items: i32,
    target: &RObject,
    concentration: &RObject,
    log_skew: &RObject,
    tail_shape: &RObject,
    buffer_c: &RObject,
    n_clusters: &RObject,
) {
    let kernel_is_tadpole = !tail_shape.is_null();
    let k = determine_k(target, concentration, log_skew, tail_shape, n_clusters, pc);
    if k == 0 {
        stop!("number of clusters must be positive.");
    }
    let mut partition_vec = parse_partition_i64(partition, pc);
    partition_vec.sort_unstable();
    if partition_vec.len() != k {
        stop!(
            "partition length ({}) does not match number of clusters ({})",
            partition_vec.len(),
            k
        );
    }

    let target_param = parse_target(target, pc);
    let concentration_param = parse_concentration(concentration, k, pc);
    let buffer_c = parse_buffer_c(buffer_c, pc);

    if !kernel_is_tadpole {
        let mut builder = TidalLorenzIpBuilder::new(n_items.into());
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };
        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }
        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();
        diagnostics_path_output(&dist, &partition_vec, kernel_is_tadpole, pc)
    } else {
        let mut builder = TadpoleLorenzIpBuilder::new(n_items.into());
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };
        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }
        if let Some(param) = parse_kernel_param(tail_shape, k, "tail_shape", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.tail_shape_values(values).stop(),
                KernelParam::Curve(spline) => builder.tail_shape_curve(spline).stop(),
            };
        }
        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();
        diagnostics_path_output(&dist, &partition_vec, kernel_is_tadpole, pc)
    }
}

#[roxido]
fn lorenzip_diagnostics_mc(
    n_samples: i32,
    n_items: i32,
    target: &RObject,
    concentration: &RObject,
    log_skew: &RObject,
    tail_shape: &RObject,
    buffer_c: &RObject,
    n_clusters: &RObject,
) {
    if n_samples <= 0 {
        stop!("n_samples must be positive.");
    }
    let n_samples_usize = n_samples as usize;
    let kernel_is_tadpole = !tail_shape.is_null();
    let k = determine_k(target, concentration, log_skew, tail_shape, n_clusters, pc);
    if k == 0 {
        stop!("number of clusters must be positive.");
    }

    let target_param = parse_target(target, pc);
    let concentration_param = parse_concentration(concentration, k, pc);
    let buffer_c = parse_buffer_c(buffer_c, pc);

    let seed = R::random_bytes::<16>();
    let mut rng = Pcg64Mcg::from_seed(seed);

    if !kernel_is_tadpole {
        let mut builder = TidalLorenzIpBuilder::new(n_items.into());
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };
        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }
        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();
        diagnostics_mc_output(&dist, kernel_is_tadpole, n_samples_usize, &mut rng, pc)
    } else {
        let mut builder = TadpoleLorenzIpBuilder::new(n_items.into());
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };
        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }
        if let Some(param) = parse_kernel_param(tail_shape, k, "tail_shape", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.tail_shape_values(values).stop(),
                KernelParam::Curve(spline) => builder.tail_shape_curve(spline).stop(),
            };
        }
        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();
        diagnostics_mc_output(&dist, kernel_is_tadpole, n_samples_usize, &mut rng, pc)
    }
}

#[roxido]
fn cocluster_mc(
    n_samples: i32,
    n_items: i32,
    target: &RObject,
    concentration: &RObject,
    log_skew: &RObject,
    tail_shape: &RObject,
    buffer_c: &RObject,
    n_clusters: &RObject,
) {
    if n_samples <= 0 {
        stop!("n_samples must be positive.");
    }
    if n_items < 2 {
        stop!("n_items must be at least 2.");
    }
    let n_samples_usize = n_samples as usize;
    let kernel_is_tadpole = !tail_shape.is_null();
    let k = determine_k(target, concentration, log_skew, tail_shape, n_clusters, pc);
    if k == 0 {
        stop!("number of clusters must be positive.");
    }

    let target_param = parse_target(target, pc);
    let concentration_param = parse_concentration(concentration, k, pc);
    let buffer_c = parse_buffer_c(buffer_c, pc);

    let seed = R::random_bytes::<16>();
    let mut rng = Pcg64Mcg::from_seed(seed);

    if !kernel_is_tadpole {
        let mut builder = TidalLorenzIpBuilder::new(n_items.into());
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };
        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }
        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();
        dist.cocluster_mc(&mut rng, n_samples_usize).stop()
    } else {
        let mut builder = TadpoleLorenzIpBuilder::new(n_items.into());
        builder = match target_param {
            TargetParam::Values(weights) => builder.lorenz_weights(weights).stop(),
            TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve).stop(),
            TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve).stop(),
        };
        builder = match concentration_param {
            ConcentrationParam::Values(values) => builder.concentration_values(values).stop(),
            ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
        };
        if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
            };
        }
        if let Some(param) = parse_kernel_param(tail_shape, k, "tail_shape", pc) {
            builder = match param {
                KernelParam::Values(values) => builder.tail_shape_values(values).stop(),
                KernelParam::Curve(spline) => builder.tail_shape_curve(spline).stop(),
            };
        }
        builder = builder.buffer_scale(buffer_c).stop();
        let dist = builder.build(k).stop();
        dist.cocluster_mc(&mut rng, n_samples_usize).stop()
    }
}

#[roxido]
fn downsample(
    n_items: i32,
    k_weights: &RObject,
    target: &RObject,
    concentration: &RObject,
    log_skew: &RObject,
    tail_shape: &RObject,
    buffer_c: &RObject,
    n_samples: i32,
    log_k_weights: bool,
) {
    if n_items <= 0 {
        stop!("n_items must be at least 1.");
    }
    if n_samples <= 0 {
        stop!("n_samples must be positive.");
    }

    let log_weights = parse_log_k_weights(k_weights, log_k_weights, pc);
    let r = log_weights.len();
    let n = n_items as usize;
    let k_max_n = n.min(r);
    let k_max_n1 = (n + 1).min(r);

    if k_max_n == 0 {
        stop!("k_weights has no support for k <= n_items.");
    }
    if k_max_n1 == 0 {
        stop!("k_weights has no support for k <= n_items + 1.");
    }

    let target_param = parse_target(target, pc);
    let buffer_c = parse_buffer_c(buffer_c, pc);

    let mut fixed_k: Option<usize> = None;
    let mut fixed_source: Option<&'static str> = None;

    if let TargetParam::Values(weights) = &target_param {
        if weights.is_empty() {
            stop!("target must be non-empty.");
        }
        update_fixed_k(&mut fixed_k, &mut fixed_source, "target", weights.len());
    }

    if concentration.is_vector() {
        let conc_vec = concentration
            .as_vector()
            .stop_str("concentration should be a vector or cubic_spline");
        let conc_f64 = conc_vec.to_f64(pc);
        let len = conc_f64.len();
        if len > 1 {
            update_fixed_k(&mut fixed_k, &mut fixed_source, "concentration", len + 1);
        }
    }

    if !log_skew.is_null() && log_skew.is_vector() {
        let log_skew_vec = log_skew
            .as_vector()
            .stop_str("log_skew should be a vector or cubic_spline");
        let log_skew_f64 = log_skew_vec.to_f64(pc);
        let len = log_skew_f64.len();
        if len > 1 {
            update_fixed_k(&mut fixed_k, &mut fixed_source, "log_skew", len + 1);
        }
    }

    if !tail_shape.is_null() && tail_shape.is_vector() {
        let tail_shape_vec = tail_shape
            .as_vector()
            .stop_str("tail_shape should be a vector or cubic_spline");
        let tail_shape_f64 = tail_shape_vec.to_f64(pc);
        let len = tail_shape_f64.len();
        if len > 1 {
            update_fixed_k(&mut fixed_k, &mut fixed_source, "tail_shape", len + 1);
        }
    }

    if let Some(k_fixed) = fixed_k {
        if k_fixed > k_max_n1 {
            stop!("k_weights has no support for required k={}.", k_fixed);
        }
        for (i, &lw) in log_weights[..k_max_n1].iter().enumerate() {
            let k = i + 1;
            if k != k_fixed && lw.is_finite() {
                stop!(
                    "k_weights must place all mass on k={} because {} fixes k.",
                    k_fixed,
                    fixed_source.unwrap_or("a parameter")
                );
            }
        }
        if !log_weights[k_fixed - 1].is_finite() {
            stop!("k_weights must assign positive mass to k={}.", k_fixed);
        }
    }

    let log_p_k_n = normalize_log_weights(&log_weights, k_max_n);
    let log_p_k_n1 = normalize_log_weights(&log_weights, k_max_n1);
    let cdf_k_n = log_probs_to_cdf(&log_p_k_n);

    let n_samples = n_samples as usize;

    let kernel_is_tadpole = !tail_shape.is_null();
    let (tv, kl) = if !kernel_is_tadpole {
        downsample_mc(
            n_items,
            n_samples,
            &log_p_k_n,
            &cdf_k_n,
            &log_p_k_n1,
            |n_items_local, k| {
                let mut builder = TidalLorenzIpBuilder::new(n_items_local.into());

                builder = match &target_param {
                    TargetParam::Values(weights) => builder.lorenz_weights(weights.clone()).stop(),
                    TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve.clone()).stop(),
                    TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve.clone()).stop(),
                };

                let concentration_param = parse_concentration(concentration, k, pc);
                builder = match concentration_param {
                    ConcentrationParam::Values(values) => {
                        builder.concentration_values(values).stop()
                    }
                    ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
                };

                if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
                    builder = match param {
                        KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                        KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
                    };
                }

                builder = builder.buffer_scale(buffer_c).stop();
                builder.build(k).stop()
            },
        )
    } else {
        downsample_mc(
            n_items,
            n_samples,
            &log_p_k_n,
            &cdf_k_n,
            &log_p_k_n1,
            |n_items_local, k| {
                let mut builder = TadpoleLorenzIpBuilder::new(n_items_local.into());

                builder = match &target_param {
                    TargetParam::Values(weights) => builder.lorenz_weights(weights.clone()).stop(),
                    TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve.clone()).stop(),
                    TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve.clone()).stop(),
                };

                let concentration_param = parse_concentration(concentration, k, pc);
                builder = match concentration_param {
                    ConcentrationParam::Values(values) => {
                        builder.concentration_values(values).stop()
                    }
                    ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
                };

                if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
                    builder = match param {
                        KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                        KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
                    };
                }
                if let Some(param) = parse_kernel_param(tail_shape, k, "tail_shape", pc) {
                    builder = match param {
                        KernelParam::Values(values) => builder.tail_shape_values(values).stop(),
                        KernelParam::Curve(spline) => builder.tail_shape_curve(spline).stop(),
                    };
                }

                builder = builder.buffer_scale(buffer_c).stop();
                builder.build(k).stop()
            },
        )
    };

    let out = RList::with_names(&["tv", "kl"], pc);
    let tv_r = tv.to_r(pc);
    let kl_r = kl.to_r(pc);
    out.set(0, tv_r).stop();
    out.set(1, kl_r).stop();
    out
}

#[allow(clippy::too_many_arguments)]
fn crp_vs_lorenz_mc<K, F>(
    n_items: i32,
    n_samples: usize,
    crp_concentration: f64,
    crp_discount: f64,
    fixed_k: Option<usize>,
    log_p_k_lorenz: Option<&[f64]>,
    log_p_k_crp: Option<&[f64]>,
    mut build_dist: F,
) -> (f64, f64)
where
    K: lorenz_ip_distr::KernelFactory,
    F: FnMut(usize) -> Option<LorenzIp<K>>,
{
    let n_items_usize: usize = n_items.try_into().unwrap();
    let mut dists: Vec<Option<LorenzIp<K>>> = Vec::with_capacity(n_items_usize + 1);
    dists.resize_with(n_items_usize + 1, || None);
    let mut attempted = vec![false; n_items_usize + 1];

    let seed = R::random_bytes::<16>();
    let mut rng = Pcg64Mcg::from_seed(seed);

    let mut counts = vec![0usize; n_items_usize + 1];
    for _ in 0..n_samples {
        let k = sample_crp_num_clusters(n_items_usize, crp_concentration, crp_discount, &mut rng);
        let k = k as usize;
        if k == 0 || k > n_items_usize {
            continue;
        }
        counts[k] += 1;
    }

    let mut tv_sum = 0.0;
    let mut kl_sum = 0.0;
    let mut total = 0usize;

    for k in 1..=n_items_usize {
        let count = counts[k];
        if count == 0 {
            continue;
        }
        let samples = sample_crp_integer_partition(
            count,
            n_items_usize,
            k,
            crp_concentration,
            crp_discount,
            &mut rng,
        );
        for i in 0..count {
            let mut partition_usize = Vec::with_capacity(k);
            let mut partition_i64 = Vec::with_capacity(k);
            for j in 0..k {
                let idx = i + count * j;
                let val = samples[idx];
                partition_usize.push(val);
                partition_i64.push(val as i64);
            }

            let log_p = log_crp_size_profile(
                &partition_usize,
                n_items_usize,
                crp_concentration,
                crp_discount,
            );

            let mut log_p_k = f64::NEG_INFINITY;
            if let Some(fixed) = fixed_k {
                if k == fixed {
                    if let Some(weights) = log_p_k_lorenz {
                        if k <= weights.len() {
                            log_p_k = weights[k - 1];
                        }
                    } else if let Some(weights) = log_p_k_crp {
                        if k <= weights.len() {
                            log_p_k = weights[k - 1];
                        }
                    } else {
                        log_p_k = 0.0;
                    }
                }
            } else if let Some(weights) = log_p_k_lorenz {
                if k <= weights.len() {
                    log_p_k = weights[k - 1];
                }
            } else if let Some(weights) = log_p_k_crp
                && k <= weights.len()
            {
                log_p_k = weights[k - 1];
            }

            let log_q = if log_p_k.is_finite() {
                let log_pmf = get_lorenz_log_pmf(
                    k,
                    &partition_i64,
                    &mut dists,
                    &mut attempted,
                    &mut build_dist,
                );
                if log_pmf.is_finite() {
                    log_p_k + log_pmf
                } else {
                    f64::NEG_INFINITY
                }
            } else {
                f64::NEG_INFINITY
            };

            let ratio = (log_q - log_p).exp();
            tv_sum += (1.0 - ratio).abs();
            kl_sum += log_p - log_q;
            total += 1;
        }
    }

    let total_f = total as f64;
    let tv = 0.5 * tv_sum / total_f;
    let kl = kl_sum / total_f;
    (tv, kl)
}

#[roxido]
fn crp_vs_lorenz(
    n_items: i32,
    target: &RObject,
    concentration: &RObject,
    log_skew: &RObject,
    tail_shape: &RObject,
    buffer_c: &RObject,
    crp_concentration: f64,
    crp_discount: f64,
    n_samples: i32,
    k_weights: &RObject,
    log_k_weights: bool,
) {
    if let Some(msg) = validate_crp_params(n_items, crp_concentration, crp_discount) {
        stop!("{}", msg);
    }
    if n_samples <= 0 {
        stop!("n_samples must be positive.");
    }

    let n_items_usize: usize = n_items.try_into().unwrap();

    let target_param = parse_target(target, pc);
    let buffer_c = parse_buffer_c(buffer_c, pc);

    let mut fixed_k: Option<usize> = None;
    let mut fixed_source: Option<&'static str> = None;

    if let TargetParam::Values(weights) = &target_param {
        if weights.is_empty() {
            stop!("target must be non-empty.");
        }
        update_fixed_k(&mut fixed_k, &mut fixed_source, "target", weights.len());
    }

    if concentration.is_vector() {
        let conc_vec = concentration
            .as_vector()
            .stop_str("concentration should be a vector or cubic_spline");
        let conc_f64 = conc_vec.to_f64(pc);
        let len = conc_f64.len();
        if len > 1 {
            update_fixed_k(&mut fixed_k, &mut fixed_source, "concentration", len + 1);
        }
    }

    if !log_skew.is_null() && log_skew.is_vector() {
        let log_skew_vec = log_skew
            .as_vector()
            .stop_str("log_skew should be a vector or cubic_spline");
        let log_skew_f64 = log_skew_vec.to_f64(pc);
        let len = log_skew_f64.len();
        if len > 1 {
            update_fixed_k(&mut fixed_k, &mut fixed_source, "log_skew", len + 1);
        }
    }

    if !tail_shape.is_null() && tail_shape.is_vector() {
        let tail_shape_vec = tail_shape
            .as_vector()
            .stop_str("tail_shape should be a vector or cubic_spline");
        let tail_shape_f64 = tail_shape_vec.to_f64(pc);
        let len = tail_shape_f64.len();
        if len > 1 {
            update_fixed_k(&mut fixed_k, &mut fixed_source, "tail_shape", len + 1);
        }
    }

    if let Some(k_fixed) = fixed_k
        && k_fixed > n_items_usize
    {
        stop!("n_clusters cannot exceed n_items.");
    }

    let log_p_k_lorenz = if k_weights.is_null() {
        None
    } else {
        let log_weights = parse_log_k_weights(k_weights, log_k_weights, pc);
        let k_max = n_items_usize.min(log_weights.len());
        Some(normalize_log_weights_allow_empty(&log_weights, k_max))
    };

    let log_p_k_crp = if k_weights.is_null() && fixed_k.is_none() {
        let log_stirling = compute_log_stirling_row(n_items_usize, n_items_usize, crp_discount);
        let log_denom = log_rising_factorial(crp_concentration, n_items_usize.saturating_sub(1));
        let mut log_p_k = vec![f64::NEG_INFINITY; n_items_usize];
        for k in 1..=n_items_usize {
            let log_s = log_stirling[k];
            let log_c = log_concentration_factor(crp_concentration, crp_discount, k);
            if log_s.is_finite() && log_c.is_finite() {
                log_p_k[k - 1] = log_s + log_c - log_denom;
            }
        }
        let log_norm = log_sum_exp(&log_p_k);
        if log_norm.is_finite() {
            for val in &mut log_p_k {
                *val -= log_norm;
            }
        }
        Some(log_p_k)
    } else {
        None
    };

    let n_samples = n_samples as usize;

    let kernel_is_tadpole = !tail_shape.is_null();
    let (tv, kl) = if !kernel_is_tadpole {
        crp_vs_lorenz_mc(
            n_items,
            n_samples,
            crp_concentration,
            crp_discount,
            fixed_k,
            log_p_k_lorenz.as_deref(),
            log_p_k_crp.as_deref(),
            |k| {
                if let Some(fixed) = fixed_k
                    && k != fixed
                {
                    return None;
                }
                let mut builder = TidalLorenzIpBuilder::new(n_items.into());

                builder = match &target_param {
                    TargetParam::Values(weights) => builder.lorenz_weights(weights.clone()).stop(),
                    TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve.clone()).stop(),
                    TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve.clone()).stop(),
                };

                let concentration_param = parse_concentration(concentration, k, pc);
                builder = match concentration_param {
                    ConcentrationParam::Values(values) => {
                        builder.concentration_values(values).stop()
                    }
                    ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
                };

                if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
                    builder = match param {
                        KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                        KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
                    };
                }

                builder = builder.buffer_scale(buffer_c).stop();
                builder.build(k).ok()
            },
        )
    } else {
        crp_vs_lorenz_mc(
            n_items,
            n_samples,
            crp_concentration,
            crp_discount,
            fixed_k,
            log_p_k_lorenz.as_deref(),
            log_p_k_crp.as_deref(),
            |k| {
                if let Some(fixed) = fixed_k
                    && k != fixed
                {
                    return None;
                }
                let mut builder = TadpoleLorenzIpBuilder::new(n_items.into());

                builder = match &target_param {
                    TargetParam::Values(weights) => builder.lorenz_weights(weights.clone()).stop(),
                    TargetParam::LorenzLinear(curve) => builder.lorenz_curve(curve.clone()).stop(),
                    TargetParam::LorenzISpline(curve) => builder.lorenz_curve(curve.clone()).stop(),
                };

                let concentration_param = parse_concentration(concentration, k, pc);
                builder = match concentration_param {
                    ConcentrationParam::Values(values) => {
                        builder.concentration_values(values).stop()
                    }
                    ConcentrationParam::Curve(spline) => builder.concentration_curve(spline).stop(),
                };

                if let Some(param) = parse_kernel_param(log_skew, k, "log_skew", pc) {
                    builder = match param {
                        KernelParam::Values(values) => builder.log_skew_values(values).stop(),
                        KernelParam::Curve(spline) => builder.log_skew_curve(spline).stop(),
                    };
                }
                if let Some(param) = parse_kernel_param(tail_shape, k, "tail_shape", pc) {
                    builder = match param {
                        KernelParam::Values(values) => builder.tail_shape_values(values).stop(),
                        KernelParam::Curve(spline) => builder.tail_shape_curve(spline).stop(),
                    };
                }

                builder = builder.buffer_scale(buffer_c).stop();
                builder.build(k).ok()
            },
        )
    };

    let out = RList::with_names(&["tv", "kl"], pc);
    let tv_r = tv.to_r(pc);
    let kl_r = kl.to_r(pc);
    out.set(0, tv_r).stop();
    out.set(1, kl_r).stop();
    out
}

#[roxido]
fn enumerate_integer_partitions(n_items: i32, n_clusters: i32) {
    if n_items < 0 {
        stop!("Number of items must be non-negative.");
    }
    if n_clusters < 0 {
        stop!("Number of clusters must be non-negative.");
    }
    if n_items == 0 {
        if n_clusters == 0 {
            return RMatrix::<i32>::from_value(0, 1, 0, pc);
        }
        stop!(
            "Number of clusters ({}) cannot exceed number of items ({}).",
            n_clusters,
            n_items
        );
    }
    if n_clusters == 0 {
        stop!("Number of clusters must be at least 1.");
    }
    if n_clusters > n_items {
        stop!(
            "Number of clusters ({}) cannot exceed number of items ({}).",
            n_clusters,
            n_items
        );
    }

    let n: i64 = n_items.into();
    let k: usize = n_clusters.try_into().unwrap();

    let partitions = rust_enumerate_integer_partitions(n, k);
    let n_partitions = partitions.len();

    if n_partitions == 0 {
        // Return empty matrix with 0 rows and k columns
        return RMatrix::<i32>::from_value(0, 0, k, pc);
    }

    let out = RMatrix::<i32>::from_value(0, n_partitions, k, pc);

    for (i, partition) in partitions.iter().enumerate() {
        for (j, &val) in partition.iter().enumerate() {
            out.set(i, j, i32::try_from(val).stop()).stop();
        }
    }

    out
}

#[roxido]
fn enumerate_partitions(n_items: i32, n_clusters: &RObject) {
    if n_items < 0 {
        stop!("Number of items must be non-negative.");
    }
    if n_items == 0 {
        if n_clusters.is_null() {
            return RMatrix::<i32>::from_value(0, 1, 0, pc);
        }
        let k_vec = n_clusters
            .as_vector()
            .stop_str("n_clusters should be a scalar");
        let k_vec = k_vec.to_i32(pc);
        let k = k_vec.get(0).stop();
        if k < 0 {
            stop!("Number of clusters must be non-negative.");
        }
        if k == 0 {
            return RMatrix::<i32>::from_value(0, 1, 0, pc);
        }
        stop!(
            "Number of clusters ({}) cannot exceed number of items ({}).",
            k,
            n_items
        );
    }

    let n: usize = n_items.try_into().unwrap();

    // Determine which k values to enumerate
    let k_values: Vec<usize> = if n_clusters.is_null() {
        // NULL means enumerate for all k from 1 to n
        (1..=n).collect()
    } else {
        // Single k value
        let k_vec = n_clusters
            .as_vector()
            .stop_str("n_clusters should be a scalar");
        let k_vec = k_vec.to_i32(pc);
        let k = k_vec.get(0).stop();
        if k <= 0 {
            stop!("Number of clusters must be at least 1.");
        }
        if k > n_items {
            stop!(
                "Number of clusters ({}) cannot exceed number of items ({}).",
                k,
                n_items
            );
        }
        vec![k as usize]
    };

    // Collect all partitions
    let mut all_data: Vec<usize> = Vec::new();
    let mut total_partitions: usize = 0;

    for k in k_values {
        let (data, n_partitions, _) = enumerate_set_partitions(n, k);
        all_data.extend(data);
        total_partitions += n_partitions;
    }

    if total_partitions == 0 {
        // Return empty matrix with 0 rows and n columns
        return RMatrix::<i32>::from_value(0, 0, n, pc);
    }

    // Create matrix: total_partitions rows x n columns
    let out = RMatrix::<i32>::from_value(0, total_partitions, n, pc);

    for i in 0..total_partitions {
        for j in 0..n {
            out.set(i, j, all_data[i * n + j] as i32).stop();
        }
    }

    out
}

#[roxido]
fn rcrpip(n: usize, n_items: i32, n_clusters: i32, discount: f64) {
    if n_items <= 0 {
        stop!("Number of items must be at least 1.");
    }
    if n_clusters <= 0 {
        stop!("Number of clusters must be at least 1.");
    }
    if n_clusters > n_items {
        stop!("Number of clusters cannot exceed number of items.");
    }
    let n_items: usize = n_items.try_into().unwrap();
    let k: usize = n_clusters.try_into().unwrap();
    if !(0.0..1.0).contains(&discount) {
        stop!("Discount must be in [0,1).");
    }
    let out = RMatrix::<i32>::from_value(0, n, k, pc);

    if k == 0 {
        return out;
    }

    let seed = R::random_bytes::<16>();
    let mut rng = Pcg64Mcg::from_seed(seed);

    let x = sample_crp_integer_partition(n, n_items, k, 1.0, discount, &mut rng);
    let slice = out.slice_mut();
    for (y, x) in slice.iter_mut().zip(x.iter()) {
        *y = i32::try_from(*x).stop();
    }
    out
}

#[roxido]
fn dcrpip(x: &RObject, n_items: i32, n_clusters: i32, discount: f64, log: bool) {
    if n_items <= 0 {
        stop!("Number of items must be at least 1.");
    }
    if n_clusters <= 0 {
        stop!("Number of clusters must be at least 1.");
    }
    if n_clusters > n_items {
        stop!("Number of clusters cannot exceed number of items.");
    }
    if !(0.0..1.0).contains(&discount) {
        stop!("Discount must be in [0,1).");
    }

    let n_items_usize: usize = n_items.try_into().unwrap();
    let k: usize = n_clusters.try_into().unwrap();

    let (n_partitions, mut x_data): (usize, Vec<i64>) = if x.is_matrix() {
        let x_mat = x.as_matrix().stop_str("x should be a matrix");
        let x_mat = x_mat.to_i32(pc);
        let nrow = x_mat.nrow();
        let ncol = x_mat.ncol();
        if ncol != k {
            stop!(
                "Number of columns ({}) does not match number of clusters ({})",
                ncol,
                k
            );
        }
        let slice = x_mat.slice();
        let mut data = Vec::with_capacity(nrow * ncol);
        for i in 0..nrow {
            for j in 0..ncol {
                data.push(slice[j * nrow + i] as i64);
            }
        }
        (nrow, data)
    } else {
        let x_vec = x.as_vector().stop_str("x should be a vector or matrix");
        let x_vec = x_vec.to_i32(pc);
        if x_vec.len() != k {
            stop!(
                "Length of x ({}) does not match number of clusters ({})",
                x_vec.len(),
                k
            );
        }
        let data: Vec<i64> = x_vec.slice().iter().map(|&v| v as i64).collect();
        (1, data)
    };

    if k > 1 {
        for i in 0..n_partitions {
            x_data[i * k..(i + 1) * k].sort_unstable();
        }
    }

    let out = RVector::from_value(f64::NAN, n_partitions, pc);
    if n_partitions == 0 {
        return out;
    }

    let log_stirling = compute_log_stirling_row(n_items_usize, k, discount);
    let log_s_k = *log_stirling.get(k).unwrap_or(&f64::NEG_INFINITY);

    for i in 0..n_partitions {
        let partition = &x_data[i * k..(i + 1) * k];
        let log_prob = if log_s_k.is_finite() {
            let mut sum = 0usize;
            let mut invalid = false;
            let mut part_usize = Vec::with_capacity(k);
            for &val in partition {
                if val <= 0 {
                    invalid = true;
                    break;
                }
                let u = val as usize;
                sum += u;
                part_usize.push(u);
            }
            if invalid || sum != n_items_usize {
                f64::NEG_INFINITY
            } else {
                log_crp_integer_partition_pmf_with_log_s(
                    &part_usize,
                    n_items_usize,
                    k,
                    discount,
                    log_s_k,
                )
            }
        } else {
            f64::NEG_INFINITY
        };
        let result = if log { log_prob } else { log_prob.exp() };
        out.set(i, result).stop();
    }

    out
}
#[roxido]
fn n_partitions(n_items: i32, n_clusters: &RObject, integer_partition: &RObject, log: bool) {
    if n_items < 0 {
        stop!("Number of items must be non-negative.");
    }
    let n = n_items as usize;

    // Parse n_clusters if provided
    let k_opt: Option<usize> = if n_clusters.is_null() {
        None
    } else {
        let k_vec = n_clusters
            .as_vector()
            .stop_str("n_clusters should be a scalar");
        let k_vec = k_vec.to_i32(pc);
        let k = k_vec.get(0).stop();
        if k < 0 {
            stop!("Number of clusters must be non-negative.");
        }
        if n > 0 && k as usize > n {
            stop!(
                "Number of clusters ({}) cannot exceed number of items ({}).",
                k,
                n
            );
        }
        Some(k as usize)
    };

    // Parse integer_partition if provided
    let partition_opt: Option<Vec<usize>> = if integer_partition.is_null() {
        None
    } else {
        let vec = integer_partition
            .as_vector()
            .stop_str("integer_partition should be a vector");
        let vec = vec.to_i32(pc);
        let partition: Vec<usize> = vec
            .slice()
            .iter()
            .map(|&v| {
                if v <= 0 {
                    stop!("All partition elements must be positive.");
                }
                v as usize
            })
            .collect();
        Some(partition)
    };

    // Compute the log result based on which arguments were provided
    let log_result = match (&k_opt, &partition_opt) {
        (None, None) => {
            // Bell number B(n)
            bell_number_log(n)
        }
        (Some(k), None) => {
            // Stirling number S(n, k)
            stirling_second_kind_log(n, *k)
        }
        (Some(k), Some(partition)) => {
            // Multinomial coefficient for set partitions
            if partition.len() != *k {
                stop!(
                    "Length of integer_partition ({}) must equal n_clusters ({}).",
                    partition.len(),
                    k
                );
            }
            multinomial_partition_log(n, partition)
        }
        (None, Some(_)) => {
            stop!("If integer_partition is provided, n_clusters must also be provided.");
        }
    };

    // Return result in appropriate format
    if log {
        log_result
    } else {
        log_result.exp().round()
    }
}

// ============================================================================
// Two-parameter CRP (Pitman-Yor) number of clusters distribution
// ============================================================================

/// Validate CRP parameters and return error message if invalid.
fn validate_crp_params(n_items: i32, concentration: f64, discount: f64) -> Option<&'static str> {
    if n_items < 1 {
        return Some("n_items must be at least 1.");
    }
    if !(0.0..1.0).contains(&discount) {
        return Some("discount must be in [0, 1).");
    }
    if concentration <= -discount {
        return Some("concentration must be greater than -discount.");
    }
    None
}

/// Sample the number of clusters from the two-parameter CRP using sequential simulation.
/// This is O(n_items) per sample.
fn sample_crp_num_clusters<R: rand::Rng + ?Sized>(
    n_items: usize,
    concentration: f64,
    discount: f64,
    rng: &mut R,
) -> i32 {
    if n_items == 0 {
        return 0;
    }

    let mut k = 1i32; // First customer always starts a table

    for i in 1..n_items {
        // Probability of new table: (concentration + k * discount) / (concentration + i)
        let p_new = (concentration + (k as f64) * discount) / (concentration + (i as f64));
        if rng.random::<f64>() < p_new {
            k += 1;
        }
    }

    k
}

/// Compute log of generalized Stirling numbers S_σ(N, k) for k = 0, 1, ..., k_max.
/// Uses two-row DP to save memory: O(k_max) space instead of O(N * k_max).
/// Returns a vector where result[k] = log S_σ(N, k).
fn compute_log_stirling_row(n_items: usize, k_max: usize, discount: f64) -> Vec<f64> {
    let cols = k_max + 1;
    let mut prev_row = vec![f64::NEG_INFINITY; cols];
    let mut curr_row = vec![f64::NEG_INFINITY; cols];

    // Base case: S_σ(0, 0) = 1
    prev_row[0] = 0.0;

    for n in 0..n_items {
        // Reset current row
        for value in curr_row.iter_mut() {
            *value = f64::NEG_INFINITY;
        }

        let max_k_for_next = (n + 1).min(k_max);
        for k in 0..=max_k_for_next {
            // Recurrence: S_σ(n+1, k) = (n - σ*k) * S_σ(n, k) + S_σ(n, k-1)
            let mut log_term1 = f64::NEG_INFINITY;
            let mut log_term2 = f64::NEG_INFINITY;

            // term1: (n - σ*k) * S_σ(n, k)
            if k <= n {
                let log_s_nk = prev_row[k];
                if log_s_nk.is_finite() {
                    let coeff = n as f64 - discount * k as f64;
                    if coeff > 0.0 {
                        log_term1 = coeff.ln() + log_s_nk;
                    }
                }
            }

            // term2: S_σ(n, k-1)
            if k > 0 && (k - 1) <= n {
                let log_s_nk1 = prev_row[k - 1];
                if log_s_nk1.is_finite() {
                    log_term2 = log_s_nk1;
                }
            }

            curr_row[k] = log_sum_exp_2(log_term1, log_term2);
        }

        std::mem::swap(&mut prev_row, &mut curr_row);
    }

    prev_row
}

/// Numerically stable log(exp(a) + exp(b)).
fn log_sum_exp_2(a: f64, b: f64) -> f64 {
    if a.is_infinite() && b.is_infinite() {
        return f64::NEG_INFINITY;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// Compute log of the concentration factor (α+σ)(α+2σ)...(α+(k-1)σ)
/// This is (α+σ)^{k-1}_σ in Pitman's notation (k-1 terms, not k terms).
/// Returns 0 (log of 1) for k <= 1 (empty product).
/// Returns NEG_INFINITY if any factor is <= 0.
fn log_concentration_factor(concentration: f64, discount: f64, k: usize) -> f64 {
    if k <= 1 {
        return 0.0; // Empty product = 1 for k=0 or k=1
    }

    let mut log_sum = 0.0;
    for i in 1..k {
        // Product from i=1 to k-1, giving factors: α+σ, α+2σ, ..., α+(k-1)σ
        let factor = concentration + (i as f64) * discount;
        if factor <= 0.0 {
            return f64::NEG_INFINITY;
        }
        log_sum += factor.ln();
    }
    log_sum
}

/// Compute log (α+1)_{N-1} = log((α+1)(α+2)...(α+N-1))
/// = log Γ(α+N) - log Γ(α+1)
fn log_rising_factorial(alpha: f64, n_minus_1: usize) -> f64 {
    if n_minus_1 == 0 {
        return 0.0; // Empty product = 1
    }

    // Use lgamma for numerical stability
    // (α+1)_{N-1} = Γ(α+N) / Γ(α+1)
    let n = n_minus_1 + 1;
    lgamma(alpha + n as f64) - lgamma(alpha + 1.0)
}

/// Compute log (a)_m = log(a(a+1)...(a+m-1)) for m >= 0.
fn log_rising_factorial_general(a: f64, m: usize) -> f64 {
    if m == 0 {
        return 0.0;
    }
    lgamma(a + m as f64) - lgamma(a)
}

/// Log gamma function using Lanczos approximation for positive values,
/// reflection formula for negative non-integer values.
fn lgamma(x: f64) -> f64 {
    // For most cases, use the standard library's implementation
    // via the gamma function relationship
    if x > 0.0 {
        // Use Stirling's approximation for large x, direct computation otherwise
        lgamma_positive(x)
    } else if x == x.floor() {
        // Negative integer: gamma is undefined (pole)
        f64::NEG_INFINITY
    } else {
        // Reflection formula: Γ(x)Γ(1-x) = π/sin(πx)
        // log|Γ(x)| = log(π) - log|sin(πx)| - log|Γ(1-x)|
        use std::f64::consts::PI;
        let sin_val = (PI * x).sin().abs();
        if sin_val == 0.0 {
            f64::NEG_INFINITY
        } else {
            PI.ln() - sin_val.ln() - lgamma_positive(1.0 - x)
        }
    }
}

/// Log gamma for positive x using Lanczos approximation.
fn lgamma_positive(x: f64) -> f64 {
    // Lanczos coefficients for g=7
    const G: f64 = 7.0;
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];

    if x < 0.5 {
        // Use reflection formula
        use std::f64::consts::PI;
        let sin_val = (PI * x).sin();
        if sin_val == 0.0 {
            return f64::INFINITY;
        }
        return PI.ln() - sin_val.abs().ln() - lgamma_positive(1.0 - x);
    }

    let x = x - 1.0;
    let mut a = COEFFICIENTS[0];
    for (i, &coeff) in COEFFICIENTS.iter().enumerate().skip(1) {
        a += coeff / (x + i as f64);
    }

    let t = x + G + 0.5;
    use std::f64::consts::PI;
    0.5 * (2.0 * PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

#[roxido]
fn rcrpk(n: usize, n_items: i32, concentration: f64, discount: f64) {
    // Validate parameters
    if let Some(msg) = validate_crp_params(n_items, concentration, discount) {
        stop!("{}", msg);
    }

    let n_items_usize: usize = n_items.try_into().unwrap();
    let result = RVector::<i32>::new(n, pc);

    // RNG seeded from R's RNG
    let seed = R::random_bytes::<16>();
    let mut rng = Pcg64Mcg::from_seed(seed);

    for i in 0..n {
        let k = sample_crp_num_clusters(n_items_usize, concentration, discount, &mut rng);
        result.set(i, k).stop();
    }

    result
}

#[roxido]
fn dcrpk(x: &RObject, n_items: i32, concentration: f64, discount: f64, log: bool) {
    // Validate parameters
    if let Some(msg) = validate_crp_params(n_items, concentration, discount) {
        stop!("{}", msg);
    }

    let n_items_usize: usize = n_items.try_into().unwrap();

    // Parse x as a vector of integers
    let x_vec = x.as_vector().stop_str("x should be a vector");
    let x_i32 = x_vec.to_i32(pc);
    let x_slice = x_i32.slice();
    let n_x = x_slice.len();

    if n_x == 0 {
        return RVector::<f64>::new(0, pc);
    }

    // Find the maximum valid k value to determine how much computation we need
    let max_k: usize = x_slice
        .iter()
        .filter(|&&k| k >= 1 && k <= n_items)
        .map(|&k| k as usize)
        .max()
        .unwrap_or(0);

    // If no valid k values, just return appropriate values
    let result = RVector::<f64>::new(n_x, pc);

    if max_k == 0 {
        // All x values are out of range
        let out_val = if log { f64::NEG_INFINITY } else { 0.0 };
        for (i, _) in x_slice.iter().enumerate() {
            result.set(i, out_val).stop();
        }
        return result;
    }

    // Compute log S_σ(N, k) for k = 0, 1, ..., max_k
    let log_stirling = compute_log_stirling_row(n_items_usize, max_k, discount);

    // Compute log denominator: log (concentration+1)_{N-1}
    let log_denom = log_rising_factorial(concentration, n_items_usize.saturating_sub(1));

    // Compute log PMF for each x value
    for (i, &k) in x_slice.iter().enumerate() {
        let log_prob = if k < 1 || k > n_items {
            // Out of range
            f64::NEG_INFINITY
        } else {
            let k_usize = k as usize;

            // log P(K = k) = log S_σ(N, k) + log C(α, σ, k) - log (α+1)_{N-1}
            let log_s = log_stirling[k_usize];
            let log_c = log_concentration_factor(concentration, discount, k_usize);

            if log_s.is_finite() && log_c.is_finite() {
                log_s + log_c - log_denom
            } else {
                f64::NEG_INFINITY
            }
        };

        let out_val = if log { log_prob } else { log_prob.exp() };
        result.set(i, out_val).stop();
    }

    result
}

#[roxido]
fn n_integer_partitions(n_items: i32, n_clusters: &RObject, log: bool) {
    if n_items < 0 {
        stop!("Number of items must be non-negative.");
    }
    let n = n_items as usize;

    // Parse n_clusters if provided
    let k_opt: Option<usize> = if n_clusters.is_null() {
        None
    } else {
        let k_vec = n_clusters
            .as_vector()
            .stop_str("n_clusters should be a scalar");
        let k_vec = k_vec.to_i32(pc);
        let k = k_vec.get(0).stop();
        if k < 0 {
            stop!("Number of clusters must be non-negative.");
        }
        if n > 0 && k as usize > n {
            stop!(
                "Number of clusters ({}) cannot exceed number of items ({}).",
                k,
                n
            );
        }
        Some(k as usize)
    };

    // Compute the log result
    let log_result = match k_opt {
        None => {
            // Total number of integer partitions p(n)
            integer_partition_count_log(n)
        }
        Some(k) => {
            // Number of integer partitions into exactly k parts p(n, k)
            integer_partition_count_with_parts_log(n, k)
        }
    };

    // Return result in appropriate format
    if log {
        log_result
    } else {
        log_result.exp().round()
    }
}
