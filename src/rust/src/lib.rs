roxido_registration!();
use roxido::*;

use dahl_randompartition::clust::Clustering;
use dahl_randompartition::cpp::CppParameters;
use dahl_randompartition::crp::CrpParameters;
use dahl_randompartition::distr::PredictiveProbabilityFunction;
use dahl_randompartition::distr::{PartitionSampler, ProbabilityMassFunction};
use dahl_randompartition::epa::{EpaParameters, SquareMatrixBorrower};
use dahl_randompartition::fixed::FixedPartitionParameters;
use dahl_randompartition::jlp::JlpParameters;
use dahl_randompartition::lsp::LspParameters;
use dahl_randompartition::mcmc::{update_neal_algorithm3, update_neal_algorithm8};
use dahl_randompartition::perm::Permutation;
use dahl_randompartition::prelude::*;
use dahl_randompartition::shrink::Shrinkage;
use dahl_randompartition::sp::SpParameters;
use dahl_randompartition::up::UpParameters;
use rand::distr::{Distribution, Uniform};
use rand::{Rng, RngExt, SeedableRng};
use rand_distr::{Beta as BetaRNG, Gamma as GammaRNG};
use rand_pcg::Pcg64Mcg;
use std::convert::{TryFrom, TryInto};
use std::ptr::NonNull;

#[roxido]
fn new_FixedPartitionParameters(anchor: &RVector) {
    let anchor = Clustering::from_slice(anchor.to_i32(pc).slice());
    let p = FixedPartitionParameters::new(anchor);
    RExternalPtr::encode(p, "fixed", pc)
}

#[roxido]
fn new_UpParameters(n_items: usize) {
    let p = UpParameters::new(n_items);
    RExternalPtr::encode(p, "up", pc)
}

#[roxido]
fn new_JlpParameters(concentration: f64, permutation: &RVector) {
    let permutation = mk_permutation(permutation, pc);
    let concentration = Concentration::new(concentration).stop_str("Invalid concentration value");
    let p = JlpParameters::new(permutation.n_items(), concentration, permutation)
        .stop_str("Invalid Jensen Liu parametrization");
    RExternalPtr::encode(p, "jlp", pc)
}

#[roxido]
fn new_CrpParameters(n_items: usize, concentration: f64, discount: f64) {
    let discount = Discount::new(discount).stop_str("Invalid discount value");
    let concentration = Concentration::new_with_discount(concentration, discount)
        .stop_str("Invalid concentration value");
    let p = CrpParameters::new_with_discount(n_items, concentration, discount)
        .stop_str("Invalid CRP parametrization");
    RExternalPtr::encode(p, "crp", pc)
}

#[roxido]
fn new_EpaParameters(
    similarity: &mut RMatrix<f64>,
    permutation: &RVector,
    concentration: f64,
    discount: f64,
) {
    let ni = similarity.nrow();
    let slice = similarity.slice_mut();
    let similarity = SquareMatrixBorrower::from_slice(slice, ni);
    let permutation = mk_permutation(permutation, pc);
    let discount = Discount::new(discount).unwrap_or_else(|| stop!("Invalid discount value"));
    let concentration = Concentration::new_with_discount(concentration, discount)
        .unwrap_or_else(|| stop!("Invalid concentration value"));
    let p = EpaParameters::new(similarity, permutation, concentration, discount);
    RExternalPtr::encode(p, "epa", pc)
}

#[roxido]
fn new_LspParameters(anchor: &RVector, shrinkage: f64, concentration: f64, permutation: &RVector) {
    let anchor = Clustering::from_slice(anchor.to_i32(pc).slice());
    let shrinkage =
        ScalarShrinkage::new(shrinkage).unwrap_or_else(|| stop!("Invalid shrinkage value"));
    let concentration =
        Concentration::new(concentration).unwrap_or_else(|| stop!("Invalid concentration value"));
    let permutation = mk_permutation(permutation, pc);
    let p =
        LspParameters::new_with_shrinkage(anchor, shrinkage, concentration, permutation).unwrap();
    RExternalPtr::encode(p, "lsp", pc)
}

#[roxido]
fn new_CppParameters(anchor: &RVector, rate: f64, baseline: &RExternalPtr, use_vi: bool, a: f64) {
    let anchor = Clustering::from_slice(anchor.to_i32(pc).slice());
    macro_rules! distr_macro {
        ($tipe:ty, $label:literal) => {{
            let p = NonNull::new(baseline.address() as *mut $tipe).unwrap();
            let baseline = unsafe { p.as_ref().clone() };
            RExternalPtr::encode(
                CppParameters::new(anchor, rate, baseline, use_vi, a).unwrap(),
                &$label,
                pc,
            )
        }};
    }
    let tag = baseline.tag().as_scalar().stop();
    let name = tag.str(pc);
    match name {
        "up" => distr_macro!(UpParameters, "cpp-up"),
        "jlp" => distr_macro!(JlpParameters, "cpp-jlp"),
        "crp" => distr_macro!(CrpParameters, "cpp-crp"),
        _ => stop!("Unsupported distribution: {}", name),
    }
}

#[roxido]
fn new_SpParameters(
    anchor: &RVector,
    shrinkage: &RVector,
    permutation: &RVector,
    grit: f64,
    baseline: &RExternalPtr,
    shortcut: bool,
) {
    let anchor = Clustering::from_slice(anchor.to_i32(pc).slice());
    let shrinkage = Shrinkage::from(shrinkage.to_f64(pc).slice()).stop_str("Invalid shrinkage");
    let permutation = mk_permutation(permutation, pc);
    let grit = match Grit::new(grit) {
        Some(grit) => grit,
        None => stop!("Grit value out of range"),
    };
    macro_rules! distr_macro {
        ($tipe:ty, $label:literal) => {{
            let p = NonNull::new(baseline.address() as *mut $tipe).unwrap();
            let baseline = unsafe { p.as_ref().clone() };
            RExternalPtr::encode(
                SpParameters::new(anchor, shrinkage, permutation, grit, baseline, shortcut)
                    .stop_str("Invalid shrinkage partition parametrization"),
                $label,
                pc,
            )
        }};
    }
    let tag = baseline.tag().as_scalar().stop();
    let name = tag.str(pc);
    match name {
        "up" => distr_macro!(UpParameters, "sp-up"),
        "jlp" => distr_macro!(JlpParameters, "sp-jlp"),
        "crp" => distr_macro!(CrpParameters, "sp-crp"),
        _ => stop!("Unsupported distribution: {}", name),
    }
}

fn mk_permutation(permutation: &RVector, pc: &Pc) -> Permutation {
    let vector = permutation
        .to_i32(pc)
        .slice()
        .iter()
        .map(|x| *x as usize)
        .collect();
    Permutation::from_vector(vector).stop_str("Invalid permutation")
}

#[derive(PartialEq, Copy, Clone)]
enum RandomizeShrinkage {
    Fixed,
    Common,
    Cluster,
    Idiosyncratic,
}

#[roxido]
fn samplePartition(
    n_samples: usize,
    n_items: usize,
    prior: &RExternalPtr,
    randomize_permutation: bool,
    randomize_shrinkage: &str,
    randomize_grit: bool,
    shrinkage_shape: f64,
    shrinkage_rate: f64,
    grit_shape1: f64,
    grit_shape2: f64,
    n_cores: usize,
) {
    let n_cores = {
        if n_cores == 0 {
            num_cpus::get()
        } else {
            n_cores
        }
    };
    let np = n_samples.max(1);
    let np_per_core = np / n_cores;
    let np_extra = np % n_cores;
    let ni = n_items.max(1);
    let chunk_size = np_per_core * ni;
    let matrix_rval = RMatrix::<i32>::new(ni, np, pc);
    let mut stick = matrix_rval.slice_mut();
    let randomize_shrinkage = match randomize_shrinkage {
        "fixed" => RandomizeShrinkage::Fixed,
        "common" => RandomizeShrinkage::Common,
        "cluster" => RandomizeShrinkage::Cluster,
        "idiosyncratic" => RandomizeShrinkage::Idiosyncratic,
        _ => stop!("Unrecognized 'randomize_shrinkage' value"),
    };
    let shrinkage_shape = Shape::new(shrinkage_shape)
        .unwrap_or_else(|| stop!("'shrinkage_shape' should be greater than 0"));
    let shrinkage_rate = Rate::new(shrinkage_rate)
        .unwrap_or_else(|| stop!("'shrinkage_rate' should be greater than 0"));
    let grit_shape1 =
        Shape::new(grit_shape1).unwrap_or_else(|| stop!("'grit_shape1' should be greater than 0"));
    let grit_shape2 =
        Shape::new(grit_shape2).unwrap_or_else(|| stop!("'grit_shape2' should be greater than 0"));
    macro_rules! distr_macro {
        ($tipe: ty, $callback: tt) => {{
            let _ = crossbeam::scope(|s| {
                let mut nonnull = NonNull::new(prior.address() as *mut $tipe).unwrap();
                let distr = unsafe { nonnull.as_mut() };
                let mut plan = Vec::with_capacity(n_cores);
                let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
                for k in 0..n_cores - 1 {
                    let (left, right) =
                        stick.split_at_mut(chunk_size + if k < np_extra { ni } else { 0 });
                    plan.push((left, distr.clone(), rng.random::<u128>()));
                    stick = right;
                }
                plan.push((stick, distr.clone(), rng.random()));
                plan.into_iter().for_each(|mut p| {
                    s.spawn(move |_| {
                        let mut rng = Pcg64Mcg::new(p.2);
                        for j in 0..p.0.len() / ni {
                            #[allow(clippy::redundant_closure_call)]
                            ($callback)(&mut p.1, &mut rng);
                            let clust = p.1.sample(&mut rng).standardize();
                            let labels = clust.allocation();
                            for i in 0..ni {
                                p.0[ni * j + i] = (labels[i] + 1).try_into().unwrap();
                            }
                        }
                    });
                });
            });
        }};
    }
    macro_rules! cpp_macro {
        ($tipe: ty) => {{
            match randomize_shrinkage {
                RandomizeShrinkage::Common | RandomizeShrinkage::Fixed => {}
                _ => stop!("Unsupported randomize_shrinkage for this distribution"),
            }
            let distr = prior.decode_ref::<$tipe>();
            let mut plan = Vec::with_capacity(n_cores);
            let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
            for k in 0..n_cores - 1 {
                let (left, right) =
                    stick.split_at_mut(chunk_size + if k < np_extra { ni } else { 0 });
                plan.push((left, distr.clone(), rng.random::<u128>()));
                stick = right;
            }
            plan.push((stick, distr.clone(), rng.random()));
            let log_like = |_i: usize, _indices: &[usize]| 0.0;
            let nup = 1;
            let _ = crossbeam::scope(|s| {
                plan.into_iter().for_each(|p| {
                    s.spawn(move |_| {
                        let perm = Permutation::natural_and_fixed(n_items);
                        let mut rng = Pcg64Mcg::new(p.2);
                        let mut current = distr.anchor.clone();
                        for i in 0..p.0.len() / ni {
                            update_neal_algorithm3(
                                nup,
                                &mut current,
                                &perm,
                                distr,
                                &log_like,
                                &mut rng,
                            );
                            current.relabel_into_slice(
                                1,
                                &mut p.0[(i * n_items)..((i + 1) * n_items)],
                            );
                        }
                    });
                });
            });
        }};
    }
    fn mk_lambda_sp<D: PredictiveProbabilityFunction + Clone>(
        randomize_permutation: bool,
        randomize_shrinkage: RandomizeShrinkage,
        randomize_grit: bool,
        shrinkage_shape: Shape,
        shrinkage_rate: Rate,
        grit_shape1: Shape,
        grit_shape2: Shape,
    ) -> impl Fn(&mut SpParameters<D>, &mut Pcg64Mcg) {
        let beta = BetaRNG::new(grit_shape1.get(), grit_shape2.get()).unwrap();
        move |distr: &mut SpParameters<D>, rng: &mut Pcg64Mcg| {
            if randomize_permutation {
                distr.permutation.shuffle(rng);
            }
            match randomize_shrinkage {
                RandomizeShrinkage::Fixed => {}
                RandomizeShrinkage::Common => {
                    distr
                        .shrinkage
                        .randomize_common(shrinkage_shape, shrinkage_rate, rng)
                }
                RandomizeShrinkage::Cluster => distr.shrinkage.randomize_common_cluster(
                    shrinkage_shape,
                    shrinkage_rate,
                    &distr.anchor,
                    rng,
                ),
                RandomizeShrinkage::Idiosyncratic => {
                    distr
                        .shrinkage
                        .randomize_idiosyncratic(shrinkage_shape, shrinkage_rate, rng)
                }
            };
            if randomize_grit {
                distr.grit = Grit::new(beta.sample(rng)).unwrap()
            }
        }
    }
    let tag = prior.tag().as_scalar().stop();
    let name = tag.str(pc);
    match name {
        "fixed" => distr_macro!(
            FixedPartitionParameters,
            (|_distr: &mut FixedPartitionParameters, _rng: &mut Pcg64Mcg| {})
        ),
        "up" => distr_macro!(
            UpParameters,
            (|_distr: &mut UpParameters, _rng: &mut Pcg64Mcg| {})
        ),
        "jlp" => distr_macro!(
            JlpParameters,
            (|_distr: &mut JlpParameters, _rng: &mut Pcg64Mcg| {})
        ),
        "crp" => distr_macro!(
            CrpParameters,
            (|_distr: &mut CrpParameters, _rng: &mut Pcg64Mcg| {})
        ),
        "epa" => {
            if randomize_permutation {
                distr_macro!(
                    EpaParameters,
                    (|distr: &mut EpaParameters, rng: &mut Pcg64Mcg| {
                        distr.permutation.shuffle(rng);
                    })
                );
            } else {
                distr_macro!(
                    EpaParameters,
                    (|_distr: &mut EpaParameters, _rng: &mut Pcg64Mcg| {})
                )
            }
        }
        "lsp" => {
            match randomize_shrinkage {
                RandomizeShrinkage::Common | RandomizeShrinkage::Fixed => {}
                _ => stop!("Unsupported randomize_shrinkage for this distribution"),
            }
            if randomize_permutation && randomize_shrinkage == RandomizeShrinkage::Common {
                distr_macro!(
                    LspParameters,
                    (|distr: &mut LspParameters, rng: &mut Pcg64Mcg| {
                        distr.permutation.shuffle(rng);
                        let gamma =
                            GammaRNG::new(shrinkage_shape.get(), 1.0 / shrinkage_rate.get())
                                .unwrap();
                        distr.shrinkage = ScalarShrinkage::new(gamma.sample(rng)).unwrap();
                    })
                );
            } else if randomize_permutation {
                distr_macro!(
                    LspParameters,
                    (|distr: &mut LspParameters, rng: &mut Pcg64Mcg| {
                        distr.permutation.shuffle(rng);
                    })
                );
            } else if randomize_shrinkage == RandomizeShrinkage::Common {
                distr_macro!(
                    LspParameters,
                    (|distr: &mut LspParameters, rng: &mut Pcg64Mcg| {
                        let gamma =
                            GammaRNG::new(shrinkage_shape.get(), 1.0 / shrinkage_rate.get())
                                .unwrap();
                        distr.shrinkage = ScalarShrinkage::new(gamma.sample(rng)).unwrap();
                    })
                );
            } else {
                distr_macro!(
                    LspParameters,
                    (|_distr: &mut LspParameters, _rng: &mut Pcg64Mcg| {})
                );
            }
        }
        "cpp-up" => {
            cpp_macro!(CppParameters<UpParameters>);
        }
        "cpp-jlp" => {
            cpp_macro!(CppParameters<JlpParameters>);
        }
        "cpp-crp" => {
            cpp_macro!(CppParameters<CrpParameters>);
        }
        "sp-up" => {
            distr_macro!(
                SpParameters<UpParameters>,
                (mk_lambda_sp::<UpParameters>(
                    randomize_permutation,
                    randomize_shrinkage,
                    randomize_grit,
                    shrinkage_shape,
                    shrinkage_rate,
                    grit_shape1,
                    grit_shape2
                ))
            );
        }
        "sp-jlp" => {
            distr_macro!(
                SpParameters<JlpParameters>,
                (mk_lambda_sp::<JlpParameters>(
                    randomize_permutation,
                    randomize_shrinkage,
                    randomize_grit,
                    shrinkage_shape,
                    shrinkage_rate,
                    grit_shape1,
                    grit_shape2
                ))
            );
        }
        "sp-crp" => {
            distr_macro!(
                SpParameters<CrpParameters>,
                (mk_lambda_sp::<CrpParameters>(
                    randomize_permutation,
                    randomize_shrinkage,
                    randomize_grit,
                    shrinkage_shape,
                    shrinkage_rate,
                    grit_shape1,
                    grit_shape2
                ))
            );
        }
        _ => stop!("Unsupported distribution: {}", name),
    }
    matrix_rval.transpose(pc)
}

pub fn sample_into_slice<S: PartitionSampler, T: Rng, F: Fn(&mut S, &mut T)>(
    n_partitions: usize,
    n_items: usize,
    matrix: &mut [i32],
    rng: &mut T,
    distr: &mut S,
    callback: F,
) {
    for i in 0..n_partitions {
        callback(distr, rng);
        let p = distr.sample(rng).standardize();
        let labels = p.allocation();
        for j in 0..n_items {
            matrix[n_partitions * j + i] = (labels[j] + 1).try_into().unwrap();
        }
    }
}

#[roxido]
fn prPartition(partition: &RMatrix, prior: &RExternalPtr) {
    let partition = partition.to_i32(pc);
    let matrix = partition.slice();
    let np = partition.nrow();
    let ni = partition.ncol();
    let log_probs_rval = RVector::new(np, pc);
    let log_probs_slice = log_probs_rval.slice_mut();
    macro_rules! distr_macro {
        ($tipe:ty) => {{
            let p = prior.address() as *mut $tipe;
            let distr = unsafe { NonNull::new(p).unwrap().as_mut() };
            for i in 0..np {
                let mut target_labels = Vec::with_capacity(ni);
                for j in 0..ni {
                    target_labels.push((matrix[np * j + i]).try_into().unwrap());
                }
                let target = Clustering::from_vector(target_labels);
                log_probs_slice[i] = distr.log_pmf(&target);
            }
        }};
    }
    let tag = prior.tag().as_scalar().stop();
    let name = tag.str(pc);
    match name {
        "fixed" => distr_macro!(FixedPartitionParameters),
        "up" => distr_macro!(UpParameters),
        "jlp" => distr_macro!(JlpParameters),
        "crp" => distr_macro!(CrpParameters),
        "epa" => distr_macro!(EpaParameters),
        "lsp" => distr_macro!(LspParameters),
        "cpp-up" => distr_macro!(CppParameters<UpParameters>),
        "cpp-jlp" => distr_macro!(CppParameters<JlpParameters>),
        "cpp-crp" => distr_macro!(CppParameters<CrpParameters>),
        "sp-up" => distr_macro!(SpParameters<UpParameters>),
        "sp-jlp" => distr_macro!(SpParameters<JlpParameters>),
        "sp-crp" => distr_macro!(SpParameters<CrpParameters>),
        _ => stop!("Unsupported distribution: {}", name),
    }
    log_probs_rval
}

#[roxido]
fn nealAlgorithm(
    partition: &RVector,
    callback: &RObject,
    is_algorithm3: bool,
    n_updates_for_partition: usize,
    prior: &RExternalPtr,
) {
    let callback = callback.as_option().map(|s| s.as_function().stop());
    let nup = u32::try_from(n_updates_for_partition).unwrap();
    let mut current = Clustering::from_slice(partition.to_i32(pc).slice());
    let perm = Permutation::natural_and_fixed(current.n_items());
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let tag = prior.tag().as_scalar().stop();
    let name = tag.str(pc);
    if is_algorithm3 {
        macro_rules! distr_macro {
            ($tipe:ty, $log_like:expr) => {{
                let p = std::ptr::NonNull::new(prior.address() as *mut $tipe).unwrap();
                let p_ref = unsafe { p.as_ref() };
                update_neal_algorithm3(nup, &mut current, &perm, p_ref, &$log_like, &mut rng)
            }};
        }
        match callback {
            None => {
                let log_like = |_i: usize, _indices: &[usize]| 0.0;
                match name {
                    "fixed" => distr_macro!(FixedPartitionParameters, log_like),
                    "up" => distr_macro!(UpParameters, log_like),
                    "jlp" => distr_macro!(JlpParameters, log_like),
                    "crp" => distr_macro!(CrpParameters, log_like),
                    "epa" => distr_macro!(EpaParameters, log_like),
                    "lsp" => distr_macro!(LspParameters, log_like),
                    "cpp-up" => distr_macro!(CppParameters<UpParameters>, log_like),
                    "cpp-jlp" => distr_macro!(CppParameters<JlpParameters>, log_like),
                    "cpp-crp" => distr_macro!(CppParameters<CrpParameters>, log_like),
                    "sp-up" => distr_macro!(SpParameters<UpParameters>, log_like),
                    "sp-jlp" => distr_macro!(SpParameters<JlpParameters>, log_like),
                    "sp-crp" => distr_macro!(SpParameters<CrpParameters>, log_like),
                    _ => stop!("Unsupported distribution: {}", name),
                }
            }
            Some(func) => {
                let log_like = |i: usize, indices: &[usize]| {
                    let index = i32::try_from(i + 1).unwrap().to_r(pc);
                    let indices_plus_one = RVector::<i32>::new(indices.len(), pc);
                    let slice = indices_plus_one.slice_mut();
                    for (x, y) in slice.iter_mut().zip(indices.iter()) {
                        *x = i32::try_from(*y).unwrap() + 1;
                    }
                    match func.call2(index, indices_plus_one, pc) {
                        Ok(x) => x.as_scalar().stop().f64(),
                        Err(e) => stop!("Error in running user-supplied function with code: {}", e),
                    }
                };
                match name {
                    "fixed" => distr_macro!(FixedPartitionParameters, log_like),
                    "up" => distr_macro!(UpParameters, log_like),
                    "jlp" => distr_macro!(JlpParameters, log_like),
                    "crp" => distr_macro!(CrpParameters, log_like),
                    "epa" => distr_macro!(EpaParameters, log_like),
                    "lsp" => distr_macro!(LspParameters, log_like),
                    "cpp-up" => distr_macro!(CppParameters<UpParameters>, log_like),
                    "cpp-jlp" => distr_macro!(CppParameters<JlpParameters>, log_like),
                    "cpp-crp" => distr_macro!(CppParameters<CrpParameters>, log_like),
                    "sp-up" => distr_macro!(SpParameters<UpParameters>, log_like),
                    "sp-jlp" => distr_macro!(SpParameters<JlpParameters>, log_like),
                    "sp-crp" => distr_macro!(SpParameters<CrpParameters>, log_like),
                    _ => stop!("Unsupported distribution: {}", name),
                }
            }
        };
    } else {
        current.exclude_label(0);
        struct Dummy {}
        let cacher = |_item: usize| Dummy {};
        macro_rules! distr_macro {
            ($tipe:ty, $log_like:expr) => {{
                let p = std::ptr::NonNull::new(prior.address() as *mut $tipe).unwrap();
                let p_ref = unsafe { p.as_ref() };
                update_neal_algorithm8(
                    nup,
                    &mut current,
                    &perm,
                    p_ref,
                    &mut $log_like,
                    cacher,
                    &mut rng,
                )
            }};
        }
        match callback {
            None => {
                let mut log_like = |_i: usize, _dummy: &Dummy, _label: usize, _is_new: bool| 0.0;
                match name {
                    "fixed" => distr_macro!(FixedPartitionParameters, log_like),
                    "up" => distr_macro!(UpParameters, log_like),
                    "jlp" => distr_macro!(JlpParameters, log_like),
                    "crp" => distr_macro!(CrpParameters, log_like),
                    "epa" => distr_macro!(EpaParameters, log_like),
                    "lsp" => distr_macro!(LspParameters, log_like),
                    "cpp-up" => distr_macro!(CppParameters<UpParameters>, log_like),
                    "cpp-jlp" => distr_macro!(CppParameters<JlpParameters>, log_like),
                    "cpp-crp" => distr_macro!(CppParameters<CrpParameters>, log_like),
                    "sp-up" => distr_macro!(SpParameters<UpParameters>, log_like),
                    "sp-jlp" => distr_macro!(SpParameters<JlpParameters>, log_like),
                    "sp-crp" => distr_macro!(SpParameters<CrpParameters>, log_like),
                    _ => stop!("Unsupported distribution: {}", name),
                }
            }
            Some(func) => {
                let mut log_like = |i: usize, _dummy: &Dummy, label: usize, is_new: bool| {
                    let index = i32::try_from(i + 1).unwrap().to_r(pc);
                    let label = i32::try_from(label).unwrap().to_r(pc);
                    let is_new = is_new.to_r(pc);
                    match func.call3(index, label, is_new, pc) {
                        Ok(x) => x.as_scalar().stop().f64(),
                        Err(e) => {
                            stop!("Error in running user-supplied function with code: {}", e)
                        }
                    }
                };
                match name {
                    "fixed" => distr_macro!(FixedPartitionParameters, log_like),
                    "up" => distr_macro!(UpParameters, log_like),
                    "jlp" => distr_macro!(JlpParameters, log_like),
                    "crp" => distr_macro!(CrpParameters, log_like),
                    "epa" => distr_macro!(EpaParameters, log_like),
                    "lsp" => distr_macro!(LspParameters, log_like),
                    "cpp-up" => distr_macro!(CppParameters<UpParameters>, log_like),
                    "cpp-jlp" => distr_macro!(CppParameters<JlpParameters>, log_like),
                    "cpp-crp" => distr_macro!(CppParameters<CrpParameters>, log_like),
                    "sp-up" => distr_macro!(SpParameters<UpParameters>, log_like),
                    "sp-jlp" => distr_macro!(SpParameters<JlpParameters>, log_like),
                    "sp-crp" => distr_macro!(SpParameters<CrpParameters>, log_like),
                    _ => stop!("Unsupported distribution: {}", name),
                }
            }
        };
    };
    let (clustering, map) = Clustering::relabel(current.allocation(), 1, None, !is_algorithm3);
    let partition = RVector::<i32>::new(clustering.n_items(), pc);
    let partition_slice = partition.slice_mut();
    for (x, y) in partition_slice.iter_mut().zip(&clustering.allocation()[..]) {
        *x = i32::try_from(*y).unwrap()
    }
    let result = RList::with_names(&["partition", "map"], pc);
    result.set(0, partition).stop();
    result.set(1, R::null()).stop();
    if !is_algorithm3 {
        let slice = &map.unwrap()[1..];
        let rval = RVector::<i32>::new(slice.len(), pc);
        let rval_slice = rval.slice_mut();
        for (x, y) in rval_slice.iter_mut().zip(slice) {
            *x = (*y).try_into().unwrap();
        }
        result.set(1, rval).stop();
    }
    result
}

#[roxido]
fn slice_sampler(
    current: f64,
    target: &RFunction,
    w: f64,
    max_number_of_steps: f64,
    on_log_scale: bool,
) {
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let uniform = Uniform::new(0.0, 1.0).unwrap();
    let mut u = || uniform.sample(&mut rng);
    let xobj = 0.0.to_r(pc);
    let mut f = |x: f64| {
        xobj.set(x);
        target.call1(xobj, pc).stop().as_scalar().stop().f64()
    };
    if w <= 0.0 {
        stop!("Width w must be strictly positive");
    }
    let x = current;
    let fx = f(x);
    // Step 1 (slicing)
    let y = if on_log_scale {
        (u() * fx.exp()).ln()
    } else {
        u() * fx
    };
    // Step 2 (stepping out)
    let mut l = x - u() * w;
    let mut r = l + w;
    let max = max_number_of_steps;
    if !R::is_finite(max) {
        while y < f(l) {
            l -= w
        }
        while y < f(r) {
            r += w
        }
    } else if max > 0.0 {
        let mut j = (u() * max).floor() as u32;
        let mut k = max as u32 - 1 - j;
        while j > 0 && y < f(l) {
            l -= w;
            j -= 1;
        }
        while k > 0 && y < f(r) {
            r += w;
            k -= 1;
        }
    }
    // Step 3 (shrinkage)
    loop {
        let x1 = l + u() * (r - l);
        let fx1 = f(x1);
        if y < fx1 {
            return x1.to_r(pc);
        }
        if x1 < x {
            l = x1;
        } else {
            r = x1;
        }
    }
}
