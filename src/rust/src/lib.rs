// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();
use roxido::*;

mod clustering;
mod timers;

//const ALLOW_DEBUGGING: bool = true;
const ALLOW_DEBUGGING: bool = false;

type LabelType = u8; // At most 255 clusters supported
type CountType = u16; // At most 65,536 clusterings per splinter
type CountSummationType = u32;

fn check_shortcut() -> bool {
    matches!(std::env::var("DBD_CHECK_SHORTCUT"), Ok(e) if e.to_lowercase() == "true")
}

use crate::clustering::Clustering;
use dahl_randompartition::clust::Clustering as Partition;
use dahl_randompartition::crp::CrpParameters;
use dahl_randompartition::distr::{FullConditional, PartitionSampler};
use dahl_randompartition::mcmc::{update_neal_algorithm3, update_neal_algorithm8};
use dahl_randompartition::perm::Permutation;
use dahl_randompartition::prelude::*;
use ndarray::prelude::*;
use rand::prelude::*;
use rand_distr::StandardNormal;
use rand_pcg::Pcg64Mcg;
use rayon::prelude::*;
use statrs::function::gamma::ln_gamma;
use std::collections::HashMap;
use std::collections::HashSet;

struct Count {
    marginals: Vec<CountType>,
    joints: Array3<CountType>,
}

struct Counts {
    x: Vec<Count>,
}

#[derive(Debug)]
struct Splinter {
    n_items: CountType,
    mapping: Vec<Option<usize>>,
    clusterings: Array2<LabelType>,
    max_n_clusters_observed: LabelType,
}

#[roxido]
fn new_splinter(n_items_in_all: usize, items: &RVector<i32>, clusterings: &RMatrix) {
    let labels = canonize_labels(clusterings, 0, pc);
    let n_items = items.len();
    let n_samples = labels.len() / n_items;
    let clusterings = Array2::from_shape_vec((n_items, n_samples).f(), labels).unwrap();
    let max_n_clusters_observed = *clusterings.iter().max().unwrap() + 1;
    let mut mapping = vec![None; n_items_in_all];
    for (index, &item) in items.slice().iter().enumerate() {
        mapping[usize::try_from(item - 1).unwrap()] = Some(index);
    }
    let splinter = Splinter {
        n_items: n_items.try_into().unwrap(),
        mapping,
        clusterings,
        max_n_clusters_observed,
    };
    RExternalPtr::encode_full(splinter, R::null(), false, pc)
}

struct Splinters {
    x: Vec<Splinter>,
    sum_of_weights: f64,
}

impl Splinters {
    fn from_r(splinters: &RList) -> Splinters {
        if !ALLOW_DEBUGGING && check_shortcut() {
            panic!("Not compiled for checking shortcut. Please change ALLOW_DEBUGGING to true.")
        }
        let n_splinters = splinters.len();
        let mut s = Vec::with_capacity(n_splinters);
        let mut sum_of_weights = 0.0;
        for i in 0..n_splinters {
            let splinter = splinters
                .get(i)
                .stop()
                .as_external_ptr()
                .stop()
                .decode_val::<Splinter>()
                .stop();
            sum_of_weights += splinter.n_items as f64;
            s.push(splinter);
        }
        Splinters {
            x: s,
            sum_of_weights,
        }
    }

    fn tally(&self, estimate: &Clustering, max_n_clusters: LabelType) -> Counts {
        let max_n_clusters = if max_n_clusters == 0 {
            self.x
                .iter()
                .fold(0, |max, x| max.max(x.max_n_clusters_observed))
        } else {
            max_n_clusters
        };
        let mut counts = Vec::with_capacity(self.x.len());
        for splinter in &self.x {
            let marginals = vec![0; estimate.max_n_clusters().into()];
            let joints = Array3::default((
                max_n_clusters as usize,
                splinter.max_n_clusters_observed as usize,
                splinter.clusterings.ncols(),
            ));
            counts.push(Count { marginals, joints })
        }
        for (item, &label) in estimate.labels().iter().enumerate() {
            let label = label as usize;
            for (splinter, count) in self.x.iter().zip(&mut counts) {
                if let Some(alias) = splinter.mapping[item] {
                    count.marginals[label] += 1;
                    let s = splinter.clusterings.slice(s!(alias, ..));
                    //let s = splinter.clusterings.index_axis(Axis(0), alias);
                    for (sample_index, &label_observed) in s.iter().enumerate() {
                        let label_observed = label_observed as usize;
                        count.joints[(label, label_observed, sample_index)] += 1;
                    }
                }
            }
        }
        Counts { x: counts }
    }

    fn update_tallies(
        &self,
        counts: &mut Counts,
        item: usize,
        label_from: LabelType,
        label_to: LabelType,
    ) {
        let label_from = label_from as usize;
        let label_to = label_to as usize;
        for (splinter, count) in self.x.iter().zip(&mut counts.x) {
            if let Some(alias) = splinter.mapping[item] {
                count.marginals[label_from] -= 1;
                count.marginals[label_to] += 1;
                let s = splinter.clusterings.slice(s!(alias, ..));
                //let s = splinter.clusterings.index_axis(Axis(0), alias);
                for (sample_index, &label_observed) in s.iter().enumerate() {
                    let label_observed = label_observed as usize;
                    count.joints[(label_from, label_observed, sample_index)] -= 1;
                    count.joints[(label_to, label_observed, sample_index)] += 1;
                }
            }
        }
    }

    fn compute_expected_loss<T: LossComputer>(&self, counts: &Counts, loss: &T) -> f64 {
        let mut sum = 0.0;
        for (splinter, count) in self.x.iter().zip(&counts.x) {
            sum += (splinter.n_items as f64)
                * loss.compute_expectation(&count.marginals, &count.joints, splinter.n_items)
        }
        sum / self.sum_of_weights
    }

    fn compute_expected_loss_change<T: LossComputer>(
        &self,
        counts: &Counts,
        item: usize,
        label_estimate_from: LabelType,
        label_estimate_to: LabelType,
        loss: &T,
    ) -> f64 {
        let le_from = label_estimate_from as usize;
        let le_to = label_estimate_to as usize;
        let mut sum = 0.0;
        for (splinter, count) in self.x.iter().zip(&counts.x) {
            if let Some(alias) = splinter.mapping[item] {
                let sum_marginal = if le_from == le_to {
                    loss.compute_change_marginal(count.marginals[le_from], true)
                } else {
                    loss.compute_change_marginal(count.marginals[le_to], false)
                };
                let sum_joint = loss.compute_change_joint(splinter, count, alias, le_from, le_to);
                sum += (splinter.n_items as f64) * (sum_marginal + sum_joint)
            }
        }
        sum
    }

    fn max_n_items(&self) -> CountType {
        self.x
            .iter()
            .map(|splinter| splinter.n_items)
            .max()
            .unwrap()
    }
}

fn core<S: LossComputer + Sync, T: Rng>(
    n_items_in_all: usize,
    max_n_clusters_estimate: LabelType,
    s: &Splinters,
    loss: &S,
    n_runs: usize,
    n_cores: usize,
    rng: &mut T,
) -> (Clustering, f64, i32, &'static str) {
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_cores)
        .build()
        .unwrap();
    let mut rngs = Vec::with_capacity(n_runs);
    for _ in 0..n_runs {
        let mut seed = [0_u8; 16];
        rng.fill_bytes(&mut seed);
        let new_rng = Pcg64Mcg::from_seed(seed);
        rngs.push(new_rng);
    }
    pool.install(|| {
        let candidates = rngs.into_par_iter().map(|mut rng| {
            let mut estimate = Clustering::new(&sample_1tok(
                n_items_in_all,
                max_n_clusters_estimate,
                &mut rng,
            ));
            estimate.set_max_n_clusters(max_n_clusters_estimate);
            let mut counts = s.tally(&estimate, max_n_clusters_estimate);
            // Optimize
            let mut n_scans = 0;
            let mut visited_states = HashSet::new();
            let mut changed = true;
            while changed {
                n_scans += 1;
                changed = false;
                let permutation = {
                    let mut x: Vec<_> = (0..n_items_in_all).collect();
                    x.shuffle(&mut rng);
                    x
                };
                for item in permutation {
                    let label_estimate = estimate.get(item) as LabelType;
                    let label_best = {
                        let labels_and_deltas =
                            estimate.available_labels(item).map(|label_candidate| {
                                let label_candidate = label_candidate as LabelType;
                                let delta = s.compute_expected_loss_change(
                                    &counts,
                                    item,
                                    label_estimate,
                                    label_candidate,
                                    loss,
                                );
                                (label_candidate, delta)
                            });
                        labels_and_deltas
                            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                            .unwrap()
                            .0
                    };
                    if ALLOW_DEBUGGING && check_shortcut() {
                        let label_best_slow =
                            loss.best_allocation_slow(item, &estimate, s, &mut counts);
                        if label_best_slow[0].0 != label_best
                            && label_best_slow[0].1 + 1.0e-10 < label_best_slow[1].1
                        {
                            println!("best = {}, but {:?}", label_best, label_best_slow);
                            panic!();
                        }
                    }
                    if label_estimate != label_best {
                        s.update_tallies(&mut counts, item, label_estimate, label_best);
                        estimate.set(item, label_best);
                        changed = true;
                    }
                }
                if changed {
                    let relabeled = relabel(estimate.labels(), 0);
                    if visited_states.contains(&relabeled) {
                        changed = false;
                    } else {
                        visited_states.insert(relabeled);
                    }
                }
            }
            let expected_loss = s.compute_expected_loss(&counts, loss);
            (estimate, expected_loss, n_scans)
        });
        let mut candidates: Vec<_> = candidates.collect();
        candidates.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        let r = candidates.swap_remove(0);
        (r.0, r.1, r.2, loss.name())
    })
}

#[roxido]
fn splintered_clustering(
    splinters: &RList,
    n_items: usize,
    max_n_clusters: usize,
    n_runs: usize,
    use_vi: bool,
    a: f64,
    n_cores: usize,
) {
    let mut timer = timers::TicToc::new();
    timer.tic();
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let n_items_in_all = n_items;
    let splinters = Splinters::from_r(splinters);
    let max_n_clusters_estimate = max_n_clusters as LabelType;
    let max_n_clusters_estimate = if max_n_clusters_estimate == 0 {
        splinters
            .x
            .iter()
            .fold(0, |max, x| max.max(x.max_n_clusters_observed))
    } else {
        max_n_clusters_estimate
    };
    let (estimate, expected_loss, n_scans, loss_name) = if use_vi {
        let loss = VILossComputer::new(a, splinters.max_n_items()).unwrap();
        core(
            n_items_in_all,
            max_n_clusters_estimate,
            &splinters,
            &loss,
            n_runs,
            n_cores,
            &mut rng,
        )
    } else {
        let loss = BinderLossComputer::new(a).unwrap();
        core(
            n_items_in_all,
            max_n_clusters_estimate,
            &splinters,
            &loss,
            n_runs,
            n_cores,
            &mut rng,
        )
    };
    let result_rval = RList::with_names(&["estimate", "info"], pc);
    let estimate_rval = RVector::<i32>::new(n_items_in_all, pc);
    relabel_into_slice(estimate.labels(), 1, estimate_rval.slice_mut());
    result_rval.set(0, estimate_rval).stop();
    let names = [
        "loss",
        "a",
        "maxNClusters",
        "expectedLoss",
        "nScans",
        "nRuns",
        "nCores",
        "seconds",
    ];
    let info_rval = RList::with_names(&names, pc);
    info_rval.set(0, loss_name.to_r(pc)).stop();
    info_rval.set(1, a.to_r(pc)).stop();
    info_rval
        .set(2, i32::try_from(max_n_clusters).unwrap().to_r(pc))
        .stop();
    info_rval.set(3, expected_loss.to_r(pc)).stop();
    info_rval.set(4, n_scans.to_r(pc)).stop();
    info_rval
        .set(5, i32::try_from(n_runs).unwrap().to_r(pc))
        .stop();
    info_rval
        .set(6, i32::try_from(n_cores).unwrap().to_r(pc))
        .stop();
    timer.toc();
    info_rval.set(7, timer.as_secs_f64().to_r(pc)).stop();
    result_rval.set(1, info_rval).stop();
    result_rval
}

#[roxido]
fn expected_loss(estimate: &RVector, splinters: &RList, use_vi: bool, a: f64) {
    let estimate = Clustering::new(estimate.to_i32(pc).slice());
    let s = Splinters::from_r(splinters);
    let counts = s.tally(&estimate, estimate.max_n_clusters());
    if use_vi {
        s.compute_expected_loss(&counts, &VILossComputer::new(a, s.max_n_items()).unwrap())
    } else {
        s.compute_expected_loss(&counts, &BinderLossComputer::new(a).unwrap())
    }
}

fn validate_mass_and_sd(mass: f64, sd: f64) {
    if mass < 0.0 {
        stop!("'mass' must be greater than 0.0");
    }
    if sd < 0.0 {
        stop!("'sd' must be greater than 0.0");
    }
}

fn mk_crp(n_items: usize, mass: f64, discount: f64) -> CrpParameters {
    let discount = Discount::new(discount).stop_str("Discount must be in (0,1).");
    let concentration = Concentration::new_with_discount(mass, discount)
        .stop_str("Concentration must be greater than 0.");
    CrpParameters::new_with_discount(n_items, concentration, discount)
        .stop_str("Invalid mass and discount parameters.")
}

#[derive(Debug, Clone)]
pub struct CrpFusedParameters {
    n_items: Vec<usize>,
    concentration: f64,
    discount: f64,
}

impl CrpFusedParameters {
    fn new(n_items: Vec<usize>, concentration: f64, discount: f64) -> Result<Self, &'static str> {
        if n_items.contains(&0) {
            return Err("All of the sizes must be greater than 0.");
        }
        if concentration <= 0.0 {
            return Err("Concentration must be greater than 0.0.");
        }
        if !((0.0..1.0).contains(&discount)) {
            return Err("Discount must be in [0.0, 1.0).");
        }
        Ok(Self {
            n_items,
            concentration,
            discount,
        })
    }
}

impl FullConditional for CrpFusedParameters {
    fn log_full_conditional(&self, item: usize, clustering: &Partition) -> Vec<(usize, f64)> {
        let discount = self.discount;
        clustering
            .available_labels_for_reallocation(item)
            .map(|label| {
                let size = clustering.size_of_without(label, item);
                let value = if size == 0 {
                    (self.concentration + (clustering.n_clusters_without(item) as f64) * discount)
                        .ln()
                        + ln_gamma(self.n_items[item] as f64 - discount)
                        - ln_gamma(1.0 - discount)
                } else {
                    let size = clustering
                        .items_of_without(label, item)
                        .iter()
                        .map(|item| self.n_items[*item])
                        .sum::<usize>();
                    ln_gamma((size + self.n_items[item]) as f64 - discount)
                        - ln_gamma(size as f64 - discount)
                };
                (label, value)
            })
            .collect()
    }
}

#[roxido]
fn dpmm_generate(n_items: usize, mass: f64, sd: f64) {
    validate_mass_and_sd(mass, sd);
    let p = mk_crp(n_items, mass, 0.0);
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let partition: Partition = p.sample(&mut rng);
    let mut parameters: Vec<f64> = vec![0.0; partition.max_label() + 1];
    for &label in partition.active_labels().iter() {
        parameters[label] = rng.sample(StandardNormal);
    }
    let data: Vec<f64> = partition
        .allocation()
        .iter()
        .map(|&label| {
            let u: f64 = rng.sample(StandardNormal);
            parameters[label] + sd * u
        })
        .collect();
    let result = RList::with_names(&["partition", "parameters", "data"], pc);
    result
        .set(
            0,
            partition
                .allocation()
                .iter()
                .map(|&label| i32::try_from(label).stop() + 1)
                .to_r(pc),
        )
        .stop();
    result.set(1, parameters.to_r(pc)).stop();
    result.set(2, data.to_r(pc)).stop();
    result
}

#[roxido]
fn neal_data() {
    [-1.48, -1.40, -1.16, -1.08, -1.02, 0.14, 0.51, 0.53, 0.78].to_r(pc)
}

#[roxido]
fn dpmm_fit(data: &RVector, n_samples: usize, burnin: usize, thin: usize, mass: f64, sd: f64) {
    let data = data.to_f64(pc);
    let data = data.slice();
    let n_items = data.len();
    let n_iterations = (n_samples - burnin) / thin;
    validate_mass_and_sd(mass, sd);
    let p = mk_crp(n_items, mass, 0.0);
    let mut state = Partition::one_cluster(n_items);
    let permutation = Permutation::natural_and_fixed(n_items);
    let variance = sd * sd;
    let precision = variance.recip();
    #[allow(clippy::excessive_precision)]
    const LN_SQRT_2PI: f64 = 0.91893853320467274178032973640561763986139747363778;
    let log_posterior_predictive = |index: usize, subset: &[usize]| {
        let sum = subset.iter().fold(0.0, |s, &i| s + data[i]);
        let variance2 = (1.0 + (subset.len() as f64) * precision).recip();
        let mean = variance2 * sum * precision;
        let std_dev = (variance2 + variance).sqrt();
        let d = (data[index] - mean) / std_dev;
        -0.5 * d * d - LN_SQRT_2PI - std_dev.ln()
    };
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    // Burnin
    update_neal_algorithm3(
        burnin as u32,
        &mut state,
        &permutation,
        &p,
        &log_posterior_predictive,
        &mut rng,
    );
    // Collect
    let samples = RMatrix::<i32>::new(n_items, n_iterations, pc);
    let mut slice = samples.slice_mut();
    for _ in 0..n_iterations {
        update_neal_algorithm3(
            thin as u32,
            &mut state,
            &permutation,
            &p,
            &log_posterior_predictive,
            &mut rng,
        );
        state.relabel_into_slice(1, &mut slice[..n_items]);
        slice = &mut slice[n_items..];
    }
    let result = RList::with_names(&["partitions"], pc);
    let _ = result.set(0, samples);
    result
}

#[roxido]
fn dpmm_fit2(data: &RVector, n_samples: usize, burnin: usize, thin: usize, mass: f64, sd: f64) {
    let data = data.to_f64(pc);
    let data = data.slice();
    let n_items = data.len();
    let n_saved_samples = (n_samples - burnin) / thin;
    validate_mass_and_sd(mass, sd);
    let p = mk_crp(n_items, mass, 0.0);
    let mut state = Partition::one_cluster(n_items);
    let mut parameters: Vec<f64> = vec![0.0; state.max_label() + 1];
    let partitions_rval = RMatrix::<i32>::new(n_items, n_saved_samples, pc);
    let partitions_slice = partitions_rval.slice_mut();
    let parameters_rval = RList::new(n_saved_samples, pc);
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let mut seed: <Pcg64Mcg as SeedableRng>::Seed = Default::default();
    rng.fill(&mut seed);
    let mut rng2 = Pcg64Mcg::from_seed(seed);
    for &label in state.active_labels().iter() {
        parameters[label] = rng.sample(StandardNormal);
    }
    let permutation = Permutation::natural_and_fixed(n_items);
    #[allow(clippy::excessive_precision)]
    const LN_SQRT_2PI: f64 = 0.91893853320467274178032973640561763986139747363778;
    let log_normalizing_constant = LN_SQRT_2PI + sd.ln();
    let precision = (sd * sd).recip();
    let cacher = |_: usize| ();
    let mut iteration_counter = 0;
    let mut in_burnin = true;
    let mut how_many = burnin;
    loop {
        for _ in 0..how_many {
            let mut log_likelihood_contribution_fn =
                |index: usize, _cache: &(), label: usize, new_cluster: bool| {
                    let p = if !new_cluster {
                        parameters[label]
                    } else {
                        let parameter: f64 = rng2.sample(StandardNormal);
                        if label >= parameters.len() {
                            parameters.resize(label + 1, 0.0);
                        }
                        parameters[label] = parameter;
                        parameter
                    };
                    let d = (data[index] - p) / sd;
                    -0.5 * d * d - log_normalizing_constant
                };
            update_neal_algorithm8(
                1,
                &mut state,
                &permutation,
                &p,
                &mut log_likelihood_contribution_fn,
                cacher,
                &mut rng,
            );
            for (label, p) in parameters
                .iter_mut()
                .take(state.max_label() + 1)
                .enumerate()
            {
                let subset = state.items_of(label);
                if !subset.is_empty() {
                    let sum = subset.iter().fold(0.0, |s, &i| s + data[i]);
                    let variance = (1.0 + (subset.len() as f64) * precision).recip();
                    let mean = variance * sum * precision;
                    let u: f64 = rng.sample(StandardNormal);
                    *p = mean + variance.sqrt() * u;
                }
            }
        }
        if in_burnin {
            in_burnin = false;
            how_many = thin;
            continue;
        }
        let mut map = HashMap::new();
        let mut next_virgin_label = 1;
        let slice_out = &mut partitions_slice
            [(iteration_counter * state.n_items())..(iteration_counter + 1) * state.n_items()];
        let mut new_parameters = Vec::new();
        for (old_label, new_label) in state.allocation().iter().zip(slice_out.iter_mut()) {
            *new_label = *map.entry(*old_label).or_insert_with(|| {
                new_parameters.push(parameters[*old_label]);
                let virgin_label = next_virgin_label;
                next_virgin_label += 1;
                virgin_label
            });
        }
        let _ = parameters_rval.set_with_pc(iteration_counter, |pc| new_parameters.to_r(pc));
        iteration_counter += 1;
        if iteration_counter == n_saved_samples {
            break;
        }
    }
    let samples = RList::with_names(&["partitions", "parameters"], pc);
    let _ = samples.set(0, partitions_rval);
    let _ = samples.set(1, parameters_rval);
    samples
}

#[roxido]
fn dpmm_fit3(data: &RList, n_samples: usize, burnin: usize, thin: usize, mass: f64, sd: f64) {
    let mut data2 = Vec::new();
    for i in 0..data.len() {
        let Ok(x) = data.get(i) else {
            stop!("Cannot access element {} in list.", i);
        };
        let Ok(x) = x.as_vector() else {
            stop!("Element {} in list is not a vector.", i);
        };
        let x = x.to_f64(pc);
        let slice = x.slice();
        if slice.is_empty() {
            stop!("Element {} in list is a vector of length 0.", i);
        }
        data2.push(slice);
    }
    let n_items = data2.len();
    let n_saved_samples = (n_samples - burnin) / thin;
    validate_mass_and_sd(mass, sd);
    let n_items_per_item = data2.iter().map(|x| x.len()).collect::<Vec<_>>();
    let p = CrpFusedParameters::new(n_items_per_item, mass, 0.0).unwrap();
    let mut state = Partition::one_cluster(n_items);
    let mut parameters: Vec<f64> = vec![0.0; state.max_label() + 1];
    let partitions_rval = RMatrix::<i32>::new(n_items, n_saved_samples, pc);
    let partitions_slice = partitions_rval.slice_mut();
    let parameters_rval = RList::new(n_saved_samples, pc);
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let mut seed: <Pcg64Mcg as SeedableRng>::Seed = Default::default();
    rng.fill(&mut seed);
    let mut rng2 = Pcg64Mcg::from_seed(seed);
    for &label in state.active_labels().iter() {
        parameters[label] = rng.sample(StandardNormal);
    }
    let permutation = Permutation::natural_and_fixed(n_items);
    #[allow(clippy::excessive_precision)]
    const LN_SQRT_2PI: f64 = 0.91893853320467274178032973640561763986139747363778;
    let log_normalizing_constant = LN_SQRT_2PI + sd.ln();
    let precision = (sd * sd).recip();
    let cacher = |_: usize| ();
    let mut iteration_counter = 0;
    let mut in_burnin = true;
    let mut how_many = burnin;
    loop {
        for _ in 0..how_many {
            let mut log_likelihood_contribution_fn =
                |index: usize, _cache: &(), label: usize, new_cluster: bool| {
                    let p = if !new_cluster {
                        parameters[label]
                    } else {
                        let parameter: f64 = rng2.sample(StandardNormal);
                        if label >= parameters.len() {
                            parameters.resize(label + 1, 0.0);
                        }
                        parameters[label] = parameter;
                        parameter
                    };
                    let mut sum = 0.0;
                    for x in data2[index] {
                        let d = (x - p) / sd;
                        sum += -0.5 * d * d - log_normalizing_constant
                    }
                    sum
                };
            update_neal_algorithm8(
                1,
                &mut state,
                &permutation,
                &p,
                &mut log_likelihood_contribution_fn,
                cacher,
                &mut rng,
            );

            for (label, parameter) in parameters
                .iter_mut()
                .take(state.max_label() + 1)
                .enumerate()
            {
                let subset = state.items_of(label);
                if !subset.is_empty() {
                    let n_items = subset.iter().map(|i| p.n_items[*i]).sum::<usize>() as f64;
                    let sum = subset
                        .iter()
                        .fold(0.0, |s, &i| s + data2[i].iter().sum::<f64>());
                    let variance = (1.0 + n_items * precision).recip();
                    let mean = variance * sum * precision;
                    let u: f64 = rng.sample(StandardNormal);
                    *parameter = mean + variance.sqrt() * u;
                }
            }
        }
        if in_burnin {
            in_burnin = false;
            how_many = thin;
            continue;
        }
        let mut map = HashMap::new();
        let mut next_virgin_label = 1;
        let slice_out = &mut partitions_slice
            [(iteration_counter * state.n_items())..(iteration_counter + 1) * state.n_items()];
        let mut new_parameters = Vec::new();
        for (old_label, new_label) in state.allocation().iter().zip(slice_out.iter_mut()) {
            *new_label = *map.entry(*old_label).or_insert_with(|| {
                new_parameters.push(parameters[*old_label]);
                let virgin_label = next_virgin_label;
                next_virgin_label += 1;
                virgin_label
            });
        }
        let _ = parameters_rval.set(iteration_counter, new_parameters.to_r(pc));
        iteration_counter += 1;
        if iteration_counter == n_saved_samples {
            break;
        }
    }
    let samples = RList::with_names(&["partitions", "parameters"], pc);
    let _ = samples.set(0, partitions_rval);
    let _ = samples.set(1, parameters_rval);
    samples
}

fn sample_1tok<T: Rng>(n_items: usize, max_clusters: LabelType, rng: &mut T) -> Vec<usize> {
    let max_clusters = max_clusters as usize;
    let mut v = Vec::with_capacity(n_items);
    v.resize_with(n_items, || rng.random_range(0..max_clusters));
    v
}

trait LossComputer {
    fn name(&self) -> &'static str;
    fn engine(&self, x: CountType, n: CountType) -> f64;
    fn compute(&self, marginal1: f64, joints: &ArrayView2<CountType>, n: CountType) -> f64;
    fn compute_change_marginal(&self, n: CountType, same: bool) -> f64;
    fn compute_change_joint(
        &self,
        splinter: &Splinter,
        count: &Count,
        alias: usize,
        label_estimate_from: usize,
        label_estimate_to: usize,
    ) -> f64;

    fn compute_expectation(
        &self,
        marginals: &[CountType],
        joints: &Array3<CountType>,
        n: CountType,
    ) -> f64 {
        let marginal1 = self.compute_part_marginal1(marginals, n);
        joints
            .axis_iter(Axis(2))
            .map(|x| self.compute(marginal1, &x, n))
            .sum::<f64>()
            / (joints.len_of(Axis(2)) as f64)
    }

    fn compute_part_marginal1(&self, marginals: &[CountType], n: CountType) -> f64 {
        marginals.iter().map(|x| self.engine(*x, n)).sum::<f64>()
    }

    fn compute_parts(&self, joints: &ArrayView2<CountType>, n: CountType) -> (f64, f64) {
        let marginal0 = joints
            .axis_iter(Axis(1))
            .map(|x| self.engine(x.iter().sum::<CountType>(), n))
            .sum::<f64>();
        let joint = joints.iter().map(|x| self.engine(*x, n)).sum::<f64>();
        (marginal0, joint)
    }

    fn best_allocation_slow(
        &self,
        item: usize,
        estimate: &Clustering,
        splinters: &Splinters,
        counts: &mut Counts,
    ) -> Vec<(LabelType, f64)>
    where
        Self: Sized,
    {
        let label_estimate = estimate.labels()[item] as LabelType;
        let iter = estimate.available_labels(item).map(|label_candidate| {
            let label_candidate = label_candidate as LabelType;
            let delta = {
                let shard_losses = splinters.compute_expected_loss(counts, self);
                splinters.update_tallies(counts, item, label_estimate, label_candidate);
                let shard_loss_candidate = splinters.compute_expected_loss(counts, self);
                splinters.update_tallies(counts, item, label_candidate, label_estimate);
                shard_loss_candidate - shard_losses
            };
            (label_candidate, delta)
        });
        let mut list = iter.collect::<Vec<_>>();
        list.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        list
    }
}

struct BinderLossComputer {
    a: f64,
    b: f64,
    a_half: f64,
    b_half: f64,
}

impl BinderLossComputer {
    fn new(a: f64) -> Option<Self> {
        if (0.0..=2.0).contains(&a) {
            Some(Self {
                a,
                b: 2.0 - a,
                a_half: a / 2.0,
                b_half: (2.0 - a) / 2.0,
            })
        } else {
            None
        }
    }
}

impl LossComputer for BinderLossComputer {
    fn name(&self) -> &'static str {
        "binder"
    }

    fn engine(&self, x: CountType, n: CountType) -> f64 {
        let y = (x as f64) / (n as f64);
        y * y
    }

    fn compute(&self, marginal1: f64, joints: &ArrayView2<CountType>, n: CountType) -> f64 {
        let (marginal0, joint) = self.compute_parts(joints, n);
        self.a * marginal0 + self.b * marginal1 - 2.0 * joint
    }

    fn compute_change_marginal(&self, n: CountType, same: bool) -> f64 {
        if same {
            self.a_half + self.b_half * (n as f64)
        } else {
            self.b_half * (n as f64)
        }
    }

    fn compute_change_joint(
        &self,
        splinter: &Splinter,
        count: &Count,
        alias: usize,
        label_estimate_from: usize,
        label_estimate_to: usize,
    ) -> f64 {
        let mut sum_joint = 0;
        for sample_index in 0..splinter.clusterings.ncols() {
            let label_observed = splinter.clusterings[(alias, sample_index)] as usize;
            sum_joint += if label_estimate_from == label_estimate_to {
                count.joints[(label_estimate_from, label_observed, sample_index)]
                    as CountSummationType
            } else {
                count.joints[(label_estimate_to, label_observed, sample_index)]
                    as CountSummationType
            }
        }
        -(sum_joint as f64) / (count.joints.len_of(Axis(2)) as f64)
    }
}

struct VILossComputer {
    a: f64,
    b: f64,
    b_half: f64,
    log2: Vec<f64>,
    nlogn_diff: Vec<f64>,
}

impl VILossComputer {
    fn new(a: f64, max_n_items: CountType) -> Option<Self> {
        if (0.0..=2.0).contains(&a) {
            let b = 2.0 - a;
            let b_half = b / 2.0;
            let mut log2 = Vec::with_capacity((max_n_items + 2).into());
            for i in 0..log2.capacity() {
                log2.push((i as f64).log2());
            }
            let mut nlogn_diff = Vec::with_capacity((max_n_items + 1).into());
            nlogn_diff.push(0.0);
            let mut previous = 0.0;
            for (i, &v) in log2.iter().enumerate().skip(2) {
                let next = (i as f64) * v;
                nlogn_diff.push(next - previous);
                previous = next;
            }
            Some(Self {
                a,
                b,
                b_half,
                log2,
                nlogn_diff,
            })
        } else {
            None
        }
    }

    fn nlogn_diff(&self, x: CountType) -> f64 {
        // return if x == 0 {
        //     0.0
        // } else {
        //     let xx = x as f64;
        //     let xx_plus_1 = (x + 1) as f64;
        //     xx_plus_1 * xx_plus_1.log2() - xx * xx.log2()
        // };
        self.nlogn_diff[x as usize]
    }

    fn log2(&self, x: CountType) -> f64 {
        // return (x as f64).log2()
        self.log2[x as usize]
    }
}

impl LossComputer for VILossComputer {
    fn name(&self) -> &'static str {
        "VI"
    }

    fn engine(&self, x: CountType, n: CountType) -> f64 {
        if x == 0 {
            0.0
        } else {
            let y = (x as f64) / (n as f64);
            y * (self.log2(x) - self.log2(n))
        }
    }

    fn compute(&self, marginal1: f64, joints: &ArrayView2<CountType>, n: CountType) -> f64 {
        let (marginal0, joint) = self.compute_parts(joints, n);
        self.a * marginal0 + self.b * marginal1 - 2.0 * joint
    }

    fn compute_change_marginal(&self, n: CountType, same: bool) -> f64 {
        if same {
            self.b_half * self.nlogn_diff(n - 1)
        } else {
            self.b_half * self.nlogn_diff(n)
        }
    }

    fn compute_change_joint(
        &self,
        splinter: &Splinter,
        count: &Count,
        alias: usize,
        label_estimate_from: usize,
        label_estimate_to: usize,
    ) -> f64 {
        let mut sum_joint = 0.0;
        for sample_index in 0..splinter.clusterings.ncols() {
            let label_observed = splinter.clusterings[(alias, sample_index)] as usize;
            sum_joint += if label_estimate_from == label_estimate_to {
                self.nlogn_diff(
                    count.joints[(label_estimate_from, label_observed, sample_index)] - 1,
                )
            } else {
                self.nlogn_diff(count.joints[(label_estimate_to, label_observed, sample_index)])
            }
        }
        -sum_joint / (count.joints.len_of(Axis(2)) as f64)
    }
}

trait UnitIncrementor {
    fn next(x: &mut Self);
}

impl UnitIncrementor for LabelType {
    #[inline]
    fn next(x: &mut LabelType) {
        if *x == LabelType::MAX {
            panic!("Overflow!");
        }
        *x += 1
    }
}

impl UnitIncrementor for i32 {
    #[inline]
    fn next(x: &mut i32) {
        if *x == i32::MAX {
            panic!("Overflow!");
        }
        *x += 1
    }
}

fn relabel_into_slice<S: Copy + Eq + std::hash::Hash, T: Copy + UnitIncrementor>(
    slice_in: &[S],
    first_label: T,
    slice_out: &mut [T],
) {
    let mut map = HashMap::new();
    let mut next_new_label = first_label;
    for (old_label, slice_item) in slice_in.iter().zip(slice_out.iter_mut()) {
        *slice_item = *map.entry(*old_label).or_insert_with(|| {
            let new_label = next_new_label;
            T::next(&mut next_new_label);
            new_label
        });
    }
}

fn relabel<S: Copy + Eq + std::hash::Hash, T: Copy + UnitIncrementor>(
    slice_in: &[S],
    first_label: T,
) -> Vec<T> {
    let mut map = HashMap::new();
    let mut out = Vec::with_capacity(slice_in.len());
    let mut next_new_label = first_label;
    for old_label in slice_in.iter() {
        out.push(*map.entry(*old_label).or_insert_with(|| {
            let new_label = next_new_label;
            T::next(&mut next_new_label);
            new_label
        }));
    }
    out
}

fn canonize_labels<T: Copy + UnitIncrementor + std::fmt::Debug>(
    rval: &RMatrix,
    first_label: T,
    pc: &Pc,
) -> Vec<T> {
    let rval = rval.to_i32(pc);
    let n_items = rval.nrow();
    let mut result = vec![first_label; rval.len()];
    let slices_out = result.chunks_exact_mut(n_items);
    let slices_in = rval.slice().chunks_exact(n_items);
    for (si, so) in slices_in.zip(slices_out) {
        relabel_into_slice(si, first_label, so)
    }
    result
}
