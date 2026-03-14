// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();
use roxido::*;

use dahl_salso::clustering::Clusterings;
use dahl_salso::optimize::{minimize_by_salso, SALSOParameters};
use dahl_salso::{LabelType, LossFunction, PartitionDistributionInformation};
use epa::epa::{sample, EpaParameters, SquareMatrixBorrower};
use epa::perm::Permutation;
use rand::prelude::SliceRandom;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Beta, Distribution};
use rand_pcg::Pcg64Mcg;
use roots::find_root_regula_falsi as find_root;
use std::convert::TryFrom;

fn sample_epa_engine<T: Rng>(
    n_samples: usize,
    n_items: usize,
    similarity: &[f64],
    mass: f64,
    n_cores: usize,
    rng: &mut T,
) -> (Vec<LabelType>, Vec<LabelType>) {
    let n_cores = if n_cores == 0 {
        std::thread::available_parallelism()
            .map(|x| x.get())
            .unwrap_or(1)
    } else {
        n_cores
    };
    let n_samples = n_samples.max(1);
    let n_samples_per_core = 1 + (n_samples - 1) / n_cores;
    let chunk_size = n_samples_per_core * n_items;
    let mut samples: Vec<LabelType> = vec![0; n_cores * chunk_size];
    let mut n_clusters: Vec<LabelType> = vec![0; n_cores * n_samples_per_core];

    let _result = crossbeam::scope(|s| {
        let mut stick1 = &mut samples[..];
        let mut stick2 = &mut n_clusters[..];
        let mut plan = Vec::with_capacity(n_cores);
        for _ in 0..n_cores - 1 {
            let (left1, right1) = stick1.split_at_mut(chunk_size);
            let (left2, right2) = stick2.split_at_mut(n_samples_per_core);
            plan.push((left1, left2, rng.random::<u128>()));
            stick1 = right1;
            stick2 = right2;
        }
        plan.push((stick1, stick2, rng.random()));
        let sim = SquareMatrixBorrower::from_slice(similarity, n_items);
        plan.into_iter().for_each(|p| {
            s.spawn(move |_| {
                let mut rng = Pcg64Mcg::new(p.2);
                let mut params =
                    EpaParameters::new(sim, Permutation::natural(n_items), mass).unwrap();
                for i in 0..n_samples_per_core {
                    params.shuffle_permutation(&mut rng);
                    let clustering = sample(&params, &mut rng);
                    let zero: LabelType = 0;
                    clustering.relabel_into_slice(zero, &mut p.0[i * n_items..(i + 1) * n_items]);
                    p.1[i] = LabelType::try_from(clustering.max_label() + 1).unwrap();
                }
            });
        });
    });
    (samples, n_clusters)
}

#[roxido]
fn sample_epa(n_samples: usize, similarity: &RMatrix<f64>, mass: f64, n_cores: usize) {
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let n_items = similarity.nrow();
    let (samples, _) = sample_epa_engine(
        n_samples,
        n_items,
        similarity.slice(),
        mass,
        n_cores,
        &mut rng,
    );
    let n_samples = samples.len() / n_items;
    let result = RMatrix::<i32>::new(n_samples, n_items, pc);
    let result_slice = result.slice_mut();
    for i in 0..n_items {
        for j in 0..n_samples {
            result_slice[i * n_samples + j] = i32::from(samples[j * n_items + i] + 1);
        }
    }
    result
}

#[roxido]
fn caviarpd_n_clusters(
    n_samples: usize,
    similarity: &RMatrix<f64>,
    mass: f64,
    use_vi: bool,
    n_runs: i32,
    max_size: i32,
    n_cores: usize,
) {
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let n_items = similarity.nrow();
    let (samples, n_clusters) = sample_epa_engine(
        n_samples,
        n_items,
        similarity.slice(),
        mass,
        n_cores,
        &mut rng,
    );
    let n_samples = samples.len() / n_items;
    let clusterings = Clusterings::unvalidated(n_samples, n_items, samples, n_clusters);
    let pdi = PartitionDistributionInformation::Draws(&clusterings);
    let a = 1.0;
    let loss_function = if use_vi {
        LossFunction::VI(a)
    } else {
        LossFunction::BinderDraws(a)
    };
    let p = SALSOParameters {
        n_items,
        max_size: LabelType::try_from(max_size).unwrap(),
        max_size_as_rf: false,
        max_scans: u32::MAX,
        max_zealous_updates: 10,
        n_runs: u32::try_from(n_runs).unwrap(),
        prob_sequential_allocation: 0.5,
        prob_singletons_initialization: 0.0,
    };
    let fit = minimize_by_salso(
        pdi,
        loss_function,
        &p,
        f64::INFINITY,
        u32::try_from(n_cores).unwrap(),
        &mut rng,
    );
    let result = fit.clustering.into_iter().max().unwrap() + 1;
    i32::try_from(result).unwrap().to_r(pc)
}

fn expected_number_of_clusters(mass: f64, n_items: usize) -> f64 {
    (0..n_items).fold(0.0, |sum, i| sum + mass / (mass + (i as f64)))
}

fn find_mass(enoc: f64, n_items: usize) -> f64 {
    let f = |mass| expected_number_of_clusters(mass, n_items) - enoc;
    match find_root(f64::EPSILON, enoc, f, &mut 1e-5_f64) {
        Ok(root) => root,
        Err(e) => {
            println!("Root finding error.... {}", e);
            1.0
        }
    }
}

#[roxido]
fn caviarpd_expected_number_of_clusters(mass: f64, n_items: usize) {
    expected_number_of_clusters(mass, n_items)
}

#[roxido]
fn caviarpd_mass(expected_number_of_clusters: f64, n_items: usize) {
    find_mass(expected_number_of_clusters, n_items)
}

// ---

#[roxido]
fn caviarpd_algorithm2(
    similarity: &RMatrix<f64>,
    min_n_clusters: f64,
    max_n_clusters: f64,
    mass: &RObject,
    n_samples: usize,
    grid_length: usize,
    n0: f64,
    tol: f64,
    use_vi: bool,
    salso_max_n_clusters: i32,
    salso_n_runs: i32,
    n_cores: usize,
) {
    let mut rng = Pcg64Mcg::from_seed(R::random_bytes::<16>());
    let n_items = similarity.nrow();
    let similarity_rval = similarity;
    let similarity = similarity_rval.slice();
    let (min_n_clusters, max_n_clusters) = {
        let x1 = min_n_clusters;
        let x2 = max_n_clusters;
        if x1 < x2 {
            (x1, x2)
        } else {
            (x2, x1)
        }
    };
    let grid_length = grid_length.max(if min_n_clusters == max_n_clusters {
        1
    } else {
        2
    });
    let salso_n_runs = salso_n_runs.max(1);
    let samples_rval = RMatrix::<i32>::new(n_samples * grid_length, n_items, pc);
    let samples_slice = samples_rval.slice_mut();
    let p = SALSOParameters {
        n_items,
        max_size: LabelType::try_from(salso_max_n_clusters).unwrap(),
        max_size_as_rf: false,
        max_scans: u32::MAX,
        max_zealous_updates: 10,
        n_runs: u32::try_from(salso_n_runs).unwrap(),
        prob_sequential_allocation: 0.5,
        prob_singletons_initialization: 0.0,
    };
    let mut previous = 1.0;
    let mut candidates_labels = Vec::with_capacity(grid_length * n_items);
    let mut candidates_n_clusters = Vec::with_capacity(grid_length);
    let masses = {
        let mut masses = if mass.is_null() {
            let step_size = (max_n_clusters - min_n_clusters) / (grid_length as f64);
            (0..grid_length)
                .map(|x| find_mass(min_n_clusters + (x as f64) * step_size, n_items))
                .collect::<Vec<_>>()
        } else {
            let mass_rval = mass.as_vector().stop().to_f64(pc);
            let mass = mass_rval.slice();
            if mass.len() == 1 {
                vec![mass[0]; grid_length]
            } else {
                mass.to_vec()
            }
        };
        masses.shuffle(&mut rng);
        masses
    };
    for (i, mass) in masses.into_iter().enumerate() {
        let (samples, n_clusters) =
            sample_epa_engine(n_samples, n_items, similarity, mass, n_cores, &mut rng);
        let clusterings =
            Clusterings::unvalidated(samples.len() / n_items, n_items, samples, n_clusters);
        for jj in 0..n_samples {
            let labels = clusterings.labels(jj);
            for (ii, value) in labels.iter().enumerate() {
                samples_slice[n_samples * (ii * grid_length + i) + jj] = i32::from(*value + 1);
            }
        }
        let pdi = PartitionDistributionInformation::Draws(&clusterings);
        let (mut lower, mut upper) = (0.0, 2.0);
        let beta = Beta::new(n0 * previous / 2.0, n0 * (1.0 - previous / 2.0)).unwrap();
        let mut a = 2.0 * beta.sample(&mut rng);
        let candidate;
        loop {
            let loss_function = if use_vi {
                LossFunction::VI(a)
            } else {
                LossFunction::BinderDraws(a)
            };
            let fit = minimize_by_salso(
                pdi,
                loss_function,
                &p,
                f64::INFINITY,
                u32::try_from(n_cores).unwrap(),
                &mut rng,
            );
            let n_clusters = fit.clustering.iter().max().unwrap() + 1;
            if upper - lower <= tol {
                candidate = fit.clustering;
                break;
            } else if (n_clusters as f64) < min_n_clusters {
                upper = a;
                a = (lower + a) / 2.0;
            } else if (n_clusters as f64) > max_n_clusters {
                lower = a;
                a = (upper + a) / 2.0;
            } else {
                candidate = fit.clustering;
                break;
            }
        }
        previous = a;
        candidates_labels.extend(candidate.iter().map(|x| LabelType::try_from(*x).unwrap()));
        candidates_n_clusters
            .push(LabelType::try_from(candidate.iter().max().unwrap() + 1).unwrap());
    }
    let candidates = Clusterings::unvalidated(
        grid_length,
        n_items,
        candidates_labels,
        candidates_n_clusters,
    );
    let pdi = PartitionDistributionInformation::Draws(&candidates);
    let loss_function = if use_vi {
        LossFunction::VI(1.0)
    } else {
        LossFunction::BinderDraws(1.0)
    };
    let fit = minimize_by_salso(
        pdi,
        loss_function,
        &p,
        f64::INFINITY,
        u32::try_from(n_cores).unwrap(),
        &mut rng,
    );
    let estimate_rval = RVector::<i32>::new(n_items, pc);
    for (src, dst) in fit.clustering.iter().zip(estimate_rval.slice_mut()) {
        *dst = i32::try_from(*src + 1).unwrap();
    }
    let result = RList::with_names(&["estimate", "samples"], pc);
    result.set(0, estimate_rval).stop();
    result.set(1, samples_rval).stop();
    result
}
