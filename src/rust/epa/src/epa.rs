// Ewens Pitman attraction partition distribution

use crate::clust::Clustering;
use crate::perm::Permutation;
// use roots::find_root_regula_falsi as find_root;

use rand::prelude::*;
use std::slice;

type SimilarityBorrower<'a> = SquareMatrixBorrower<'a>;

#[derive(Debug, Clone)]
pub struct EpaParameters<'a> {
    similarity: SimilarityBorrower<'a>,
    permutation: Permutation,
    mass: f64,
}

impl<'a> EpaParameters<'a> {
    pub fn new(
        similarity: SimilarityBorrower<'a>,
        permutation: Permutation,
        mass: f64,
    ) -> Option<Self> {
        if similarity.n_items() != permutation.n_items() {
            None
        } else {
            Some(Self {
                similarity,
                permutation,
                mass,
            })
        }
    }

    pub fn shuffle_permutation<T: Rng>(&mut self, rng: &mut T) {
        self.permutation.shuffle(rng);
        /*
        match std::env::var("DBD_METHOD").as_deref() {
            Ok("jumps" | "biased") => {
                self.permutation = {
                    let mut permutation = Vec::with_capacity(self.permutation.n_items());
                    let mut available: Vec<_> = (0..self.permutation.n_items()).collect();
                    let start = rng.gen_range(0..available.len());
                    let mut current_index = available.swap_remove(start);
                    permutation.push(current_index);
                    while !available.is_empty() {
                        let index_and_weights = available
                            .iter()
                            .map(|i| self.similarity[(current_index, *i)])
                            .enumerate();
                        let (index, _) =
                            Clustering::select(index_and_weights, false, 0, Some(rng), false);
                        current_index = available.swap_remove(index);
                        permutation.push(current_index);
                    }
                    Permutation::from_vector(permutation).unwrap()
                }
            }
            _ => self.permutation.shuffle(rng),
        }
        */
    }
}

/// A data structure representing a square matrix.
///
#[derive(Debug)]
pub struct SquareMatrix {
    data: Vec<f64>,
    n_items: usize,
}

impl SquareMatrix {
    pub fn zeros(n_items: usize) -> Self {
        Self {
            data: vec![0.0; n_items * n_items],
            n_items,
        }
    }

    pub fn ones(n_items: usize) -> Self {
        Self {
            data: vec![1.0; n_items * n_items],
            n_items,
        }
    }

    pub fn identity(n_items: usize) -> Self {
        let ni1 = n_items + 1;
        let n2 = n_items * n_items;
        let mut data = vec![0.0; n2];
        let mut i = 0;
        while i < n2 {
            data[i] = 1.0;
            i += ni1
        }
        Self { data, n_items }
    }

    pub fn data(&self) -> &[f64] {
        &self.data[..]
    }

    pub fn data_mut(&mut self) -> &mut [f64] {
        &mut self.data[..]
    }

    pub fn view(&mut self) -> SquareMatrixBorrower<'_> {
        SquareMatrixBorrower::from_slice(&self.data[..], self.n_items)
    }

    pub fn n_items(&self) -> usize {
        self.n_items
    }
}

#[derive(Debug, Copy, Clone)]
pub struct SquareMatrixBorrower<'a> {
    data: &'a [f64],
    n_items: usize,
}

impl std::ops::Index<(usize, usize)> for SquareMatrixBorrower<'_> {
    type Output = f64;
    fn index(&self, (i, j): (usize, usize)) -> &Self::Output {
        &self.data[self.n_items * j + i]
    }
}

impl<'a> SquareMatrixBorrower<'a> {
    pub fn from_slice(data: &'a [f64], n_items: usize) -> Self {
        assert_eq!(data.len(), n_items * n_items);
        Self { data, n_items }
    }

    /// # Safety
    ///
    /// You're on your own.
    pub unsafe fn from_ptr(data: *const f64, n_items: usize) -> Self {
        let data = slice::from_raw_parts(data, n_items * n_items);
        Self { data, n_items }
    }

    pub fn n_items(&self) -> usize {
        self.n_items
    }

    /// # Safety
    ///
    /// You're on your own.
    pub unsafe fn get_unchecked(&self, (i, j): (usize, usize)) -> &f64 {
        self.data.get_unchecked(self.n_items * j + i)
    }

    pub fn data(&self) -> &[f64] {
        self.data
    }

    pub fn sum_of_triangle(&self) -> f64 {
        let mut sum = 0.0;
        for i in 0..self.n_items {
            for j in 0..i {
                sum += unsafe { *self.get_unchecked((i, j)) };
            }
        }
        sum
    }

    pub fn sum_of_row_subset(&self, row: usize, columns: &[usize]) -> f64 {
        let mut sum = 0.0;
        for j in columns {
            sum += unsafe { *self.get_unchecked((row, *j)) };
        }
        sum
    }
}

pub fn sample<T: Rng>(parameters: &EpaParameters, rng: &mut T) -> Clustering {
    let ni = parameters.similarity.n_items();
    let (mass, path): (f64, Option<Vec<f64>>) = (parameters.mass, None);
    /*
    let (mass, path) = match std::env::var("DBD_METHOD").as_deref() {
        Ok("jumps") => {
            let mut path: Vec<_> = std::iter::once(1.0)
                .chain((1..ni).map(|i| {
                    parameters.similarity[(
                        parameters.permutation.get(i - 1),
                        parameters.permutation.get(i),
                    )]
                }))
                .collect();
            let avg = path.iter().sum::<f64>() / (ni as f64);
            for p in &mut path {
                *p = (avg / *p).min(100.0);
            }
            let mass = parameters.mass;
            let expected_number_of_clusters =
                (0..ni).fold(0.0, |sum, i| sum + mass / (mass + (i as f64)));
            let f = |m| {
                (0..ni).fold(0.0, |sum, i| {
                    let p = m * path[i];
                    sum + p / (p + (i as f64))
                }) - expected_number_of_clusters
            };
            let root = find_root(f64::EPSILON, 10.0 * mass, &f, &mut 1e-5_f64);
            (
                root.unwrap_or_else(|e| {
                    println!("Root finding error.... {e}");
                    mass
                }),
                Some(path),
            )
        }
        _ => (parameters.mass, None),
    };
    */
    let mut clustering = Clustering::unallocated(ni);
    for i in 0..ni {
        let ii = parameters.permutation.get(i);
        let jump_density = match path {
            Some(ref path) => path[i],
            None => 1.0,
        };
        let kt = (i as f64)
            / parameters
                .similarity
                .sum_of_row_subset(ii, parameters.permutation.slice_until(i));
        let labels_and_weights = clustering
            .available_labels_for_allocation_with_target(None, ii)
            .map(|label| {
                let n_items_in_cluster = clustering.size_of(label);
                let weight = if n_items_in_cluster == 0 {
                    mass * jump_density
                } else {
                    kt * parameters
                        .similarity
                        .sum_of_row_subset(ii, &clustering.items_of(label)[..])
                };
                (label, weight)
            });
        let subset_index = Clustering::select(labels_and_weights, false, 0, Some(rng), false).0;
        clustering.allocate(ii, subset_index);
    }
    clustering
}
