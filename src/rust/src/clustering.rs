use crate::CountType;
use crate::LabelType;

use std::collections::HashMap;

pub struct Clustering {
    labels: Vec<LabelType>,
    sizes: Vec<CountType>,
    n_occupied: LabelType,
}

impl Clustering {
    pub fn new<T: Copy + Eq + std::hash::Hash>(labels: &[T]) -> Self {
        let n_items = labels.len();
        let mut new_labels = Vec::with_capacity(n_items);
        let mut sizes = Vec::new();
        let mut map = HashMap::new();
        let mut next_new_label = 0;
        for &item in labels.iter() {
            let c = *map.entry(item).or_insert_with(|| {
                let c = next_new_label;
                next_new_label += 1;
                c
            });
            new_labels.push(c);
            let cc = c as usize;
            if cc < sizes.len() {
                sizes[cc] += 1;
            } else {
                sizes.push(1);
            }
        }
        let n_occupied = sizes.len().try_into().unwrap();
        Self {
            labels: new_labels,
            sizes,
            n_occupied,
        }
    }

    pub fn set_max_n_clusters(&mut self, max: LabelType) {
        let len: LabelType = self.sizes.len().try_into().unwrap();
        use std::cmp::Ordering;
        match len.cmp(&max) {
            Ordering::Less => self.sizes.resize(max.into(), 0),
            Ordering::Equal => {}
            Ordering::Greater => {
                panic!(
                    "'labels' implies {} clusters, which is greater than max_n_clusters = {}.",
                    len, max
                )
            }
        }
    }

    #[allow(dead_code)]
    pub fn n_items(&self) -> usize {
        self.labels.len()
    }

    pub fn labels(&self) -> &Vec<LabelType> {
        &self.labels
    }

    #[allow(dead_code)]
    pub fn sizes(&self) -> &Vec<CountType> {
        &self.sizes
    }

    pub fn set(&mut self, index: usize, label: LabelType) {
        let old = &mut self.labels[index];
        self.sizes[*old as usize] -= 1;
        if self.sizes[*old as usize] == 0 {
            self.n_occupied -= 1;
        }
        if self.sizes[label as usize] == 0 {
            self.n_occupied += 1;
        }
        self.sizes[label as usize] += 1;
        *old = label;
    }

    pub fn get(&mut self, index: usize) -> LabelType {
        self.labels[index]
    }

    pub fn max_n_clusters(&self) -> LabelType {
        self.sizes.len() as LabelType
    }

    pub fn available_labels(&self, index: usize) -> LabelIterator<'_> {
        let found_empty = self.sizes[self.labels[index] as usize] == 1;
        LabelIterator {
            label: 0,
            iter: self.sizes.iter(),
            n_found: 0,
            n_occupied: self.n_occupied,
            found_empty,
        }
    }
}

impl std::ops::Index<usize> for Clustering {
    type Output = LabelType;
    fn index(&self, index: usize) -> &Self::Output {
        &self.labels[index]
    }
}

pub struct LabelIterator<'a> {
    label: LabelType,
    iter: std::slice::Iter<'a, CountType>,
    n_found: LabelType,
    n_occupied: LabelType,
    found_empty: bool,
}

impl Iterator for LabelIterator<'_> {
    type Item = LabelType;
    fn next(&mut self) -> Option<LabelType> {
        if self.n_found == self.n_occupied && self.found_empty {
            return None;
        }
        loop {
            let label = self.label;
            self.label += 1;
            match self.iter.next() {
                Some(&size) => {
                    if size == 0 {
                        if !self.found_empty {
                            self.found_empty = true;
                            return Some(label);
                        } else {
                            continue;
                        }
                    } else {
                        self.n_found += 1;
                        return Some(label);
                    }
                }
                None => {
                    return None;
                }
            }
        }
    }
}
