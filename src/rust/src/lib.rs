// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();

use fastrand::Rng;
use rayon::ThreadPool;
use rayon::prelude::*;
use roxido::*;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use walltime::TicToc;

// Integer sqrt root using a binary search algorithm
fn isqrt(x: usize) -> usize {
    let mut lower = 0;
    let mut upper = x + 1;
    while lower != upper - 1 {
        let m = (lower + upper) / 2;
        if m * m <= x {
            lower = m;
        } else {
            upper = m;
        }
    }
    lower
}

#[allow(dead_code)]
trait Dag {
    // Returns the number of nodes
    fn n_nodes(&self) -> usize;

    // Returns true if the graph is a DAG (i.e., no cycles)
    fn is_dag(&self) -> bool {
        !self.contains_cycles()
    }

    // Returns true if the graph has cycles (i.e., not a DAG)
    fn contains_cycles(&self) -> bool;

    // Adds a directed edge from parent to child
    fn add_edge(&mut self, parent: usize, child: usize) -> bool;

    // Deletes a directed edge parent to child
    fn delete_edge(&mut self, parent: usize, child: usize) -> bool;

    // Returns true if there is a directed edge from parent to child
    fn has_edge(&self, parent: usize, child: usize) -> bool;

    fn index(&self, parent: usize, child: usize) -> usize {
        self.n_nodes() * child + parent
    }

    fn unindex(&self, index: usize) -> (usize, usize) {
        (index % self.n_nodes(), index / self.n_nodes())
    }
}

#[derive(Debug, Clone)]
pub struct DagAdjacencyVector {
    adjacency_list: Vec<Vec<usize>>,
}

impl Dag for DagAdjacencyVector {
    fn n_nodes(&self) -> usize {
        self.adjacency_list.len()
    }

    fn contains_cycles(&self) -> bool {
        let mut visited = vec![false; self.adjacency_list.len()];
        let mut stack = vec![false; self.adjacency_list.len()];
        for node in 0..self.adjacency_list.len() {
            if !visited[node] && self.dfs(node, &mut visited, &mut stack) {
                return true;
            }
        }
        false
    }

    fn add_edge(&mut self, parent: usize, child: usize) -> bool {
        let mut visited = vec![false; self.adjacency_list.len()];
        let mut stack = vec![false; self.adjacency_list.len()];
        visited[parent] = true;
        stack[parent] = true;
        if self.dfs(child, &mut visited, &mut stack) {
            return false;
        }
        let list = &mut self.adjacency_list[parent];
        if !list.contains(&child) {
            list.push(child);
        }
        true
    }

    fn delete_edge(&mut self, parent: usize, child: usize) -> bool {
        let list = &mut self.adjacency_list[parent];
        let mut result = false;
        if let Some(index) = list.iter().position(|x| *x == child) {
            list.swap_remove(index);
            result = true;
        }
        result
    }

    fn has_edge(&self, parent: usize, child: usize) -> bool {
        let list = &self.adjacency_list[parent];
        list.contains(&child)
    }
}

impl DagAdjacencyVector {
    // Creates a new graph
    fn new(n_nodes: usize) -> Self {
        Self {
            adjacency_list: vec![Vec::new(); n_nodes],
        }
    }

    pub fn from_slice(x: &[i32]) -> Result<Self, &'static str> {
        Self::from_slice_with_n(x, isqrt(x.len()))
    }

    pub fn from_slice_with_n<T: HasEqualsZero>(x: &[T], n: usize) -> Result<Self, &'static str> {
        if n * n != x.len() {
            return Err("Expected a square matrix");
        }
        let mut result = Self::new(n);
        for index in x
            .iter()
            .enumerate()
            .filter(|y| !(*y.1).equals_zero())
            .map(|z| z.0)
        {
            let (parent, child) = result.unindex(index);
            if !result.add_edge(parent, child) {
                return Err("Expected a DAG");
            }
        }
        Ok(result)
    }

    pub fn from_r(x: &RMatrix, pc: &Pc) -> Result<Self, &'static str> {
        let x = x.to_i32(pc);
        let [n, n2] = x.dim();
        if n != n2 {
            return Err("Expected a square matrix");
        }
        Self::from_slice_with_n(x.slice(), n)
    }

    // Depth-first search with cycle detection
    fn dfs(&self, node: usize, visited: &mut Vec<bool>, stack: &mut Vec<bool>) -> bool {
        if stack[node] {
            return true; // Cycle found
        }
        if visited[node] {
            return false; // Already visited
        }
        visited[node] = true;
        stack[node] = true;
        let neighbors = &self.adjacency_list[node];
        for &neighbor in neighbors {
            if self.dfs(neighbor, visited, stack) {
                return true;
            }
        }
        stack[node] = false;
        false
    }

    fn vec(&self) -> Vec<u8> {
        let n_nodes = self.n_nodes();
        let mut result = vec![0_u8; n_nodes * n_nodes];
        for (parent, children) in self.adjacency_list.iter().enumerate() {
            for &child in children.iter() {
                let index = self.index(parent, child);
                result[index] = 1;
            }
        }
        result
    }
}

struct Gsh {
    proportions: Vec<f64>,
    a: f64,
}

impl Gsh {
    pub fn from_r(array: &RArray, a: f64, pc: &Pc) -> Self {
        if !(0.0..=2.0).contains(&a) {
            stop!("'a' must be between 0.0 and 2.0");
        }
        let dim = array.dim();
        if !(2..=3).contains(&dim.len()) {
            stop!("'x' must be two or three dimensional");
        }
        let n_items = dim[0];
        if n_items != dim[1] {
            stop!("The first two elements of 'x' should give a square matrix")
        }
        let len = n_items.pow(2);
        let mut proportions;
        match dim.get(2) {
            Some(n_samples) => {
                proportions = Vec::with_capacity(len);
                let array2 = array.to_i32(pc);
                let slice = array2.slice();
                let n_samples_f64 = *n_samples as f64;
                for j in 0..proportions.capacity() {
                    let mut count = 0;
                    for i in 0..*n_samples {
                        count += slice[i * len + j]
                    }
                    proportions.push(f64::from(count) / n_samples_f64)
                }
            }
            None => {
                proportions = vec![0.0; len];
                proportions.copy_from_slice(array.to_f64(pc).slice())
            }
        }
        Self { proportions, a }
    }

    pub fn expected_loss<T: HasEqualsZero>(&self, x: &[T]) -> Result<f64, &'static str> {
        if x.len() != self.proportions.len() {
            return Err("Unexpected length");
        }
        let mut ac = 0.0;
        let mut bc = 0.0;
        for (x, y) in x.iter().zip(self.proportions.iter()) {
            if x.equals_zero() {
                bc += *y
            } else {
                ac += 1.0 - *y
            }
        }
        Ok(self.a * ac + (2.0 - self.a) * bc)
    }

    pub fn recommendation(&self, index: usize) -> bool {
        2.0 * self.proportions[index] > self.a
    }
}

pub trait HasEqualsZero {
    fn equals_zero(&self) -> bool;
}

impl HasEqualsZero for i32 {
    fn equals_zero(&self) -> bool {
        *self == 0
    }
}

impl HasEqualsZero for u8 {
    fn equals_zero(&self) -> bool {
        *self == 0
    }
}

struct GshBuilder {
    n_items: usize,
    counts: Vec<usize>,
    a: f64,
    n_samples_total: usize,
    n_samples_processed: usize,
    index_of_candidates: Vec<usize>,
    candidates: Vec<DagAdjacencyVector>,
    gsh_option: Option<Gsh>,
}

impl GshBuilder {
    fn new(n_items: usize, n_samples: usize, n_starts: usize, a: f64) -> Self {
        let counts = vec![0; n_items * n_items];
        let seed = u64::from_le_bytes(R::random_bytes::<8>());
        let mut rng = Rng::with_seed(seed);
        let n_candidates = if n_starts > 0 {
            n_starts.min(n_samples)
        } else {
            n_samples
        };
        let mut index_of_candidates = rng.choose_multiple(0..n_samples, n_candidates);
        index_of_candidates.sort_unstable_by(|a, b| b.cmp(a));
        Self {
            n_items,
            counts,
            a,
            n_samples_total: n_samples,
            n_samples_processed: 0,
            index_of_candidates,
            candidates: Vec::with_capacity(n_candidates),
            gsh_option: None,
        }
    }

    fn process<T: HasEqualsZero>(&mut self, x: &[T]) -> Result<(), &'static str> {
        if x.len() != self.counts.len() {
            return Err("Unexpected length.");
        }
        for (y, z) in self.counts.iter_mut().zip(x) {
            if !z.equals_zero() {
                *y += 1;
            }
        }
        if let Some(&index) = self.index_of_candidates.last()
            && self.n_samples_processed == index
        {
            let Ok(dag) = DagAdjacencyVector::from_slice_with_n(x, self.n_items) else {
                return Err("Network is not a DAG.");
            };
            self.index_of_candidates.pop();
            self.candidates.push(dag);
        }
        self.n_samples_processed += 1;
        Ok(())
    }

    fn sample(
        &mut self,
        x: &[i32],
        n_samples: usize,
        prob: f64,
        rng: &mut Rng,
    ) -> Result<(), &'static str> {
        if n_samples == 0 {
            return Ok(());
        }
        if x.len() != self.counts.len() {
            return Err("Unexpected length.");
        }
        let Ok(mut dag) = DagAdjacencyVector::from_slice_with_n(x, self.n_items) else {
            return Err("Network is not a DAG.");
        };
        let mut counter = 0;
        loop {
            for (parent, children) in dag.adjacency_list.iter().enumerate() {
                for &child in children.iter() {
                    let index = dag.index(parent, child);
                    self.counts[index] += 1;
                }
            }
            if let Some(&index) = self.index_of_candidates.last()
                && self.n_samples_processed == index
            {
                self.index_of_candidates.pop();
                self.candidates.push(dag.clone());
            }
            counter += 1;
            self.n_samples_processed += 1;
            if counter >= n_samples {
                break;
            }
            engine_random(&mut dag, prob, rng);
        }
        Ok(())
    }

    fn finalize(&mut self) -> Result<(&Gsh, Vec<DagAdjacencyVector>), &'static str> {
        if self.gsh_option.is_none() {
            if self.n_samples_processed != self.n_samples_total {
                return Err(
                    "Number of processed networks does not match the originally declared number.",
                );
            }
            let mut proportions = Vec::with_capacity(self.counts.len());
            let n_samples = self.n_samples_total as f64;
            self.counts
                .iter()
                .for_each(|x| proportions.push(*x as f64 / n_samples));
            let gsh = Gsh {
                proportions,
                a: self.a,
            };
            self.gsh_option = Some(gsh);
        }
        Ok((self.gsh_option.as_ref().unwrap(), self.candidates.clone()))
    }
}

fn engine_optimize<A: Dag>(dag: &mut A, gsh: &Gsh, rng: &mut Rng) -> u32 {
    let mut leftovers = Vec::new();
    let mut permutation = (0..dag.n_nodes().pow(2)).collect::<Vec<_>>();
    rng.shuffle(&mut permutation);
    for index in permutation {
        let (parent, child) = dag.unindex(index);
        if gsh.recommendation(index) {
            if !dag.add_edge(parent, child) {
                leftovers.push((parent, child));
            }
        } else {
            dag.delete_edge(parent, child);
        }
    }
    let mut not_optimal_counter = 0;
    rng.shuffle(&mut leftovers);
    for (parent, child) in leftovers.into_iter() {
        if !dag.add_edge(parent, child) {
            not_optimal_counter += 1
        };
    }
    not_optimal_counter
}

fn engine_random<A: Dag>(dag: &mut A, prob: f64, rng: &mut Rng) -> u32 {
    let len = dag.n_nodes().pow(2);
    let u1 = rng.f64();
    let u2 = rng.f64();
    let which = rng.choose_multiple(0..len, sample_binomial(len, prob, u1, u2));
    let mut not_optimal_counter = 0;
    for index in which {
        let (parent, child) = dag.unindex(index);
        if !dag.delete_edge(parent, child) && !dag.add_edge(parent, child) {
            not_optimal_counter += 1
        }
    }
    not_optimal_counter
}

fn sample_binomial(n: usize, p: f64, u1: f64, u2: f64) -> usize {
    debug_assert!((0.0..=1.0).contains(&p), "p must be in [0, 1]");
    debug_assert!((0.0..=1.0).contains(&u1), "u1 must be in [0, 1]");
    debug_assert!((0.0..=1.0).contains(&u2), "u2 must be in [0, 1]");
    if n == 0 || p == 0.0 {
        return 0;
    }
    if p == 1.0 {
        return n;
    }
    let mean = n as f64 * p;
    let std_dev = (n as f64 * p * (1.0 - p)).sqrt();
    // Use the normal approximation for large n with Box-Muller transformation
    let z0 = (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos();
    let normal = mean + std_dev * z0;
    normal.round().clamp(0.0, n as f64) as usize
}

#[roxido]
fn get_external_ptr_tag(x: &RExternalPtr) {
    unsafe { crate::rbindings::R_ExternalPtrTag(x.sexp()) }
}

#[roxido]
fn is_dag(dag: &RMatrix) {
    DagAdjacencyVector::from_r(dag, pc).is_ok()
}

#[roxido]
fn initialize_expected_gsh_loss(samples: &RArray, a: f64) {
    let gsh = Gsh::from_r(samples, a, pc);
    RExternalPtr::encode(gsh, "gsh", pc)
}

#[roxido]
fn gsh_loss_builder_new(n_items: usize, n_samples: usize, n_starts: usize, a: f64) {
    let gsh_builder = GshBuilder::new(n_items, n_samples, n_starts, a);
    RExternalPtr::encode_full(gsh_builder, "gsh_builder".to_r(pc), false, pc)
}

#[roxido]
fn gsh_loss_builder_process(network: &RVector, gsh_builder: &mut RExternalPtr) {
    let gsh_builder = RExternalPtr::decode_mut::<GshBuilder>(gsh_builder);
    gsh_builder.process(network.to_i32(pc).slice()).stop();
}

#[roxido]
fn gsh_loss_builder_sample(
    network: &RVector,
    gsh_builder: &mut RExternalPtr,
    n_samples: usize,
    prob: f64,
) {
    let gsh_builder = RExternalPtr::decode_mut::<GshBuilder>(gsh_builder);
    let seed = u64::from_le_bytes(R::random_bytes::<8>());
    let mut rng = Rng::with_seed(seed);
    gsh_builder
        .sample(network.to_i32(pc).slice(), n_samples, prob, &mut rng)
        .stop();
}

#[roxido]
fn compute_expected_gsh_loss(network: &RMatrix, gsh: &RExternalPtr) {
    let summary = RExternalPtr::decode_ref::<Gsh>(gsh);
    summary.expected_loss(network.to_i32(pc).slice()).stop()
}

struct Timers {
    overall: TicToc,
    transform: TicToc,
    proportions: TicToc,
    search: TicToc,
}

impl Timers {
    pub fn new() -> Self {
        Self {
            overall: TicToc::new(),
            transform: TicToc::new(),
            proportions: TicToc::new(),
            search: TicToc::new(),
        }
    }
}

fn bantha_core<'a>(
    gsh: &Gsh,
    candidates: Vec<DagAdjacencyVector>,
    mut timers: Timers,
    mut rng: Rng,
    pool: ThreadPool,
    n_cores: usize,
    pc: &'a Pc,
) -> &'a RMatrix<i32> {
    timers.overall.tic();
    let found = Arc::new(AtomicBool::new(false));
    let n_candidates = candidates.len();
    // Compute
    timers.search.tic();
    let (dag, expected_loss, not_optimal_counter) = if n_cores == 1 {
        let mapped = candidates.into_iter().map(|mut dag| {
            let not_optimal_counter = engine_optimize(&mut dag, gsh, &mut rng);
            let expected_loss = gsh.expected_loss(&dag.vec()).unwrap();
            (dag, expected_loss, not_optimal_counter)
        });
        mapped
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap()
    } else {
        pool.install(|| {
            let mapped = candidates
                .into_par_iter()
                .map_init(Rng::new, |rng, mut dag| {
                    if found.load(Ordering::Relaxed) {
                        return (dag, f64::INFINITY, u32::MAX);
                    }
                    let not_optimal_counter = engine_optimize(&mut dag, gsh, rng);
                    let expected_loss = gsh.expected_loss(&dag.vec()).unwrap();
                    if not_optimal_counter == 0 {
                        found.store(true, Ordering::Relaxed);
                    }
                    (dag, expected_loss, not_optimal_counter)
                });
            mapped
                .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
                .unwrap()
        })
    };
    let result = RMatrix::<i32>::from_value(0, dag.n_nodes(), dag.n_nodes(), pc);
    let slice = result.slice_mut();
    for (parent, children) in dag.adjacency_list.iter().enumerate() {
        for &child in children.iter() {
            let index = dag.index(parent, child);
            slice[index] = 1;
        }
    }
    timers.search.toc();
    timers.overall.toc();
    let settings = RList::with_names(&["a", "n_candidates", "n_cores"], pc);
    settings.set(0, gsh.a.to_r(pc)).unwrap();
    settings
        .set(1, i32::try_from(n_candidates).unwrap().to_r(pc))
        .unwrap();
    settings
        .set(2, i32::try_from(n_cores).unwrap().to_r(pc))
        .unwrap();
    result.set_attribute(RSymbol::from("settings").unwrap(), settings);
    result.set_attribute(
        RSymbol::from("n_not_optimal").unwrap(),
        i32::try_from(not_optimal_counter).unwrap().to_r(pc),
    );
    result.set_attribute(
        RSymbol::from("expected_loss").unwrap(),
        expected_loss.to_r(pc),
    );
    let elapsed_seconds_vector = RVector::<f64>::new(4, pc);
    elapsed_seconds_vector
        .set_names(["overall", "transform", "proportions", "search"].to_r(pc))
        .unwrap();
    let slice = elapsed_seconds_vector.slice_mut();
    slice[0] = timers.overall.as_secs_f64();
    slice[1] = timers.transform.as_secs_f64();
    slice[2] = timers.proportions.as_secs_f64();
    slice[3] = timers.search.as_secs_f64();
    result.set_attribute(
        RSymbol::from("elapsed_seconds").unwrap(),
        elapsed_seconds_vector,
    );
    result
}

#[roxido]
fn bantha_big_data(gsh_builder: &mut RExternalPtr, n_cores: usize) {
    let gsh_builder = RExternalPtr::decode_mut::<GshBuilder>(gsh_builder);
    let (gsh, candidates) = gsh_builder.finalize().stop();
    let timers = Timers::new();
    //  Prepare multithreadhing
    let seed = u64::from_le_bytes(R::random_bytes::<8>());
    let rng = Rng::with_seed(seed);
    let n_cores = if n_cores != 0 {
        n_cores
    } else {
        match std::thread::available_parallelism() {
            Ok(x) => x.get(),
            Err(_) => 0,
        }
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_cores)
        .build()
        .unwrap();
    bantha_core(gsh, candidates, timers, rng, pool, n_cores, pc)
}

#[roxido]
fn bantha_psm(gsh: &mut RExternalPtr, candidates: &RArray, n_cores: usize) {
    let gsh = RExternalPtr::decode_ref::<Gsh>(gsh);
    let candidates = candidates.to_i32(pc);
    let slice = candidates.slice();
    let dim = candidates.dim();
    if dim.len() != 3 {
        stop!("'candidates' must be three dimensional");
    }
    let n_items = dim[0];
    if n_items != dim[1] {
        stop!("The first two elements of 'x' should give a square matrix")
    }
    let n_samples = dim[2];
    let mut candidates2: Vec<DagAdjacencyVector> = Vec::with_capacity(n_samples);
    let len = n_items.pow(2);
    for i in 0..n_samples {
        let Ok(dag) =
            DagAdjacencyVector::from_slice_with_n(&slice[(i * len)..((i + 1) * len)], n_items)
        else {
            stop!("Element {} in 'candidates' is not a DAG", i);
        };
        candidates2.push(dag);
    }
    let timers = Timers::new();
    //  Prepare multithreadhing
    let seed = u64::from_le_bytes(R::random_bytes::<8>());
    let rng = Rng::with_seed(seed);
    let n_cores = if n_cores != 0 {
        n_cores
    } else {
        match std::thread::available_parallelism() {
            Ok(x) => x.get(),
            Err(_) => 0,
        }
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_cores)
        .build()
        .unwrap();
    bantha_core(gsh, candidates2, timers, rng, pool, n_cores, pc)
}

#[roxido]
fn bantha(samples: &RArray, a: f64, n_starts: usize, n_cores: usize) {
    let mut timers = Timers::new();
    timers.overall.tic();
    // Validate input
    if !(0.0..=2.0).contains(&a) {
        stop!("'a' must be between 0.0 and 2.0");
    }
    let dim = samples.dim();
    if dim.len() != 3 {
        stop!("'x' must be three dimensional");
    }
    let n_items = dim[0];
    if n_items != dim[1] {
        stop!("The first two elements of 'x' should give a square matrix")
    }
    let n_samples = dim[2];
    let n_candidates = if n_starts > 0 {
        n_starts.min(n_samples)
    } else {
        n_samples
    };
    //  Prepare multithreadhing
    let seed = u64::from_le_bytes(R::random_bytes::<8>());
    let mut rng = Rng::with_seed(seed);
    let n_cores = if n_cores != 0 {
        n_cores
    } else {
        match std::thread::available_parallelism() {
            Ok(x) => x.get(),
            Err(_) => 0,
        }
    };
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(n_cores)
        .build()
        .unwrap();
    // Compute
    let tasks = rng.choose_multiple(0..n_samples, n_candidates);
    let len = n_items.pow(2);
    timers.transform.tic();
    let slice = samples.to_i32(pc).slice();
    timers.transform.toc();
    let n_samples_f64 = n_samples as f64;
    timers.proportions.tic();
    let proportions = if n_cores == 1 {
        let mut proportions = Vec::with_capacity(len);
        for j in 0..proportions.capacity() {
            let mut count = 0;
            for i in 0..n_samples {
                count += slice[i * len + j];
            }
            proportions.push(count as f64 / n_samples_f64);
        }
        proportions
    } else {
        let mut proportions = vec![0.0; len];
        let chunck_size = if len % n_cores != 0 {
            len / n_cores + 1
        } else {
            len / n_cores
        };
        let chuncks = proportions.par_chunks_mut(chunck_size).enumerate();
        pool.install(|| {
            chuncks.into_par_iter().for_each(|(index1, x)| {
                for (index2, y) in x.iter_mut().enumerate() {
                    let j = chunck_size * index1 + index2;
                    let mut count = 0;
                    for i in 0..n_samples {
                        count += slice[i * len + j];
                    }
                    *y = count as f64 / n_samples_f64;
                }
            })
        });
        proportions
    };
    let gsh = Gsh { proportions, a };
    let mut candidates = Vec::with_capacity(tasks.len());
    timers.proportions.toc();
    if n_cores == 1 {
        tasks.iter().for_each(|&i| {
            let slice = &slice[i * len..(i + 1) * len];
            let Ok(dag) = DagAdjacencyVector::from_slice_with_n(slice, n_items) else {
                stop!("Found an element of 'samples' which is not a DAG.");
            };
            candidates.push(dag)
        });
    } else {
        let mapped: Vec<_> = pool
            .install(|| {
                tasks.par_iter().map(|&i| {
                    let slice = &slice[i * len..(i + 1) * len];
                    let Ok(dag) = DagAdjacencyVector::from_slice_with_n(slice, n_items) else {
                        return Err("Found an element of 'samples' which is not a DAG.");
                    };
                    Ok(dag)
                })
            })
            .collect();
        for x in mapped.into_iter() {
            let dag = x.stop();
            candidates.push(dag);
        }
    };
    timers.overall.toc();
    bantha_core(&gsh, candidates, timers, rng, pool, n_cores, pc)
}
