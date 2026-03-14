// The 'roxido_registration' macro is called at the start of the 'lib.rs' file.
roxido_registration!();
use num_traits::PrimInt;
use rand::{RngExt, rng};
use rand_distr::weighted::WeightedAliasIndex;
use rand_distr::{Beta, Distribution, Normal};
use roxido::*;
use std::collections::HashMap;
use std::fmt::{self, Display, Formatter};

const TOL: f64 = 0.0;

#[derive(Debug)]
pub enum BalancedTree<T> {
    Leaf {
        index: usize, // breadth-first index *within the leaves*
        weight: f64,
        atom: T,
    },
    Node {
        index: usize, // breadth-first index of the internal node
        break_proportion: f64,
        weight: f64,
        left: Box<BalancedTree<T>>,
        right: Box<BalancedTree<T>>,
    },
}

impl<T> BalancedTree<T> {
    pub fn new<FW, FA>(depth: u16, mut break_fn: FW, mut atom_fn: FA) -> Self
    where
        FW: FnMut(usize) -> f64,
        FA: FnMut() -> T,
    {
        // We still want leaves numbered 0,1,2,… in creation order
        let mut leaf_indexer = 0;
        Self::build(
            depth,
            1.0, // root weight
            0,   // root index in breadth-first order
            &mut break_fn,
            &mut atom_fn,
            &mut leaf_indexer,
        )
    }

    fn build<FW, FA>(
        levels_remaining: u16,
        weight: f64,
        index: usize, // the node’s *desired* breadth-first index
        break_fn: &mut FW,
        atom_fn: &mut FA,
        leaf_indexer: &mut usize,
    ) -> Self
    where
        FW: FnMut(usize) -> f64,
        FA: FnMut() -> T,
    {
        if levels_remaining == 0 || weight == 0.0 {
            // leaf: give it a sequential index among leaves
            let index = *leaf_indexer;
            *leaf_indexer += 1;
            BalancedTree::Leaf {
                index,
                weight,
                atom: atom_fn(),
            }
        } else {
            // internal node: index already known
            let break_proportion = break_fn(index);
            BalancedTree::Node {
                index,
                break_proportion,
                weight,
                left: Box::new(Self::build(
                    levels_remaining - 1,
                    weight * break_proportion,
                    2 * index + 1, // breadth-first index for the left child
                    break_fn,
                    atom_fn,
                    leaf_indexer,
                )),
                right: Box::new(Self::build(
                    levels_remaining - 1,
                    weight * (1.0 - break_proportion),
                    2 * index + 2, // breadth-first index for the right child
                    break_fn,
                    atom_fn,
                    leaf_indexer,
                )),
            }
        }
    }

    /// Compute **all** leaf weights in one pass.
    /// Returns a Vec<f64> of length 2^depth, indexed by 0-based leaf label.
    pub fn breaks_atomfn_to_discrete_measure<FA>(
        breaks: &[f64],
        mut atom_fn: FA,
    ) -> Result<DiscreteMeasure<T>, &'static str>
    where
        FA: FnMut() -> T,
    {
        let leaf_count = breaks.len() + 1;
        let mut atoms = Vec::with_capacity(leaf_count);
        for _ in 0..leaf_count {
            atoms.push(atom_fn());
        }
        Self::breaks_atoms_to_discrete_measure(breaks, atoms)
    }

    pub fn breaks_atoms_to_discrete_measure(
        breaks: &[f64],
        atoms: Vec<T>,
    ) -> Result<DiscreteMeasure<T>, &'static str> {
        // Verify that len + 1 is a power of two, i.e. breaks = 2^d − 1
        let n_plus_one = breaks.len() + 1;
        if !n_plus_one.is_power_of_two() {
            return Err("Length of 'breaks' is not 2^depth - 1.");
        };
        // Because len + 1 = 2^depth
        let mut depth = 0;
        let mut x = breaks.len() + 1;
        while x > 1 {
            x >>= 1;
            depth += 1;
        }
        let leaf_count = 1usize << depth; // 2^depth
        if atoms.len() != leaf_count {
            return Err("Length of 'atoms' is not 2^depth.");
        }
        if depth == 0 {
            // Single node tree: weight is one because there are no breaks.
            return Ok(DiscreteMeasure {
                weights: vec![1.0],
                atoms,
            });
        }
        let first_leaf = breaks.len(); // heap index of leftmost leaf
        let mut weights = vec![0.0; leaf_count]; // output
        let mut stack = Vec::<(usize, f64)>::with_capacity(depth + 1);
        // (heap_index, cumulative_weight_so_far)
        stack.push((0, 1.0)); // start at root
        while let Some((idx, w)) = stack.pop() {
            if idx < breaks.len() {
                // internal node
                let b = breaks[idx];
                let left = idx * 2 + 1;
                let right = left + 1;
                // push right first so left is processed first (keeps order nice)
                stack.push((right, w * (1.0 - b))); // go right
                stack.push((left, w * b)); // go left
            } else {
                // leaf: store its final weight
                weights[idx - first_leaf] = w;
            }
        }
        Ok(DiscreteMeasure { weights, atoms })
    }

    /// Count leaves whose weight is strictly > 0.0.
    pub fn nonzero_leaves(&self) -> usize {
        // use an explicit stack to avoid deep recursion limits
        let mut stack = vec![self];
        let mut count = 0;
        while let Some(node) = stack.pop() {
            match node {
                BalancedTree::Leaf { weight, .. } => {
                    if *weight > TOL {
                        count += 1;
                    }
                }
                BalancedTree::Node {
                    weight,
                    left,
                    right,
                    ..
                } => {
                    if *weight > TOL {
                        // push children only when the parent’s weight is non-zero
                        stack.push(left);
                        stack.push(right);
                    }
                }
            }
        }
        count
    }
}

impl<T: Display> Display for BalancedTree<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        fn fmt_inner<T: Display>(
            tree: &BalancedTree<T>,
            f: &mut Formatter<'_>,
            indent: usize,
        ) -> fmt::Result {
            let pad = "  ".repeat(indent);
            // helper: convert index according to the chosen base
            let idx = |i: usize| if f.alternate() { i + 1 } else { i };
            match tree {
                BalancedTree::Leaf {
                    index,
                    weight,
                    atom,
                } => writeln!(
                    f,
                    "{pad}Leaf(index = {}, weight = {:.6}, atom = {atom})",
                    idx(*index),
                    weight,
                ),
                BalancedTree::Node {
                    index,
                    break_proportion,
                    weight,
                    left,
                    right,
                } => {
                    writeln!(
                        f,
                        "{pad}Node(index = {}, break_proportion = {break_proportion}, weight = {:.6})",
                        idx(*index),
                        weight,
                    )?;
                    fmt_inner(left, f, indent + 1)?;
                    fmt_inner(right, f, indent + 1)
                }
            }
        }

        fmt_inner(self, f, 0)
    }
}

struct CompactTreeDisplay<'a, T>(&'a BalancedTree<T>);

impl<'a, T: Display> Display for CompactTreeDisplay<'a, T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        fn fmt_inner<T: Display>(
            tree: &CompactTreeDisplay<T>,
            f: &mut Formatter<'_>,
            indent: usize,
        ) -> fmt::Result {
            let pad = "  ".repeat(indent);
            // helper: convert index according to the chosen base
            let idx = |i: usize| if f.alternate() { i + 1 } else { i };
            match tree.0 {
                BalancedTree::Leaf {
                    index,
                    weight,
                    atom,
                } => {
                    if *weight > TOL {
                        writeln!(
                            f,
                            "{pad}Leaf(index = {}, weight = {:.6}, atom = {atom})",
                            idx(*index),
                            weight,
                        )
                    } else {
                        Ok(())
                    }
                }
                BalancedTree::Node {
                    index,
                    break_proportion,
                    weight,
                    left,
                    right,
                } => {
                    if *weight > TOL {
                        writeln!(
                            f,
                            "{pad}Node(index = {}, break_proportion = {break_proportion}, weight = {:.6})",
                            idx(*index),
                            weight,
                        )?;
                        fmt_inner(&CompactTreeDisplay(left), f, indent + 1)?;
                        fmt_inner(&CompactTreeDisplay(right), f, indent + 1)
                    } else {
                        Ok(())
                    }
                }
            }
        }

        fmt_inner(self, f, 0)
    }
}

impl BalancedTree<f64> {
    pub fn new_beta_gaussian(
        depth: u16,
        alpha: f64,
        beta: f64,
        mean: f64,
        sd: f64,
    ) -> Result<Self, &'static str> {
        let mut rng1 = rng();
        let mut rng2 = rng();

        let beta_dist = Beta::new(alpha, beta).map_err(|_| "'alpha' or 'beta' are not valid.")?;
        let normal_dist = Normal::new(mean, sd).map_err(|_| "'sd' must be greater than 0.0")?;

        let mut break_fn = |_| beta_dist.sample(&mut rng1);
        let mut atom_fn = || normal_dist.sample(&mut rng2);

        Ok(Self::new(depth, &mut break_fn, &mut atom_fn))
    }
}

impl<T: Clone> BalancedTree<T> {
    pub fn to_discrete_measure(&self) -> DiscreteMeasure<T> {
        let mut weights = Vec::new();
        let mut atoms = Vec::new();
        self.to_discrete_measure_recursion(&mut weights, &mut atoms);
        DiscreteMeasure { weights, atoms }
    }

    fn to_discrete_measure_recursion(&self, weights: &mut Vec<f64>, atoms: &mut Vec<T>) {
        match self {
            BalancedTree::Leaf { weight, atom, .. } => {
                if *weight > 0.0 {
                    weights.push(*weight);
                    atoms.push(atom.clone());
                }
            }
            BalancedTree::Node { left, right, .. } => {
                left.to_discrete_measure_recursion(weights, atoms);
                right.to_discrete_measure_recursion(weights, atoms);
            }
        }
    }
}

#[derive(Debug)]
pub struct DiscreteMeasure<T> {
    weights: Vec<f64>,
    atoms: Vec<T>,
}

impl<T: Clone> DiscreteMeasure<T> {
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }

    pub fn atoms(&self) -> &[T] {
        &self.atoms
    }

    /// Count leaves whose weight is strictly > 0.0.
    pub fn nonzero_leaves(&self) -> usize {
        self.weights.iter().filter(|x| **x > 0.0).count()
    }

    pub fn sample_partition<I>(&self, n_items: usize) -> Vec<I>
    where
        I: PrimInt,
    {
        let mut rng = rng();
        let alias = WeightedAliasIndex::new(self.weights.clone()).expect("Invalid weights.");
        let mut map = HashMap::new();
        let mut result = Vec::new();
        let mut next_label = I::zero();
        for _ in 0..n_items {
            let index = alias.sample(&mut rng);
            result.push(*map.entry(index).or_insert_with(|| {
                let label = next_label;
                next_label = next_label + I::one();
                label
            }));
        }
        result
    }

    pub fn sample_partition_into<I>(&self, n_items: usize, slice: &mut [I])
    where
        I: PrimInt,
    {
        let mut rng = rng();
        let alias = WeightedAliasIndex::new(self.weights.clone()).expect("Invalid weights.");
        let mut map = HashMap::new();
        let mut next_label = I::zero();
        for (counter, y) in slice.iter_mut().enumerate() {
            if counter % n_items == 0 {
                map.clear();
                next_label = I::zero();
            }
            let index = alias.sample(&mut rng);
            *y = *map.entry(index).or_insert_with(|| {
                let label = next_label;
                next_label = next_label + I::one();
                label
            });
        }
    }

    pub fn sample(&self) -> T {
        let mut rng = rng();
        let cum: Vec<f64> = self
            .weights
            .iter()
            .scan(0.0, |s, &w| {
                *s += w;
                Some(*s)
            })
            .collect();
        let u = rng.random::<f64>() * cum.last().unwrap();
        let index = cum.partition_point(|&c| c < u);
        self.atoms[index].clone()
    }

    pub fn sample_into(&self, slice: &mut [T]) {
        let mut rng = rng();
        let alias = WeightedAliasIndex::new(self.weights.clone()).expect("Invalid weights.");
        for y in slice.iter_mut() {
            let index = alias.sample(&mut rng);
            *y = self.atoms[index].clone();
        }
    }
}

#[roxido]
fn new_tree_beta_gaussian(depth: usize, alpha: f64, beta: f64, mean: f64, sd: f64) {
    let tree =
        BalancedTree::new_beta_gaussian(depth.try_into().stop(), alpha, beta, mean, sd).stop();
    let ptr = RExternalPtr::encode(tree, "treesb_tree", pc);
    ptr.set_class(RVector::<char>::from_value("treesb", 1, pc));
    ptr
}

fn safe_usize_plus1_to_f64(x: usize) -> Result<f64, &'static str> {
    const MAX_SAFE_F64_INT: usize = 9_007_199_254_740_992;
    let y = x.checked_add(1).ok_or("Overflow.")?;
    if y > MAX_SAFE_F64_INT {
        return Err("Not exactly representable.");
    }
    Ok(y as f64)
}

#[roxido]
fn new_tree_from_funcs(depth: usize, break_fn: &RFunction, atom_fn: &RFunction) {
    let break_fn = |index: usize| {
        let pc2 = &mut Pc::default(); // Avoid blowing the protection stack
        let index = RScalar::from_value(safe_usize_plus1_to_f64(index).stop(), pc2);
        let result = break_fn.call1(index, pc2).stop();
        result.as_scalar().stop().as_f64().stop().get()
    };
    let atom_fn = || {
        let pc2 = &mut Pc::default(); // Avoid blowing the protection stack
        let result = atom_fn.call0(pc2).stop();
        result.as_scalar().stop().as_f64().stop().get()
    };
    let tree = BalancedTree::new(depth.try_into().stop(), break_fn, atom_fn);
    let ptr = RExternalPtr::encode(tree, "treesb_tree", pc);
    ptr.set_class(RVector::<char>::from_value("treesb", 1, pc));
    ptr
}

#[roxido]
fn new_dm_beta_gaussian(depth: usize, alpha: f64, beta: f64, mean: f64, sd: f64, method: &str) {
    let dm = match method {
        "via_tree" => {
            let tree =
                BalancedTree::new_beta_gaussian(depth.try_into().stop(), alpha, beta, mean, sd)
                    .stop();
            tree.to_discrete_measure()
        }
        "direct" => {
            let mut rng = rng();
            let beta_dist = Beta::new(alpha, beta).stop();
            let count = (1usize << depth) - 1;
            let mut breaks = Vec::with_capacity(count);
            for _ in 0..count {
                breaks.push(beta_dist.sample(&mut rng));
            }
            let normal_dist = Normal::new(mean, sd).stop();
            let atom_fn = || normal_dist.sample(&mut rng);
            BalancedTree::breaks_atomfn_to_discrete_measure(&breaks, atom_fn).stop()
        }
        _ => {
            stop!("'method' not recognized.")
        }
    };
    let ptr = RExternalPtr::encode(dm, "treesb_dm", pc);
    ptr.set_class(RVector::<char>::from_value("treesb", 1, pc));
    ptr
}

#[roxido]
fn new_dm_from_breaks_and_atoms(breaks: &[f64], atoms: &[f64]) {
    let dm = BalancedTree::breaks_atoms_to_discrete_measure(breaks, atoms.to_vec()).stop();
    let ptr = RExternalPtr::encode(dm, "treesb_dm", pc);
    ptr.set_class(RVector::<char>::from_value("treesb", 1, pc));
    ptr
}

#[roxido]
fn nonzero_leaves(x: &RExternalPtr) {
    let count = match x.tag_str() {
        "treesb_tree" => {
            let tree = x.decode_ref::<BalancedTree<f64>>();
            tree.nonzero_leaves()
        }
        "treesb_dm" => {
            let dm = x.decode_ref::<DiscreteMeasure<f64>>();
            dm.nonzero_leaves()
        }
        _ => {
            stop!("Unsupported type.")
        }
    };
    i32::try_from(count).stop()
}

#[roxido]
fn dm_sample_partition(x: &RExternalPtr, n_items: usize, n_samples: usize) {
    match x.tag_str() {
        "treesb_dm" => {
            let dm = x.decode_ref::<DiscreteMeasure<f64>>();
            let result = RMatrix::<i32>::new(n_items, n_samples, pc);
            let slice = result.slice_mut();
            dm.sample_partition_into(n_items, slice);
            result
        }
        _ => stop!("Unsupported type."),
    }
}

#[roxido]
fn dm_sample(x: &RExternalPtr, n_samples: usize) {
    match x.tag_str() {
        "treesb_dm" => {
            let dm = x.decode_ref::<DiscreteMeasure<f64>>();
            let result = RVector::from_value(0.0, n_samples, pc);
            let slice = result.slice_mut();
            dm.sample_into(slice);
            result
        }
        _ => stop!("Unsupported type."),
    }
}

#[roxido]
fn dm_weights(x: &RExternalPtr) {
    match x.tag_str() {
        "treesb_dm" => {
            let dm = x.decode_ref::<DiscreteMeasure<f64>>();
            (&dm.weights[..]).to_r(pc)
        }
        _ => stop!("Unsupported type."),
    }
}

#[roxido]
fn dm_atoms(x: &RExternalPtr) {
    match x.tag_str() {
        "treesb_dm" => {
            let dm = x.decode_ref::<DiscreteMeasure<f64>>();
            (&dm.atoms[..]).to_r(pc)
        }
        _ => stop!("Unsupported type."),
    }
}

#[roxido]
fn print_treesb(x: &RExternalPtr, include_zeros: bool) {
    match x.tag_str() {
        "treesb_tree" => {
            let tree = x.decode_ref::<BalancedTree<f64>>();
            if include_zeros {
                rprint!("{tree:#}");
            } else {
                rprint!("{:#}", CompactTreeDisplay(tree));
            }
        }
        "treesb_dm" => {
            let dm = x.decode_ref::<DiscreteMeasure<f64>>();
            rprintln!("{dm:?}");
        }
        _ => {
            stop!("Unsupported type.")
        }
    }
}

#[roxido]
fn to_discrete_measure(x: &RExternalPtr) {
    let ptr = match x.tag_str() {
        "treesb_tree" => {
            let tree = x.decode_ref::<BalancedTree<f64>>();
            let dm = tree.to_discrete_measure();
            RExternalPtr::encode(dm, "treesb_dm", pc)
        }
        _ => {
            stop!("Unsupported type.")
        }
    };
    ptr.set_class(RVector::<char>::from_value("treesb", 1, pc));
    ptr
}

/// Path information from the root to a given leaf (0-based heap indexing).
#[derive(Debug)]
pub struct PathInfo {
    nodes: Vec<usize>,     // zero-based node indexing
    directions: Vec<bool>, // 'false' is "go left", 'true' is "go right"
}

impl PathInfo {
    pub fn weight(&self, breaks: &[f64]) -> f64 {
        let mut sum = 0.0;
        for (&node, &direction) in self.nodes.iter().zip(self.directions.iter()) {
            let b = breaks[node];
            sum += if direction { 1.0 - b } else { b };
        }
        sum
    }
}

/// Compute both the internal-node indices and left/right directions in one pass.
/// * `depth` may be 0‥= `usize::BITS - 1`  (≤ 63 on a 64-bit build).
/// * `leaf_index` must satisfy `0 ≤ leaf_index < 2^depth` (0-based labels).
pub fn path_to_leaf(leaf_index: usize, depth: usize) -> Result<PathInfo, &'static str> {
    if depth > usize::BITS as usize - 1 {
        return Err("'depth' too large for usize.");
    };
    if leaf_index >= (1usize << depth) {
        return Err("'leaf' index is out of range");
    };

    // Heap index of the leaf (0-based): first leaf = 2^depth − 1
    let mut idx = leaf_index + ((1usize << depth) - 1);

    // Pre-allocate final size and fill from the back toward the front
    let mut nodes = vec![0usize; depth];
    let mut directions = vec![false; depth];

    for lvl in (0..depth).rev() {
        let parent = (idx - 1) >> 1; // parent’s heap index
        nodes[lvl] = parent; // store node
        directions[lvl] = (idx & 1) == 0; // true = right, false = left
        idx = parent; // climb up for next iteration
    }

    Ok(PathInfo { nodes, directions })
}

/// Like `path_to_leaf`, but for a *target internal node*.
///
/// * `depth`  – height of the full tree (0-based leaves start at 2^depth−1).
/// * `node_index` must satisfy 0 ≤ node_index < 2^depth − 1  (i.e. it is an
///   internal node; leaves are not allowed).
///
/// Returns a `PathInfo` whose
///   * `nodes[i]`      is the i-th internal ancestor (root first);
///   * `directions[i]` is `false` if you go **left** at `nodes[i]`,
///     `true`  if you go **right** on the way to `node_index`.
///
/// Special cases
/// -------------
/// * `depth == 0` ⇒ the tree has no internal nodes at all – the function
///   returns `nodes = []`, `directions = []`.
/// * `node_index == 0` (the root itself) ⇒ there are no ancestors, so both
///   vectors are empty.
pub fn path_to_internal(node_index: usize, depth: usize) -> Result<PathInfo, &'static str> {
    // ---------- argument checks -------------------------------------------
    if depth > usize::BITS as usize - 1 {
        return Err("'depth' too large for usize.");
    }
    let first_leaf = (1usize << depth) - 1; // 2^d − 1
    if node_index >= first_leaf {
        return Err("index is not an internal node for this depth.");
    }

    // ---------- determine path length (number of ancestors) ---------------
    let mut len = 0;
    let mut tmp = node_index;
    while tmp > 0 {
        tmp = (tmp - 1) >> 1; // move to parent
        len += 1;
    }

    // ---------- allocate and fill *from the back* -------------------------
    let mut nodes = vec![0usize; len];
    let mut directions = vec![false; len];

    let mut child = node_index;
    for lvl in (0..len).rev() {
        let parent = (child - 1) >> 1;
        nodes[lvl] = parent;
        directions[lvl] = (child & 1) == 0; // true = right, false = left
        child = parent;
    }

    Ok(PathInfo { nodes, directions })
}

#[roxido]
fn path_to_leaf_r(leaf_index: usize, depth: usize) {
    let path = path_to_leaf(leaf_index - 1, depth).stop();
    path_to_r(path, pc).sexp()
}

#[roxido]
fn path_to_internal_r(node_index: usize, depth: usize) {
    let path = path_to_internal(node_index - 1, depth).stop();
    path_to_r(path, pc).sexp()
}

fn path_to_r(path: PathInfo, pc: &mut Pc) -> &RList {
    let nodes = RVector::<i32>::new(path.nodes.len(), pc);
    let nodes_slice = nodes.slice_mut();
    for (x, y) in path.nodes.iter().zip(nodes_slice) {
        *y = i32::try_from(*x + 1).stop();
    }
    let directions = RVector::<bool>::new(path.directions.len(), pc);
    let directions_slice = directions.slice_mut();
    for (x, y) in path.directions.iter().zip(directions_slice) {
        *y = if *x { R::TRUE() } else { R::FALSE() };
    }
    let list = RList::with_names(&["nodes", "directions"], pc);
    list.set(0, nodes).stop();
    list.set(1, directions).stop();
    list
}

/// A perfect binary tree whose *internal* nodes carry “break” values `b`.
/// For a given leaf, each internal node contributes
///
///   dir = left  →  weight +=  b
///   dir = right →  weight += (1 - b)
///
/// where the direction is determined by the leaf’s bit pattern
/// in the usual 0-based heap layout.
///
/// Example (`depth = 2`):
///
///                  (0)
///               /      \
///           (1)         (2)        <- internal node labels
///          /   \       /   \
///        [0]  [1]    [2]  [3]      <- leaf node labels
///
/// Every internal node has a break given by [b0, b1, b2]
///
/// leaf-0 weight = b0 * b1  
/// leaf-1 weight = b0 * (1-b1)  
/// leaf-2 weight = (1-b0) * b2  
/// leaf-3 weight = (1-b0) * (1-b2)
///
#[derive(Debug)]
pub struct TreeSb {
    /// Break value for every internal node using 0-based heap indices.
    /// Length must be 2^depth − 1.
    breaks: Vec<f64>,
}

impl TreeSb {
    /// Construct and sanity-check.
    pub fn new(breaks: Vec<f64>) -> Result<Self, &'static str> {
        // Verify that len + 1 is a power of two, i.e. breaks = 2^d − 1
        let n_plus_one = breaks.len() + 1;
        if !n_plus_one.is_power_of_two() {
            return Err("Length of 'breaks' is not 2^depth - 1");
        };
        Ok(Self { breaks })
    }

    /// Depth of the tree (number of edges from root to any leaf).
    #[inline]
    fn depth(&self) -> usize {
        // Because len + 1 = 2^depth
        let mut depth = 0;
        let mut x = self.breaks.len() + 1;
        while x > 1 {
            x >>= 1;
            depth += 1;
        }
        depth
    }

    /// Compute **all** leaf weights in one pass.
    /// Returns a Vec<f64> of length 2^depth, indexed by 0-based leaf label.
    pub fn weights(&self) -> Vec<f64> {
        let depth = self.depth();
        let leaf_count = 1usize << depth; // 2^depth
        if depth == 0 {
            // Single node tree: weight is one because there are no breaks.
            return vec![1.0];
        }

        let first_leaf = self.breaks.len(); // heap index of leftmost leaf
        let mut weights = vec![0.0; leaf_count]; // output
        let mut stack = Vec::<(usize, f64)>::with_capacity(depth + 1);

        // (heap_index, cumulative_weight_so_far)
        stack.push((0, 1.0)); // start at root

        while let Some((idx, w)) = stack.pop() {
            if idx < self.breaks.len() {
                // internal node
                let b = self.breaks[idx];
                let left = idx * 2 + 1;
                let right = left + 1;

                // push right first so left is processed first (keeps order nice)
                stack.push((right, w * (1.0 - b))); // go right
                stack.push((left, w * b)); // go left
            } else {
                // leaf: store its final weight
                weights[idx - first_leaf] = w;
            }
        }

        weights
    }
}

pub fn dm_weights_from_breaks_unchecked(breaks: &[f64], weights: &mut [f64], depth: usize) {
    if depth == 0 {
        // Single node tree: weight is one because there are no breaks.
        weights[0] = 1.0;
        return;
    }
    let first_leaf = breaks.len(); // heap index of leftmost leaf
    let mut stack = Vec::<(usize, f64)>::with_capacity(depth + 1);
    // (heap_index, cumulative_weight_so_far)
    stack.push((0, 1.0)); // start at root
    while let Some((idx, w)) = stack.pop() {
        if idx < breaks.len() {
            // internal node
            let b = breaks[idx];
            let left = idx * 2 + 1;
            let right = left + 1;
            // push right first so left is processed first (keeps order nice)
            stack.push((right, w * (1.0 - b))); // go right
            stack.push((left, w * b)); // go left
        } else {
            // leaf: store its final weight
            weights[idx - first_leaf] = w;
        }
    }
}

pub fn dm_weights_from_breaks(breaks: &[f64], weights: &mut [f64]) -> Result<(), &'static str> {
    let mut n_plus_one = breaks.len() + 1;
    if !n_plus_one.is_power_of_two() {
        return Err("Length of 'breaks' is not 2^depth - 1");
    };
    let mut depth = 0;
    while n_plus_one > 1 {
        n_plus_one >>= 1;
        depth += 1;
    }
    let leaf_count = 1usize << depth; // 2^depth
    if weights.len() != leaf_count {
        return Err("Length of 'weights' should be 2^depth.");
    }
    dm_weights_from_breaks_unchecked(breaks, weights, depth);
    Ok(())
}

#[roxido]
fn dm_weights_from_breaks_r(breaks: &RVector<f64>) {
    let mut n_plus_one = breaks.len() + 1;
    if !n_plus_one.is_power_of_two() {
        stop!("Length of 'breaks' is not 2^depth - 1");
    };
    let mut depth = 0;
    while n_plus_one > 1 {
        n_plus_one >>= 1;
        depth += 1;
    }
    let leaf_count = 1usize << depth; // 2^depth
    let weights = RVector::<f64>::new(leaf_count, pc);
    let weights_slice = weights.slice_mut();
    dm_weights_from_breaks_unchecked(breaks.slice(), weights_slice, depth);
    weights
}
