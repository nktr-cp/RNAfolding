# RNA Secondary Structure Prediction

This repository implements multiple algorithms for RNA secondary structure prediction, focusing on dynamic programming approaches.

- Nussinov Algorithm: Predicts the optimal secondary structure by maximizing the number of base pairs.
- Zuker Algorithm: Finds a thermodynamically stable structure by minimizing free energy.

---

## Nussinov Algorithm Overview

The Nussinov algorithm was first introduced by Nussinov et al. (1978) in *SIAM Journal of Applied Mathematics*. This algorithm uses dynamic programming to predict RNA secondary structures by maximizing the number of complementary base pairs.

Given a nucleotide sequence of length L, the algorithm assumes:
- δ(i, j) = 1 if the i-th and j-th nucleotides can form a base pair (A-U, G-C, or G-U). 
- δ(i, j) = 0 if they cannot form a pair or if j - i - 1 < 3 (a loop cannot be shorter than 3 unpaired bases).

### Initialization
```
γ(i, i) = 0  for i = 0 to L-1
γ(i, i - 1) = 0  for i = 1 to L-1
```

### Recursion
```
γ(i, j) = max {
    γ(i + 1, j),
    γ(i, j - 1),
    γ(i + 1, j - 1) + δ(i, j),
    max(i ≤ k < j) [γ(i, k) + γ(k + 1, j)]
}
```
- γ(i, j) represents the maximum number of base pairs for the subsequence from position i to j.

### Traceback
Use stack for traceback to reconstruct the optimal base-pairing structure:

#### Initialization
Push the pair `(0, L-1)` onto a stack, where L is the length of the RNA sequence.

#### Recursion
Repeat the following until the stack is empty:
   
   - Pop `(i, j)` from the stack.
   
   - If `i >= j`, skip this iteration and continue with the next.
   
   - Otherwise, check which of the following cases holds:
     - If `γ(i + 1, j) == γ(i, j)`, push `(i + 1, j)` onto the stack.
     - If `γ(i, j - 1) == γ(i, j)`, push `(i, j - 1)` onto the stack.
     - If `γ(i + 1, j - 1) + δ(i, j) == γ(i, j)`:
       - Record the base pair `(i, j)`.
       - Push `(i + 1, j - 1)` onto the stack.
     
   - If none of the above cases hold, loop over all `k` from `i` to `j - 1` and check if:
     - `γ(i, k) + γ(k + 1, j) == γ(i, j)`.
     - If true, push `(k + 1, j)` and `(i, k)` onto the stack, then break from the loop.
   
