
// Product factorization algorithm:
// - fn1: factorize_product: input: t_1 t_2 ... t_n, tensor registry, and (optional) resource constraints
//   - for each i in 1..n
//     - for each j in i+1..n
//       - initialize path (will record best contraction path + its cost)
//       - contract_all(t_1 t_2 ... t_n, i, j, registry, path, nflops=0)
//
// - fn2: contract_all: product = t_1 ... t_n, i, j, registry, path, nflops
//     - record tt_{ij} = t_i t_j: add it to the registry (need to be able to hash a binary tensor expr), evaluate its "cost" (# of operations, size, etc.) and add to the running total of the current cost
//     - if length(product) > 2:
//       - red_product = t_1 ... t_{i-1} tt_{ij} t_{i+1} ... t_{j-1} t_{j+1} ... t_n
//       - for each k in 1..n
//         - for each l in k+1..n
//           - contract_all(red_product, k, l, registry, path, nflops)
//     - else:
//       compare the cost to path, update path if needed
//