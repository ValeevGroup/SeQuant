//
// Created by Bimal Gaudel on 12/20/19.
//

// A tree data structure of positions of SeQuant Expr objects
// has to be implemented. Choosing positions to be ordinals makes it
// easier for eliminating redundant path computations.
// 
// Such a tree for a product of four tensors will have paths that look like
// this:
// (0, 1, 2, 3)
// (0, 1, 3, 2)
// (0, 2, 1, 3)
// (0, 2, 3, 1)
// (0, 3, 1, 2)
// (0, 3, 2, 1)
// (1, 2, 0, 3)
// (1, 2, 3, 0)
// ((0, 3), (1, 2))
// (1, 3, 0, 2)
// (1, 3, 2, 0)
// ((0, 2), (1, 3))
// (2, 3, 0, 1)
// (2, 3, 1, 0)
// ((0, 1), (2, 3))
//
// Once we have this data structure handy, it'll be easy to compute the
// computational cost of each path and we can reform the optimal product
// expression of the given Product.
