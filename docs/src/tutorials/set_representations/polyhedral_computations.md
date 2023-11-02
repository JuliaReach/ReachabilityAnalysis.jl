```@meta
DocTestSetup  = quote
    using ReachabilityAnalysis
end
CurrentModule = ReachabilityAnalysis
```

## Concrete polyhedra

## Polyhedra backends

## Limitations

Computing with concrete polyhedra in high dimensions is generally expensive. In particular, converting between vertex and constraint representations (so-called dual representations) should be used only if it is strictly necessary. However, there are some operations that are cheap:

- Intersecting two (or more) sets in constraint representation, or whose `constraints_list` can be computed efficiently. Such computation only requires concatenating the constraints and removing redundant inequalities (operation that requires the solution of linear programs).

- Taking linear maps of sets in vertex representation, $MX$. This operation requires to map each vertex of $X$ under the transformation $M$. Linear transformations can also be done efficiently in constraint representation provided that the matrix $M$ is invertible. LazySets handles other cases ($M$ not invertible, and the sets either in constraint or in vertex representation), but they are generally expensive in high dimensions. However, using specific classes of sets (e.g. zonotopes).
