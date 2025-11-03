# Endomorphism Rings

Compute the projective modules of the endomorphism ring of a quiver algebra and its dual.

---

## Overview

This Python project takes a quiver algebra A and computes the projective modules of its endomorphism ring:

End_A(A + D(A))

where D(A) = Hom_k(A, k) is the dual of A as a vector space over a field k.

---

## Quiver Algebras

A quiver algebra is a finite-dimensional algebra over a field, constructed from a directed graph (a quiver).

1. Define a directed graph Q.  
2. Let kQ be the vector space (over k) whose basis consists of all paths in Q.  
   - Each vertex v_i contributes a "stationary path" e_i.  
3. Define multiplication in kQ by concatenating paths:  
   - For paths p and q, define p * q = pq if concatenation exists, and 0 otherwise.  
   - Extend this linearly over kQ.  
4. The space kQ is typically infinite-dimensional.  
   Example: a graph with one vertex and one loop corresponds to the polynomial ring k[x].  
5. A quiver algebra is obtained by taking a quotient:  

   A = kQ / I  

   where I is an ideal making A finite-dimensional.  
   Example: A = k[x] / <x^2> is 2-dimensional.

---

## Duality and Endomorphism Rings

Since A is a finite-dimensional vector space over k, we can take its dual:

D(A) = Hom_k(_A A, k)

We then compute the endomorphism ring:

End_A(A + D(A))

---

## Representing Modules as Graphs

Some A-modules M can be represented as directed graphs:

- The initial vertices correspond to elements in a generating set of M.  
- From each vertex, propagate by adding new vertices for all elements h such that g = h * a for some ring element a, and continue recursively.  
- For each ring element a, add an arrow labeled a from vertex v_g to vertex v_h if g = h * a (that is, multiplication by a maps h to g).  
- Vertices may be identified with each other, and in general these graphs do not have a simple interpretation because the generating set is not necessarily (and usually is not) a basis.  
- All elements of M are given by linear combinations of the vertices in the graph, though not necessarily in a unique way.

Example:  
The polynomial ring k[x] can be represented by an infinite graph:

    0 -> 1 -> 2 -> 3 -> ...

Here, vertex i corresponds to the element x^i, and the arrows represent multiplication by x.  
A polynomial is then a linear combination of the vertices in this graph.

---

## Projective Modules as Graphs

Projective modules are direct summands of free modules.  
In the quiver algebra setting, the elements of an indecomposable projective module correspond to the paths in the quiver that start at a particular vertex.  
Therefore, the diagram of a projective module can be generated from a single stationary path.  
The code automatically builds these diagrams from the quiver algebra data.

Injective modules are projective modules over the opposite quiver (where all arrows are reversed), so we can create diagrams for them using the same approach.

---

## Computing the Endomorphism Ring

To compute End_A(A + D(A)), we find all module homomorphisms between these graph representations.

A module homomorphism M -> N is uniquely determined by a pair (Q, S), where:

- Q is a quotient submodule of M.  
- S is a submodule of N.  
- Q and S are isomorphic, representing the image of the homomorphism.

In graph terms:

- Q is a quotient graph, closed under inverse module actions (moving against arrows).  
- S is a subgraph, closed under module actions (moving with arrows).

Each valid (Q, S) pair defines a module homomorphism and thus an element of the endomorphism ring.

The endomorphism ring itself can also be represented as a graph:
- Vertices represent homomorphisms.  
- Arrows represent ring multiplication.

---

## Summary

This project:
- Constructs quiver algebras.  
- Models projective and injective modules as directed graphs.  
- Computes the endomorphism ring End_A(A + D(A)) through graph-based homomorphism analysis.
