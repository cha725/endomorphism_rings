# Endomorphism Rings

Compute the projective modules of the endomorphism ring of a quiver algebra and all its radical powers.

---

## Overview

This Python project takes a quiver algebra A and computes the projective modules of its endomorphism ring:

$$\mathrm{End}_A(A + \mathrm{rad}^i(A)  \mid  i > 0)$$

where $\mathrm{rad}(A)$ is the radical of $A$, $\mathrm{rad}^2(A) = \mathrm{rad}(\mathrm{rad}(A))$ is the second radical power of $A$ and so on.
Note: so long as the quiver algebra $A$ is admissible there exists an $n$ such that $\mathrm{rad}^n(A) = 0$.

---

## Quiver Algebras

A quiver algebra is a finite-dimensional algebra over a field, constructed from a directed graph (a quiver).

1. Define a directed graph $Q$.  
2. Let $kQ$ be the vector space (over a field $k$) whose basis consists of all paths in $Q$.  
   - Each vertex $v_i$ contributes a "stationary path" $e_i$.  
3. Define multiplication in $kQ$ by concatenation of paths:  
   - For paths $p$ and $q$, define $p * q = pq$ if concatenation exists, and $0$ otherwise.  
   - Extend this linearly over $kQ$.  
4. The space $kQ$ is typically infinite-dimensional.  
   Example: a graph with one vertex and one loop corresponds to the polynomial ring $k[x]$.  
5. A quiver algebra is obtained by taking a quotient:  

   $A = kQ / I$  

   where $I$ is an ideal making $A$ finite-dimensional.  
   Example: $A = k[x] / <x^2>$ is 2-dimensional.

---

## Radicals

We consider the Jacobsen radical of a module.
For a right module $M$, $\mathrm{rad}(M)$ is the intersection of all the maximal right ideals of $M$.

So, $\mathrm{rad}(A)$ is the intersection of all the maximal right ideals of $A$ when treated as a right $A$-module.

---

## Representing Modules as Graphs

Some $A$-modules $M$ can be represented as directed graphs:

- Fix an additive generating set $G$ of $M$.
- For each generator in $g \in G$, add a vertex $v_g$ to the graph.
- For each algebra element $a$, add an arrow labeled $a$ from vertex $v_g$ to vertex $v_h$ if $h = g * a$ (that is, multiplication by $a$ maps $g$ to $h$).  
- In general these graphs do not have a simple interpretation because the generating set is not necessarily (and usually is not) a basis.  
- All elements of $M$ are given by linear combinations of the vertices in the graph, though not necessarily in a unique way.

Example:  
The polynomial ring $k[x]$ has an additive generator set $\{ x^i  \mid  i > 0\}$ which can be represented by an infinite graph:
$$ \mathrm{id}_{k[x]} \xrightarrow{x} x \xrightarrow{x} x^2 \xrightarrow{x} \cdots $$
A polynomial is a linear combination of the vertices in this graph.

---

## Projective Modules as Graphs

Projective modules are direct summands of free modules.  
In the quiver algebra setting, the elements of an indecomposable projective module correspond to the paths in the quiver that start at a particular vertex.  
Therefore, the diagram of a projective module can be generated from a single stationary path.  
The code automatically builds these diagrams from the quiver algebra data.

Injective modules are projective modules over the opposite quiver (where all arrows are reversed), so we can create diagrams for them using the same approach.

---

## Computing the Endomorphism Ring

To compute $$\mathrm{End}_A(A + \mathrm{rad}^i(A)  \mid  i > 0)$$, we find all module homomorphisms between these graph representations.

A module homomorphism $M \to N$ is uniquely determined by a pair $(Q, S)$, where:

- $Q$ is a quotient submodule of $M$.  
- $S$ is a submodule of $N$.  
- $Q$ and $S$ are isomorphic, representing the image of the homomorphism.

In graph terms:

- $Q$ is a quotient graph, closed under inverse module actions (moving against arrows).  
- $S$ is a subgraph, closed under module actions (moving with arrows).

Each valid $(Q, S)$ pair defines a module homomorphism and thus an element of the endomorphism ring.

The endomorphism ring itself can also be represented as a graph:
- Vertices represent homomorphisms.  
- Arrows represent ring multiplication.

---

## Summary

This project:
- Constructs quiver algebras.  
- Models projective and injective modules as directed graphs.  
- Computes the endomorphism ring $$\mathrm{End}_A(A + \mathrm{rad}^i(A)  \mid  i > 0)$$ through graph-based homomorphisms.
