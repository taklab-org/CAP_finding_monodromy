# Codes of "Computer-assisted proofs for finding the monodromy of Picard--Fuchs differential equations for a family of K3 toric hypersurfaces"

This repository contains the C++ codes associated with the paper:
"Computer-assisted proofs for finding the monodromy of Picard--Fuchs differential equations for a family of K3 toric hypersurfaces"
by T Ishige and A Takayasu.

**Abstract**  In this paper, we present a numerical method for rigorously finding the monodromy of linear differential equations. Beginning at a base point where certain particular solutions are explicitly given by series expansions, we first compute the value of fundamental system of solutions using interval arithmetic to rigorously control truncation and rounding errors. The solutions are then analytically continued along a prescribed contour encircling the singular points of the differential equation via a rigorous integrator. From these computations, the monodromy matrices are derived, generating the monodromy group of the differential equation. This method establishes a mathematically rigorous framework for addressing the monodromy problem in differential equations. For a notable example, we apply our computer-assisted proof method to resolve the monodromy problem for a Picard-Fuchs differential equation associated with a family of K3 toric hypersurfaces.

These codes require the [kv library](https://github.com/mskashi/kv) (a C++ Library for rigorous numerics) version 0.4.57. 

---

### Explicit position of singular points (sec 3.1)

```
c++ -I.. -O3 -DNDEBUG -DKV_FASTROUND verify_singular_points.cc
```

This code provides the explicit position of sigular points that are intersection of singular locus and genelic line on C^2.

### Rigorous inclusion of the fundamental solutions (sec 3.2)

Using interval arithmetic and the truncation error bounds provided in Theorems 3.4, 3.5, 3.6, and 3.7, the rigorous inclusion of each fundamental solution value is obtained via the code:

```
c++ -I.. -O3 -DNDEBUG -DKV_FASTROUND verify_fundamental_sol.cc
```

This code provides the values of fundametal solutions at the base point $(\lambda_0, \mu_0) = (2^{-10},2^{-10})$. It takes more that 5 mins.

### Analytic continuation using *double* precision (sec 4)

The fundamental solution matrix is analytically continued along the loop $\Sigma_1$. Note that this code requires OpenMP API for parallelizing integration of each fundamental solution.

```
c++ -I.. -O3 -DNDEBUG -DKV_FASTROUND find_monodromy_path1.cc -fopenmp
```

This code provides rigorous inclusion of $(\Sigma_1)_\ast\mathrm{Id}$. It took approximately 319 mins ($\approx$ 5 hours) and failed in including the solution under the target magin of error.

### Analytic continuation using *DD* precision (sec 4)

Using DD arithmatic, the fundamental solution matrix is analytically continued precisely.

```
c++ -I.. -O3 -DNDEBUG -DKV_FASTROUND find_monodromy_path1_dd.cc -fopenmp
```

This code provides tight enclosure of $(\Sigma_1)_\ast\mathrm{Id}$. On our workstation (AMD(R) EPYC(TM) 9754 @ 2.25 GHz processors), the computation took approximately 1207 minutes ($\approx$ 20 hours) and successfully enclosed the solution under the target margin of error.


For the other paths, one can execute the files `find_monodromy_pathi_dd.cc`, where i corresponds the index of each path. Note that it took 1316 mins ($\approx$ 22 hours) for path2_alt, 1095 mins ($\approx$ 18.3 hours) for path3, 2613 mins ($\approx$ 43.5 hours) for path4, 801 mins ($\approx$ 13.4 hours) for path5, and 825 mins ($\approx$ 13.8 hours) for path6.


Copyright (C) 2024  T Ishige and A Takayasu.
