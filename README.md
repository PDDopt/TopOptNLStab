# TopOptNLStab
Matlab code for geometrically nonlinear topology optimization, including stability constraints.

This code solves 2D topology optimization problems using a three-field density method, including geometrically nonlinear effects and stability (buckling) constraints, and can be used to run the examples published in [1].

Geometrically noninear modelling is achieved using the co-rotational method [2], and some formulations utilise the arc-length method suggested in [3]. The code utilises some of the framework from the top88 matlab code [4], and some ideas from papers exploring geometrically nonlinear topology optimization [5 - 10].

#Running the code
To run the code, you will also need a copy of the Matlab implmentation of MMA, which can be downloaded from: https://www.smoptit.se

The main file is: topnlstab.m. Note that various options are hard-coded and you need to mannually change this file to run different optimization problems.

To change from linear to nonlinear end compliance - comment/uncomment lines 129-130

To change the noninear stability constraint formulation - comment/uncomment lines 140-143

To run the various example in [1], use the following commands:
Cantilever problem in sections 5: topnlstab(0.025,160,40,0.4,3,3,1,1) & set objective to linear compliance (line 129)
Snap-through problem in section 5: topnlstab(0.0125,240,80,0.1,3,2,2,1) & set objective to linear compliance (line 129)
Cantilever problem in section 6.1: topnlstab(0.025,160,40,0.4,3,3,1,2) & set objective to nonlinear end compliance (line 130), change load on line 43
Column problem in section 6.2: topnlstab(0.00833333,240,120,0.35,3,4,3,2) & set objective to nonlinear end compliance (line 130)

# References
1. Dunning, P.D. Stability constraints for geometrically nonlinear topology optimization. Structural and Multidisciplinary Optimization 66, 253 (2023). https://doi.org/10.1007/s00158-023-03712-8

2. Crisfield, M., Moita, G. A co-rotational formulation for 2-d continua including incompatible modes. International Journal for Numerical Methods in Engineering 39(15), 2619–2633 (1996)

3. Lam, W., Morley, C. Arc-length method for passing limit points in structural calculation. Journal of Structural Engineering 118(1), 169–185 (1992)

4. Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B.S., & Sigmund, O. Efficient topology optimization in MATLAB using 88 lines of code. Structural and Multidisciplinary Optimization 43, 1-16 (2011)
   
6. Wang, F., Lazarov, B.S., Sigmund, O., Jensen, J.S. Interpolation scheme for fictitious domain techniques and topology optimization of finite strain elastic problems. Computer Methods in Applied Mechanics and Engineering 276, 453–472 (2014)

7. Dalklint, A.,Wallin, M., Tortorelli, D.A. Structural stability and artificial buckling modes in topology optimization. Structural and Multidisciplinary Optimization 64(4), 1751–1763 (2021)

8. Pedersen, N.L., Pedersen, P. Buckling load optimization for 2d continuum models, with alternative formulation for buckling load estimation. Structural and Multidisciplinary Optimization 58, 2163–2172 (2018)

9. Kemmler, R., Lipka, A., Ramm, E. Large deformations and stability in topology optimization. Structural and Multidisciplinary Optimization 30(6), 459–476 (2005)

10. Lindgaard, E., Dahl, J. On compliance and buckling objective functions in topology optimization of snap-through problems. Structural and Multidisciplinary Optimization 47(3), 409–421 (2013)
