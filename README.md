# TopOptNLStab
Matlab code for geometrically nonlinear topology optimization, including stability constraints.

This code solves 2D topology optimization problems using a three-field density method, including geometrically nonlinear effects and stability (buckling) constraints, and can be used to run the examples published in [1].

Geometrically nonlinear modelling is achieved using the co-rotational method [2], and some formulations utilise the arc-length method suggested in [3]. The code utilises some of the framework from the top88 matlab code [4], and some ideas from papers exploring geometrically nonlinear topology optimization [5 - 9].

# Disclaimer
Although the code has been extensively tested, I cannot guarantee that it is free from errors. No claim is made about the efficiency of this code, or its compactness (number of lines). Some files use parfor for speedup, but the code should run without this, by replacing the 'parfor', with 'for'.

# Running the code
To run the code, you will also need a copy of the Matlab implementation of MMA, which can be downloaded from: https://www.smoptit.se

The main file is: topnlstab.m. Note that various options are hard-coded and you need to manually change this file to run different optimization problems.

* To change from linear to nonlinear end compliance - comment/uncomment lines 129-130.

* To change the nonlinear stability constraint formulation - comment/uncomment lines 140-143.

To run the various example in [1], use the following commands:

* Cantilever in section 5: ```topnlstab(0.025,160,40,0.4,3,3,1,1)``` & set objective to linear compliance (line 129).

* Snap-through in section 5: ```topnlstab(0.0125,240,80,0.1,3,2,2,1)``` & set objective to linear compliance (line 129).

* Cantilever in section 6.1: ```topnlstab(0.025,160,40,0.4,3,3,1,2)``` & set objective to nonlinear end compliance (line 130), change load on line 43.

* Column in section 6.2: ```topnlstab(0.00833333,240,120,0.35,3,4,3,2)``` & set objective to nonlinear end compliance (line 130).

There is also a script (cshape.m) that runs the c-shape benchmark problem [10], as seen in section 3.3 of [1].

# References
1. Dunning, P.D. Stability constraints for geometrically nonlinear topology optimization. Structural and Multidisciplinary Optimization 66, 253 (2023). https://doi.org/10.1007/s00158-023-03712-8

2. Crisfield, M., Moita, G. A co-rotational formulation for 2-d continua including incompatible modes. International Journal for Numerical Methods in Engineering 39(15), 2619–2633 (1996). https://doi.org/10.1002/(SICI)1097-0207(19960815)39:15%3C2619::AID-NME969%3E3.0.CO;2-N

3. Lam, W., Morley, C. Arc-length method for passing limit points in structural calculation. Journal of Structural Engineering 118(1), 169–185 (1992). https://doi.org/10.1061/(ASCE)0733-9445(1992)118:1(169)

4. Andreassen, E., Clausen, A., Schevenels, M., Lazarov, B.S., & Sigmund, O. Efficient topology optimization in MATLAB using 88 lines of code. Structural and Multidisciplinary Optimization 43, 1-16 (2011). https://doi.org/10.1007/s00158-010-0594-7
   
5. Wang, F., Lazarov, B.S., Sigmund, O., Jensen, J.S. Interpolation scheme for fictitious domain techniques and topology optimization of finite strain elastic problems. Computer Methods in Applied Mechanics and Engineering 276, 453–472 (2014). https://doi.org/10.1016/j.cma.2014.03.021

6. Dalklint, A.,Wallin, M., Tortorelli, D.A. Structural stability and artificial buckling modes in topology optimization. Structural and Multidisciplinary Optimization 64(4), 1751–1763 (2021). https://doi.org/10.1007/s00158-021-03012-z

7. Pedersen, N.L., Pedersen, P. Buckling load optimization for 2d continuum models, with alternative formulation for buckling load estimation. Structural and Multidisciplinary Optimization 58, 2163–2172 (2018). https://doi.org/10.1007/s00158-018-2030-3

8. Kemmler, R., Lipka, A., Ramm, E. Large deformations and stability in topology optimization. Structural and Multidisciplinary Optimization 30(6), 459–476 (2005). https://doi.org/10.1007/s00158-005-0534-0

9. Lindgaard, E., Dahl, J. On compliance and buckling objective functions in topology optimization of snap-through problems. Structural and Multidisciplinary Optimization 47(3), 409–421 (2013). https://doi.org/10.1007/s00158-012-0832-2
    
10. Yoon, G.H., Kim, Y.Y. Element connectivity parameterization for topology optimization of geometrically nonlinear structures. International Journal of Solids and Structures 42(7), 1983–2009 (2005). https://doi.org/10.1016/j.ijsolstr.2004.09.005
