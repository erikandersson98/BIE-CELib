# BIE-CELib
A boundary integral equation library implementing several integral operators arising in different problems of the Laplace and Helmholtz equations.

BIE-CELib was developed as part of a master's degree project in numerical analysis at LTH, Sweden. The underlying quadrature is Gauss-Legendre quadrature, and close evaluation (the CE in BIE-CELib) was implemented using product integration and regularization.

The paper "Implementation and study of boundary integral operators related to PDE:s in the plane." summarizes the theory used for the implementation, in sections 2-6. The tests used in the included examples are described in section 7.

The examples which entail solving a boundary value problem are prefaced 'PDE_' in the 'Examples' folder. The tests which involve more direct testing of the integral operators are prefaced 'Op_'. The boundary conditions, wavenumber, screening parameters, etc. can easily be changed, resulting in slightly different boundary value problems. To change the boundary, replace the 'zinit' method in the examples.

If the user wants to use the included operators to implement a different problem, such as an interior Helmholtz problem, then the structure of the examples can be reused. Simply change the construction of the system matrix to reflect the new problem, as well as the formula for field evaluation in the 'FieldComp' method.  

The implementations of the integral operators for target points on the boundary can be found in the 'Operators' folder.

The 'FieldPotentials' folder contains the implementation of the operators for target points not on the boundary, these functions are used in the 'FieldComp' methods in the BVPs examples.

The 'Misc' folder includes the methods for Gauss-Legendre quadrature, a function for determining if a point is to the left or right of the boundary, two functions for determining which panels a point is close to, and one function for the solving of linear systems of equations using a GMRES method.

From the results of the degree project, the HypC_LocRegP and LogC_ProdIntP methods are recommended for use in the solving of BVPs.

Remark: The output of the LogC methods need to be multiplied by awzp, but the same does not apply for the output of the CauC or HypC methods. Internally, the naming convention differentiates between wCcmp, wHcmp and wLcorr, and serves to remind of this difference. The logarithmic term is a 'correction' weight, which needs to be multiplied by awzp, while the Cauchy and Hypersingular terms are 'compensation' weights, which do not need to be multiplied by awzp.
