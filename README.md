# Sparse unmixing via variable splitting and augmented Lagrangian methods (SUNSAL)

## Description

SUNSAL solves the following l2-l1 optimization  problem
[size(M) = (L,p); size(X) = (p,N)]; size(Y) = (L,N)]

	min_X  (1/2) ||M X-y||^2_F + lambda ||X||_1

	where ||X||\_1 = sum_i sum_j |X_{i,j}|.

CONSTRAINTS ACCEPTED:

1) POSITIVITY:  X >= 0;
2) ADDONE:  X_j = 1 for all j;

NOTES:
  1) The optimization w.r.t each column of X is decoupled. Thus, SUNSAL solves N simultaneous problems.

  2) SUNSAL solves the following  problems:

     a) BPDN - Basis pursuit denoising l2-l1
               (lambda > 0, POSITIVITY = False, ADDONE, False)

     b) CBPDN - Constrained basis pursuit denoising l2-l1
               (lambda > 0, POSITIVITY = True, ADDONE, False)

     c) CLS   - Constrained least squares
                (lambda = 0, POSITIVITY = True, ADDONE, False)

     c) FCLS   - Fully constrained least squares
                (lambda >=0 , POSITIVITY = True, ADDONE, True)
                 In this case, the regularizer ||X||_1  plays no role,
                 as it is constant.

## Line of Attack

SUNSAL solves the above optimization problem by introducing a variable splitting and then solving the resulting constrained optimization with the augmented Lagrangian method of multipliers (ADMM).

	min_{X,Z}  (1/2) ||M X-y||^2_F + lambda ||Z||_1

	subject to: X_j = 1 for all j; Z >= 0; X = Z

Augmented Lagrangian (scaled version):

    L(X,Z,D) = (1/2) ||M X-y||^2_F + lambda ||Z||_1 + mu/2||X-Z-D||^2_F

where D are the scale Lagrange multipliers

ADMM:

    do
        X  <-- arg min L(X,Z,D)
                   X, s.t: sum(X) = ones(1,N));
        Z  <-- arg min L(X,Z,D)
                   Z, s.t: Z >= 0;
        D  <-- D - (X-Z);
    while ~stop_rule

More details on the method:

J. Bioucas-Dias and M. Figueiredo, "Alternating direction algorithms for constrained sparse regression: Application to hyperspectral unmixing",
in 2nd IEEE GRSS Workshop on Hyperspectral Image and SignalProcessing-WHISPERS'2010, Raykjavik, Iceland, 2010.

# Usage
x,res_p,res_d,i = sunsal(M,y,AL_iters=1000,lambda_0=0.,positivity=False,addone=False,tol=1e-4,x0 = None,verbose=False)

## Required inputs

M - [L(channels) x p(endmembers)] endmembers matrix

y - pixels matrix with  L(channels) x N(pixels).
    each pixel is a linear mixture of p endmembers signatures y = M*x + noise,

## Optional inputs

AL_ITERS - Minimum number of augmented Lagrangian iterations - Default: 1000

lambda_0 - regularization parameter. lambda is either a scalar or a vector with N components (one per column of x) - Default: 0.

positivity  = {True, False}; Enforces the positivity constraint: X >= 0 - Default: False

addone  = {True, False}; Enforces the positivity constraint: X >= 0 - Default: False

tol    - tolerance for the primal and  dual residuals - Default: 1e-4;

verbose   = {True, False}; Default: False

## Output variables

x      - estimated mixing matrix [pxN]
res_p  - primal residual
res_d  - dual residual
i      - number of iteration until convergence 

# Requirements

Scipy needs to be installed.

# Authors

Software translated from matlab to python by Adrien Lagrange (ad.lagrange@gmail.com), 2018.

Initial matlab author: Jose Bioucas-Dias, 2009

# License

SUNSAL is distributed under the terms of the GNU General Public License 2.0.
