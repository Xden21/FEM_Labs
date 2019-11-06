Function{
  // BH data EIcore
  MatEIcore_h = {0.,
    2.1827e+02,3.8272e+02,5.2036e+02,7.1167e+02,
    8.4921e+02,1.0405e+03,1.4741e+03,2.0409e+03,
    3.0109e+03,5.0572e+03,7.5335e+03,1.0037e+04,
    1.2486e+04,1.5015e+04,1.7464e+04,2.0021e+04,
    2.2040e+04,2.5000e+04} ;

  MatEIcore_b = {0.,
    2.0329e-01,4.0287e-01,6.0986e-01,8.0575e-01,
    1.0053e+00,1.1975e+00,1.4008e+00,1.5154e+00,
    1.5967e+00,1.6706e+00,1.7076e+00,1.7335e+00,
    1.7446e+00,1.7593e+00,1.7667e+00,1.7815e+00,
    1.7889e+00,1.7963e+00} ;

  MatEIcore_b2 = List[MatEIcore_b]^2 ;
  MatEIcore_nu = List[MatEIcore_h]/List[MatEIcore_b] ;
  MatEIcore_nu(0) = MatEIcore_nu(1);

  MatEIcore_nu_b2  = ListAlt[MatEIcore_b2, MatEIcore_nu] ;
  nu_EIcore[] = InterpolationLinear[SquNorm[$1]]{List[MatEIcore_nu_b2]} ;
  dnudb2_EIcore[] = dInterpolationLinear[SquNorm[$1]]{List[MatEIcore_nu_b2]} ;
  h_EIcore[] = nu_EIcore[$1] * $1 ;
  // diff. reluctivity tensor
  dhdb_EIcore[] = TensorDiag[1,1,1]*nu_EIcore[$1#1] + 2*dnudb2_EIcore[#1] * SquDyadicProduct[#1] ;
  // non-linear part of diff. reluctivity tensor
  dhdb_EIcore_NL[] = 2*dnudb2_EIcore[$1] * SquDyadicProduct[$1] ;

  // Alternatively, we could use an analytical law.
  // e.g. analytical brauer law
  k1= 940;
  k2= 0.023;
  k3= 4.1;
  nu_brauer[] = k1 + k2 * Exp[k3*SquNorm[$1]] ;  // scalar reluctivity
  dnudb2_brauer[] = k2*k3 * Exp[k3*SquNorm[$1]] ; // derivative w.r.t. to the square norm of the real b-vector

  // the 2 functions used for the Newton-Raphson scheme, real b-vector as input
  h_brauer[] = nu_brauer[$1]*$1 ;
  dhdb_brauer[] = TensorDiag[1,1,1] * nu_brauer[$1] + 2*dnudb2_brauer[$1] * SquDyadicProduct[$1]  ; // diff. reluctivity tensor
  dhdb_brauer_NL[] = 2*dnudb2_brauer[$1] * SquDyadicProduct[$1]  ; // diff. reluctivity tensor

}
