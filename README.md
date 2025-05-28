# oblique-shock-calculator
MATLAB code for computing oblique shock properties based on upstream Mach number and deflection angle.

% OBLIQUE_SHOCK Compute oblique shock properties
 sol = oblique_shock(M1, theta_deg) returns a struct with flow
 properties after an oblique shock.


% Inputs:
 M1    - Upstream Mach number (>1)
 theta - Flow deflection angle (deg)


% Output:
  sol - Struct with shock properties
