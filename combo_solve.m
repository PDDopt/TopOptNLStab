% Function to find converged equilibrium at target load factor
% using combo of N-R and arc-length
function [U,load_fac] = combo_solve(load_fac_in,load_fac_tar,xPhys,U,data)
n_solve = 0; n_ass = 0; % counters
exit = 0;
load_fac = load_fac_in;
%U = zeros(length(data.F),1);
while abs(load_fac - load_fac_tar) > 1e-6 && exit < 3
    % try N-R
    [U,load_fac,ns,na,dU,dlf] = nr_up(U,load_fac,load_fac_tar-load_fac,xPhys,data);
    n_solve = n_solve + ns; n_ass = n_ass + na; 
    % check if arc-length needed
    if load_fac < load_fac_tar && exit < 2
       fprintf('\nSwitching to arc-length at lf %f',load_fac);
       [U,load_fac,ns,na] = arc_up(U,load_fac,load_fac_tar-load_fac,dU,dlf,xPhys,data);
       n_solve = n_solve + ns; n_ass = n_ass + na;
    end
    exit = exit+1;
end
fprintf('\nConverged to lf %f, with %i solves, %i assembles\n',load_fac,n_solve,n_ass);
end