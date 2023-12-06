% Function to perfrom a N-R update step
% to reach a target load factor
function [U,load_fac,n_solve,n_ass,dU,dlf] = nr_up(U,load_fac_in,load_inc_in,xPhys,data)
% read some data
exact = data.exact;
load_target = load_fac_in + load_inc_in; % target load factor
load_fac = load_fac_in;
load_inc = load_inc_in;
active = data.freedofs;
% controls
tol = 5e-9;
ittd = 5; % target no. iteration per load increment
relax=1; % under-relaxation factor
sn = sign(load_inc); % check which direction we are going
dU=0; dlf=0;
n_solve = 0; n_ass = 0; % counters
% start the loop
while (sn*load_fac < sn*load_target) && (abs(load_inc) > 1e-2*abs(load_inc_in))
      load_inc = sn*min(abs(load_target-load_fac),abs(load_inc));
      load_fac = load_fac + load_inc;
      U_trial = U; % starting disp for this increment
      nl_loop = 0;
      nl_res = 1;
      nl_res_old = 1;
      while nl_res > 0
          [K,~,Fres] = corot_cris2(U_trial,load_fac,xPhys,data,exact,1);
          n_ass = n_ass + 1;
          % check convergence, or divergence
          if(nl_loop > 0),  nl_res = norm(Fres(active),Inf) / (load_fac*norm(data.F)); end
          if nl_res < tol
              fprintf('\nlam %12.6e , disp %12.6e, DEL_lam %12.6e, n_ass %i, n_slv %i',...
                  load_fac,(U_trial'*data.F)/(norm(data.F)),load_inc,n_ass,n_solve);
              dU = U_trial(active) - U(active); % converged disp increment
              dlf = load_inc; % converged load increment
              load_inc = load_inc * (ittd / nl_loop); % adjust load increment
              U = U_trial;
              break;
          elseif(nl_res_old < nl_res && nl_loop > 1)
              relax = 0.25*relax;
              % slow or non-convergence, reduce load increment
              if relax < 0.05
                  load_fac = load_fac - load_inc; % reset
                  load_inc = load_inc*0.25;
                  relax = 1;
                  break;
              end
              % Otherwise, reduce under relaxation factor and try again
              U_trial(active) = U_trial(active) - 0.75*delU;
              delU = 0.25*delU;
          else
              % if not converged and not failing, then solve
              n_solve = n_solve + 1;
              delU = K(active,active)\Fres(active);
              U_trial(active) = U_trial(active) + delU; % displacement increment
              nl_res_old = nl_res;
              relax = 1;
          end
          nl_loop = nl_loop + 1;
      end
end
% end loop
end