%% NONLINEAR GEOMETRY LOAD FACTOR - USING ARC-LENGTH METHOD
function [val, dlam] = lf_arc(x,data,loop)
  % READ SOME DATA
  outname = data.outname;
  nelx = data.nelx;
  nely = data.nely;
  exact = data.exact; % 0 = simplified, 1 = exact, 2 = with strain energy interpolation
  % Threshold params
  beta = data.beta;
  eta = 0.5;
  th_nm = tanh(beta*eta);
  th_dm = th_nm + tanh(beta*(1-eta));
  % UPDATE xPhys (filter & threshold)
  xFilt = reshape(x,nely,nelx);
  xFilt(:) = data.Hnew*xFilt(:); % filter
  xPhys = zeros(nely,nelx);
  xPhys(:) = (th_nm + tanh(beta*(xFilt(:)-eta)))/th_dm; % threshold
  %% FE-ANALYSIS
  U = zeros(2*(nely+1)*(nelx+1),1); % initialize disp = 0
  % Arc-length loop
  lam = 0; % starting load factor
  DEL_lam = 1; % starting load increment
  f0_fact = -1; % load scaling factor
  del_len = -1; % arc length (squared) - auto set in 1st iteration
  len_min = 0; % auto set in 1st iteration
  n_solve = 0; % count number of linear solves
  n_ass = 0; % count number of tangent assemblies / residual evals
  tol = 5e-9; % convergence tolerance for residual
  lim_tol = 1e-8; % tolerance for identificaiton of limit point
  steps = 0; % counter for number of steps in the arc-length method
  exit = 0; % flag for exit when limit point is found
  ittd = 3; % desired number of iterations for convergence (smaller = smaller arc-length steps)
  DEL_UC = zeros(length(data.freedofs),1); %last converged DEL_U
  lim_passed=0;
  while exit == 0
      % Start of this increment
      U_trial = U; % start from last converged disp
      nl_loop = 0;
      nl_res = 1;
      nl_res_old = 1e10;
      n_div = 0; % number of divergent increments
      upmode = 1;
      if steps > 0
          mu = sqrt( del_len / (DEL_UC'*DEL_UC + f0_fact*DEL_lam*DEL_lam) );
          DEL_U = mu*DEL_UC;
          U_trial(data.freedofs) = U(data.freedofs) + DEL_U;
          DEL_lam = mu*DEL_lam;
      else
          DEL_lam = 0; % total del load factor
          DEL_U = zeros(length(data.freedofs),1); % total del disp (only freedofs)
      end
      while nl_res > 0
          [K,~,Fres] = corot_cris2(U_trial,lam+DEL_lam,xPhys,data,exact,1);
          n_ass = n_ass + 1;
          % check convergence, or divergence
          if(nl_loop > 0)
              nl_res = norm(Fres(data.freedofs),Inf) / ((lam+DEL_lam)*norm(data.F));
              if (nl_res_old-tol < nl_res && nl_loop > 1), n_div = n_div + 1; end
          end
          if nl_res < tol
              % output
              if data.out > 1
                  name = sprintf('%s_arc%i_%i.vtk',outname,loop,steps);
                  top88line_outvtk(name,data.h,nelx,nely,xPhys,0,U_trial);
              end
              DK = eig_lim4(1,U_trial,lam+DEL_lam,xPhys,data,loop,0); % get 1st eigenvalue
              steps = steps+1;
              if data.out > 0
              fprintf('\nlam %12.6e , disp %12.6e, DEL_lam %12.6e, del_len %12.6e, eig %12.4e, n_ass %i, n_slv %i',...
                  lam+DEL_lam,(U_trial'*data.F)/(norm(data.F)),DEL_lam,del_len,DK(1),n_ass,n_solve);
              end
              % Update everything (if not past limit / birfurcation point)
              if DEL_lam > 0 && DK(1) > 1
                  U = U_trial;
                  DEL_UC = DEL_U;
                  lam = lam + DEL_lam; % save converged point
                  % adaptive arc-length
                  del_len = del_len * (ittd / nl_loop)^2;
                  if del_len < len_min, del_len = 2*len_min; end
                  if DEL_lam < lim_tol, lim_passed = 1; end % if change in lam small, assume we are close to a limit point
              else
                  % if DEL_lam < lim_tol - assume we are near limit point
                  if DEL_lam < lim_tol
                      % if too far past - then go back (get closer to limit point)
                      if DEL_lam < -1e-3, del_len = del_len*0.0625;
                      else, lim_passed = 1;
                      end
                  else
                      lim_passed = 2; % if DEL_lam > lim_tol, assume birfircation point
                  end
              end
              % Check if we have passed (or are near) a limit / birfircation point
              if lam > 0 && lim_passed > 0
                  [D,V] = eig_lim4(1,U,lam,xPhys,data,loop,data.out);
                  exit = 1; % exit  
              end
              break;
          elseif(n_div > 2)
               % slow or non-convergence, reduce load increment
               if data.out > 1, fprintf('\nCUT with nl_loop=%i, nlres=%12.4e, nlres_old=%12.4e, upmode=%i',...
                    nl_loop,nl_res,nl_res_old,upmode); end
                    del_len = del_len*0.25;
               break;
          else
              % if not converged and not failing, then solve for delU_t & delU_bar
              n_solve = n_solve + 1;
              Udel = K(data.freedofs,data.freedofs)\[data.F(data.freedofs),Fres(data.freedofs)];
              delU_t = Udel(:,1); delU_bar = Udel(:,2);
              if f0_fact < 1, f0_fact = full(delU_t'*delU_t); end % axis scaling
              if del_len < 0, del_len = (data.lf_start.value)^2*f0_fact; len_min = 1e-6*del_len; end % first iteration
              % solve for update
              [del_lam,DEL_U_new,quad] = arc_up1(f0_fact,DEL_U,DEL_UC,DEL_lam,del_len,delU_bar,delU_t,data.F(data.freedofs),Fres(data.freedofs));
              if quad == 1
                [del_lam, DEL_U_new, quad] = arc_up2(f0_fact,DEL_U,DEL_UC,DEL_lam,del_len,delU_bar,delU_t,data.F(data.freedofs),Fres(data.freedofs));
                upmode=2;
                  if quad == 1
                     g = full((Fres'*data.F) / (data.F'*data.F));
                     delU_bar = delU_bar - g*delU_t;
                     delU_cr = DEL_U + delU_bar;
                     Ucr = U; Ucr(data.freedofs) = Ucr(data.freedofs) + delU_cr;
                     [~,~,F_cr] = corot_cris2(Ucr,0,xPhys,data,exact,1);
                     lam_cr = (-F_cr'*data.F) / (data.F'*data.F);
                     del_lam_cr = lam_cr - lam;
                     mu = sqrt(del_len)/sqrt(delU_cr'*delU_cr + f0_fact*del_lam_cr*del_lam_cr);
                     del_lam = mu*del_lam_cr - DEL_lam;
                     DEL_U_new = mu*delU_cr;
                     upmode=3;
                  end
              end
              DEL_U = DEL_U_new;
              DEL_lam = DEL_lam + del_lam;
              U_trial(data.freedofs) = U(data.freedofs) + DEL_U; % update
              nl_res_old = nl_res;
          end
          nl_loop = nl_loop + 1;
      end
  end
lam = lam*D(1,1);
data.lf_start.value = 0.25*lam; % adative starting load factor (for next opt iter)
val = data.lf_min - lam; % constraint value
data.lf_cr.value = lam; % save critical value
%% COMPUTE GRADIENTS
if nargout > 1
    if lim_passed == 1, dlam = nlgeom_sens(xPhys, xFilt, U, 1, 0, V(:,1), [0,-1,0], 2, 0, 0, data); % for limit point
    else, dlam = nlgeom_sens(xPhys, xFilt, U, 1, 2, V(:,1), [0,-1,1], 2, 0, 0, data); % for birfurcation
    end
end
 %% OUTPUT
 name = sprintf('%s_lf%i.vtk',outname,loop);
 if nargout > 1 && data.out > 0
    top88line_outvtk(name,data.h,nelx,nely,xPhys,reshape(dlam,nely,nelx),real(U));
 else
    top88line_outvtk(name,data.h,nelx,nely,xPhys,0,real(U));
 end
 if lim_passed == 1, type = 'limit';
 else, type = 'bifur';
 end
 fprintf('\n It.:%5i Vol.:%7.6f lf.:%12.6e con: %12.12e type: %s, n_slv.:%i n_ass.:%i res.:%12.4e\n',loop, ...
     mean(xPhys(:)),lam,val,type,n_solve,n_ass,full(nl_res));
end
%% Arc-updated functions
function [del_lam, DEL_U, quad] = arc_up1(f0_fact,DEL_U,DEL_UC,DEL_lam,del_len,delU_bar,delU_t,f0,res)
    % Solve quadratic for update
    g = full((res'*f0) / (f0'*f0));
    delU_bar = delU_bar - g*delU_t;
    DEL_lamg = DEL_lam - g;
    delUU = DEL_U + delU_bar;
    aa = delU_t'*delU_t + f0_fact;
    bb = 2*(delU_t'*delUU + DEL_lamg*f0_fact);
    cc = delUU'*delUU + DEL_lamg*DEL_lamg*f0_fact - del_len;
    % check for complex roots
    check = bb*bb - 4*aa*cc;
    if check < 0
        quad = 1; del_lam = 0; DEL_U = DEL_U;
        return;
    else
        quad = 0;
    end
    % 2 solutions
    qf1 = -bb/(2*aa);
    qf2 = sqrt(check)/(2*aa);
    del_lam1 = qf1 + real(qf2);
    del_lam2 = qf1 - real(qf2);
    % choose best solution
    DEL_U1 = DEL_U + delU_bar + del_lam1*delU_t;
    DEL_U2 = DEL_U + delU_bar + del_lam2*delU_t;
    DOT1 = (DEL_U1'*DEL_UC);
    DOT2 = (DEL_U2'*DEL_UC); 
    if DOT2 > DOT1
      DEL_U = DEL_U2;
      del_lam = del_lam2 - g;
    else
      DEL_U = DEL_U1;
      del_lam = del_lam1 - g;
    end
end
% modified
function [del_lam, DEL_U, quad] = arc_up2(f0_fact,DEL_U,DEL_UC,DEL_lam,del_len,delU_bar,delU_t,f0,res)
    % factors and update
    g = full((res'*f0) / (f0'*f0));
    delU_bar = delU_bar - g*delU_t;
    % solve for eta
    utt = full(delU_t'*delU_t);
    eta_fact = f0_fact + utt;
    utb = full(delU_t'*delU_bar);
    ubb = full(delU_bar'*delU_bar);
    dut = full(DEL_U'*delU_t);
    dub = full(DEL_U'*delU_bar);
    duu = full(DEL_U'*DEL_U);
    eta_fact2 = DEL_lam - g;
    % eta equn
    aa = eta_fact*ubb - utb*utb;
    bb = 2*(eta_fact*dub - (f0_fact*eta_fact2 + dut)*utb);
    cc =  eta_fact*(duu - del_len) + f0_fact*utt*eta_fact2*eta_fact2 ...
            - 2*dut*f0_fact*eta_fact2 - dut*dut;
    % eta solve
    if abs(aa) > 0 && abs(bb) > 0
        ef1 = -bb/(2*aa);
        ef2 = sqrt(bb*bb - 4*aa*cc)/(2*aa);
        % make eta2 >= eta1
        if ef2 > 0
            eta1 = ef1 - ef2; eta2 = ef1 + ef2;
        else
            eta1 = ef1 + ef2; eta2 = ef1 - ef2;
        end
        % choose appropriate eta
        zeta = 0.01*abs(eta2 - eta1);
        if eta2 < 1, eta = eta2 - zeta;
        elseif -bb/aa < 1, eta = eta2 + zeta;
        elseif eta1 < 1, eta = eta1 - zeta;
        else, eta = eta1 + zeta;
        end
    elseif abs(bb) > 0, eta = -cc/bb;
    else, eta = 1;
    end
    % eta check (should be between 0 & 1)
    if eta > 1, eta = 1;
    elseif eta < 0, eta = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % solve quadratic for del_lam
    aa = utt + f0_fact;
    bb = 2*(dut + eta*utb + eta_fact2*f0_fact);
    cc = duu + 2*eta*dub + eta*eta*ubb + eta_fact2*eta_fact2*f0_fact - del_len;
    % 2 solutions
    qf1 = -bb/(2*aa);
    qf2 = sqrt(bb*bb - 4*aa*cc)/(2*aa);
    % check for complex roots
    if (bb*bb - 4*aa*cc) < 0
        quad = 1; del_lam = 0; DEL_U = DEL_U;
        return;
    else
        quad = 0;
    end
    del_lam1 = qf1 + real(qf2);
    del_lam2 = qf1 - real(qf2);
    % choose best solution
    DEL_U1 = DEL_U + eta*delU_bar + del_lam1*delU_t;
    DEL_U2 = DEL_U + eta*delU_bar + del_lam2*delU_t;
    DOT1 = (DEL_U1'*DEL_UC); 
    DOT2 = (DEL_U2'*DEL_UC); 
    if DOT2 > DOT1
      DEL_U = DEL_U2;
      del_lam = del_lam2 - g;
    else
      DEL_U = DEL_U1;
      del_lam = del_lam1 - g;
    end
end