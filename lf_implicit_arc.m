% Function to estimate limit load using drop tangent stiffness
% Path traced using arc-length method
function [val,dlam] = lf_implicit_arc(x,data,loop,alpha_target)
% READ SOME DATA
outname = data.outname;
nelx = data.nelx;
nely = data.nely;
exact = data.exact; % 0 = simplified, 1 = exact, 2 = strain interpolation
active = data.freedofs;
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
% data for iterations
TS0 = -1; % initial tangent stiffness
TS = -1e10; % current tangent stiffness
U = zeros(2*(data.nely+1)*(data.nelx+1),1); % initialize disp = 0
U0 = U;
DEL_UC = zeros(length(active),1); % last converged DEL_U
% data for arc-length
lam = 0; % load factor
f0_fact = -1; % scaling factor
del_len = -1; % arc length (squared)
len_min = 0;
% Newton-Raphson loop
n_solve = 0; % count number of linear solves
n_ass = 0; % count number of tangent assemblies / residual evals
tol = 5e-9; % convergence tolerance for residual
steps = 0; % counter for number of steps in the arc-length method
exit = 0; % flag for exit when limit point is found
ittd = 3; % desired number of iterations for convergence (smaller = smaller step in arc-length)
alpha = 1; % drop in stiffness factor
DEL_lamC = 1; % last converged DEL_lam;
U_trial = U; % start from zeros
while exit == 0
  % Start of this increment
  nl_loop = 0;
  nl_res = 1;
  nl_res_old = 1e10;
  n_div = 0; % number of divergent increments
  upmode = 1;
  % use previous step to make first guess for this step
  if steps > 0
      mu = sqrt( del_len / (DEL_UC'*DEL_UC + f0_fact*DEL_lamC*DEL_lamC) );
      DEL_U = mu*DEL_UC;
      U_trial(active) = U(active) + DEL_U;
      DEL_lam = mu*DEL_lamC;
  else
      DEL_lam = 0; % total del load factor
      DEL_U = zeros(length(active),1); % total del disp (only freedofs)
  end
  while nl_res > 0
      [K,~,Fres] = corot_cris2(U_trial,lam+DEL_lam,xPhys,data,exact,1);
      n_ass = n_ass + 1;
      % check convergence, or divergence
      if(nl_loop > 0)
          nl_res = norm(Fres(active),Inf) / ((lam+DEL_lam)*norm(data.F));
          if (nl_res_old-tol < nl_res && nl_res > 0 && nl_loop > 1), n_div = n_div + 1; end % check for divergence
      end
      if nl_res < tol
        % compute current tangent stiffness and alpha
        delU = K(active,active)\data.F(active); % for tangent stiffness
        TS = data.F(active)'*delU; % current tangent stiffness measure
        alpha = TS0/TS;
        % output
        if data.out > 1
            name = sprintf('%s_arc%i_%i.vtk',outname,loop,steps);
            top88line_outvtk(name,data.h,nelx,nely,xPhys,0,U_trial);
        end
        steps = steps+1;
        if data.out > 0
          fprintf('\nlam %12.6e , disp %12.6e, DEL_lam %12.6e, del_len %12.6e, alpha %12.4e, n_ass %i, n_slv %i',...
              lam+DEL_lam,(U_trial'*data.F)/(norm(data.F)),DEL_lam,del_len,alpha,n_ass,n_solve);
        end
        if data.out > 1, eig_lim4(3,U_trial,lam,xPhys,data,loop); end % compute eigs for extra info
        % Update everything (if not past target point)
        if DEL_lam > 0 && alpha > alpha_target
              U = U_trial;
              DEL_UC = DEL_U;
              DEL_lamC = DEL_lam;
              lam = lam + DEL_lam; % save converged point
              % adaptive arc-length
              del_len = del_len * (ittd / nl_loop)^2;
              if del_len < len_min, del_len = 2*len_min; end
        end
        % if we are past a limit point point - go back
        if DEL_lam < -1e-8
            del_len = del_len*0.25; % need to be +ve
        % otherwise, check if we are past target alpha value - and refine
        elseif alpha < alpha_target || DEL_lam < 0
            lam_u = lam; % upper bound
            lam_l = lam+DEL_lam; % lower bound
            lam_t = 0.5*(lam_u+lam_l); % initial guess
            DEL_UC = zeros(length(U),1); DEL_UC(active) = delU;
            itr = 0; % refinement iteration counter
            % refine estimate using bisection method
            while ((abs(lam_u-lam_l) > 1e-4) ...
                    && abs(alpha-alpha_target) > 1e-2) && itr < 10
                DEL_lam = lam_t - lam;
                [U,lam,ns,na] = nr_up(U,lam,DEL_lam,xPhys,data); % N-R to find equilibrium for new trial point
                n_ass = n_ass + na; n_solve = n_solve + ns;
                K = corot_cris2(U,lam,xPhys,data,exact,1);
                DEL_UC(active) = K(active,active)\data.F(active); % for tangent stiffness
                TS = data.F(active)'*DEL_UC(active); % tangent stiffness measure of trial
                alpha = TS0/TS; % alpha for trial
                if data.out > 0
                    fprintf('\nlam %12.6e , disp %12.6e, DEL_lam %12.6e, alpha %12.4e, n_ass %i, n_slv %i',...
                        lam,(U'*data.F)/(norm(data.F)),DEL_lam,alpha,n_ass,n_solve);
                end
                % update bounds
                if alpha < alpha_target, lam_l = lam;
                else, lam_u = lam;
                end
                lam_t = 0.5*(lam_u+lam_l); % new trial point (bi-section method)
                itr = itr+1;
            end
            if data.out > 1, eig_lim4(3,U,lam,xPhys,data,loop,1); end % compute eigs for extra info
            exit = 1; % stop here
        end
        break;
      elseif(n_div > 2)
           % slow or non-convergence, reduce load increment
           if data.out > 1, fprintf('\nCUT with nl_loop=%i, nlres=%12.4e, nlres_old=%12.4e, upmode=%i',...
                   nl_loop,nl_res,nl_res_old,upmode); end
            del_len = del_len*0.25;
            break;
      else
          % if not past target point and not failing, then solve for delU_t & delU_bar
          n_solve = n_solve + 1;
          Udel = K(active,active)\[data.F(active),Fres(active)];
          delU_t = Udel(:,1); delU_bar = Udel(:,2);
          if TS0 < 0, U0(active) = delU_t; TS0 = data.F'*U0; end % initial tangent stiffness
          if f0_fact < 1, f0_fact = full(delU_t'*delU_t); end % axis scaling
          if del_len < 0, del_len = (data.lf_start.value)^2*f0_fact; len_min = 1e-6*del_len; end % first iteration
          % solve for update
          [del_lam,DEL_U_new,quad] = arc_up1(f0_fact,DEL_U,DEL_UC,DEL_lam,del_len,delU_bar,delU_t,data.F(active),Fres(active));
          if quad == 1
            [del_lam, DEL_U_new, quad] = arc_up2(f0_fact,DEL_U,DEL_UC,DEL_lam,del_len,delU_bar,delU_t,data.F(active),Fres(active));
            upmode=2;
              if quad == 1
                 g = full((Fres'*data.F) / (data.F'*data.F));
                 delU_bar = delU_bar - g*delU_t;
                 delU_cr = DEL_U + delU_bar;
                 Ucr = U; Ucr(active) = Ucr(active) + delU_cr;
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
          U_trial(active) = U(active) + DEL_U; % update
          nl_res_old = nl_res;
      end
      nl_loop = nl_loop + 1;
  end
end
data.lf_cr.value = lam;
data.lf_start.value = 0.25*lam; % adative starting load factor (for next opt iter)
val = data.lf_min - lam; % constraint function value
%% COMPUTE GRADIENTS
if nargout > 1
    dlam = nlgeom_sens(xPhys, xFilt, U, 1, 2, DEL_UC, [0,-1,1], 2, -alpha, U0, data);
end
%% OUTPUT
name = sprintf('%s_lf%i.vtk',outname,loop);
if nargout > 1 && data.out > 0
    top88line_outvtk(name,data.h,nelx,nely,xPhys,reshape(dlam,nely,nelx),real(U));
else
    top88line_outvtk(name,data.h,nelx,nely,xPhys,0,real(U));
end
    fprintf('\n It.:%5i Vol.:%7.6f lf.:%12.6e con: %12.12e n_slv.:%i n_ass.:%i res.:%12.4e\n',loop, ...
                mean(xPhys(:)),lam,val,n_solve,n_ass,full(nl_res));
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