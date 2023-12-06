% Function to perfrom a arc-length update step
% for a prescribed arc-length
function [U,lam,n_solve,n_ass,DEL_UC,DEL_lamC] = arc_up(U,load_fac_in,load_inc_in,DEL_UC,DEL_lamC,xPhys,data)
% read some data
exact = data.exact;
lam = load_fac_in;
lam_max = load_fac_in + load_inc_in;
% starting info (depening on if we are at zero, or not)
if norm(U) > 1e-10
    f0_fact = U'*U / (load_fac_in)^2; % scaling factor
    target_len = 2*load_inc_in*sqrt(f0_fact); % target arc-length depends on load_inc
    del_len = (0.1*target_len)^2;
    len_min = 1e-8*del_len;
    steps = 1;
% we are at zero, so determine values in 1st iteration
else
    f0_fact = -1;
    target_len = 1;
    del_len = 1;
    len_min = 0;
    steps = 0;
    DEL_lam = load_inc_in;
end
U_trial = zeros(length(U),1);
ittd = 3; % target no. iteration per load increment
n_solve = 0; n_ass = 0; % counters
total_len = 0;
tol = 5e-9; % convergence tolerance
% start iterations
while (del_len > len_min && lam < lam_max && abs(DEL_lamC) > tol && total_len < target_len)
  % Start of this increment
  nl_loop = 0;
  nl_res = 1;
  nl_res_old = 1e10;
  n_div = 0; % number of divergent increments
  upmode = 1;
   if steps > 0
      mu = sqrt( del_len / (DEL_UC'*DEL_UC + f0_fact*DEL_lamC*DEL_lamC) );
      DEL_U = mu*DEL_UC;
      U_trial(data.freedofs) = U(data.freedofs) + DEL_U;
      DEL_lam = mu*DEL_lamC;
   end
   % start iterations for this increment
   while nl_res > 0
      [K,~,Fres] = corot_cris2(U_trial,lam+DEL_lam,xPhys,data,exact,1);
      n_ass = n_ass + 1;
      % check convergence, or divergence
      if(nl_loop > 0)
          nl_res = norm(Fres(data.freedofs),Inf) / ((lam+DEL_lam)*norm(data.F));
          if (nl_res_old-tol < nl_res && nl_loop > 1), n_div = n_div + 1; end
      end
      if nl_res < tol
          if data.out > 0
              fprintf('\nlam %12.6e, disp %12.6e, DEL_lam %12.6e, del_len %12.6e, n_ass %i, n_slv %i',...
                    lam+DEL_lam,(U_trial'*data.F)/(norm(data.F)),DEL_lam,sqrt(del_len),n_ass,n_solve);
          end
          % Update everything
          steps = steps+1;
          U = U_trial;
          DEL_UC = DEL_U;
          lam = lam + DEL_lam; % save converged point
          DEL_lamC = DEL_lam;
          total_len = total_len + sqrt(del_len);
          % adaptive arc-length
          del_len = del_len * (ittd / nl_loop)^2;
          break;
      % if diverged too many times, reduce del_len and try again
      elseif(n_div > 2)
           if data.out > 1, fprintf('\nCUT with nl_loop=%i, nlres=%12.4e, nlres_old=%12.4e, upmode=%i',...
                   nl_loop,nl_res,nl_res_old,upmode); end
            del_len = del_len*0.25;
            break;
      else
          % if not converged and not failing, then solve for delU_t & delU_bar
          n_solve = n_solve + 1;
          Udel = K(data.freedofs,data.freedofs)\[data.F(data.freedofs),Fres(data.freedofs)];
          delU_t = Udel(:,1); delU_bar = Udel(:,2);
          % set arc-length parameters, if needed
          if f0_fact < 0
               f0_fact = U'*U / (load_fac_in)^2; % scaling factor
               target_len = 2*load_inc_in*sqrt(f0_fact); % target arc-length depends on load_inc
               del_len = (0.1*target_len)^2;
               len_min = 1e-8*del_len;
          end
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
    DOT1 = (DEL_U1'*DEL_UC); % + f0_fact*DEL_lam*(DEL_lam + del_lam1 - g);
    DOT2 = (DEL_U2'*DEL_UC); % + f0_fact*DEL_lam*(DEL_lam + del_lam2 - g);
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
    DOT1 = (DEL_U1'*DEL_UC); %+ f0_fact*DEL_lam*(DEL_lam + del_lam1 - g);
    DOT2 = (DEL_U2'*DEL_UC); % + f0_fact*DEL_lam*(DEL_lam + del_lam2 - g);
    if DOT2 > DOT1
      DEL_U = DEL_U2;
      del_lam = del_lam2 - g;
    else
      DEL_U = DEL_U1;
      del_lam = del_lam1 - g;
    end
end