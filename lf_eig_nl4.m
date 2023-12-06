% Nonlinear Limit load - using eig analysis at a converged point (lamT)
% Equilibrium using combined N-R + arc-length
function [val, dlam] = lf_eig_nl4(x,data,loop,lamT)
%% READ SOME DATA
outname = data.outname;
nelx = data.nelx;
nely = data.nely;
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
% Solve Equilibrium, to load factor lamT (if needed)
U = data.Usave.value(:,1); lam = 1; % start point
if lamT > lam, [U,lam] = combo_solve(lam,lamT,xPhys,U,data); end
%% Solve eignevalue problem
num_eig = 6;
[D,phi] = eig_lim4(num_eig,U,lam,xPhys,data,loop,data.out);
% processing for multiple eigs
if num_eig > 1
    pn = 50; % p for KS-function
    muk = zeros(num_eig,1);
    for i=1:num_eig, muk(i) = 1/(D(i)); end
    [ks_val,ksens] = KS_fun(muk,pn); % K-S
    for i=1:num_eig, ksens(i) = ksens(i)*(muk(i)/ks_val)^2; end
    eig_cr = 1/ks_val;    
else
    eig_cr = D(1); ksens(1) = 1;
end
val = (data.lf_min / lam) - eig_cr;
data.lf_cr.value = eig_cr*lam; % approx limit load factor
%% Sensitivity analysis
if nargout > 1
    dlam = zeros(nelx*nely,1);
    % sum for each eigenvalue
    for i=1:num_eig
        dlam = dlam + nlgeom_sens(xPhys, xFilt, U, D(i), 2, phi(:,i), ksens(i)*[0,-1,1], 3, 0, 0, data);
    end
end
%% OUTPUT
name = sprintf('%s_lam_%i.vtk',outname,loop);
if nargout > 1 && data.out > 0
    top88line_outvtk(name,data.h,nelx,nely,xPhys,reshape(dlam,nely,nelx),U);
else
    top88line_outvtk(name,data.h,nelx,nely,xPhys,0,U);
end
colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow;
fprintf('\n It.:%5i Con.:%12.4e Vol.:%7.6f lf.:%7.2f \n',loop,val,mean(xPhys(:)),lam);
end
%K-S aggregation function
function [val,sens] = KS_fun(x,p)
    val = 1;
    sens = zeros(length(x),1); sens(1) = 1;
    for i=2:length(x)
        sens(i) = exp(p*(x(i)-x(1)));
        val = val + sens(i);
    end
    % sensitivity factors
    for i=1:length(x)
        sens(i) = sens(i) / val;
    end
    val = x(1) + log(val)/p;
end
