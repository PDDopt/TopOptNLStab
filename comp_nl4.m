% Nonlinear End Compliance obj fucntion
% Equilibrium using combined N-R + arc-length
function [val, dc] = comp_nl4(x,data,loop)
% READ SOME DATA
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
% Solve Equilibrium, to load factor 1
if loop > 20, U = data.Usave.value(:,1); % use U from last iteration as starting guess
else, U = zeros(2*(nely+1)*(nelx+1),1);
end
lam = min(1,data.lf_start.value);
if norm(x - data.xeval.value,inf) > 1e-10
    [U,lam] = combo_solve(0,1,xPhys,U,data);
    data.Usave.value(:,1) = U; % save converged U
    data.xeval.value = x;
end
val = lam*(U'*data.F); % objective
%% Sensitivity analysis
if nargout > 1
    dc = lam * nlgeom_sens(xPhys, xFilt, U, 1, 1, U, [0,-1,0], 1, 0, 0, data);
end
%% OUTPUT
name = sprintf('%s_%i.vtk',outname,loop);
if nargout > 1 && data.out > 0
    top88line_outvtk(name,data.h,nelx,nely,xPhys,reshape(dc,nely,nelx),U);
else
    top88line_outvtk(name,data.h,nelx,nely,xPhys,0,U);
end
if data.out > 0, colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow; end
fprintf(' It.:%5i Obj.:%12.12e Vol.:%7.6f lf.:%7.2f \n',loop,val,mean(xPhys(:)),lam);
end