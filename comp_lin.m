%% LINEAR COMPLIANCE
function [val, dc] = comp_lin(x,data,loop)
% READ SOME DATA
nelx = data.nelx;
nely = data.nely;
h = data.h;
E0 = data.E0;
Emin = data.Emin;
penal = data.penal;
outname = data.outname;
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
sK = reshape(data.KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1); % modified SIMP
K = sparse(data.iK,data.jK,sK); K = (K+K')/2;
U(data.freedofs) = K(data.freedofs,data.freedofs)\data.F(data.freedofs);
% Compute objective (end compliance)
val = U'*data.F;
data.comp.value = val;
%% COMPUTE GRADIENTS
if nargout > 1
    ce = reshape(sum((U(data.edofMat)*data.KE).*U(data.edofMat),2),nely,nelx);
    dc=zeros(nely,nelx);
    fact = -penal*(E0-Emin);
    parfor i=1:nely
      for j=1:nelx
          tanhf = tanh(beta*(xFilt(i,j)-eta)); sech2f = 1-tanhf*tanhf;
          dc(i,j) = (beta/th_dm)*sech2f*fact*xPhys(i,j)^(penal-1)*ce(i,j);
      end
    end
    dc = data.Hsens*dc(:); % filter
end
 %% OUTPUT
 name = sprintf('%s_%i.vtk',outname,loop);
 if nargout > 1 && data.out > 0
     top88line_outvtk(name,h,nelx,nely,xPhys,reshape(dc,nely,nelx),U);
 else
     top88line_outvtk(name,h,nelx,nely,xPhys,0,U);
 end
 if data.out > 0, colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal; axis off; drawnow; end
 fprintf('\n It.:%5i Comp lin.:%12.6e Vol.:%7.12f\n',loop,val,mean(xPhys(:)));
end