%  Generic function to calculate sensitivies for nonlinear geometric
%  stiffness and stability functions
function sens = nlgeom_sens(xPhys, xFilt, U, D, adjoint, phi, mu, dm, lm, U0, data)
%% read some data  
nelx = data.nelx;
nely = data.nely;
E0 = data.E0;
Emin = data.Emin;
penal = data.penal;
KE = data.KE;
exact = data.exact;
% exact = 1; % hybrid sens
edofMat = data.edofMat;
h = data.h;
H = 0.5*[-h;-h;h;-h;h;h;-h;h]; % centre
% Threshold params
beta = data.beta;
eta = 0.5;
th_nm = tanh(beta*eta);
th_dm = th_nm + tanh(beta*(1-eta));
%% Adjoint solve
[KU,KS,~,Fint,Sen] = corot_cris2(U,1,xPhys,data,exact,2); K = KU+KS;
P = phi; % default
% adjoint based on external forces
if adjoint == 1
    G = data.F;
% adjoint based on dK/dU
elseif adjoint == 2
    np = norm(U)/norm(phi); % ensure U and phi are of same magnitude
    eps = 1e-7*np; spe = 1/eps;
    [dKU,dKS] = corot_cris2(U-eps*phi,1,xPhys,data,exact,2); % FD on dK/dU
    G = spe*(phi'*((KU-dKU) + D*(KS-dKS)))'; % assuming all K's are symmetric
end
if adjoint > 0
    P(data.freedofs) = K(data.freedofs,data.freedofs)'\(G(data.freedofs));
end
%% Sens calculation
sens = zeros(nely*nelx,1);
% dE/dx factors
fact = penal*(E0-Emin)*(beta/th_dm); % constant factor
dE_dx(:) = fact*((sech(beta*(xFilt(:)-eta))).^2).*xPhys(:).^(penal-1);
% Strain energy part
if abs(mu(1)) > 0
    sens = sens + mu(1)*dE_dx(:).*Sen(:);
end
% Internal force part
if abs(mu(2)) > 0
   sens = sens + mu(2)*dE_dx'.*sum(Fint.*P(edofMat),2);
end
% K tangent part
if abs(mu(3)) > 0
   for i=1:nelx*nely
      xdef = U(edofMat(i,:)) + H; % deformed coords in global
      [KTU,KTS,~,~] = Ktan_Cris2(h,xdef,KE);
      sens(i) = sens(i) + mu(3)*dE_dx(i)*(phi(edofMat(i,:))'*(KTU + D*KTS)*phi(edofMat(i,:)));
   end
end
% Linear part
if abs(lm) > 0
    ce = sum((U0(data.edofMat)*data.KE).*U0(data.edofMat),2);
    sens = sens + (lm*dE_dx'.*ce);
end
%% Denominator
if dm == 2
    d = P'*(data.F);
elseif dm == 3
    d = phi(data.freedofs)'*KS(data.freedofs,data.freedofs)*phi(data.freedofs);
end
if dm > 1
    sens = sens/d;
end
%% Filtering
sens = data.Hsens*sens;
end