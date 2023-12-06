%% C-shape benchmark
clc; clear all;
%% MESH
h = 1;
nelx = 10;
nely = 10;
nthk = 1;
%% MATERIAL PROPERTIES
data.E0 = 1; % E0 x thk
data.Emin = data.E0*1e-9;
nu = 0.3;
data.penal = 3;
%% PREPARE FINITE ELEMENT ANALYSIS
data.KE = KE_Q6(1,nu);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
data.edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
data.edofMat = repmat(data.edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
data.iK = reshape(kron(data.edofMat,ones(8,1))',64*nelx*nely,1);
data.jK = reshape(kron(data.edofMat,ones(1,8))',64*nelx*nely,1);
data.ass = zeros(nelx*nely,2);
for i=1:(nelx*nely)
    data.ass(i,1) = ((i-1)*64)+1;
    data.ass(i,2) = i*64;
end
%% LOADING & BCs
data.F = sparse([2*(nely)*(nelx+1)+2,2*(nely+1)*(nelx+1)-1],[1,1],9*[-0.002,0.003],2*(nely+1)*(nelx+1),1);
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
data.freedofs = setdiff(alldofs,fixeddofs);
%% ELEMENT SETS
solid = 1:(nely*nthk); % left edge
for i=1:nthk, solid=union(solid,union(i:nely:nely*nelx,nely-i+1:nely:nely*nelx)); end % top & bot
%% PHYSICAL DENSITIES
x = zeros(nelx*nely,1);
x(solid) = 1;
%% RUN ARC-LENGTH
data.nelx = nelx;
data.nely = nely;
data.h = h;
data.exact = 1; % 1 = without strain energy interpolation, 2 = with strain energy interpolation
data.sbeta = 500;
data.seta = 1e-2; % beta and eta for strain energy interpolation
data.out = 2;
data.lf_start = libpointer('doublePtr',1); % starting guess for load factor
%% SOLVE & OUTPUT
[U,~] = combo_solve(0,1,x,zeros(length(alldofs),1),data);
top88line_outvtk('c_shape.vtk',h,nelx,nely,x,0,U);
% END
%% Stiffness matrix for square 4-node bilinear plane stress elements with incompatible modes
function K = KE_Q6(E,v)
g = (1-v)/2;
Efact = E/(1-v*v);
K=Efact*...
[ - v^2/12 + g/4 + 1/3,            g/4 + v/4,   v^2/12 + g/4 - 1/3,            v/4 - g/4, - v^2/12 - g/4 - 1/6,          - g/4 - v/4,   v^2/12 - g/4 + 1/6,            g/4 - v/4;...
            g/4 + v/4, - v^2/12 + g/4 + 1/3,            g/4 - v/4,   v^2/12 - g/4 + 1/6,          - g/4 - v/4, - v^2/12 - g/4 - 1/6,            v/4 - g/4,   v^2/12 + g/4 - 1/3;...
   v^2/12 + g/4 - 1/3,            g/4 - v/4, - v^2/12 + g/4 + 1/3,          - g/4 - v/4,   v^2/12 - g/4 + 1/6,            v/4 - g/4, - v^2/12 - g/4 - 1/6,            g/4 + v/4;...
            v/4 - g/4,   v^2/12 - g/4 + 1/6,          - g/4 - v/4, - v^2/12 + g/4 + 1/3,            g/4 - v/4,   v^2/12 + g/4 - 1/3,            g/4 + v/4, - v^2/12 - g/4 - 1/6;...
 - v^2/12 - g/4 - 1/6,          - g/4 - v/4,   v^2/12 - g/4 + 1/6,            g/4 - v/4, - v^2/12 + g/4 + 1/3,            g/4 + v/4,   v^2/12 + g/4 - 1/3,            v/4 - g/4;...
          - g/4 - v/4, - v^2/12 - g/4 - 1/6,            v/4 - g/4,   v^2/12 + g/4 - 1/3,            g/4 + v/4, - v^2/12 + g/4 + 1/3,            g/4 - v/4,   v^2/12 - g/4 + 1/6;...
   v^2/12 - g/4 + 1/6,            v/4 - g/4, - v^2/12 - g/4 - 1/6,            g/4 + v/4,   v^2/12 + g/4 - 1/3,            g/4 - v/4, - v^2/12 + g/4 + 1/3,          - g/4 - v/4;...
            g/4 - v/4,   v^2/12 + g/4 - 1/3,            g/4 + v/4, - v^2/12 - g/4 - 1/6,            v/4 - g/4,   v^2/12 - g/4 + 1/6,          - g/4 - v/4, - v^2/12 + g/4 + 1/3];
 end