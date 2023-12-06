function topnlstab(h,nelx,nely,volfrac,penal,rmin,problem,numcon)
% hard-coded benchmark problems
% 1 = cantilever (also the default)
% 2 = snap-though problem
% 3 = compressed column
%% MATERIAL PROPERTIES
if problem == 3
    data.E0 = 1; % E0 x thk
    nu = 0.3;
else
    data.E0 = 0.3e9; % E0 x thk
    nu = 0.4;
end
data.Emin = data.E0*1e-9;
%% PREPARE FINITE ELEMENT ANALYSIS
data.KE = KE_Q6(1,nu); % square bilinear element with incompatible modes (scaled with E=1)
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
data.edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
data.edofMat = repmat(data.edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
data.iK = reshape(kron(data.edofMat,ones(8,1))',64*nelx*nely,1);
data.jK = reshape(kron(data.edofMat,ones(1,8))',64*nelx*nely,1);
% data for parallel assembly of global matrices 
data.ass = zeros(nelx*nely,2);
for i=1:(nelx*nely)
    data.ass(i,1) = ((i-1)*64)+1;
    data.ass(i,2) = i*64;
end
% loading and supports
if problem == 2
    % SNAP-THROUGH PROBLEM
    data.F = sparse((nelx/2)*2*(nely+1)+2,1,-30000,2*(nely+1)*(nelx+1),1);
    fixeddofs = [1:2*(nely+1),2*(nelx)*(nely+1)+1:2*(nelx+1)*(nely+1)];
elseif problem == 3
    % COMPRESSED COLUMN (rotated 90deg clockwise)
    felem = ceil(nelx/60);
    lcDof = 2*nodenrs(nely/2+1+[-felem:felem],end)-1;
    modF = 1e-2/(nelx*h)/(length(lcDof)-1);                         
    data.F = sparse(lcDof,1,-modF,2*(nely+1)*(nelx+1),1);
    [data.F(lcDof(1)),data.F(lcDof(end))] = deal(data.F(lcDof(1))/2,data.F(lcDof(end))/2); 
    fixeddofs = (1:2*(nely+1));
else
    % CANTILEVER LOADS AND SUPPORTS (right-centre load)
    data.F = sparse(2*(nely+1)*(nelx+1)-nely,1,-100000,2*(nely+1)*(nelx+1),1);
    fixeddofs = (1:2*(nely+1));
end
alldofs = (1:2*(nely+1)*(nelx+1));
data.freedofs = setdiff(alldofs,fixeddofs);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
  for j1 = 1:nely
    e1 = (i1-1)*nely+j1;
    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
      for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
        e2 = (i2-1)*nely+j2;
        k = k+1;
        iH(k) = e1;
        jH(k) = e2;
        sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
      end
    end
  end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
data.Hnew = H;
data.Hsens = H;
for i=1:k
    data.Hnew(iH(i),jH(i)) = data.Hnew(iH(i),jH(i)) / Hs(iH(i));
    data.Hsens(iH(i),jH(i)) = data.Hsens(iH(i),jH(i)) / Hs(jH(i));
end
%% PREPARE DATA STRUCTURE
% mesh parameters
data.h = h; % element edge length
data.nelx = nelx; % no. elements in x
data.nely = nely; % no. elements in y
% SIMP parameters
data.penal = penal; % starting penalisation factor
data.beta = 2; % starting density thresold sharpness parameter 
data.exact = 1; % use exact co-rotational method (without strain energy interpolation)
data.sbeta = 500;
data.seta = 1e-2; % beta and eta for strain energy interpolation
% output controls
data.out = 1; % output level, 0=min, 1=normal, 2=full
if problem == 2, data.outname = 'snap'; elseif problem == 3, data.outname = 'column'; else, data.outname = 'cant'; end
% constraint limits
data.lf_min = 2.0; % minimum critical load factor (if constraint used)
if problem==3, volfrac=0; end
data.volfrac = volfrac; % volume fraction constraint
% data recorders
data.comp = libpointer('doublePtr',1); % compliance
data.lf_start = libpointer('doublePtr',1); % starting guess for load factor
data.lf_cr = libpointer('doublePtr',1); % critical load factor
data.Usave = libpointer('doublePtr',zeros(2*(nely+1)*(nelx+1),2)); % save converged U (displacement vector)
data.xeval = libpointer('doublePtr',zeros(nelx*nely,1)); % last design vector evaluated
%% SETUP AND RUN
% set initial design to be just feasible for volume constraint
if problem == 3
    th_fact = 1; % start at fully solid design domain for column problem
else
    th_nm = tanh(data.beta*0.5);
    th_dm = th_nm + tanh(data.beta*(1-0.5));
    th_fact = (1/data.beta)*atanh(volfrac*th_dm - th_nm)+ 0.5;
end
x = th_fact*ones(nely*nelx,1);
% upper and lower design variable limits
xmin = zeros(length(x),1);
xmax = ones(length(x),1);
% MMA settings
if numcon ~= 2, numcon = 1; end
a0_mma = 1; a_mma = zeros(numcon,1);
c_mma = 1e3*ones(numcon,1); d_mma = ones(numcon,1);
xOld = x; xOld1 = xOld ;
low = zeros(length(x),1); upp = zeros(length(x),1);
% scaling factors
Sobj = 1;
Svol = 1/(nelx*nely);
Slf = 1;
% stopping criteria
ftol_rel = 1e-4;
maxeval = 500;
%% MMA loop
loop = 1; odiff = 1; xdiff = 1;
while (loop < maxeval) && (odiff > ftol_rel || conmax > 1e-4 || data.beta < 12)
   % Objective
   % [obj,dg0] = comp_lin(x,data,loop); % linear compliance
   [obj,dg0] = comp_nl4(x,data,loop); % nonlinear end compliance
   if loop == 1
       Sobj = 10/obj; dg0=10*dg0/obj; obj=10; % scale obj value so it is 10 in 1st iteration
   else
       obj=obj*Sobj; dg0=Sobj*dg0;
   end
   % constraints
   [vf,dvf] = volume(x,data); vf=vf*Svol; dvf=dvf*Svol; % volume (with scaling)
   % nonlinear stability constraint
   if numcon == 2
       %[lf,dlf] = lf_arc(x,data,loop); % direct method
       [lf,dlf] = lf_eig_nl4(x,data,loop,1); % eigenvalue anlaysis at gamma=1
       %[lf,dlf] = lf_eig_nl4(x,data,loop,data.lf_min); % eigenvalue anlaysis at gamma=lf_min
       %[lf,dlf] = lf_implicit_arc(x,data,loop,0.5); % drop in tangent stiffness (alpha is last input)
       lf = lf*Slf; dlf=dlf*Slf; % apply scaling
   else
       lf = 0; % no stability constraint
   end
   % monitor predicted vs actual function values
   if loop > 1
       fprintf("Predicted comp = %f, Actual = %f",obj_pred,obj/Sobj);
       fprintf("\nPredicted vol frac = %f, Actual = %f\n",vf_pred+volfrac,vf+volfrac);
       if numcon>1, fprintf("Predicted lf = %f, Actual = %f\n",lf_pred,data.lf_cr.value); end
   end
   % update
   if problem==3
      % swap comp & vol if compressed column problem (obj=vol)
      if numcon == 2
         [x0,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(2,length(x),loop,x,xmin,xmax,xOld,xOld1,vf,dvf,[obj-20;lf],[dg0,dlf]',low,upp,a0_mma,a_mma,c_mma,d_mma);
      else
         [x0,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(1,length(x),loop,x,xmin,xmax,xOld,xOld1,vf,dvf,obj-20,dg0',low,upp,a0_mma,a_mma,c_mma,d_mma);
      end
      conmax=max(obj-20,lf); fprintf("Conmax = %12.4e\n",conmax);
   else
      if numcon == 2
          [x0,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(2,length(x),loop,x,xmin,xmax,xOld,xOld1,obj,dg0,[vf;lf],[dvf,dlf]',low,upp,a0_mma,a_mma,c_mma,d_mma);
      else
          [x0,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(1,length(x),loop,x,xmin,xmax,xOld,xOld1,obj,dg0,vf,dvf',low,upp,a0_mma,a_mma,c_mma,d_mma);
      end
      conmax=max(vf,lf); fprintf("Conmax = %12.4e\n",conmax);
   end
   xOld1 = xOld ; xOld = x; x = x0; % record history
   % make predictions (by linear extrapolation)
   obj_pred = (obj + dg0'*(x-xOld))/Sobj;
   vf_pred = (vf + dvf'*(x-xOld));
   if numcon==2, lf_pred = data.lf_cr.value - (dlf'*(x-xOld))/Slf; end
   % check stopping criteria
   if problem==3, obj=vf; end
   if loop > 1
       odiff = abs((obj-obj_old)/obj_old);
       xdiff = norm(x-xOld,inf);
   end
   fprintf("odiff = %12.4e, xdiff = %12.4e\n",odiff,xdiff);
   loop = loop+1; % update
   obj_old = obj;
   % continuation on threshold projection sharpness & penalty
   %if inner_loop > 149 && odiff < 10*ftol_rel && data.beta < 12          
   if loop > 39 && mod(loop,20)==0 && data.beta < 12
        data.beta = data.beta + 2; % increase threshold sharpness
        fprintf("\n__Increasing beta to %f",data.beta);
   end
   if loop > 39 && mod(loop,20)==0 && data.penal < 3
        data.penal = data.penal + 0.25; % increase penalisation
   end
end
fprintf('\n\n! Optimisation complete - entering post-proc\n');
%% POST PROCESSING
data.beta = 12; % reset for post-proc
loop = 2000;
val = comp_nl4(xOld,data,loop); % solve to save U
loop = loop + 1;
fprintf('\n\n__Load factor using path tracking (lf_arc)\n');
data.lf_start.value = 1;
val = lf_arc(xOld,data,loop);
for i=0:1
    gamT = 1 + i*1;
    fprintf('\n\n__Load factor using eigen-analysis at lf=%f\n',gamT);
    loop = loop + 1;
    val = lf_eig_nl4(xOld,data,loop,gamT);
end
for i=0:2
    alpha = 0.1 + i*0.4;
    data.lf_start.value = 1;
    fprintf('\n\n__Load factor using drop in stiffness at alpha=%f\n',alpha);
    loop = loop + 1;
    val = lf_implicit_arc(xOld,data,loop,alpha);
end
% END
end
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