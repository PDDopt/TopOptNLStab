%% Crisfield Corot 2D quad element
%  return KT is split into 2 parts
function [KTU,KTS,fint,sen] = Ktan_Cris2(h,xdef,KE)
H = 0.5*[-h;-h;h;-h;h;h;-h;h]; % centre
% work out rotation angle
c = (0.5/h)*[-1,-1,1,-1,1,1,-1,1];
d = (0.5/h)*[-1,1,-1,-1,1,-1,1,1];
a = c*xdef; b = d*xdef;
if abs(b) < 1e-10, phi = 0;
elseif abs(a) < 1e-10, phi = pi()/ 2;
else phi = atan(-b/a); end
% adjust for rotations beyond +/- Pi/2
if a < -1e-10
    if b < -1e-10, phi = phi + pi();
    else, phi = phi - pi(); end
end
% rotation matrix
cp = cos(phi); sp = sin(phi);
R = [cp,sp;-sp,cp]; R = mybd(R);
dR = [-sp,cp;-cp,-sp]; dR = mybd(dR);
uloc = R*xdef - H;
% internal forces & strain energy
floc = KE*uloc;
if nargout > 2, fint = R'*floc; end
if nargout == 4, sen = 0.5*floc'*uloc; end
%% tangent stiffness matrix
% simple part (dependent on U)
KTU = R'*KE*R;
% stress stiffness part
KTS = zeros(8,8);
denom = (a*a + b*b);
if denom > 1e-10
    w = ( b*c - a*d ) / denom; % derivative of phi
    N = dR'*floc + R'*KE*dR*xdef;
    for i=1:8, KTS(:,i) = N*w(i); end
end
end
% custom blkdiag (more efficient)
function R = mybd(R1)
    R=zeros(8,8);
    R(1:2,1:2) = R1;
    R(3:4,3:4) = R1;
    R(5:6,5:6) = R1;
    R(7:8,7:8) = R1;
end