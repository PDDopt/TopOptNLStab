%% VOLUME (in units of no. elements)
function [val, dv] = volume(x,data)
 % READ SOME DATA
  nelx = data.nelx;
  nely = data.nely;
  % Threshold params
  beta = data.beta;
  eta = 0.5;
  th_nm = tanh(beta*eta);
  th_dm = th_nm + tanh(beta*(1-eta));
  % COMPUTE VOLUME CONSTRAINT VIOLATION
  % filter and threshold
  xFilt = reshape(x,nely,nelx);
  xFilt(:) = data.Hnew*xFilt(:); % filter
  xPhys = zeros(nely,nelx);
  xPhys(:) = (th_nm + tanh(beta*(xFilt(:)-eta)))/th_dm; % threshold
  val = sum(xPhys(:)) - nelx*nely*data.volfrac;
  % COMPUTE GRADIENTS
  if nargout > 1
      dv=zeros(nely,nelx);
      parfor i=1:nely
          for j=1:nelx
              tanhf = tanh(beta*(xFilt(i,j)-eta)); sech2f = 1-tanhf*tanhf;
              dv(i,j) = (beta/th_dm)*sech2f; % including threshold
          end
      end
      % filter sens
      dv = data.Hsens*dv(:);
  end
end