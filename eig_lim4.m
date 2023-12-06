% Eigenvalue solve at a converged point
function [D,phi] = eig_lim4(num_eig,U,lam,xPhys,data,loop,out)
    % read data
    nelx = data.nelx;
    nely = data.nely;
    active = data.freedofs;
    num_eig_slv = 2*num_eig + 1; % solve for more, in case we gte -ve values
    [KU,KS] = corot_cris2(U,1,xPhys,data,2,2); % always use strain energy interpolation to avoid false modes
    dK = decomposition(KU(active,active),'chol','lower');
    matFun = @(x) dK\(KS(active,active)*x); % matrix action function
    phi = zeros(2*(nely+1)*(nelx+1),num_eig_slv); % eigenvectors
    goteigs=0;
    while goteigs == 0
      [V,D,fg] = eigs(matFun,length(active),num_eig_slv,'sa');  
      [D,I] = sort(-1./diag(D)); % ensure eigs are in ascending order
      % screen out negative eigenvalues (if necessary)
      cnt=0;
      for ee=1:num_eig_slv
          if real(D(ee)) > 0
              cnt = cnt+1;
              D(cnt) = real(D(ee));
              if cnt <= num_eig
                  sn = sign(data.F(active)'*V(:,I(ee))); % ensure consistent sign
                  phi(active,cnt) = sn*(V(:,I(ee)));
              end
          end
      end
      % check we have enough positive eigenvalues
      if cnt < num_eig, num_eig_slv = num_eig_slv + 2*(num_eig + 1);
      else, num_eig_slv = cnt; goteigs = 1;
      end
    end
    % output
    if out > 0
        fprintf('\n%i Non-neg eigenvalues at LF %f (%i) = ',num_eig_slv,lam,fg);
        for ee=1:num_eig_slv
          fprintf(' %12.4e',D(ee));
          if out > 0 && ee <= min(num_eig,3)
             name = sprintf('%s_eig%i_%i.vtk',data.outname,ee,loop);
             top88line_outvtk(name,data.h,nelx,nely,xPhys,0,real(phi(:,ee)));
          end
        end   
    end
end