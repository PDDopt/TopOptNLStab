%% Co-rotational global stiffness matrix and residual force
%  return KT is split into 2 parts
function [KU,KS,Fres,Fint,Sen] = corot_cris2(U,lam,xPhys,data,exact,mode)
    % read data
    nelx = data.nelx;
    nely = data.nely;
    h = data.h;
    E0 = data.E0;
    Emin = data.Emin;
    penal = data.penal;
    edofMat = data.edofMat;
    KE = data.KE;
    % strain interpolation method
    sbeta = data.sbeta; seta = data.seta;
    numer =  tanh(sbeta*seta);
    denom =  numer + tanh(sbeta*(1-seta));
    % setup useful data structures
    H = 0.5*[-h;-h;h;-h;h;h;-h;h]; % centre
    Fres = full(lam*data.F); % initialise residual force vector
    sKU_ass = zeros(64,nelx*nely);
    sKS_ass = zeros(64,nelx*nely);
    compF = 0; compS = 0;
    Fint = zeros(nelx*nely,8,1);
    Sen = zeros(nelx*nely,1);
    if nargout > 2, compF = 1; end % check if we need residual
    if nargout == 5, compS = 1; end % check if we need strain energy
    Ee = (Emin+xPhys.^penal*(E0-Emin)); % Young's modulus
    Ue = U(edofMat(:,:))'; % element global displacement vector
    %% assemble global tangent stiffness and internal forces (in parallel)
    parfor i=1:(nelx*nely)
      % strain interpolation
      gam = 1;
      if exact==2
          gam = (numer + tanh(sbeta*(xPhys(i)^penal-seta)))/denom;
      end
      xdef = gam*Ue(:,i) + H; % deformed coords in global
      % Assemble element tangent stiffness into global
      if compF
         if compS, [KTU,KTS,Fint(i,:),Sen(i)]= Ktan_Cris2(h,xdef,KE);
         else, [KTU,KTS,Fint(i,:)]= Ktan_Cris2(h,xdef,KE); end
         if exact==2
            FLin = KE*Ue(:,i);
            Fint(i,:) = gam*Fint(i,:)' + (1-gam*gam)*FLin;
            if compS, Sen(i) = Sen(i) + 0.5*(1-gam*gam)*(FLin'*Ue(:,i)); end
        end
      else
          [KTU,KTS] = Ktan_Cris2(h,xdef,KE);
      end
      if exact==2, KTU = gam*gam*KTU + (1-gam*gam)*KE; KTS = gam*gam*KTS; end
      sKU_ass(:,i) = Ee(i)*reshape(KTU,64,1);
      sKS_ass(:,i) = Ee(i)*reshape(KTS,64,1);
    end
    %% post parallel summing
    if mode == 2 % return split matrix (i.e. for eigenvalue calculations)
        KU = sparse(data.iK,data.jK,reshape(sKU_ass,nelx*nely*64,1)); KU = (KU+KU')/2; % last command ensures K is symmetric
        KS = sparse(data.iK,data.jK,reshape(sKS_ass,nelx*nely*64,1)); KS = (KS+KS')/2; % last command ensures K is symmetric
    else % else, just return KT = KU + KS
        KU = sparse(data.iK,data.jK,reshape(sKU_ass+sKS_ass,nelx*nely*64,1)); KU = (KU+KU')/2; % last command ensures K is symmetric
        KS = 0;
    end
    % for fint
    if compF
      for i=1:(nelx*nely)
          Fres(edofMat(i,:)) = Fres(edofMat(i,:)) - Ee(i)*Fint(i,:)'; 
      end
    end
end