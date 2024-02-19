function J = func_getJ_AVM_C(FEM)

    diagx        = sparse(1:FEM.N,1:FEM.N,FEM.x);
    Khat         = FEM.Kinihat + FEM.A1hat*diagx*FEM.A1hat' + FEM.A2hat*diagx*FEM.A2hat';

    J = zeros(length(FEM.y),FEM.N);
    
    V     =  Khat\FEM.F;
    gamma =  Khat'\FEM.M'; % this additional solve is required for the AVM
    
    for ii = 1:FEM.N
       % temp = -gamma'*((FEM.A1hat(:,ii)*FEM.A1hat(:,ii)' + FEM.A2hat(:,ii)*FEM.A2hat(:,ii)')*V);
        
        temp = -(gamma'*FEM.A1hat(:,ii))*(FEM.A1hat(:,ii)'*V) + -(gamma'*FEM.A2hat(:,ii))*(FEM.A2hat(:,ii)'*V);
        
        J(:,ii) = temp(FEM.idxM);
    end