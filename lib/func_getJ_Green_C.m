function J = func_getJ_Green_C(FEM)


Nelec = length(FEM.Y);

%idx_row = FEM.idxM-kron((0:Nelec:Nelec*(Nelec-1))',ones(Nelec-1,1));
%idx_col = kron((1:Nelec)',ones(Nelec-1,1));

idx_row = FEM.idx_row;
idx_col = FEM.idx_col;

GtA1=(FEM.Gtilde.'*FEM.A1hat); 
GtA2=(FEM.Gtilde.'*FEM.A2hat);

J = (GtA1(idx_row,:).*GtA1(idx_col,:) + GtA2(idx_row,:).*GtA2(idx_col,:));
