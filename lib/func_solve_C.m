function [FEM] = func_solve_C(FEM)


C          = sparse(1:FEM.N,1:FEM.N,FEM.x);
Khat       = FEM.Kinihat + FEM.A1hat*C*FEM.A1hat' + FEM.A2hat*C*FEM.A2hat';
[L U]      = lu(Khat(FEM.perm,FEM.perm));
FEM.Gtilde = U \ (L\FEM.Ftilde(FEM.perm,:));
FEM.Gtilde = FEM.Gtilde(FEM.iperm,:);
FEM.Y      = FEM.Gtilde.'*FEM.F;
FEM.y      = FEM.Y(FEM.idxM);

   