function [DY,Dy] = func_Jop_C(FEM,Deltax)

   diagDx = sparse(1:FEM.N,1:FEM.N,Deltax);
   DK = FEM.A1hat*(diagDx*FEM.A1hat') + FEM.A2hat*(diagDx*FEM.A2hat');      
   DY = FEM.Gtilde.'*DK*FEM.Gtilde;
   Dy = DY(FEM.idxM);