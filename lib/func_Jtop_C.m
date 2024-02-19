function [Jtr] = func_Jtop_C(FEM,r)


    R = zeros(size(FEM.Y));
    R(FEM.idxM) = r;
    R = sparse(R);
   
    GtA1=(FEM.Gtilde.'*FEM.A1hat);
    GtA2=(FEM.Gtilde.'*FEM.A2hat);
    
    Jtr = sum((GtA1.*(R*GtA1) + GtA2.*(R*GtA2)));
    
    Jtr = Jtr(:);
    
%     %%%%
%     Jtr1 = Jtr(:);
%     
%     
%     
%     A1tG=(FEM.A1hat'*FEM.Gtilde);
%     GtA2=(FEM.Gtilde'*FEM.A2hat);
%     
%     Jtr2 = sum((GtA1.*(R*GtA1) + GtA2.*(R*GtA2)));
%     
%     Jtr2 = Jtr2(:);