function   func_FEMplot(FEM)

Node =  FEM.Node;
Inzid =  FEM.Inzid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Material = real(FEM.x);
FEM.r_pipe = 0.052;
temp=linspace(0,2*pi,10000);
circx = FEM.r_pipe*cos(temp);
circy = FEM.r_pipe*sin(temp);


figure; hold on, set(gca,'FontSize',26);set(gcf,'Color','White');
set(gcf,'Renderer','OpenGL')

for ii = 1:FEM.N
    jj = FEM.idxROI(ii);
    idx = FEM.Inzid(jj,:);
    X= FEM.Node(idx,1);
    Y= FEM.Node(idx,2);
    fill(X,Y,(Material(ii)),'EdgeColor','none')   
end 
plot(circx,circy,'-k','linewidth',3)
axis off
axis equal
drawnow
title('Conductivity in S/m')
colormap gray
oldcmap = colormap;
colormap( flipud(oldcmap) );
colorbar
caxis([0 0.01]);

Material = imag(FEM.x)/(FEM.eps0*abs(FEM.s));
%temp=linspace(0,2*pi,10000);
% circx = FEM.r_pipe*cos(temp);
% circy = FEM.r_pipe*sin(temp);


figure; hold on, set(gca,'FontSize',26);set(gcf,'Color','White');
set(gcf,'Renderer','OpenGL')

for ii = 1:FEM.N
    jj = FEM.idxROI(ii);
    idx = FEM.Inzid(jj,:);
    X= FEM.Node(idx,1);
    Y= FEM.Node(idx,2);
    fill(X,Y,(Material(ii)),'EdgeColor','none')   
end 
plot(circx,circy,'-k','linewidth',3)
axis off
axis equal
drawnow
title('Relative Permittivity')
colormap gray
oldcmap = colormap;
colormap( flipud(oldcmap) );
colorbar
caxis([0 12]);



