% Code for the OIPE Paper
% 
% Precomputation steps for an exemplary Complex ECT sensor
%
% EMS 2024
% Contact: neumayer@tugraz.at
clear all, close all, clc

addpath .\lib

load MESH

EPS_tube   = 3;         % relative permittivity of tube
EPS_air    = 1;         % air; also used for backside
omega      = 2*pi*10E6; % curcular frequency, e.g 10 MHz
s          = 1i*omega; 
eps0       = 8.854E-12;

Node =  MESH.Node;
Inzid =  MESH.Inzid;
Material = MESH.Material;
Nelec = MESH.NrElec;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of FE  matrices by a Gaussian quadrature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Evaluation points and weights for Gaussian quadrature
    xis  = [0 1 0 1/2 1/2 0 1/3]; 
    etas = [0 0 1 0 1/2 1/2 1/3];
    ws  = [1/40 1/40 1/40 1/15 1/15 1/15 27/120];

    % Evaluate FE basis functions and derivatives for Gaussian quadrature points
    Nall = zeros(3,length(ws));
    dNxiTall = zeros(3,length(ws));
    dNetaTall = zeros(3,length(ws));
    for ii = 1:length(xis)
        Nall(:,ii) = func_Nilin(xis(ii),etas(ii));           % Linear basis functions
        [ dNxiT, dNetaT ] = func_dNilin(xis(ii),etas(ii));   % Derivatives
        dNxiTall(:,ii) = dNxiT;
        dNetaTall(:,ii) = dNetaT; 
    end

    %Allocate storage
    Ke = zeros(3,3,size(Inzid,1));                    % FE element matrices
    A1 = sparse(1,1,0,size(Node,1),  size(Inzid,1));  % Eigenvectors of element matrices
    A2 = sparse(1,1,0,size(Node,1),  size(Inzid,1));

    %Loop over finite elements
    for jj = 1:size(Inzid,1)
        Kel = zeros(3);  
        inzid = Inzid(jj,:);
        xnode = Node(inzid,1);
        ynode = Node(inzid,2);  
        %Gaussian quadrature
        for ii = 1:length(xis)              
            dNxi=dNxiTall(:,ii);
            dNeta=dNetaTall(:,ii);      
            J = [dNxi'*xnode dNxi'*ynode; dNeta'*xnode dNeta'*ynode];  % Jacobian
            dNxy =   J\[dNxi,dNeta]'; 
            detJw = det(J)*ws(ii);     
         
            Kel    = Kel +  dNxy' * dNxy * detJw;     % Laplace Equation
        end
  
        Ke(:,:,jj) = Kel;
  
        [V,D] = eig(Kel);
        D = diag(D);
        [~,idx]= min(abs(D));
        idxv = 1:3; idxv(idx)=[];
        a1 = sqrt(D(idxv(1)))*V(:,idxv(1));
        a2 = sqrt(D(idxv(2)))*V(:,idxv(2));
        A1(inzid,jj) = a1;
        A2(inzid,jj) = a2;  
      
    end

      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of FE  matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    idx1 = find(Material == 50);  % Tube
    idx2 = find(Material == 150); % Backside

    % Material vector
    mat        = s*eps0*ones(size(Inzid,1),1); 
    mat(idx1)  = s*eps0*EPS_tube; 
    mat(idx2)  = s*eps0*EPS_air;
    
    % Material vector to form Kinihat
    mat0       = zeros(size(Inzid,1),1); 
    mat0(idx1) = s*eps0*EPS_tube; 
    mat0(idx2) = s*eps0*EPS_air;

    % Assembly FE-Stiffness matrix without Boundary conditions, e.g. Eq.(2.12)
    C = sparse(1:length(mat),1:length(mat),mat);
    Kfull = A1*C*A1' + A2*C*A2';

    F = sparse(1,1,0,size(Node,1),Nelec);
    M = F';
    
    idxShield =  find(MESH.BndMrk==100);            

    idxElec = [];
    for ii = 1 : Nelec  
        idx = find(MESH.BndMrk==ii);    
        F(:,ii) = -sum(Kfull(:,idx),2);
        F(idx,ii) = 1;
        idxElec = [idxElec;idx];
        M(ii,:) =  sum(Kfull(idx,:),1);   
    end

    [~,idxGreen] = find(M);            %Non zero elements of M
    idxGreen = unique(idxGreen);       %For these indizes the Green's functions have to be computed
    Mtilde = M(:,idxGreen);
    
    idxDirich = [idxShield; idxElec];     % Nodes with Dirichlet boundary conditions
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly E_Y and Ftilde
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    E_Y       = sparse(idxGreen, 1:length(idxGreen), ones(length(idxGreen),1), size(Kfull,1), length(idxGreen));
    Ftilde    = E_Y * Mtilde.';    % Here we neet transpose (.'); not the hermitian 
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assembly Kinihat and Aihat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    diagx = sparse(1:length(mat0),1:length(mat0),mat0);
   
    Kinihat = A1*diagx*A1' + A2*diagx*A2';  % FE stiffness matrix w/o boundary conditions

    Kinihat(idxDirich,:) = 0;  Kinihat(:,idxDirich) = 0;   
    Kinihat = Kinihat + sparse(idxDirich, idxDirich, ones(size(idxDirich)), size(Kinihat,1), size(Kinihat,1));

    
    idxROI = find(mat0 == 0);
    
    A1hat = A1(:,idxROI);    
    A2hat = A2(:,idxROI); 
    
    idxROI = find(mat0 == 0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Permutation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
   
    diagx = speye(length(idxROI),length(idxROI));

    Khat = Kinihat + A1hat*diagx*A1hat' + A2hat*diagx*A2hat';
    %perm = amd(Khat);
    perm = symamd(Khat);
    [~,iperm] = sort(perm);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2nd Order Reg Matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    index_in_matval = find(MESH.Material == 0);
    L = fem_make_reg_matrix_1st_order(MESH.Inzid,index_in_matval,MESH.Node);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Store all relevant matrices in FEM Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    FEM.N       = length(idxROI);
    
    FEM.Kinihat = Kinihat;
    FEM.A1hat   = A1hat;
    FEM.A2hat   = A2hat;
    FEM.M       = M; 
    FEM.Mtilde  = Mtilde; 
    FEM.F       = F; 
    FEM.Ftilde  = Ftilde; 
    FEM.E_Y     = E_Y; 
    
    FEM.perm    = perm;
    FEM.iperm   = iperm;
    FEM.idxROI  = idxROI;    
    
    % Indices of elements in Y used for measurements
    FEM.idxM =  find(ones(Nelec)-eye(Nelec)); % Measurements between electrodes
    
    
    FEM.idx_row = FEM.idxM-kron((0:Nelec:Nelec*(Nelec-1))',ones(Nelec-1,1));
    FEM.idx_col = kron((1:Nelec)',ones(Nelec-1,1));
    
    
    % Additional element, e.g. for visualization
    FEM.Node  = Node;
    FEM.Inzid = Inzid;
    FEM.L     = L;
    
    FEM.s     = s;
    FEM.eps0  = eps0;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward map implementation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
   sigma = 0.05;
   epsr  = 1;


   FEM.x  =   sigma*ones(FEM.N,1) +   FEM.s*FEM.eps0*epsr*ones(FEM.N,1);  % Material vector
   
   
   % Solution with scalar potential
   diagx        = sparse(1:FEM.N,1:FEM.N,FEM.x);
   Khat         = FEM.Kinihat + FEM.A1hat*diagx*FEM.A1hat' + FEM.A2hat*diagx*FEM.A2hat';
   [L U]        = lu(Khat(FEM.perm,FEM.perm));
   V            =  U \ (L \FEM.F(FEM.perm,:));
   V            = V(FEM.iperm,:);    
   Y1           = FEM.M*V; 
   FEM.y        = Y1(FEM.idxM) ;
    
   % Solution with Green's functions
   diagx        = sparse(1:FEM.N,1:FEM.N,FEM.x);
   Khat         = FEM.Kinihat + FEM.A1hat*diagx*FEM.A1hat' + FEM.A2hat*diagx*FEM.A2hat';
   [L U]        = lu(Khat(FEM.perm,FEM.perm));
   Gtilde       =  U \ (L \FEM.Ftilde(FEM.perm,:));
   Gtilde       = Gtilde(FEM.iperm,:);        
   Y2           = Gtilde.'*FEM.F;     % Here we neet transpose (.'); not the hermitian    
   FEM.Gtilde   = Gtilde;
   FEM.y        = Y2(FEM.idxM) ;
   
   

   save SENSOR_FEM_C FEM
   
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

figure; hold on, set(gcf,'Color','White'), set(gca,'FontSize',16);   
    spy(ones(Nelec)-eye(Nelec))
    title('Elements of Y used as measurements')
    xlabel('Active electrode')
    ylabel('Measurement')
  
    
figure; hold on, set(gcf,'Color','White'), set(gca,'FontSize',16);   
    subplot(2,1,1), hold on, set(gca,'FontSize',16);   
        stem(FEM.idxM, real(Y1(FEM.idxM)),'LineWidth',2)
        stem(FEM.idxM, real(Y2(FEM.idxM)),'--','LineWidth',2)
        xlabel('Index of measurement')
        ylabel('Real y')
        legend('Standard computation','Greens functions')
        
    subplot(2,1,2), hold on, set(gca,'FontSize',16);   
        stem(FEM.idxM, real(Y1(FEM.idxM)- Y2(FEM.idxM)),'LineWidth',2)
        xlabel('Index of measurement')
        legend('Difference')
    
figure; hold on, set(gcf,'Color','White'), set(gca,'FontSize',16);   
    subplot(2,1,1), hold on, set(gca,'FontSize',16);   
        stem(FEM.idxM, imag(Y1(FEM.idxM)),'LineWidth',2)
        stem(FEM.idxM, imag(Y2(FEM.idxM)),'--','LineWidth',2)
        xlabel('Index of measurement')
        ylabel('Imag y')
        legend('Standard computation','Greens functions')
        
    subplot(2,1,2), hold on, set(gca,'FontSize',16);   
        stem(FEM.idxM, imag(Y1(FEM.idxM)- Y2(FEM.idxM)),'LineWidth',2)
        xlabel('Index of measurement')
        legend('Difference')        
       
    
figure; hold on, set(gcf,'Color','White'), set(gca,'FontSize',16);
     for ii=1:size(Inzid,1)
        XY = Node(Inzid(ii,:),:)  ;
        XY = [XY;XY(1,:)];
        Ve = V(Inzid(ii,:),1);  % Electrode 1
        Ve = abs([Ve;Ve(1)]);
        fill3(XY(:,1),XY(:,2),Ve,Ve);
     end
    axis square
    axis off
    title('Scalar potential')

figure; hold on, set(gcf,'Color','White'), set(gca,'FontSize',16);
     for ii=1:size(Inzid,1)
        XY = Node(Inzid(ii,:),:)  ;
        XY = [XY;XY(1,:)];
        Ve = Gtilde(Inzid(ii,:),1);  % Electrode 1
        Ve = abs([Ve;Ve(1)]);
        fill3(XY(:,1),XY(:,2),Ve,Ve);
     end
    axis square
    axis off
    title('Greens function')

