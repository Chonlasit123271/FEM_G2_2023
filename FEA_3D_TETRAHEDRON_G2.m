%this script was developed by the first year masters students from Asian Institute of Technology as a
%part of their FEM class project. the purpose of this code is to find the
%displacement at the top of a cantilever column and also the elemental
%stresses.

close all; clear all;  clc
tic
%-------------------Input Element Properties--------------------------------
E=2.1e10;            %Young Modulus 
nu=0.25;             %Poisson ratio
%-------------------Input support matrix and Load matrix--------------------
%[ Node ux uy uz]          ux=1 Fixed ;ux=0 Free
support=[1 1 1 1;
         2 1 1 1;
         3 1 1 1;
         4 1 1 1];

nSupport=size(support,1);

%Uses pre-created mesh from gmsh
%Run mesh file name without .m:
structuremesh025             % Name of the file without .m exported from Gmsh
XYZCoord = msh.POS;          % Node Coordinates
ELEMCon = msh.TETS(:,1:4);   % Element Connectivity 
NE=height(ELEMCon);          % Number of Elements
nNode=height(XYZCoord);      % Number of Nodes

% Plot the tetrahedron using tetramesh
%tetramesh(ELEMCon, XYZCoord);
%------(optional): Uncomment to see the plot
% Customize the plot    
% xlabel('X-axis');
% ylabel('Y-axis');
% zlabel('Z-axis');
% title('Tetrahedron');
%

% Node and Element Container 
for p=1:nNode
    NODE(p).X=XYZCoord(p,1);NODE(p).Y=XYZCoord(p,2);NODE(p).Z=XYZCoord(p,3);
end
for eNo=1:NE
    ELEMENT(eNo).con=ELEMCon(eNo,:);
end

%Shape Function
syms xi eta zeta
N(1)=1-xi-eta-zeta;
N(2)=zeta;               % x= xi ; y=eta ; z=zeta
N(3)=xi;                
N(4)=eta;
nN=length(N);

%----------------------------FORMING Global Stiffness matrix--------------------------------
for eNo=1:NE
    %Element Coordinates
    for i=1:nN
        xp(i)=XYZCoord(ELEMCon(eNo,i),1);
        yp(i)=XYZCoord(ELEMCon(eNo,i),2);
        zp(i)=XYZCoord(ELEMCon(eNo,i),3);
        ELEMENT(eNo).X(i)=NODE(ELEMENT(eNo).con(i)).X;
        ELEMENT(eNo).Y(i)=NODE(ELEMENT(eNo).con(i)).Y;
        ELEMENT(eNo).Z(i)=NODE(ELEMENT(eNo).con(i)).Z;
    end
    %Coordinate Mapping
    %Mapping
    x=0;y=0;z=0;
    for i=1:nN
    x=x+xp(i)*N(i);
    y=y+yp(i)*N(i);
    z=z+zp(i)*N(i);
    end
    %Jacobian
    J=[diff(x,xi) diff(y,xi) diff(z,xi);
      diff(x,eta) diff(y,eta) diff(z,eta);
      diff(x,zeta) diff(y,zeta) diff(z,zeta)];
    %Differentiating the N matrix to get doMat, dN
    for i=1:nN
        dN(:,i)=inv(J)*[diff(N(i),xi); diff(N(i),eta); diff(N(i),zeta)];
    end
    %Arranging the terms to get the B matrix
    %B=zeros(6,3*length(N));
    B = zeros(6, 3*nN);
    for i = 1:nN
        j = 3*i - 2;
        k = 3*i - 1;
        l = 3*i;
        B(:, [j, k, l]) = [dN(1,i), 0, 0; 
                     0, dN(2,i), 0;
                     0, 0, dN(3,i);
                     dN(2,i), dN(1,i), 0;
                     0, dN(3,i), dN(2,i);
                     dN(3,i), 0, dN(1,i)];
    end
    
    ELEMENT(eNo).B=B;
     % Define constitutive matrix (D)
    D = E / ((1+nu)*(1-2*nu)) * [1-nu, nu, nu, 0, 0, 0;
                                nu, 1-nu, nu, 0, 0, 0;
                                nu, nu, 1-nu, 0, 0, 0;
                                0, 0, 0, (1 - 2*nu) / 2, 0, 0;
                                0, 0, 0, 0, (1 - 2*nu) / 2, 0;
                                0, 0, 0, 0, 0, (1 - 2*nu) / 2];
   
    %Integrand is
    f= double(B'*D*B*det(J));
    %Direct Integration option 1: Direct Stiffness method
    ELEMENT(eNo).stiffness=double(int(int(int(f, zeta, 0, 1 - xi - eta), eta, 0, 1 - xi), xi, 0, 1));
    
    %Gauss Points option 2: Gauss Points mehtod
    % [xg,wg]=lgwt(2,-1,1);
    % %Gaussian evaluation of integral
    % ELEMENT(eNo).stiffness=wg(2)*wg(2)*subs(f,[xi,eta,zeta],[xg(2),xg(2)])+wg(1)*wg(2)*subs(f,[xi,eta,zeta],[xg(1),xg(2)])+wg(1)*wg(1)*subs(f,[xi,eta,zeta],[xg(1),xg(1)])+wg(2)*wg(1)*subs(f,[xi,eta,zeta],[xg(2),xg(1)]);
end
%-------------------------------------------------------------------------%
%                  Global Stiffness Matrix Calculation                    %
%-------------------------------------------------------------------------%
% This part calculates the global stiffness matrix. Basically; for each
% element it takes 12x12 part from the element stiffness matrix and puts to
% the correct spot on the global stiffness matrix. This process loops until
% all elements all parts placed in to the global stiffness matrix.
nDof=nNode*3;
KG = zeros(nDof,nDof);
for eNo=1:NE
    kElem=  ELEMENT(eNo).stiffness;
    for j=1:nN
        for i=1:nN
            n = ELEMCon(eNo,i);
            m = ELEMCon(eNo,j);
            KG(3*n-2,3*m-2) = KG(3*n-2,3*m-2)+kElem(3*i-2,3*j-2);              %K11
            KG(3*n-2,3*m-1) = KG(3*n-2,3*m-1)+kElem(3*i-2,3*j-1);              %K12
            KG(3*n-2,3*m) = KG(3*n-2,3*m)+kElem(3*i-2,3*j);                    %K13
            KG(3*n-1,3*m-2) = KG(3*n-1,3*m-2)+kElem(3*i-1,3*j-2);              %K21
            KG(3*n-1,3*m-1) = KG(3*n-1,3*m-1)+kElem(3*i-1,3*j-1);              %K22
            KG(3*n-1,3*m) = KG(3*n-1,3*m)+kElem(3*i-1,3*j);                    %K23
            KG(3*n,3*m-2) = KG(3*n,3*m-2)+kElem(3*i,3*j-2);                    %K13
            KG(3*n,3*m-1) = KG(3*n,3*m-1)+kElem(3*i,3*j-1);                    %K23
            KG(3*n,3*m) = KG(3*n,3*m)+kElem(3*i,3*j);                          %K33
        end
    end
end
        % Checking Global stiffness Matrix Correct or not
    % if any(diag(KG) < 0)          
    %     disp('The Stiffness is not Correct')
    % else
    %     disp('The Stiffness maybe Correct')
    % end
    %-------------------------------------------------------------------------%
    %           Apply Support Conditions to the Global Stiffnes Matrix        %
    %-------------------------------------------------------------------------%
    % This part makes zeros all columns and rows where the supports are except
    % the diagonal element of the matrix. Diagonal element set to 1. I choose
    % this method because with this way sort of displacement evaulated are not
    % changes. And later it is easier to use evaluated values. Only negative
    % side of this approach is we have to be careful not to put force where the
    % support is fixed.
for i=1:nSupport
    n = support(i,1);
    if (support(i,2) == 1)      %%%ux =0
        KG(3*n-2,:) = 0;
        KG(:,3*n-2) = 0;
        KG(3*n-2,3*n-2) = 1;
    end
    if (support(i,3) == 1)      %%%uy=0
        KG(3*n-1,:) = 0;
        KG(:,3*n-1) = 0;
        KG(3*n-1,3*n-1) = 1;
    end
    if   (support(i,4) == 1)    %%%uz=0
        KG(3*n,:) = 0;
        KG(:,3*n) = 0;
        KG(3*n,3*n) = 1;
    end
end

%-------------------------------------------------------------------------%
%                       Load Vector Forming                               %
%-------------------------------------------------------------------------%
% In this part load vector created. If there is a load vector get the value
% from load matrix. If not not load value set to zero. 

f = zeros(nDof,1);

%Load conditions [Node Fx Fy Fz]

load=[6 1 0 0                  
      7 1 0 0];                 % Manual Input for node of acting force
nLoad=size(load,1);

for i=1:nLoad
   n = load(i,1);
   f(3*n-2) = load(i,2);
   f(3*n-1) = load(i,3);
   f(3*n) = load(i,4);
end

u_Disp_load = KG\f;


% Have to mannually prepare the incremental displacement load in seperate file. Then Call for it
Data=xlsread('Displacement_data.xlsx')  
Increment=Data(:,1);
u_disp=Data(:,2);          %Value of each step of displacement
N_INC=height(Increment);   % Number of increment
Load_step=zeros(N_INC,1); 
% Forming Value of Load step 
for i=1:N_INC
    Load_step(i,:)=u_disp(i)/u_Disp_load(19,1); 
end

%---------------------------------SOLVER-----------------------------------%
for j=1:N_INC   
    fvec=zeros(nDof,1);
    load=[6 Load_step(j)  0 0           % Manual Input for node of acting force
          7 Load_step(j) 0 0];            % %Load conditions [Node Fx Fy Fz]
    nLoad=size(load,1);
    for i=1:nLoad
       n = load(i,1);
       fvec(3*n-2) = load(i,2);
       fvec(3*n-1) = load(i,3);
       fvec(3*n) = load(i,4);
    end    
    udisp(:,j)=KG\fvec;  %Global nodal displacement
end


figure()
% Plot the input displacement 
plot(Increment,u_disp);
xlabel('Increment');
ylabel('Displacement');
title('Input displacement');
grid on
figure()
% Plot the input displacement 
plot(Increment,Load_step);
xlabel('Increment');
ylabel('Force');
title('Force from input displacement');
grid on


for i = 1:nNode
    u(:,i)=[udisp(3*i-2,N_INC);udisp(3*i-1,N_INC); udisp(3*i,N_INC)];
    NODE(i).u = [udisp(3*i-2,N_INC);udisp(3*i-1,N_INC); udisp(3*i,N_INC)];   
end

% 3D Element Displacements Calculation
for i = 1:NE
    temp = [];
    for j = 1:nN
        temp = [temp; NODE(ELEMCon(i,j)).u];
    end
    ELEMENT(i).u = temp;
end

% 3D Element Stresses Calculation
SCP = [-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1]; % Stress calculation points
for i = 1:NE
    uelem = ELEMENT(i).u;
    Belem = ELEMENT(i).B; % Assuming B matrix is defined for 3D elements
    sigma = D * Belem * uelem;
    for j = 1:4
        sigma2(:,j)=double(subs(sigma, [xi ,eta ,zeta], SCP(j,:)));
        ELEMENT(i).sigma(:,j) = double(subs(sigma, [xi ,eta ,zeta], SCP(j,:))); 
    end
end
%-- For showing the stress in each element of the strcuture-------
% for i=1:NE
%     disp(ELEMENT(i).sigma(1,1)) %
% end


%-------------------3D Displacement Deformation Visual of structure----------
Cor_dis=u'+XYZCoord;
%Deformation scaling factor (adjust as needed)
scale_factor = 0.1;
% Plot the tetrahedron using tetramesh for undeformed and deformed meshes
figure;
tetramesh(ELEMCon, XYZCoord, 'FaceColor', 'b', 'EdgeColor', 'k');  % Undeformed in blue
hold on;
tetramesh(ELEMCon, Cor_dis, 'FaceColor', 'r', 'EdgeColor', 'none'); % Deformed in red

% Customize the plot
xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');
title('Tetrahedron Deformation');

% Add reference lines for deformation
for i = 1:nNode
    plot3([XYZCoord(i,1), Cor_dis(i,1)], [XYZCoord(i,2), Cor_dis(i,2)], [XYZCoord(i,3), Cor_dis(i,3)], 'k--');
end



%---------------------------------Paraview Visualization of Displacement and Stress-----------------------------------%
%------This part of the code will generate the .vtk to open for viewing
%------the Displacement and Stress of the element-------------------
filename = 'model.vtk';        % Exported File name to open in paraview ("model.vtk")
fid = fopen(filename, 'w');
fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, 'Model Data\n');
fprintf(fid, 'ASCII\n');
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');

% Write node coordinates
fprintf(fid, 'POINTS %d float\n', nNode);
for i = 1:nNode
    fprintf(fid, '%f %f %f\n', XYZCoord(i, :));
end

% Write element connectivity
fprintf(fid, 'CELLS %d %d\n', NE, NE*(nN+1));
for i = 1:NE
    fprintf(fid, '%d', nN);
    for j = 1:nN
        fprintf(fid, ' %d', ELEMCon(i, j)-1); % Subtract 1 to convert to zero-based indexing
    end
    fprintf(fid, '\n');
end

% Write element types
fprintf(fid, 'CELL_TYPES %d\n', NE);
for i = 1:NE
    fprintf(fid, '10\n'); % Assuming tetrahedral elements (VTK cell type 10)
end

% Write displacement data
fprintf(fid, 'POINT_DATA %d\n', nNode);
fprintf(fid, 'VECTORS Displacement float\n');
for i = 1:nNode
    fprintf(fid, '%f %f %f\n', NODE(i).u);
end

% Write stress data
fprintf(fid, 'CELL_DATA %d\n', NE);
fprintf(fid, 'SCALARS Stress float\n');
fprintf(fid, 'LOOKUP_TABLE default\n');
for i = 1:NE
    for j = 1:4 % Assuming you have 4 stress components
        fprintf(fid, '%f\n', ELEMENT(i).sigma(:,j)');
    end
end

fclose(fid);
toc

%---------------------------------Thank you-----------------------------------------------%

