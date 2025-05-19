close all; clear all; clc;
tic

%% --------------------------Inputs----------------------------------------
%-------------------Inputs for our structure (column)
height = 3;          %height of column in m
breadth = 0.25;      %breadth of column in m
width = 0.25;        %width of column in m

%-------------------Inputs for Element Properties--------------------------------
E=3.163*10^10;            %Young Modulus in N/m^2
nu=0.2;                  %Poisson ratio
rho = 2400;          %Density of concrete in kg/m^3
g = 9.81;            %Acceleration due to gravity, m/s^2
SigmaMax_C = 35*10^6; % maximum compressive stress in N/m^2 (assuming 25 MPa cocrete)
SigmaMax_T = 10^6*0.62*sqrt(SigmaMax_C/10^6);  % maximum tensile stress in N/m^2 = 0.62*sqrt(SigmaMax_C) from ACI 318-19 


%% -----------Run mesh file to get the node coordinates and nodes:
MeshG2            
nodeCoordinates = msh.POS;          % Node Coordinates in m
ELEMCon = msh.TETS(:,1:4);          % Element Connectivity 
NE=length(ELEMCon);                 % Number of Elements
nNode=length(nodeCoordinates);      % Number of Nodes
nDOF = nNode*3;

%% ---------------------Inputs for Newmark Beta----------------------------
%We are using Newmark Beta with constant acceleration (trapezoidal rule).
%The parameters used in this method are given below:
beta = 1/4;
gamma = 1/2;

% Node and Element Container 
for p=1:nNode
    NODE(p).X=nodeCoordinates(p,1);NODE(p).Y=nodeCoordinates(p,2);NODE(p).Z=nodeCoordinates(p,3);
end
for eNo=1:NE
    ELEMENT(eNo).con=ELEMCon(eNo,:);
end

%% ---------------Shape Function  % x= xi ; y=eta ; z=zeta-----------------
syms xi eta zeta
N(1)=1-xi-eta-zeta;
N(2)= xi;              
N(3)= eta;                
N(4)= zeta;
nN=length(N);

%for TET4
N_matrix = [N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0 0;
            0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4) 0;
            0 0 N(1) 0 0 N(2) 0 0 N(3) 0 0 N(4)];
 
%for mass matrix
int_N = double(int(int(int(N_matrix'*N_matrix, zeta, 0, 1 - xi - eta), eta, 0, 1 - xi), xi, 0, 1));

% for body force
b = [0; 0; -rho*g];
int_Nb = double(int(int(int(N_matrix'*b, zeta, 0, 1 - xi - eta), eta, 0, 1 - xi), xi, 0, 1));
%% ----FORMING Element Mass Matrix, Stiffness matrix, Body Force Vector----
for eNo=1:NE
    %Element Coordinates
    for i=1:nN
        xp(i)=nodeCoordinates(ELEMCon(eNo,i),1);
        yp(i)=nodeCoordinates(ELEMCon(eNo,i),2);
        zp(i)=nodeCoordinates(ELEMCon(eNo,i),3);
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
     
    % Assuming Constant Stiffness 
    ELEMENT(eNo).E = E;  %Young's Modulus of Elasticity (we are using this after failure criteria check)

    %% -------------------------Jacobian-----------------------------------

    % J=[diff(x,xi) diff(y,xi) diff(z,xi);
    %   diff(x,eta) diff(y,eta) diff(z,eta);
    %   diff(x,zeta) diff(y,zeta) diff(z,zeta)];

    %Constant for an element(for TET4 only)
    J = [ELEMENT(eNo).X(2)-ELEMENT(eNo).X(1) ELEMENT(eNo).Y(2)-ELEMENT(eNo).Y(1) ELEMENT(eNo).Z(2)-ELEMENT(eNo).Z(1);
         ELEMENT(eNo).X(3)-ELEMENT(eNo).X(1) ELEMENT(eNo).Y(3)-ELEMENT(eNo).Y(1) ELEMENT(eNo).Z(3)-ELEMENT(eNo).Z(1);
         ELEMENT(eNo).X(4)-ELEMENT(eNo).X(1) ELEMENT(eNo).Y(4)-ELEMENT(eNo).Y(1) ELEMENT(eNo).Z(4)-ELEMENT(eNo).Z(1)];
      ELEMENT(eNo).det_J = det(J); 
      ELEMENT(eNo).J = J;

    %-------------Differentiating the N matrix to get doMat, dN---------
    % Option 1:
    % for i=1:nN
    %     dN(:,i)=inv(J)*[diff(N(i),xi); diff(N(i),eta); diff(N(i),zeta)];
    % end

    % Option 2: The derivative of N with respect to parent axes is constant
    %(for TET4 only)
    dN = J\[-1 1 0 0; -1 0 1 0; -1 0 0 1];
    % ELEMENT(eNo).dN = dN;

%% -------------------------Mass Matrix------------------------------------
 ELEMENT(eNo).Mass = rho*int_N*double(det(J));

%% --------------------------B matrix and D matrix-------------------------
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
                                0, 0, 0, (1 - 2*nu)/2, 0, 0;
                                0, 0, 0, 0, (1-2*nu)/2, 0;
                                0, 0, 0, 0, 0, (1-2*nu)/2];

   ELEMENT(eNo).D = D;
   
   %% ---------------------Element Stiffness Matrix------------------------
    %Integrand is
    f_k = double(B'*D*B*det(J));

    % %Direct Integration option 1: Direct Stiffness method
    % ELEMENT(eNo).stiffness=double(int(int(int(f_k, zeta, 0, 1 - xi - eta), eta, 0, 1 - xi), xi, 0, 1));

    %Gauss Points option 2: Gauss Points mehtod
    
     %Option 3: only fo TET4 elements
        %For Tetrahedral element with 4 nodes (linear element), there are no
        %variables in integral. So, for faster computation, we are using the
        %constant value i.e. 
        % int(int(int(1, zeta, 0, 1 - xi - eta), eta, 0, 1 -xi), xi, 0, 1) = 1/6
        ELEMENT(eNo).stiffness = 1/6*f_k;

%% --------Element body force (Constant for each element)------------------

ELEMENT(eNo).BodyForce = int_Nb*double(det(J));

end

% %% --------------------Global Mass Matrix Calculation----------------------
MG = zeros(nDOF,nDOF);
for eNo=1:NE
    mElem=  ELEMENT(eNo).Mass;
    for j=1:nN
        for i=1:nN
            n = ELEMCon(eNo,i);
            m = ELEMCon(eNo,j);
            MG(3*n-2,3*m-2) = MG(3*n-2,3*m-2)+mElem(3*i-2,3*j-2);              %m11
            MG(3*n-2,3*m-1) = MG(3*n-2,3*m-1)+mElem(3*i-2,3*j-1);              %m12
            MG(3*n-2,3*m) = MG(3*n-2,3*m)+mElem(3*i-2,3*j);                    %m13
            MG(3*n-1,3*m-2) = MG(3*n-1,3*m-2)+mElem(3*i-1,3*j-2);              %m21
            MG(3*n-1,3*m-1) = MG(3*n-1,3*m-1)+mElem(3*i-1,3*j-1);              %m22
            MG(3*n-1,3*m) = MG(3*n-1,3*m)+mElem(3*i-1,3*j);                    %m23
            MG(3*n,3*m-2) = MG(3*n,3*m-2)+mElem(3*i,3*j-2);                    %m13
            MG(3*n,3*m-1) = MG(3*n,3*m-1)+mElem(3*i,3*j-1);                    %m23
            MG(3*n,3*m) = MG(3*n,3*m)+mElem(3*i,3*j);                          %m33
        end
    end
end

%% -----------------Global Stiffness Matrix Calculation--------------------

% Initialize the global stiffness matrix
KG = zeros(nDOF, nDOF);

% Loop over all elements
for eNo = 1:NE
    % Retrieve element stiffness matrix
    kElem = ELEMENT(eNo).stiffness;
    
    % Loop over local nodes in the element
    for j = 1:nN
        for i = 1:nN
          
            n = ELEMCon(eNo, i);  % Global node index for local node i
            m = ELEMCon(eNo, j);  % Global node index for local node j
            
            % Assemble element stiffness into global stiffness matrix
            KG(3*n-2, 3*m-2) = KG(3*n-2, 3*m-2) + kElem(3*i-2, 3*j-2);  % K11
            KG(3*n-2, 3*m-1) = KG(3*n-2, 3*m-1) + kElem(3*i-2, 3*j-1);  % K12
            KG(3*n-2, 3*m)   = KG(3*n-2, 3*m)   + kElem(3*i-2, 3*j);    % K13
            
            KG(3*n-1, 3*m-2) = KG(3*n-1, 3*m-2) + kElem(3*i-1, 3*j-2);  % K21
            KG(3*n-1, 3*m-1) = KG(3*n-1, 3*m-1) + kElem(3*i-1, 3*j-1);  % K22
            KG(3*n-1, 3*m)   = KG(3*n-1, 3*m)   + kElem(3*i-1, 3*j);    % K23
            
            KG(3*n,   3*m-2) = KG(3*n,   3*m-2) + kElem(3*i,   3*j-2);  % K31
            KG(3*n,   3*m-1) = KG(3*n,   3*m-1) + kElem(3*i,   3*j-1);  % K32
            KG(3*n,   3*m)   = KG(3*n,   3*m)   + kElem(3*i,   3*j);    % K33
        end
    end
end


KG_initial = KG;


%% ---------------------Global Body force----------------------------------
Body_ForceG = zeros(nDOF,1);
for eNo=1:NE
    BF_Elem=  ELEMENT(eNo).BodyForce;
    for i=1:nN
            n = ELEMCon(eNo,i);
            Body_ForceG(3*n-2) = Body_ForceG(3*n-2)+BF_Elem(3*i-2);             
            Body_ForceG(3*n-1) = Body_ForceG(3*n-1)+BF_Elem(3*i-1);                          
            Body_ForceG(3*n) = Body_ForceG(3*n)+BF_Elem(3*i);                         
    end
end

%% ------------------------Apply Support Conditions------------------------

%Identifying nodes for boundary conditions

% Nodes where z =0m (zero displacement in all directions)
nodes_z0 = find(nodeCoordinates(:,3) == 0);

% Nodes where x == 0.25m and z == 3m (displacement u in x-direction only)
nodes_x025_z3 = find(nodeCoordinates(:,1) == breadth & nodeCoordinates(:,3) == height);

bdof_s = [];
bddis_s = [];
% Zero displacement at nodes_z0 (x,y,z)---> same for every iteration
for i = 1:length(nodes_z0)
    node = nodes_z0(i);
    dofs = 3*(node-1) + (1:3); % DOFs for x,y,z
    bdof_s = [bdof_s, dofs];
    bddis_s = [bddis_s, zeros(1,3)];
end


%% -------------------Incremental Displacement Vector---------------------- 
Totaltime = 100;             %Gradual application of load during total time, seconds
N_INC = 1000;
del_t = Totaltime/N_INC;                 %time interval for each step, seconds (del_t should be small enough to capture nonlinear behavior; such as Totaltime = 10, and N_INC=1000)
disp_top = 4/100*height;                 %Given, displacement at top 4% of height 
disp_step = disp_top/N_INC;              %Linear Increment of displacment at top of column
applied_d = disp_step*ones(N_INC,1);
R_ForceG = zeros(nDOF,1);                %initially, no residual force
displ = zeros(nDOF,1);
vel = zeros(nDOF,1);
acc = zeros(nDOF,1);
output_disp = [zeros(nDOF,1)];
fail = [];
 
%% -----------------------Iteration starts---------------------------------
for INC = 1:32
   
    %Displacement Bounary condition
       
        bdof_d = [];
        bddis_d = [];
        % Displacement u at nodes_x025_z3 (only x axis)
        for i = 1:length(nodes_x025_z3)
            node = nodes_x025_z3(i);
            dof_x = 3*(node-1) + 1; % x-direction DOF only
            bdof_d = [bdof_d, dof_x];
            bddis_d = [bddis_d, applied_d(INC)];
        end
        
        bdof = [bdof_d bdof_s];
        bddis =[ bddis_d bddis_s];

        % Modify force vector to account for prescribed displacement
         % F_eff = Body_ForceG + R_ForceG + MG*vel/(beta*del_t) + MG*acc*(1-2*beta)/(2*beta);  
         if INC ==1
         F_eff = Body_ForceG + MG*vel/(beta*del_t) + MG*acc*(1-2*beta)/(2*beta);
         else
          F_eff = R_ForceG + MG*vel/(beta*del_t) + MG*acc*(1-2*beta)/(2*beta);
         end
          KG_eff = KG + MG/(beta*del_t*del_t);
          KG1 = KG_eff;          
          
          %% --------------Applying Boundary conditions--------------------
          for i = 1:length(bdof)
            dof = bdof(i);
            known_disp = bddis(i);
            KG_dof = KG1(:, dof);
            KG_dof(dof) = -1; 
            F_eff(dof) = 0;
            F_eff = F_eff - KG_dof * known_disp;

            % Zero out row and column in stiffness matrix
            KG1(:, dof) = 0;
            KG1(dof, :) = 0;
        
            % Set diagonal to 1 to enforce displacement
            KG1(dof, dof) = 1;
        
         end
        del_u = KG1\F_eff;        % in each iteration   
        
        %store the total displacement at each iteration
        output_disp = [output_disp output_disp(:,INC)+del_u];


        %Element Displacement, Strains and Dispalcement
        displ = displ + del_u;      %displacement increases at time step
        
        %% ------------------------Newmark Beta----------------------------
        % It is used to get the velocity and acceleration for next time step

        acc_next = 1/(beta*del_t*del_t)*del_u - 1/(beta*del_t)*vel - (1-2*beta)/(2*beta)*acc;
        vel_next = gamma/(beta*del_t)*del_u - (gamma-beta)/(beta)*vel - del_t*(gamma-2*beta)/(2*beta)*acc;
        acc = acc_next;
        vel = vel_next;
        % Dispalcement vector
        for i = 1:nNode
            u(:,i)=[displ(3*i-2);displ(3*i-1); displ(3*i)];
            NODE(i).u = [displ(3*i-2);displ(3*i-1); displ(3*i)];   
        end
   
 %% --------------3D Element Displacements Calculation---------------------
     for eNo = 1:NE
            temp = [];
            for j = 1:nN
                temp = [temp; NODE(ELEMCon(eNo,j)).u];
            end

            ELEMENT(eNo).u = temp;
            Strain = ELEMENT(eNo).B*ELEMENT(eNo).u;    
           
            
            %Tensor form:
            Strain_Tensor = [Strain(1) Strain(4)/2 Strain(6)/2;    
                      Strain(4)/2 Strain(2) Strain(5)/2;
                      Strain(6)/2 Strain(5)/2 Strain(3)];
            ELEMENT(eNo).Strain_Tensor = Strain_Tensor;
            
            %% ---------Eigen analysis for principal strains---------------
        
            % Calculating Principal strains and Rotation vector by solving Eigen Equations
            % R = Rotation Vector and P_Strain = Principal Strains
            [V, P_S] = eig(Strain_Tensor);
            R = V;
            P_Strain = [P_S(1,1);
                        P_S(2,2);
                        P_S(3,3)];           
             ELEMENT(eNo).P_Strain = P_Strain;
         ELEMENT(eNo).Cal_P_Stress = ELEMENT(eNo).E*P_Strain;  
         %% ---Actual Principal strains, stresses for each element-------
           for j = 1:3    
                if P_Strain(j)< 0  
                [P_Stress(j)] = GetStressHDCompression(ELEMENT(eNo).E, SigmaMax_C, abs(P_Strain(j)));  % Compression stress
                else
                 [P_Stress(j)] = ELEMENT(eNo).E*P_Strain(j);    %tension
                end
                Act_P_Stress(j) = P_Stress(j) * sign(P_Strain(j));  % Apply sign convention for stress
           
               %% ------------------Failure criteria-----------------------
                
               ELEMENT(eNo).Act_P_Stress =  Act_P_Stress';
               
           end
       %for HD Stress Strain model, the concrete stress will not
       %reach maximumm compressive stress level so, only checking for tension
                  
           if any(Act_P_Stress >= SigmaMax_T)       
              
               fail = [fail; eNo INC];
                  disp("Element "+num2str(eNo)+" failed in tension at prescribed displacment of "+num2str(applied_d(INC)*INC)+" at step "+num2str(INC)+ ".");       
                       %1. Remove the element
                       %2. Reduce the stiffness (we use this approach)
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       ELEMENT(eNo).E = ELEMENT(eNo).E/10000000;
                       ELEMENT(eNo).stiffness = ELEMENT(eNo).stiffness/10000000; 
                       % ELEMENT(eNo).E = 1;
                       % ELEMENT(eNo).stiffness = ones(12,12);
           end

                 ELEMENT(eNo).Residual_P_Stress = ELEMENT(eNo).Cal_P_Stress - ELEMENT(eNo).Act_P_Stress;

              
        %% -------------Residual Stress in Cartesian Coordinates-----------
        % diag---> diagonal matrix
        R_tensor = R*diag(ELEMENT(eNo).Residual_P_Stress)*R';               %Residual stress tensor in Cartesian planes
        ELEMENT(eNo).R_Stress = [R_tensor(1,1);
                                      R_tensor(2,2);
                                      R_tensor(3,3);
                                      R_tensor(1,2);
                                      R_tensor(2,3);
                                      R_tensor(1,3)];
        
        %% -----Residual force for an element in Cartesian Coordinates-----
        
        f= (ELEMENT(eNo).B)'*ELEMENT(eNo).R_Stress*ELEMENT(eNo).det_J;
        
        %Direct Integration option 1:
        %     ELEMENT(eNo).R_Force=double(int(int(int(f, zeta, 0, 1 - xi - eta), eta, 0, 1 - xi), xi, 0, 1));
        
        %Gauss Points option 2: Gauss Points mehtod
            
        %Option 3: only fo TET4 elements
        %For Tetrahedral element with 4 nodes (linear element), there are no
        %variables in integral. So, for faster computation, we are using the
        %constant value i.e. 
        % int(int(int(1, zeta, 0, 1 - xi - eta), eta, 0, 1 -xi), xi, 0, 1) = 1/6
        ELEMENT(eNo).R_Force = 1/6*f;
     end

    %% ----------------Re_assembling Global Stiffness matrix---------------
    KG = zeros(nDOF, nDOF);
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

        %% ------------------Assembly of Residual force at each step--------------------
         R_ForceG = zeros(nDOF,1);
        for eNo=1:NE
            f_Elem=  ELEMENT(eNo).R_Force;
            for i=1:nN
                    n = ELEMCon(eNo,i);
                    R_ForceG(3*n-2) = R_ForceG(3*n-2)+f_Elem(3*i-2);             
                    R_ForceG(3*n-1) = R_ForceG(3*n-1)+f_Elem(3*i-1);                          
                    R_ForceG(3*n) = R_ForceG(3*n)+f_Elem(3*i);                         
            end
        end
end

COORD_deformed = zeros(size(nodeCoordinates));
for i = 1:size(nodeCoordinates,1)
    COORD_deformed(i,:) = nodeCoordinates(i,:) + (NODE(i).u)';
end

    % Plot the mesh
    figure
    tetramesh(ELEMCon,COORD_deformed,'FaceAlpha',0.5);
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    title('Deformated Mesh of the Portal Frame'); 
    view(3); % Show the plot directly in 3D
toc

keyboard 

%% ============== Post-processing ====================
% ----------------------- 1. Deflected Displacement -----------------------------
% Export Displacment, Node coodinate, Element connectivity to Excel files
Displacement_inc = output_disp; 
writematrix(Displacement_inc, 'Displacement with increments.xlsx'); 

Node_coordinates = nodeCoordinates; 
writematrix(Node_coordinates, 'nodecoordinates.xlsx'); 

Element_connectivity = ELEMCon; 
writematrix(Element_connectivity, 'ElemCon.xlsx'); 

% User parameters
amplification = 20;           % Displacement amplification factor
filename_disp = 'Displacement with increments.xlsx';
filename_elem = 'ElemCon.xlsx';
filename_node = 'nodecoordinates.xlsx';

gif_filename = 'mesh_deformation.gif';  % Output GIF file
video_filename = 'mesh_deformation.mp4'; % Output video file
frame_delay = 0.5;           % Delay between frames (seconds)

% Read data using xlsread for compatibility
[disp_data, ~, ~] = xlsread(filename_disp);
[elemcon, ~, ~]   = xlsread(filename_elem);
[nodecoord, ~, ~] = xlsread(filename_node);

% Read data using xlsread for compatibility
[disp_data, ~, ~] = xlsread(filename_disp);
[elemcon, ~, ~]   = xlsread(filename_elem);
[nodecoord, ~, ~] = xlsread(filename_node);

% Validate sizes
[num_dofs, n_inc] = size(disp_data);
num_nodes = size(nodecoord,1);
if num_dofs ~= 3*num_nodes
    error('Displacement rows must equal 3*N_nodes');
end

% Reshape displacement: N_nodes x 3 x N_increments
U = reshape(disp_data, 3, num_nodes, n_inc);
U = permute(U, [2 1 3]);

% Prepare figure
fig = figure('Color','w');
view(3); grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Mesh Deformation');
hold on;

% Define local edges of a tetrahedral element (4-node elements)
edges = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

% Loop through increments and capture frames
gif_created = false;
for k = 1:n_inc
    cla;
    deformed = nodecoord + amplification * squeeze(U(:,:,k));

    % Plot elements
    for e = 1:size(elemcon,1)
        coords = deformed(elemcon(e,:), :);
        for ed = 1:size(edges,1)
            p1 = coords(edges(ed,1), :);
            p2 = coords(edges(ed,2), :);
            plot3([p1(1),p2(1)], [p1(2),p2(2)], [p1(3),p2(3)], 'b-');
        end
    end

    % Plot nodes
    scatter3(deformed(:,1), deformed(:,2), deformed(:,3), 20, 'r', 'filled');

    % Title and increment counter
    title(sprintf('Increment %d of %d', k, n_inc));

    % Compute displacement magnitudes
    displacement_magnitudes = sqrt(sum((amplification * squeeze(U(:,:,k))).^2, 2));
    [max_disp, max_idx] = max(displacement_magnitudes);

    % Display max displacement at top of figure
    xlims = xlim; ylims = ylim; zlims = zlim;
    disp_text = sprintf('                                                                                                                              Maximum displacement = %.4f (Node %d)', max_disp/amplification, max_idx); % unscaled
    text(mean(xlims), ylims(2) + 0.05*(ylims(2)-ylims(1)), zlims(2), disp_text, ...
        'HorizontalAlignment', 'center', 'FontSize', 6, 'FontWeight', 'bold', 'Color', 'k');

    drawnow;

    % GIF capture
    frame = getframe(fig);
    [imind, cm] = rgb2ind(frame2im(frame), 256);
    if ~gif_created
        imwrite(imind, cm, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', frame_delay);
        gif_created = true;
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end
end

fprintf('Animated GIF saved as %s with %.1f s delay per frame.\n', gif_filename, frame_delay);
drawnow;
% -------------------------------------------------------------------------------------------- %

% ----------------------- 2. Display failed elements (animation) --------------------------
fprintf('\n\n========= Tension Failure Report =========\n');
fprintf('Element\t|\tStep\t|\tApplied Displacement (m)\n');
fprintf('------------------------------------------\n');
for i = 1:size(fail,1)
    eID = fail(i,1);
    stepID = fail(i,2);
    disp_val = applied_d(stepID) * stepID;
    fprintf('%d\t|\t%d\t|\t%.6f\n', eID, stepID, disp_val);
end
fprintf('------------------------------------------\n');

% 2.1 Visualization of Failed Elements Over Time 
unique_steps = unique(fail(:,2));  % List of steps where failure occurred

for s = 1:length(unique_steps)
    stepID = unique_steps(s);
    failed_in_step = fail(fail(:,2) == stepID, 1);  % Elements failed in this step

    figure;
    hold on;

    % Show full undeformed mesh (optional)
    tetramesh(ELEMCon, nodeCoordinates, ...
              'FaceAlpha', 0.1, ...
              'EdgeAlpha', 0.1, ...
              'FaceColor', [0.7 0.7 0.7]);

    % Show deformed mesh
    tetramesh(ELEMCon, COORD_deformed, ...
              'FaceAlpha', 0.2, ...
              'EdgeAlpha', 0.2, ...
              'FaceColor', [0.4 0.4 0.4]);

    % Plot failed elements in green for this increment
    for i = 1:length(failed_in_step)
        eNo = failed_in_step(i);
        nodes = ELEMCon(eNo, :);
        coords = COORD_deformed(nodes, :);
        patch('Faces', [1 2 3; 1 2 4; 1 3 4; 2 3 4], ...
              'Vertices', coords, ...
              'FaceColor', [1 0 0], ...
              'FaceAlpha', 1.0, ...
              'EdgeColor', 'none');
    end

    xlabel('X (m)', 'FontWeight', 'bold');
    ylabel('Y (m)', 'FontWeight', 'bold');
    zlabel('Z (m)', 'FontWeight', 'bold');
    title(sprintf('Tension Failure at Increment %d', stepID), 'FontWeight', 'bold');
    view(3);
    axis equal;
    grid on;
    legend({'Original', 'Deformed', 'Failed Elements (Green)'});
end

% 2.2  Visualization of All Failed Elements Together 
figure;
hold on;

% Show full undeformed mesh (optional)
tetramesh(ELEMCon, nodeCoordinates, ...
          'FaceAlpha', 0.1, ...
          'EdgeAlpha', 0.1, ...
          'FaceColor', [0.7 0.7 0.7]);

% Show deformed mesh
tetramesh(ELEMCon, COORD_deformed, ...
          'FaceAlpha', 0.2, ...
          'EdgeAlpha', 0.2, ...
          'FaceColor', [0.4 0.4 0.4]);

% Plot all failed elements
for i = 1:size(fail,1)
    eNo = fail(i,1);
    nodes = ELEMCon(eNo, :);
    coords = COORD_deformed(nodes, :);
    patch('Faces', [1 2 3; 1 2 4; 1 3 4; 2 3 4], ...
          'Vertices', coords, ...
          'FaceColor', [1 0 0], ...  % Red for combined
          'FaceAlpha', 1.0, ...
          'EdgeColor', 'none');
end

xlabel('X (m)', 'FontWeight', 'bold');
ylabel('Y (m)', 'FontWeight', 'bold');
zlabel('Z (m)', 'FontWeight', 'bold');
title('All Tension Failures Across All Increments', 'FontWeight', 'bold');
view(3);
axis equal;
grid on;
legend({'Original', 'Deformed', 'Failed Elements (Red)'});

% 2.3 Animation of Failed Elements Over Time
unique_steps = unique(fail(:,2));  % List of time steps with failure

figure;
hold on;

% Plot original and deformed mesh (once)
tetramesh(ELEMCon, nodeCoordinates, ...
          'FaceAlpha', 0.1, 'EdgeAlpha', 0.1, ...
          'FaceColor', [0.7 0.7 0.7]);
tetramesh(ELEMCon, COORD_deformed, ...
          'FaceAlpha', 0.2, 'EdgeAlpha', 0.2, ...
          'FaceColor', [0.4 0.4 0.4]);

xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
view(3); axis equal; grid on;

h_failed = [];  % Handle for failed patches

for s = 1:length(unique_steps)
    stepID = unique_steps(s);
    failed_in_step = fail(fail(:,2) == stepID, 1);  % Failed elements at this step

    % Remove previous failed patches
    if ~isempty(h_failed)
        delete(h_failed);
    end

    % Plot new failed elements
    h_failed = gobjects(length(failed_in_step), 1);
    for i = 1:length(failed_in_step)
        eNo = failed_in_step(i);
        nodes = ELEMCon(eNo, :);
        coords = COORD_deformed(nodes, :);
        h_failed(i) = patch('Faces', [1 2 3; 1 2 4; 1 3 4; 2 3 4], ...
                            'Vertices', coords, ...
                            'FaceColor', [1 0 0], 'FaceAlpha', 1.0, ...
                            'EdgeColor', 'none');
    end

    title(sprintf('Tension Failure at Increment %d', stepID));
    drawnow;
    pause(1.0);  % Adjust speed here (in seconds)
end
% 2.4  Cumulative Animation of Failed Elements Over Time with Video Export

unique_steps = unique(fail(:,2));  % Unique load increments
cumulative_failed = [];            % Store cumulative failed elements

% Set up video writer 
video_filename = 'FailedElementsAnimation.mp4';  % Change to .avi if needed
v = VideoWriter(video_filename, 'MPEG-4');  % Use 'Motion JPEG AVI' for AVI
v.FrameRate = 1;  % 1 frame per second (adjust as needed)
open(v);

% Create figure 
fig = figure;
hold on;

% Plot undeformed mesh
tetramesh(ELEMCon, nodeCoordinates, ...
          'FaceAlpha', 0.1, 'EdgeAlpha', 0.1, ...
          'FaceColor', [0.7 0.7 0.7]);

% Plot deformed mesh
tetramesh(ELEMCon, COORD_deformed, ...
          'FaceAlpha', 0.2, 'EdgeAlpha', 0.2, ...
          'FaceColor', [0.4 0.4 0.4]);

xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
view(3); axis equal; grid on;

for s = 1:length(unique_steps)
    stepID = unique_steps(s);

    % Get elements failed at current step
    failed_in_step = fail(fail(:,2) == stepID, 1);  

    % Add to cumulative list
    cumulative_failed = unique([cumulative_failed; failed_in_step]);

    % Plot new failed elements only
    for i = 1:length(failed_in_step)
        eNo = failed_in_step(i);
        nodes = ELEMCon(eNo, :);
        coords = COORD_deformed(nodes, :);
        patch('Faces', [1 2 3; 1 2 4; 1 3 4; 2 3 4], ...
              'Vertices', coords, ...
              'FaceColor', [1 0 0], ...
              'FaceAlpha', 1.0, ...
              'EdgeColor', 'none');
    end

    title(sprintf('Cumulative Failures up to Increment %d', stepID));
    drawnow;

    %  Capture and write video frame 
    frame = getframe(fig);
    writeVideo(v, frame);

    pause(0.5);  % Optional: slows down animation in real-time view
end

% Finish and close video 
close(v);
fprintf('Video saved to: %s\n', video_filename);
% -------------------------------------------------------------------------------------------%
%%
% ------------------------------- 3. Principal Stresses ------------------------------------- %
% To plot the actual principal stresses, we use ontained 'Act_P_Stress'
% To see the data of  'Act_P_Stress', Workspace >> ELEMENT >> Act_P_Stress
% Export Act_P_Stress to 'MeshG2_Postprocessing.m' to plot the contour stresses.
