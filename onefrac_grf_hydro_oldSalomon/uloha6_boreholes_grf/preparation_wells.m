function [ A,b,freeNode,u0,flux_obtain_matrix,pressure_obtain_matrix,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac,fractures_positions ] = preparation_wells()
%% FEM with fractures
addpath(genpath('files'))
addpath(genpath('model_problems'))

%% GEOMETRY PARAMETERS
h_elem=100;
frac_start_end={[0.1 0.2], [0.8 0.9]};
            
%% GEOMETRY ASSEMBLING
[node,elem,bdFlag]=rect_mesh(10,10,h_elem,h_elem); % triangulace
bdFlag(bdFlag>0)=10;
wells=[3 3;3 5;3 7;
    5 3;5 5;5 7;
    7 3;7 5;7 7;];

[fractures, fractures_positions, no_fractures] = create_fractures( frac_start_end, node, h_elem );
[ node,elem,bdFlag,fractures] = create_wells( node,elem,bdFlag,fractures,wells,0.1);
[fractures_cell,fracture_matrice,intersections,lengths] = fracture2cells_geometry( fractures );
no_intersections = size(intersections,1);
[ node,elem ,bdFlag,fractures_cell,fracture_matrice] = multi_fracture_tear( node,elem,fractures_cell ,bdFlag,fracture_matrice);
[ node ] = stretch_domain( node );

%% MATERIAL AND BOUNDARY PARAMETERS
p=1e6; % pressure on Dir. b. c.
mat_omega = 1e-15; % material - matrice
k = mat_omega*ones(length(elem),1);
mat_frac = 1e-9*ones(no_fractures,1); % material - fractures
alfa_inter = 1e-8*ones(no_intersections,1);

no_types = 4;
u0 = cell(no_types,1);
b = cell(no_types,1);
freeNode = cell(no_types,1);
wells_value_neumann=[ 
    -2; 
    -4;
    -5;
    -6;
    -8;];
for i=1:size(wells_value_neumann,1)
    bdFlag(bdFlag==-wells_value_neumann(i))=wells_value_neumann(i);
end

for cislo_ulohy = 1:no_types
    switch cislo_ulohy
        case 1
    wells_value=[
        1 p; 
        3 0;
        7 0;
        9 0;
        10 p*0.01];
        case 2
    wells_value=[
        1 0; 
        3 p;
        7 0;
        9 0;
        10 p*0.01];
        case 3
    wells_value=[
        1 0; 
        3 0;
        7 p;
        9 0;
        10 p*0.01];
        case 4
    wells_value=[
        1 0; 
        3 0;
        7 0;
        9 p;
        10 p*0.01];       
    end
    
    [flux_obtain_matrix,pressure_obtain_matrix,u0_, A, b_, freeNode_] = FEM_wells( node, elem, bdFlag, k,wells_value,wells_value_neumann);

    u0{cislo_ulohy} = u0_;
    b{cislo_ulohy} = b_;
    freeNode{cislo_ulohy} = freeNode_;
end

end

