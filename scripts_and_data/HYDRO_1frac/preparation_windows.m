function [ A,b,freeNode,u0,downEdge,rightEdge,upEdge,leftEdge,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac,fractures_positions ] = preparation_windows(n_windows)
%% FEM with fractures
addpath(genpath('files'))
addpath(genpath('model_problems'))

%% GEOMETRY PARAMETERS
h_elem=100;
frac_start_end={[0.2 0.05],[0.2 0.95]};%{[0.05 0.15], [0.8 0.9]};
            
%% GEOMETRY ASSEMBLING
[node,elem,bdFlag]=rect_mesh(10,10,h_elem,h_elem); % triangulace
[fractures, fractures_positions, no_fractures] = create_fractures( frac_start_end, node, h_elem );
[fractures_cell,fracture_matrice,intersections,lengths] = fracture2cells_geometry( fractures );
no_intersections = size(intersections,1);
[ node,elem ,bdFlag,fractures_cell,fracture_matrice] = multi_fracture_tear( node,elem,fractures_cell ,bdFlag,fracture_matrice);

%% MATERIAL AND BOUNDARY PARAMETERS
p=1e6; % pressure on Dir. b. c.
f = @(x)(0+0*x(:,1)+0*x(:,2)); % zatizeni
g_N=@(x)(0+0*x(:,1)+0*x(:,2)); % Neumann
mat_omega = 1e-15; % material - matrice
k = mat_omega*ones(length(elem),1);

mat_frac = 1e-9*ones(no_fractures,1); % material - fractures
alfa_inter = 1e-8*ones(no_intersections,1);


no_types = 4;
u0 = cell(no_types,1);
b = cell(no_types,1);
freeNode = cell(no_types,1);
downEdge = cell(no_types,1);
rightEdge = cell(no_types,1);
upEdge = cell(no_types,1);
leftEdge = cell(no_types,1);
ee = eps;
for cislo_ulohy = 1:no_types
    Dirichlet_windows = D_windows( cislo_ulohy,p,n_windows,ee );
    Dirichlet_windows(:,2:3)=Dirichlet_windows(:,2:3)*10;
    [ u0_, A, b_, freeNode_, downEdge_, rightEdge_, upEdge_, leftEdge_ ] = FEM_windows( node, elem, h_elem, bdFlag, k, f, Dirichlet_windows, g_N );
    u0{cislo_ulohy} = u0_;
    b{cislo_ulohy} = b_;
    freeNode{cislo_ulohy} = freeNode_;
    downEdge{cislo_ulohy} = downEdge_;
    rightEdge{cislo_ulohy} = rightEdge_;
    upEdge{cislo_ulohy} = upEdge_;
    leftEdge{cislo_ulohy} = leftEdge_;
end

end

