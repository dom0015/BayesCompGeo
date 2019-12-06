function [ Q_all ] = observation_windows( A,b_,freeNode_,u0_,downEdge_,rightEdge_,upEdge_,leftEdge_,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac,chosen_windows,ln_d,n_windows )

%% MATERIALS DEPENDING ON d
d = exp(ln_d);
no_fractures = length(lengths);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = d(i)*ones(lengths(i)-1,1);
end

ALFA = cell(no_fractures,1);
MAT_FRAC = cell(no_fractures,1);
for i=1:no_fractures
    MAT_FRAC{i} = mat_frac(i).*D{i};
    ALFA{i} = mat_frac(i)./D{i};
end

[fracture_matrice] = fracture2cells_parameters( fracture_matrice,ALFA,MAT_FRAC );

%% RESULTS
Q_all = [];
for i=chosen_windows
    u0 = u0_{i};
    b = b_{i};
    freeNode= freeNode_{i};
    downEdge = downEdge_{i};
    rightEdge = rightEdge_{i};
    upEdge = upEdge_{i};
    leftEdge = leftEdge_{i};
    
    [I,M]=interaction_matrix( fracture_matrice,node );
    [F,b_f,u_f,freenode_f] = fractures_matrix( node,fracture_matrice,intersections,alfa_inter,lengths);
    [B, freeNode, b, u0 ] = matrices_assembling( A, I, M, F, freeNode, freenode_f, b, b_f, u0, u_f );
    u0(freeNode)=B(freeNode,freeNode)\b(freeNode);

    %% extract flow
    ee=eps;
    Dirichlet_windows = D_windows( i,p,n_windows,ee );
    Dirichlet_windows(:,2:3)=Dirichlet_windows(:,2:3)*10;
    Q = extract_flow( A, u0(1:length(A)), node, elem, h_elem, Dirichlet_windows, downEdge, rightEdge, upEdge, leftEdge );
    Q_all = [Q_all Q];
end
Q_all = Q_all';
end

