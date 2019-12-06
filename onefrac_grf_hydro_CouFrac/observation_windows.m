function [ Q_all ] = observation_windows( A,b_,freeNode_,u0_,downEdge_,rightEdge_,upEdge_,leftEdge_,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac,chosen_windows,ln_d,data_generator,n_windows,visual )

%% MATERIALS DEPENDING ON d
no_fractures = length(lengths);
D = cell(no_fractures,1);
idx=0;
for i=1:no_fractures
    si=length(data_generator{i}.random_coef);
    D{i} = generate_frac_fin( ln_d(idx+1:idx+si), data_generator{i} );
    idx=idx+si;
    %figure(8); plot(D{i}); set(gca,'YScale','log'); hold on; figure(5); hold on
end

ALFA = cell(no_fractures,1);
MAT_FRAC = cell(no_fractures,1);
for i=1:no_fractures
    MAT_FRAC{i} = mat_frac(i).*D{i}.*D{i};
    ALFA{i} = mat_frac(i)+0*D{i};
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

    if visual
        figure;
        h = trisurf(elem,node(:,1),node(:,2),u0(1:length(A)));
        h.EdgeColor = 'none';
        colormap jet(1000)
        colorbar
        axis equal
        view(0,90);
        figure(10+i); hold on
        for j=1:no_fractures
            plot(D{j}); 
        end
    end
    %% extract flow
    ee=eps;
    Dirichlet_windows = D_windows( i,p,n_windows,ee );
    Dirichlet_windows(:,2:3)=Dirichlet_windows(:,2:3)*10;
    Q = extract_flow( A, u0(1:length(A)), node, elem, h_elem, Dirichlet_windows, downEdge, rightEdge, upEdge, leftEdge );
    Q_all = [Q_all Q];
end
Q_all = Q_all';
end

