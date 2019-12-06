function [ Q_all ] = observation_wells( A,b_,freeNode_,u0_,flux_obtain_matrix,pressure_obtain_matrix,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac,chosen_confs,ln_d,data_generator )

%% MATERIALS DEPENDING ON d
no_fractures = length(lengths);
D = cell(no_fractures,1);
for i=1:no_fractures
    D{i} = generate_frac_fin( ln_d, data_generator );
    %figure(8); plot(D{i}); set(gca,'YScale','log'); hold on; figure(5); hold on
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
for i=chosen_confs
    u0 = u0_{i};
    b = b_{i};
    freeNode= freeNode_{i};
    
    [I,M]=interaction_matrix( fracture_matrice,node );
    [F,b_f,u_f,freenode_f] = fractures_matrix( node,fracture_matrice,intersections,alfa_inter,lengths);
    [B, freeNode, b, u0 ] = matrices_assembling( A, I, M, F, freeNode, freenode_f, b, b_f, u0, u_f );
    u0(freeNode)=B(freeNode,freeNode)\b(freeNode);

    %% extract flow
    FLUX=flux_obtain_matrix* u0(1:length(A));
    PRES=pressure_obtain_matrix* u0(1:length(A));
    Q_all = [Q_all; FLUX(1:end-1); PRES];
end
end

