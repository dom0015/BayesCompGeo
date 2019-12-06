
chosen_windows = [1 2 3 4];
%% observation operator G
[ A,b,freeNode,u0,downEdge,rightEdge,upEdge,leftEdge,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac ] = preparation_windows();
G=@(x)observation_windows( A,b,freeNode,u0,downEdge,rightEdge,upEdge,leftEdge,fracture_matrice,intersections,alfa_inter,lengths,node,elem,h_elem,p,mat_frac,chosen_windows,x );

u=[-4 -9];
disp(G(u)')