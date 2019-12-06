function [ data_generator, realny_material] = set_fracture( n,fracture_positions,intersections )
%SET_FRACTURE Summary of this function goes here
%   Detailed explanation goes here
random_coef=[
   -2.3234
   -0.9789
    1.5875
   -0.3528
    0.8644
   -0.4754
    0.8674
   -0.2757
    0.5039
   -0.1675];
random_coef=random_coef(1:n);
%random_coef=[random_coef; 0.1; -0.2];

sigma=1;
lambda=0.5;
[ eig_value,eig_func ] = KL_1D( n,sigma,lambda );

f=@(x)10.^(-2-2*(abs(2*(x-0.5))).^2); % -3

tmp=fracture_positions(1:end-1,:)-fracture_positions(2:end,:);
seg_len=sqrt(tmp(:,1).^2+tmp(:,2).^2);
seg_positions=cumsum(seg_len);
scale_factor=seg_positions(end);
seg_positions=[0;seg_positions/seg_positions(end)];
frac_midpoints=(seg_positions(1:end-1)+seg_positions(2:end))/2;

prior_mean=f(frac_midpoints);

no_midpoints=length(frac_midpoints);
gen_basis=zeros(no_midpoints,n);

for i=1:n
    gen_basis(:,i)=sqrt(eig_value(i))*eig_func{i}(frac_midpoints);
end

no_inter=length(intersections);
start=1;
for i=1:no_inter
    fin=intersections(i);
    if fin>1
        temp=zeros(no_midpoints,1);
        temp(start:fin)=1;
        start=fin+1;
        gen_basis=[gen_basis temp];
    end
    if i==no_inter
        temp=zeros(no_midpoints,1);
        temp(start:no_midpoints)=1;
        gen_basis=[gen_basis temp];
    end
end

data_generator.gen_basis=gen_basis;
data_generator.prior_mean=prior_mean;
data_generator.frac_midpoints=frac_midpoints*scale_factor;
data_generator.random_coef=random_coef;
realny_material= sample_fracture( random_coef,frac_midpoints );
realny_material=prior_mean.*realny_material;
end

