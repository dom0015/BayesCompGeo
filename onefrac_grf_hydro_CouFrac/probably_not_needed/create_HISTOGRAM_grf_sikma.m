%% Bayesian solution of an inverse problem
addpath(genpath('files_MH'))
addpath(genpath('files'))

s=10;                % number of random parameters (length of u)

%% GEOMETRY PARAMETERS
h_elem=100;
frac_start_end={[0.05 0.15], [0.8 0.9]};
            
%% GEOMETRY ASSEMBLING
[node,elem,bdFlag]=rect_mesh(10,10,h_elem,h_elem); % triangulace
[fractures, fracture_positions, no_fractures] = create_fractures( frac_start_end, node, h_elem );


[data_generator, ~] = set_fracture(s,fracture_positions{1} );

lS=length(SAMPLES);
lF=-1+length(fracture_positions{1});
FW=zeros(lS,lF);
for i=1:lS
	params=SAMPLES(i,:)';
	FW(i,:) = generate_frac_fin( params, data_generator );
end


SAMPLES_size=size(SAMPLES);
midpoints=data_generator.frac_midpoints;



logFW=log10(FW);

min_logFW=min(logFW(:));
max_logFW=max(logFW(:));

lF=size(logFW,2);


steps=100;
HISTS=zeros(steps,lF);
Quantiles=zeros(lF,7);
bins=linspace(min_logFW,max_logFW,steps);
for i=1:lF
	HH = hist(logFW(:,i),bins);
	HISTS(:,i)=HH';
	for j=-3:3
	Quantiles(i,j+4)=quantile(logFW(:,i),normcdf(j));
	end
end

figure; imagesc(midpoints,bins,HISTS)
set(gca,'YDir','normal')

CO=corr(logFW);
figure; imagesc(CO)

%label='PRIOR_windows_grf_sikma'
save(['GRFHIST_' label '.mat'],'HISTS','SAMPLES_size','bins','Quantiles','label','midpoints','CO')
