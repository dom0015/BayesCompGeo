[data_generator, ~] = set_fracture(10,[linspace(0,1,81)' zeros(81,1)]);
basis=data_generator.gen_basis;
temp=basis*MH_SAMPLES(:,2:end)'+data_generator.prior_mean;
temp2=(basis*randn(10,1e6)/2+data_generator.prior_mean)*2+randn(1,1e6)-6;
k_cela=2*temp'+MH_SAMPLES(:,1);
hi=[];
hi2=[];
hi3=[];
for i=1:80
    hi(:,i)=my_histw(temp(i,:)',MH_MULTIPLICITY,[-6,0],60);
    hi2(:,i)=my_histw(k_cela(:,i),MH_MULTIPLICITY,[-17,-9],60);
    hi3(:,i)=my_histw(temp2(i,:)',ones(1e6,1),[-17,-9],60);
end
figure; imagesc([0 1],[-6 0],hi)
set(gca,'YDir','normal')

hold on

% load('uloha2.mat')
plot(linspace(0,1,80),log10(D{1}),'r');
plot(linspace(0,1,80),basis*prior_mean(2:end)+data_generator.prior_mean,'k')

figure;
plot(linspace(-8,-4,100),my_histw(MH_SAMPLES(:,1),MH_MULTIPLICITY,[-8,-4],100))

figure; imagesc([0 1],[-17,-9],hi2)
set(gca,'YDir','normal')

hold on

% load('uloha2.mat')
plot(linspace(0,1,80),log10(D{1})*2-7,'r');
plot(linspace(0,1,80),data_generator.prior_mean*2-6,'k')

figure;
plot(linspace(-8,-4,100),my_histw(MH_SAMPLES(:,1),MH_MULTIPLICITY,[-8,-4],100))

figure; imagesc([0 1],[-17,-9],hi3)
set(gca,'YDir','normal')

hold on

% load('uloha2.mat')
plot(linspace(0,1,80),log10(D{1})*2-7,'r');
plot(linspace(0,1,80),data_generator.prior_mean*2-6,'k')
