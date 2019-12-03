function [ HIST ] = my_hist2d(SAMPLES,MULTIPLICITY)

n_x=60;
n_y=60;

x_min=min(SAMPLES(:,1));
x_max=max(SAMPLES(:,1));
y_min=min(SAMPLES(:,2));
y_max=max(SAMPLES(:,2));

xx=linspace(x_min,x_max,n_x);
yy=linspace(y_min,y_max,n_y);

HIST=zeros(n_x,n_y);

for i=1:size(SAMPLES,1)
    sample_x=SAMPLES(i,1);
    ind_x=max(find(xx<=sample_x));
    sample_y=SAMPLES(i,2);
    ind_y=max(find(yy<=sample_y));
    HIST(ind_x,ind_y)=HIST(ind_x,ind_y)+MULTIPLICITY(i);
end

imagesc([x_min x_max],[y_min y_max],HIST)
set(gca,'YDir','normal')