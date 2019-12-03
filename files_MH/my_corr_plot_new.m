function [] = my_corr_plot_new(SAMPLES,WEIGHTS,u_real,prior_mu,prior_std,type,nibs_hist,nibs_2dhist)
P_quantiles=[0.001,0.999];


if nargin < 8
    nibs_2dhist=[];
end
if nargin < 7
    nibs_hist=[];
end
if nargin < 6
    type=[];
end
if nargin < 5
    prior_mu=[];
    prior_std=[];
end
if nargin < 3
    u_real=[];
end

if isempty(nibs_hist)
    nibs_hist=30;
end
if isempty(nibs_2dhist)
    nibs_2dhist=50;
end
if isempty(type)
    type=2;%-'common-size';% 1-'single';%3-'common',
end

n=size(SAMPLES,2);
histw=cell(n,1);
vinterval=cell(n,1);
maxmin=zeros(n,2);

for i=1:n
    [maxmin(i,:)]=weighted_quantiles(SAMPLES(:,i),WEIGHTS,P_quantiles);
end

if type == 1
    maxmin_scaled=maxmin;
end

if type == 2
    max_size=max(maxmin(:,2)-maxmin(:,1))/2;
    maxmin_scaled= [(maxmin(:,1)+maxmin(:,2))/2 - max_size (maxmin(:,1)+maxmin(:,2))/2 + max_size];
end

if type == 3
    maxmin_scaled= repmat([min(maxmin(:)) max(maxmin(:))],4,1);
end

for i=1:n
    [histw{i}, vinterval{i}] = my_histw_cut(SAMPLES(:,i),WEIGHTS,maxmin(i,:),nibs_hist);
end


figure('Renderer', 'painters', 'Position', [10 10 960 960])
axes_all=cell(n);
for i=1:n
    for j=1:n
        if i==j
            continue
        end
        [H,x_axis,y_axis] = my_hist2w_cut(SAMPLES(:,[i j]),WEIGHTS,maxmin_scaled(i,:),maxmin_scaled(j,:),nibs_2dhist,nibs_2dhist);
        axes_all{i,j}=subplot(n,n,(i-1)*n+j);
        hold on;
        imagesc(y_axis,x_axis,H);
        if ~isempty(prior_mu)
            prior_viz([prior_mu(j),prior_mu(i)],[prior_std(j),prior_std(i)]);
        end
        imagesc(y_axis,x_axis,H)
        alpha(0.9)
        
        if ~isempty(u_real)
            plot(u_real(j),u_real(i),'r.','MarkerSize',15)
        end
        xlim(maxmin_scaled(j,:))
        ylim(maxmin_scaled(i,:))
        grid on
        box on
    end
end

hist_max_y=zeros(n,1);
for i=1:n
    hist_max_y(i)=max(histw{i});
end
hist_max_y=max(hist_max_y);

for i=1:n
    j=i;
    axes_all{i,j}=subplot(n,n,(i-1)*n+j);
    bar(vinterval{i}, histw{i},'FaceColor',[0 .2 .6],'EdgeColor',[0 .25 .7],'LineWidth',0.5);
    hold on;
    
    if ~isempty(prior_mu)
        x_pri=linspace(prior_mu(i)-4*prior_std(i),prior_mu(i)+4*prior_std(i),1e3);
        y_pri=normpdf(x_pri,prior_mu(i),prior_std(i));
        plot(x_pri,y_pri,'m-','LineWidth',2)
    end
    
    if ~isempty(u_real)
        plot(u_real(i),0,'r.','MarkerSize',15)
    end
    xlim(maxmin_scaled(i,:))
    ylim([0 hist_max_y])
    grid on
    box on
end

so=0.07;
eo=0.07;
mo=0.005;
ol=(1-so-eo-(n-1)*mo)/n;
startpos=so:(ol+mo):(1-eo-ol);
poss=[zeros(n^2,2) ones(n^2,2)*ol];
for i=1:n
    for j=1:n
        poss((i-1)*n+j,1)=startpos(j);
        poss((i-1)*n+j,2)=startpos(n-i+1);
    end
end


for i=1:n
    for j=1:n
        if i==n
            set(axes_all{i,j},'XAxisLocation','bottom');
        end
        if i==1
            set(axes_all{i,j},'XAxisLocation','top');
        end
        if j==1
            set(axes_all{i,j},'YAxisLocation','left');
        end
        if j==n
            set(axes_all{i,j},'YAxisLocation','right');
        end
        
        set(axes_all{i,j},'YDir','normal')
        if i>1&&i<n
            set(axes_all{i,j},'xticklabel',[])
        end
        if j<n&&j>1
            set(axes_all{i,j},'yticklabel',[])
        end
        set(axes_all{i,j}, 'Position', poss((i-1)*n+j,:))
        ax = axes_all{i,j};
        ax.LineWidth = 1;
        ax.GridLineStyle = '-';
        ax.GridColor = 'k';
        if i~=j
            ax.GridAlpha = 1; % maximum line opacity
        end
    end
end


set(findall(gcf,'-property','FontSize'),'FontSize',13)
end

function [Q]=weighted_quantiles(X,W,P)
[X,idx]=sort(X);
W=W(idx);
W=W/sum(W);
W_F=[0; cumsum(W)];

Q_i=P;
for i=1:length(P)
    Q_i(i)=find(W_F>=P(i),1);
end
Q=X(Q_i-1);
end

function [  ] = prior_viz( mu0, gamma )
n_std=3; % pocet smerodatnych odchylek od prioru

x_bounds=[mu0(1)-n_std*gamma(1),mu0(1)+n_std*gamma(1)];
y_bounds=[mu0(2)-n_std*gamma(2),mu0(2)+n_std*gamma(2)];

line(x_bounds,mu0(2)*[1 1],'Color','white')
line(mu0(1)*[1 1],y_bounds,'Color','white')

t=linspace(0,2*pi,1000);
x_0=cos(t);
y_0=sin(t);
plot(mu0(1)+1*gamma(1)*x_0,mu0(2)+1*gamma(2)*y_0,'w','LineWidth',4)
plot(mu0(1)+2*gamma(1)*x_0,mu0(2)+2*gamma(2)*y_0,'w','LineWidth',3)
plot(mu0(1)+3*gamma(1)*x_0,mu0(2)+3*gamma(2)*y_0,'w','LineWidth',2)
end

function [H,x_axis] = my_histw_cut(x,w,x_lim,x_bins)
xx1=(x-x_lim(1))/(x_lim(2)-x_lim(1))*x_bins;
w(xx1<0|xx1>x_bins)=0;
xx1=min(max(ceil(xx1),1),x_bins);
H=full(sparse(xx1,xx1*0+1,w,x_bins,1));
H=H/(sum(H(:)*(x_lim(2)-x_lim(1))/x_bins));
x_axis=linspace(x_lim(1),x_lim(2),x_bins+1);
x_axis=(x_axis(1:end-1)+x_axis(2:end))/2;
end

function [H,x_axis,y_axis] = my_hist2w_cut(x,w,x_lim,y_lim,x_bins,y_bins)
xx1=(x(:,1)-x_lim(1))/(x_lim(2)-x_lim(1))*x_bins;
xx2=(x(:,2)-y_lim(1))/(y_lim(2)-y_lim(1))*y_bins;
w(xx1<0|xx1>x_bins)=0;
w(xx2<0|xx2>y_bins)=0;
xx1=min(max(ceil(xx1),1),x_bins);
xx2=min(max(ceil(xx2),1),y_bins);
H=full(sparse(xx1,xx2,w,x_bins,y_bins));
H=H/(sum(H(:)*(x_lim(2)-x_lim(1))/x_bins*(y_lim(2)-y_lim(1))/y_bins));
x_axis=linspace(x_lim(1),x_lim(2),x_bins+1);
x_axis=(x_axis(1:end-1)+x_axis(2:end))/2;
y_axis=linspace(y_lim(1),y_lim(2),y_bins+1);
y_axis=(y_axis(1:end-1)+y_axis(2:end))/2;
end