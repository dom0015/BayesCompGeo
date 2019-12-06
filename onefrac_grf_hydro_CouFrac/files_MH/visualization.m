function [  ] = visualization( mu0, gamma, u_real, G, n )
%VISUALIZATION Summary of this function goes here
%   Detailed explanation goes here
n_std=2; % pocet smerodatnych odchylek od prioru
%figure;
x_bounds=[mu0(1)-n_std*gamma,mu0(1)+n_std*gamma];
y_bounds=[mu0(2)-n_std*gamma,mu0(2)+n_std*gamma];

if nargin>3
    xx=linspace(x_bounds(1),x_bounds(2),n);
    yy=linspace(y_bounds(1),y_bounds(2),n);
    X=zeros(n);
    for i=1:n
        for j=1:n
            X(i,j)=G([xx(i) yy(j)]');
        end
    end
    imagesc(xx,yy,X);
end
hold on
lw=1.5;
line([u_real(1) u_real(1)],y_bounds,'Color','c','LineWidth',lw)
line(x_bounds,[u_real(2) u_real(2)],'Color','c','LineWidth',lw)
line(x_bounds,mu0(2)*[1 1],'Color','m','LineWidth',lw)
line(mu0(1)*[1 1],y_bounds,'Color','m','LineWidth',lw)
r=gamma; t = linspace(0,2*pi,1000);plot(mu0(1)+r*cos(t),mu0(2)+r*sin(t),'m','LineWidth',lw)

end

