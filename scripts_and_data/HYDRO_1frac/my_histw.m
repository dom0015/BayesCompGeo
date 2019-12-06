function [H,x_axis] = my_histw(x,w,x_lim,x_bins)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
x1=x(:,1);


xx1=(x1-x_lim(1))/(x_lim(2)-x_lim(1))*x_bins;


xx1=min(max(ceil(xx1),1),x_bins);


H=full(sparse(xx1,xx1*0+1,w,x_bins,1));
H=H/(sum(H(:)*(x_lim(2)-x_lim(1))/x_bins));
x_axis=linspace(x_lim(1),x_lim(2),x_bins+1);
x_axis=(x_axis(1:end-1)+x_axis(2:end))/2;
end