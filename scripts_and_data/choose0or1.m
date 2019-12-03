function [ filename ] = choose0or1( name, i, folder_name, label )
%CHOOSE0OR1 Summary of this function goes here
%   Detailed explanation goes here
    cd([folder_name '/res'])
    t0=dir([name '_' label '_' num2str(i) '_0.mat']);
    t1=dir([name '_' label '_' num2str(i) '_1.mat']);
    cd('..');
    cd('..');
    if isempty(t0)
        filename=[folder_name '/res/' name '_' label '_' num2str(i) '_1.mat'];
        disp([name '1'])
    elseif t0.bytes > t1.bytes 
        filename=[folder_name '/res/' name '_' label '_' num2str(i) '_0.mat'];
        disp([name '0'])
    else
        filename=[folder_name '/res/' name '_' label '_' num2str(i) '_1.mat'];
        disp([name '1'])
    end

end
