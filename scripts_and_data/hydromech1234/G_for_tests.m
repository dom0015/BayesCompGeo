function [outputArg1,outputArg2] = G_for_tests(x)
%G_FOR_TESTS Summary of this function goes here
%   Detailed explanation goes here
o = x(1)+x(2)*x(2)/(-6);
o = sum(x);
outputArg1 = [o o o o]';
outputArg2 = [sum(x) sum(x) sum(x) sum(x)]';
end

