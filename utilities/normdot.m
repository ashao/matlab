function [ dot ] = normdot(x,y)

x = x/norm(x);
y = y/norm(y);
dot = x'*y;