function y = logistic(b,x)
    y = b(1) ./ (1 + b(2)* exp(-b(3)*(x+b(4) )))+b(5);   
end