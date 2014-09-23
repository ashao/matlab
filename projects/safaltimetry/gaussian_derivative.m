function [ out J ] = gaussian_derivative(params,x)
   
    x0 = params(1);
    sigma = params(2);
    amp = params(3);
    exponential = exp(-(x-x0).^2/sigma^2*0.5);
    out =  -(x-x0)/sigma^2 .* exponential*amp;

    if nargout > 1
       J(:,3) = exponential;
       J(:,2) = (x-x0).^2/sigma^3 .* exponential;
       J(:,1) = (x-x0)./sigma^2 .* exponential;
    end

    
end