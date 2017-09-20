function [ area, b ] = myfunc( x, method )

% x: is a list of lengths
% method: a string that specifies the shape

% area: the results
% b: a second output here for demo purposes

area = nan;

b = 2;

if strcmpi(method,'circle')
    
    area = pi*x.^2;
    
elseif strcmpi(method, 'square')
    
    area = x.^2;
    
else
    
    error('The requested method <%s> is not yet implemented\n', method);
    
end

end

