function x = choose_distribution(caso,n)
switch caso      
    case 0                  % Unifform(0,1)
        x = 2*(rand(1,n)-.5)*sqrt(3);
    case 1                  % Normal(0,1)
        x = randn(1,n);
    case 2                  % Weibull
%         x =
    case 3                  % Gamma(3,1)
        x = gamrnd(3,1,1,n);
    case 4                  % Double Exponential
%         x =
    case 5                  % Lognormal
%         x =
    case 6                  % t(5)             
        x = trnd(5,1,n);
end

x = x(1:n);