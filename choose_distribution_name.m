function Ldist = choose_distribution_name(caso)
switch caso    
    case 0                  % Unifform(0,1)
        Ldist = 'U(0,1)';
    case 1                  % Normal(0,1)
        Ldist = 'N(0,1)';
    case 2                  % Weibull
        Ldist = 'W';
    case 3                  % Gamma(3,1)
        Ldist = 'G(3,1)';
    case 4                  % Double Exponential
        Ldist = '2E';
    case 5                  % Lognormal
        Ldist = 'LgN';
    case 6                  % t(5)             
        Ldist = 't(5)';
end
end