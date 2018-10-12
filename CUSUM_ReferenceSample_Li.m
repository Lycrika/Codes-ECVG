%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function simulates the CUSUM Reference Sample control chart by Li et al. (2010).
% Input variables and parameters:
% X = Reference Sample;
% Y = Monitoring Sample
% n = size(X);
% m = size(Y);
% H = Upper Control Limit;
% k = variable of CUSUM
% ARL = Performance desirable
% caso = [0 1 2 3 4 5 6 ] = ''UNWGPLt'' = distribution to use (Uniform).
% sigma2 = 1; There is a change in variance? - Yes, then sigma2 !=1.
% mu2 =  [0 .25 .5 1 1.5 2 3]; Change in mean.
% [RL, var] are counters to calculate the ARL and SDRL.
%
% Other variables:
% k = Cusum parameter;
% % % Par_For = {0, 1} = {NO, YES}
% 
% Output variables:
% [RL, var] are counters to calculate the ARL and SDRL.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function  [RL, var, count] = CUSUM_ReferenceSample_Li(X,n,m,H,k,ARL,caso,sigma2,mu2,RL,var)
% close all, hold on%,clear all
%%% BORRAR
% caso = 1; n = 100; m = 5; mu2 = 0; sigma2 = 1; h = 5.071; ARL = 500; RL = 0;     var = 0;
% X = choose_distribution(caso,n);
% H = h*sqrt(m*n*(n+m+1)/12);
% k = .5*sqrt(m*n*(n+m+1)/12);
count = 0;
Cplus = 0; Cminus = 0;     
% subplot(2,1,1); hold on
% plot(-n+1:0,X,'.b')
while (Cplus < H && Cminus < H )  && count <= ARL*10
    y = choose_distribution(caso,m)*(sigma2)+mu2;
    Y = [X(1:n), y(1:m)];         
    [Cplus,Cminus] = CUSUM_wrs(Y(1:n+m),k,n,m,Cplus,Cminus);
    count = count+1;
%     subplot(2,1,1); hold on
%     plot(count*ones(1,m),y,'-.k')
%     subplot(2,1,2); hold on
%     plot(count,Cplus,'*b')
%     plot(count,Cminus,'*g')
%     pause(.01)
end % while             % [count, RL, R, r]
%     plot(1:count, ones(1,count)*H, '-r')
% % % disp(count)
RL = count + RL;
var = count^2 + var;
end

% % % % parfor r=1:rep  %%% Parallel toolbox
% % % % [RL, var] = CUSUM_ReferenceSample_Li(X,n,m,H,k,ARL,caso,sigma2,mu2,RL,var)
% % % % end % parfor(rep) ARL;
% % % % ARL(R) = RL/rep;% sum(RL)/rep;        % Son importantes todos los datos
% % % % SDRL(R) = sqrt((var - rep*(ARL(R)^2))/(rep-1));
