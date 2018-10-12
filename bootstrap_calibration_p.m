%%%%%%%%%%% H_calibration is the original function. 

function h_arl = bootstrap_calibration_p(X, n, m, ARL, rep, s2, mu2, k)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Experimental data %%%%%%%%%%%%%%%%%%%%%%%%%
% clear all; clc
% caso = 0; n = 100; m = 5; X = choose_distribution(caso,n); 
% ARL = 100; rep = 1000; s2 = 1; mu2 = 0; k = .5*sqrt(m*n*(n+m+1)/12);
%%%%%%%%%%%%%%%%% Refreshed Parameters of the simulation %%%%%%%%%%%%%%%%%%
h = ones(1,10)*CalibratedCL(n,m); 
hh = 3.27; %(CalibratedCL(n,m)+3.27)/2;
H = h(1)*sqrt(m*n*(n+m+1)/12); 
%%%%%%%%%%%%%%%%%%%%%%%%%%% arl: initial value 
RL = 0;                 var = 0;       
arl = zeros(1,10);      sdrl = zeros(1,10);
%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the arl with rep, replicates. %%%%%%%%%
for r =1:rep    %%% Calculate the first arl to compare vs ARL(1);
    [RL, var] = bootstrap_CUSUM_ReferenceSample_Li(X,n,m,H,k,s2,mu2,RL,var,ARL);
end
arl(1) = RL/rep;          
sdrl(1) = sqrt((var - rep*(arl(1)^2))/(rep-1));
% fprintf('\nh = %.2f,\t arl = %.2f(%d),\t sdrl = %.2f, \t n = %d, \t m = %d.\n',h(1), arl(1),ARL(1), sdrl(1),n,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Verify |ARL-arl| < error 
aux = 1;        error = .05*ARL(1);       b = .05;       a = .05; 
while abs(ARL(1)-arl(aux))>= error  % error = .05*ARL(1);
    p = ARL/arl(aux);
    aux = aux +1;
    if p <= 1-a
        p = 1-a;
    elseif p >= 1+a
        p = 1+a;
    end
    h(aux) = p*h(aux-1);
    if h(aux) <= (1-b)*min(h(1),hh)
        h(aux) = (1-b)*min(h(1),hh);
    elseif h(aux) >= (1+b)*max(h(1),hh)
        h(aux) = (1+b)*max(h(1),hh);
    end
%%%% WHILE: difference between ARLs is more significant than error. THEN: Calibrate
%     if h(aux) == h(aux+1)       %%% In case, the h(1) = h(2) 
%         if ARL < arl(aux)
%             h(aux+1) = h(aux)*.95;
%         else
%             h(aux+1) = h(aux)*1.05;
%         end
%     end
    RL = 0;     var = 0;                %%% arl: initial value                        %%% Next iteration
    H = h(aux)*sqrt(m*n*(n+m+1)/12);
    for r =1:rep 
        [RL, var] = bootstrap_CUSUM_ReferenceSample_Li(X,n,m,H,k,s2,mu2,RL,var,ARL);
    end 
    arl(aux) = RL/rep; 
    sdrl(aux) = sqrt((var - rep*(arl(aux)^2))/(rep-1));
%     fprintf('h = %.2f,\t arl = %.2f(%d),\t sdrl = %.2f, \t n = %d, \t m = %d. \n',h(aux), arl(aux),ARL(1), sdrl(aux),n,m)
%     if COUNT < 100
%         h(aux+1) = inter_extra_polation([h(aux-1) h(aux)],[arl(aux-1) arl(aux)],ARL(1));
%     end
    if aux > 10 || h(aux)<= 0 || h(aux)>= 2*h(1)
        aux = 1; fprintf('Too many calibarion iterations.\n')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Conclude with the correspondent h %%%%%%%%%%%%
h_arl = h(aux);
end