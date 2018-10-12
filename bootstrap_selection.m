function y = bootstrap_selection(x, n, m)
    y = x (randi(n,1,m));
end


% % % % 
% % % % 
% % % % %%%%%%%%%%%%%%%%% Refreshed Parameters of the simulation %%%%%%%%%%%%%%%%%%
% % % % h(1) = CalibratedCL(n,m); 
% % % % h(2) = 3.27; % h = 5.071; ARL = 500;
% % % % H = h(1)*sqrt(m*n*(n+m+1)/12); k = .5*sqrt(m*n*(n+m+1)/12);
% % % %           % arl: initial value
% % % % RL = 0;     var = 0;    arl = 0;  
% % % % 
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the arl %%%%%%%%%%%%%%%%%%%%%%%%%
% % % % for r=1:rep  %%% Parallel toolbox
% % % % % % % %     if NRS == 1
% % % % % % % %         X = choose_distribution(caso,n);
% % % % % % % %     end
% % % %     [RL, var] = CUSUM_ReferenceSample_Li(X,n,m,H,k,ARL(1),caso,s2,mu2,RL,var);
% % % % %         rl(r) = count;                    %Verificar que arl = mean(rl)
% % % % end % parfor(rep) ARL;
% % % % %     rl' %     mean(rl)%     sdr= std(rl)  %Verificar que arl = mean(rl)
% % % % arl(1) = RL/rep; % sum(RL)/rep;        % Son importantes todos los datos
% % % % sdrl(1) = sqrt((var - rep*(arl^2))/(rep-1));
% % % % fprintf('h = %.2f,\t arl = %.2f(%d),\t sdrl = %.2f, \t n = %d, \t m = %d, \t Dist = %s \n',h(1), arl(1),ARL(1), sdrl(1),n,m,Ldist)
% % % % 
% % % % aux = 1;
% % % % 
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Verify |ARL-arl| < error 
% % % % while abs(ARL(1)-arl(aux))>= error  
% % % % % %     WHILE: difference between ARLs is more significant than error.
% % % % % %     THEN: Calibrate
% % % %     RL = 0;     var = 0;
% % % %     X = choose_distribution(caso,n);
% % % %     aux = aux+1;
% % % %     H = h(aux)*sqrt(m*n*(n+m+1)/12);
% % % %     for r=1:rep  %%% Parallel toolbox
% % % %         if NRS == 1
% % % %             X = choose_distribution(caso,n);
% % % %         end
% % % %         [RL, var] = CUSUM_ReferenceSample_Li(X,n,m,H,k,1000,caso,s2,mu2,RL,var);
% % % % %         rl(r) = count;                    %Verificar que arl = mean(rl)
% % % %     end % parfor(rep) ARL;
% % % % %     rl' %     mean(rl)%     sdr= std(rl)  %Verificar que arl = mean(rl)
% % % %     arl(aux) = RL/rep; % sum(RL)/rep;        % Son importantes todos los datos
% % % %     sdrl(aux) = sqrt((var - rep*(arl(aux)^2))/(rep-1));
% % % %     fprintf('h = %.2f,\t arl = %.2f(.95*%d),\t sdrl = %.2f, \t n = %d, \t m = %d, \t Dist = %s \n',h(aux), arl(aux),ARL(1), sdrl(aux),n,m,Ldist)
% % % % 
% % % % %     aux = 2;
% % % %     h(aux+1) = inter_extra_polation([h(aux-1) h(aux)],[arl(aux-1) arl(aux)],ARL(1));
% % % % end
% % % % 
