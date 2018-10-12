%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calibrates the control limits to achieve a prespecified ARL_0 (in control) performance. 
% The simulation to use is a Monte-Carlo simulation with N  = number of replications (N = 10^4). 
% The control limits are calibrated for a CUSUM Reference Sample control chart [Li et al (2010)]
% 
% Parameters to define a priori:
% ARL = [10, 200, 500, 1000] = Average Run-Length
% n = [50, 100, 200, 500, 1000] = size of reference sample.
% m = [1, 3, 5, 7] = size of the monitored sample.
% caso = [0 1 2 3 4 5 6 ] = ''UNWGPLt'' = distribution to use (Uniform).
% NRS = New Reference Sample = {No = 0, Yes = 1}
%
% New parameters:
% arl = ARL performance for every H.
% error = 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function H_calibration(ARL, n, m, caso)
clc, clear all, tic; Notas = {...
    'Esta es una simulación para automatizar la calibración de los límites de control'}
%%%%%%%%%%%%%%%%%%% To save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stat,struc] = fileattrib;
Carpeta = struc.Name; disp(Carpeta)
Fecha = clock; Fecha = num2str(Fecha(1:5),'%02d');
%%%%%%%%%%%%%%%%% FIXED Parameters of the simulation %%%%%%%%%%%%%%%%%%%%%%
CASOS = [0]; dist = 'UNWGPLt';  %% Selects the data distribution
no = [3000 2000];                     % Reference sample size         n = #(X)
mo = [10 7 5 3 1];                       % Batch sample size to monitor  m = #(Y)
mu2 = 0; s2 = 1;                    % \mu_2 and \sigma_2 => No change. 
NRS = 1; % New Reference Sample = {No = 0, Yes = 1}
ARL(1) = 100; error = .05*ARL(1);    % error = 5 percento of the ARL
rep = 100000;

for caso = CASOS        %%%%%%%%% DISTRIBUTION 
mm = numel(mo);
for im = 1:mm              %%%%%%%%% MONITORING SAMPLE SIZE 
    m = mo(im);     nn = numel(no);
for in = 1:nn                %%%%%%%%% REFERENCE SAMPLE SIZE
    n = no(in);

%%%%%%%%%%%%%%%%% Refreshed Parameters of the simulation %%%%%%%%%%%%%%%%%%
Ldist = dist(caso+1);
h(1) = 3.27; 
h(2) = 3.502; % h = 5.071; ARL = 500;
H = h(1)*sqrt(m*n*(n+m+1)/12); k = .5*sqrt(m*n*(n+m+1)/12);
          % arl: initial value
RL = 0;     var = 0;    arl = 0;  
X = choose_distribution(caso,n);
%%%%%%%%%%%%%%%%%%% To save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Archivo = strcat('P1_H_calibration_ARL(',num2str(ARL(1)),')rep(',num2str(rep),')_',Fecha(1:4),...
        Fecha(7:8),Fecha(11:12),'(',Fecha(15:16),Fecha(19:20),')');                     %             Archivo = strcat('P1_Li(CUSUM)',Ldist,'_',num2str(Fecha(1)),num2str(Fecha(2),'%02d'),num2str(Fecha(3),'%02d'),'(',num2str(Fecha(4),'%02d'),num2str(Fecha(5),'%02d'),')');
Direccion = strcat(Carpeta,'\',Archivo,'.xls');
Hoja = strcat(Ldist,'(',num2str(n),',',num2str(m),')');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Compute the arl %%%%%%%%%%%%%%%%%%%%%%%%%

for r=1:rep  %%% Parallel toolbox
    if NRS == 1
        X = choose_distribution(caso,n);
    end
    [RL, var] = CUSUM_ReferenceSample_Li(X,n,m,H,k,ARL(1),caso,s2,mu2,RL,var);
%         rl(r) = count;                    %Verificar que arl = mean(rl)
end % parfor(rep) ARL;
%     rl' %     mean(rl)%     sdr= std(rl)  %Verificar que arl = mean(rl)
arl(1) = RL/rep; % sum(RL)/rep;        % Son importantes todos los datos
sdrl(1) = sqrt((var - rep*(arl^2))/(rep-1));
fprintf('h = %.2f,\t arl = %.2f(%d),\t sdrl = %.2f, \t n = %d, \t m = %d, \t Dist = %s \n',h(1), arl(1),ARL(1), sdrl(1),n,m,Ldist)

aux = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Verify |ARL-arl| < error 
while abs(ARL(1)-arl(aux))>= error  
% %     WHILE: difference between ARLs is more significant than error.
% %     THEN: Calibrate
    RL = 0;     var = 0;
    X = choose_distribution(caso,n);
    aux = aux+1;
    H = h(aux)*sqrt(m*n*(n+m+1)/12);
    for r=1:rep  %%% Parallel toolbox
        if NRS == 1
            X = choose_distribution(caso,n);
        end
        [RL, var] = CUSUM_ReferenceSample_Li(X,n,m,H,k,1000,caso,s2,mu2,RL,var);
%         rl(r) = count;                    %Verificar que arl = mean(rl)
    end % parfor(rep) ARL;
%     rl' %     mean(rl)%     sdr= std(rl)  %Verificar que arl = mean(rl)
    arl(aux) = RL/rep; % sum(RL)/rep;        % Son importantes todos los datos
    sdrl(aux) = sqrt((var - rep*(arl(aux)^2))/(rep-1));
    fprintf('h = %.2f,\t arl = %.2f(.95*%d),\t sdrl = %.2f, \t n = %d, \t m = %d, \t Dist = %s \n',h(aux), arl(aux),ARL(1), sdrl(aux),n,m,Ldist)

%     aux = 2;
    h(aux+1) = inter_extra_polation([h(aux-1) h(aux)],[arl(aux-1) arl(aux)],ARL(1));
end

Time = toc;


Fecha2 = clock; Fecha2 = num2str(Fecha2(1:5));
datosList = {h(aux),arl(aux),sdrl(aux),ARL, n, m, rep,Ldist, Time, Fecha2};
DATOSList = {'h','arl', 'sdrl', 'ARL', 'n','m','Replicates','Distribution', 'TimeAc',Fecha};
%%%%%%%%%%%%%%%%% GUARDAR RESUMEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlswrite(Direccion,DATOSList,'RESUMEN','A1'); 
xlswrite(Direccion,datosList,'RESUMEN',strcat('A',num2str((im-1)*length(no)+in+2)));
%%%%%%%%%%%%%%%%% GUARDAR DATOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlswrite(Direccion,DATOSList,Hoja, 'A1' )
xlswrite(Direccion,datosList,Hoja, 'A3' )
%%%%%%%%%%%%%%%%% GUARDAR INFORMACIÓN %%%%%%%%%%%%%%%%%%%%%%%%%
xlswrite(Direccion,{'h','ARL'},Hoja, strcat('A',num2str(1)) )
xlswrite(Direccion,h',Hoja, strcat('A',num2str(1+2)) )
xlswrite(Direccion,arl',Hoja, strcat('B',num2str(1+2)) )
xlswrite(Direccion,sdrl',Hoja, strcat('C',num2str(1+2)) )
%%%%%%%%%%%%%%%%% GUARDAR NOTAS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlswrite(Direccion,{'Notas:'},Hoja, num2str((im-1)*length(no)+in+5)) )
xlswrite(Direccion,Notas,Hoja, num2str((im-1)*length(no)+in+5)) )
% size(ARL)
clear h H arl sdrl 
% fprintf('h = \t %f\n', h)
% fprintf('arl = \t %f.\n',arl)
end %n
end %m
end %caso

disp(Direccion)
disp(Fecha)
disp(num2str(clock))


% end % function