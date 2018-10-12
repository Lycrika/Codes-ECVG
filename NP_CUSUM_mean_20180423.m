%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Li et al. (2010)-CUSUM control chart
% % Measuring the performance on different scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%%%%%%%%%%%%% parpool('local') clear all; clearvars
clear all
clc
%%%%%%%%%%%%%%%%%%% To save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stat,struc] = fileattrib;
Carpeta = struc.Name; disp(Carpeta)
Fecha = clock; Fecha = num2str(Fecha(1:5),'%02d');
%%%%%%%%%%%%%%%%% Parameters of the simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
CASOS = [0]; dist = 'UNWGPLt';  %% Selects the data distribution
REP = 200;          %500;       % Compute the AARL, SDARL, ...  REP =  1,000
rep = 10000;        %10000;     % Compute the ARL  & SDRL       rep = 10,000
no = [100];                     % Reference sample size         n = #(X)
 mo = [1 3 7 10];                       % Batch sample size to monitor  m = #(Y) 
%%%%%%%%%%%%%%%%% PARAMETERS OUT OF CONTROL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sigma2 = 1;
MU2 =  [0 .25 .5 1 1.5 2 3]; % [3 2 1.5 1 .5 .25 0];  %

for caso = CASOS        %%%%%%%%% DISTRIBUTION 
for n = no              %%%%%%%%% REFERENCE SAMPLE SIZE 
for m = mo              %%%%%%%%% MONITORING SAMPLE SIZE
    Ldist = dist(caso+1);    %Ldist = choose_distribution_name(caso);    
    k = .5*sqrt(m*n*(n+m+1)/12);
    H = CalibratedCL(n,m)*sqrt(m*n*(n+m+1)/12);
    tic                 %%%%%%%%% STARTING TO COUNTING TIME
    % for ia = 1:kk    %     alpha = ALPHA(ia);    
    uu = numel(MU2);
    for iu = 1:uu
        mu2 = MU2(iu);  %%%%%%%%% IN or OUT OF CONTROL MEAN
        ARL = zeros(1,REP);
        SDRL = zeros(1,REP);
        Fecha2 = clock; Fecha2 = num2str(Fecha2(1:5));
        
        % AARL; SDARL; Each one of the ARL is going to be repeat 'REP' times.
        for R =1:REP  % I can't use parfor here, because I'm saving ARL(j).
            % ARL; Each one of the RL is going to be repeat 'rep' times.
            % The values have to be restarted
            fprintf('\nREPLICATE = %d,\t replicates = %d, \t n=%d,  \t m= %d, \t mu2 = %f, \t Time=%d\n',[R, rep, n, m, mu2, toc])

                RL = 0;     var =0; %H = 353
            x = choose_distribution(caso,n);
%                 size(x)
            for r=1:rep  %%% Parallel toolbox 
                count = 0;
                Cplus = 0; Cminus = 0;
                while (Cplus < H && Cminus < H )   
                    y = [x(1:n), choose_distribution(caso,m)*(sigma2)+mu2];         
                    [Cplus,Cminus] = CUSUM_wrs(y,k,n,m,Cplus,Cminus);
                    count = count+1;
                end %while
%                     [count, RL, R, r]
                RL = count + RL;        % Para obtener el ARL, no es necesario guardar todos los RLs.
                var = count^2 + var;
            end % parfor(rep) ARL;
            ARL(R) = RL/rep;% sum(RL)/rep;        % Son importantes todos los datos
            SDRL(R) = sqrt((var - rep*(ARL(R)^2))/(rep-1));
%             fprintf('\nDistibution = %s\t', Ldist)
            fprintf('ARL(%d) = %.2f,\tSDRL(%d) = %.2f\t',  R, ARL(R), R, SDRL(R));  disp(num2str((fix(clock)),'%02d'))


        end %for(REP) AARL; SDARL;
        %%%%%%%%%%%%%%%%% OBTENER DATOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        AARL = mean(ARL);
        SDARL = std(ARL);
        ASDRL = mean(SDRL);
        SDSDRL = std(SDRL);
        Time=toc;  %%%%%%%%%%%%%%%%%%%%%%%%% DETENER CONTADOR DE TIEMPO
        %%%%%%%%%%%%%%%%% IMPRIMIR DATOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\nAARL = %.2f,\tSDARL = %.2f\tASDRL = %.2f,\tSDSDRL = %.2f\t\n',  AARL, SDARL, ASDRL, SDSDRL);  disp(num2str((fix(clock)),'%02d'))

        %%%%%%%%%%%%%%%%% GUARDAR INFORMACIÓN %%%%%%%%%%%%%%%%%%%%%%%%%
        Archivo = strcat('P1_Li(CUSUM)',Ldist,'(',num2str(n),',',num2str(m),')_REP(',num2str(REP),')',Fecha(1:4),...
                Fecha(7:8),Fecha(11:12),'(',Fecha(15:16),Fecha(19:20),')');                     %             Archivo = strcat('P1_Li(CUSUM)',Ldist,'_',num2str(Fecha(1)),num2str(Fecha(2),'%02d'),num2str(Fecha(3),'%02d'),'(',num2str(Fecha(4),'%02d'),num2str(Fecha(5),'%02d'),')');
        Direccion = strcat(Carpeta,'\',Archivo,'.xls');
        Hoja = strcat('Mu2 = ',num2str(mu2));
        datosList = {AARL,SDARL,ASDRL,SDSDRL,REP, rep, n, m, mu2,Ldist, Time, Fecha2};
        DATOSList = {'ARL(R)','SDRL(R)', 'AARL', 'SDARL','ASDRL', 'SDSDRL',...
                'REPLICATES', 'replicates', 'n','m', 'mu2','Distribution', 'TimeAc',Fecha};
        %%%%%%%%%%%%%%%%% GUARDAR RESUMEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xlswrite(Direccion,DATOSList,'RESUMEN','A1'); 
        xlswrite(Direccion,datosList,'RESUMEN',strcat('C',num2str(iu+2)));
        %%%%%%%%%%%%%%%%% GUARDAR DATOS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xlswrite(Direccion,DATOSList,Hoja, 'A1' )
        xlswrite(Direccion,datosList,Hoja, 'C3' )
        xlswrite(Direccion,ARL',Hoja, 'A3' )
        xlswrite(Direccion,SDRL',Hoja, 'B3' )
        %%%%%%%%%%%%%%%%% BORRAR INFORMACIÓN %%%%%%%%%%%%%%%%%%%%%%%%%%
        clear ARL SDRL
    end %MU2
    % end %ALPHA
    disp(Direccion)
    disp(Fecha)
    disp(num2str(clock))
end %m
end %n
end %caso

