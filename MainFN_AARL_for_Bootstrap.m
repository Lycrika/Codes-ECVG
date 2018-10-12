function MainFN_AARL_for_Bootstrap

tic; Notas = {'Esta es la función principal para calcular las medidas basadas en ARL y SDRL, donde antes se ajustan los límites de control para cada muestra de referencia.'};
%%%%%%%%%%%%%%%%%%% To save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,struc] = fileattrib;
Carpeta = struc.Name; disp(Carpeta)
Fecha = clock; Fecha = num2str(Fecha(1:5),'%02d');
%%%%%%%%%%%%%%%%% FIXED Parameters of the simulation %%%%%%%%%%%%%%%%%%%%%%

CASOS = [0]; dist = 'UNWGPLt';  %% Selects the data distribution
no = [30 50 100 200 500 1000 3000 5000 10000];                     % Reference sample size         n = #(X)
mo = [1 3 5 7 10];                       % Batch sample size to monitor  m = #(Y)
mu2 = 0;       s2 = 1;                   % \mu_2 and \sigma_2 => No change. 
% % % NRS = 0; % New Reference Sample = {No = 0, Yes = 1}
ARL(1) = 100; error = .05*ARL(1);    % error = 5 percento of the ARL
rep = 100000;


for caso = CASOS        %%%%%%%%% DISTRIBUTION 
mm = numel(mo);
for im = 1:mm              %%%%%%%%% MONITORING SAMPLE SIZE 
    m = mo(im);     nn = numel(no);
for in = 1:nn                %%%%%%%%% REFERENCE SAMPLE SIZE
    n = no(in);
%%%%%%%%%%%%%%%%%%% To save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ldist = dist(caso+1);
Archivo = strcat('P1_H_calibration_ARL(',num2str(ARL(1)),')rep(',num2str(rep),')_',Fecha(1:4),...
        Fecha(7:8),Fecha(11:12),'(',Fecha(15:16),Fecha(19:20),')');                     %             Archivo = strcat('P1_Li(CUSUM)',Ldist,'_',num2str(Fecha(1)),num2str(Fecha(2),'%02d'),num2str(Fecha(3),'%02d'),'(',num2str(Fecha(4),'%02d'),num2str(Fecha(5),'%02d'),')');
Direccion = strcat(Carpeta,'\',Archivo,'.xls');
Hoja = strcat(Ldist,'(',num2str(n),',',num2str(m),')');


%%%%%%%%%%%%%%%%%%% Select the reference sample %%%%%%%%%%%%%%%%%%%%%%%%%%%
X = choose_distribution(caso,n);
%%%%%%%%%%%%%%%%%%% Bootstrap calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = bootstrap_calibration(X, n, m, h, ARL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Calculate the AARL, etc. for the same reference sample
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






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
xlswrite(Direccion,{'Notas:'},Hoja, num2str((im-1)*length(no)+in+5) )
xlswrite(Direccion,Notas,Hoja, num2str((im-1)*length(no)+in+5)  ) 
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




end

