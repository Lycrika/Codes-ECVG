% % % Crsitina hizo sus desmadres
% % % This code reproduces the paper "Keefe, Woodall & Jones (2015) The conditional In-Control Performance of Self-Starting Control Charts"
% % % Obtains the A(ARL) y SD(ARL) of SSCUSUM

clear all
% clc;
H = 353; k = .5;
rep1=100; %parameter i 10000
rep2=10; %parameter l 1000
n=[2 3 5 10 20 30]; % parameter p
tic
%for p=2:3
p=2; n(p)
% parpool('local')
parfor l=1:rep2
%     l
    % step 1
     y=randn(1,n(p));
     %j=n(1)
     
    for i=1:rep1 % step 4
%         i
        j=n(p);  % % n(p) <=j  ; M <= m
%         y=randn(1,n(p));
        Cplus1=0;
        Cminus1=0; %i
        %step 2 and 3            
        while Cplus1 <= H && Cminus1 <= H
            j=j+1       % n(p) <=j  ; M <= m
            y(j)=randn;
            [Cplus1,Cminus1] = wr_fn(k,n(p),y,j,Cplus1,Cminus1) ; %(k,M,x,m,Cplus1,Cminus1)
        end
        RL(i) = j-n(p);
        j=n(p);
        pause(1)
    end
    ARL(l)=mean(RL); 
%     SDRL(l)=std(RL); 
end
AARL=mean(ARL)
SDARL=std(ARL)
toc

fprintf('\nAARL = \n')
fprintf('%f\n',AARL')
fprintf('\nSDARL = \n')
fprintf('%f\n',SDARL')






