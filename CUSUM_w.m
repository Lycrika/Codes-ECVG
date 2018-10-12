% % % Function 'wr_fn' of 'NPsscusum_sdarl_fn'
% % % This code reproduces the paper "Keefe, Woodall & Jones (2015) The conditional In-Control Performance of Self-Starting Control Charts"
% % % Obtains the A(ARL) y SD(ARL) of SSCUSUM

function [Cplus1,Cminus1] = CUSUM_w(x,k,n,m,Cplus1,Cminus1) 

Ranks = tiedrank(x);
wr = sum(Ranks(n+1:m+n));
Cplus2  = max(0, wr - m*(m+n+1)/2 - k + Cplus1);
Cminus2 = max(0,-wr + m*(m+n+1)/2 - k + Cminus1); 
Cplus1=Cplus2;
Cminus1=Cminus2;

end

% function [Cplus1,Cminus1] = wr_fn(k,M,x,m,Cplus1,Cminus1) 
% 
% rangos = tiedrank(x);
% wr = sum(rangos(M+1:m,1));
% Cplus2 = max(0, wr - m*(m+M+1)/2 - k + Cplus1);
% Cminus2 = max(0,-wr + m*(m+M+1)/2 - k + Cminus1); 
% Cplus1=Cplus2;
% Cminus1=Cminus2;
% 
% end

% T=sqrt((n-1)/n)*(y-ybar1)/(stdy1);
% norminv(tcdf(T, n-2));

% ybar2=(((n-1)*ybar1)+y)/n;
% vary2=((n-2)/(n-1))*(stdy1^2)+((1/(n-1))*(y-ybar1)^2);            
% stdy2=sqrt(vary2);
% ybar1=ybar2;
% stdy1=stdy2;