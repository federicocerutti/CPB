function [wout,A] = forwardprop_credal(wcond,wmarg)

pcmin = wcond(:,1);
pcmax = wcond*[1;0;1];

pmmin = wmarg(:,1);
pmmax = wmarg*[1;0;1];

pmmin = [1-pmmin pmmin];
pmmax = [1-pmmax pmmax];


N = size(wmarg,1);
Ns = 2^N;

A = ones(1,Ns);
for i=1:N,
    NN = 2^(i-1);
    ix = bitand(bitshift(0:Ns-1,i-N),1);
    A = [A.*(ones(NN,1)*pmmin(i,ix+1)); A.*(ones(NN,1)*pmmax(i,ix+1))];
end
pomin = min(A*pcmin);
pomax = max(A*pcmax);

wout = [pomin 1-pomax pomax-pomin];





    
    