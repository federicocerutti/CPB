function [wout,A,P] = backprop_credal(wcond,wmarg,parent,wlike)

pcmin = wcond(:,1);
pcmax = wcond*[1;0;1];

pmmin = wmarg(:,1);
pmmax = wmarg*[1;0;1];


if nargin<4,
    wlike = [0.5 0.5 0];
end

plmin = wlike(1);
plmax = wlike(1)+wlike(3);



N = size(wmarg,1)+1;
Ns = 2^N;
parent = parent+1;

pmmin = [plmin; pmmin];
pmmax = [plmax; pmmax];

pmmin = [1-pmmin pmmin];
pmmax = [1-pmmax pmmax];



%instantiate all possible value of other parent and central node values

A = ones(1,Ns);
sz = 1;
for i=[1:parent-1 parent+1:N],
    ix = bitand(bitshift(0:Ns-1,i-N),1);
    A = [A.*(ones(sz,1)*pmmin(i,ix+1)); A.*(ones(sz,1)*pmmax(i,ix+1))];
    sz = sz*2;
end


P = [1-pcmin 1-pcmax; pcmin pcmax];

%index for when parent node is zero and one

ix = bitand(bitshift(0:Ns-1,parent-N),1);

%when index is of parent node is one, swithch max and min value of
%conditionals to maximixze/minimize overall message

P(ix==1,:) = P(ix==1,[2 1]);


B = A(:,ix==1)*P(ix==1,:);
C = A(:,ix==0)*P(ix==0,:);

pp = B./(B+C);

pomin = min(pp(:));
pomax = max(pp(:));

wout = [pomin 1-pomax pomax-pomin];



