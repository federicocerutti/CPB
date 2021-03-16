
function wout = backprop_gbt(wcond,wmarg,parent,wlike)

swap = [2 1 3];

N = size(wmarg,1);
Ns = 2^N;
ix = bitand(bitshift(0:Ns-1,parent-N),1);

wcondp = wcond(ix==1,:);
wcondn = wcond(ix==0,:);

plyp = wcondp*[1;0;1];
plnyp = wcondp*[0;1;1];
plyn = wcondn*[1;0;1];
plnyn = wcondn*[0;1;1];

wmarg = wmarg([1:parent-1 parent+1:N],:);

N = N-1;
NN = max(N,1);
Ns = Ns/2;

ix = double(dec2bin(0:Ns-1,NN))-48;
Ns = 3^N;


wout = zeros(1,3);
for i=0:Ns-1,
    st = double(dec2base(i,3,N))-48;
    belief = 1;
    ix2 = ix;
    for j=1:N,
        if st(j)<2,
            ix2 = ix2(ix2(:,j)==st(j),:);
        end
        belief = belief*wmarg(j,swap(st(j)+1));
    end
    loc = bin2dec(char(ix2+48))+1;
    
    plx = 1-prod(1-plyp(loc));
    plnx = 1-prod(1-plyn(loc));
    plu = 1-(1-plx)*(1-plnx);
    plx = plx/plu;
    plnx = plnx/plu;
    wout = wout + [1-plnx 1-plx plx+plnx-1]*belief*wlike(1);
    
    plx = 1-prod(1-plnyp(loc));
    plnx = 1-prod(1-plnyn(loc));
    plu = 1-(1-plx)*(1-plnx);
    plx = plx/plu;
    plnx = plnx/plu;
    wout = wout + [1-plnx 1-plx plx+plnx-1]*belief*wlike(2);
    
    wout = wout+[0 0 1]*belief*wlike(3);
    
end





