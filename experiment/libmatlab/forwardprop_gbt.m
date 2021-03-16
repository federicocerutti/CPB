function wout = forwardprop_gbt(wcond,wmarg)

swap = [2 1 3];

N = size(wmarg,1);
Ns = 2^N;

ix = double(dec2bin(0:Ns-1,N))-48;
Ns = 3^N;

wc = ones(1,3);
wout = zeros(1,3);

for i=0:Ns-1, 
    st = double(dec2base(i,3,N))-48; 
    ix2 = ix;
    belief = 1;
    for j=1:N,
        if st(j)<2,
            ix2 = ix2(ix2(:,j)==st(j),:);
        end
        belief = belief*wmarg(j,swap(st(j)+1));
    end
    loc = bin2dec(char(ix2+48))+1;
    wc(1:2) = prod(wcond(loc,1:2),1);
    wc(3) = 1-wc(1)-wc(2);
    wout = wout + wc*belief;
end





    
    