function wv = balloon_vbs_full(wmarg)

global fusionfunction

Ns = size(wmarg,1);
Nv = floor(log(Ns)/log(2));

Ns2 = 2*Ns;
Ne = 2^Ns2;
wout = zeros(Ne-1,1);
wout(end) = 1;

for i=1:Ns,
    X = zeros(size(wout));
    ix = double(dec2bin(i-1,Nv))-48;
    ix = ix(end:-1:1);
    
    mask = ones(1,Ns2);
    ix = [0 ix];
    mask(bin2dec(char(ix+48))+1) = 0;
    k = bin2dec(char(mask+48));
    X(k) = wmarg(i,1);
    
    mask = ones(1,Ns2);
    ix(1) = 1;
    mask(bin2dec(char(ix+48))+1) = 0;
    k = bin2dec(char(mask+48));
    X(k) = wmarg(i,2);
    
    X(end) = wmarg(i,3);

    
    wout = feval(fusionfunction,wout,X);
    
end

wv = wout;