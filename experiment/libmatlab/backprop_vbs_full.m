function wout = backprop_vbs_full(wcond,wmarg,parent,wlike)

global fusionfunction


Nv = size(wmarg,1)+1;
index = double(dec2bin(0:2^Nv-1,Nv))-48;

Ns = 2^(2^Nv);


wout = wcond;

for i=1:Nv-1,
    if i~=parent,
        X = zeros(size(wcond));

        ix = zeros(1,2^Nv);
        ix(index(:,Nv+1-i)==1)= 1;
        k = bin2dec(char(ix+48));
        X(k) = wmarg(i,1);

        ix = zeros(1,2^Nv);
        ix(index(:,Nv+1-i)==0)= 1;
        k = bin2dec(char(ix+48));
        X(k) = wmarg(i,2);

        X(end) = wmarg(i,3);

        wout = feval(fusionfunction,wout,X);
    end
end

X = zeros(size(wcond));

ix = zeros(1,2^Nv);
ix(index(:,1)==1)= 1;
k = bin2dec(char(ix+48));
X(k) = wlike(1);

ix = zeros(1,2^Nv);
ix(index(:,1)==0)= 1;
k = bin2dec(char(ix+48));
X(k) = wlike(2);

X(end) = wlike(3);

wout = feval(fusionfunction,wout,X);

wz = zeros(1,3);

for i=1:Ns-1,
    ix = double(dec2bin(i,2^Nv))-48;
    val = unique(index(ix == 1,Nv+1-parent));
    if length(val)==2,
        wz(3) = wz(3)+wout(i);
    elseif val==0,
        wz(2) = wz(2)+wout(i);
    else
        wz(1) = wz(1)+wout(i);
    end
end

wout = wz;

