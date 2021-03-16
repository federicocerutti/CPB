function wout = approx_backprop(wcond,wmarg,parent,wlike)

mc = wcond*[1;0;0.5];
m_in = wmarg*[1;0;0.5];

sc = 2./wcond(:,3);
s = 2./wmarg(:,3);

if nargin<4,
    wlike = [0.5 0.5 0];
end

ml = wlike*[1;0;0.5];
sl = 2./wlike(:,3);

N = length(m_in);
Ns = 2^(N+1);

mo = ones(1,Ns);
M = ones(Ns,Ns);
for i=[1:parent-1 parent+1:N],
    ix = bitand(bitshift(0:Ns-1,i-N),1);
    mo = mo.*((1-m_in(i)).^(1-ix)).*(m_in(i).^ix);
    ix = ones(Ns,1)*(ix+1);
    m = [1-m_in(i); m_in(i)];
    iz = ix';
    M2 = m(iz).*m(ix)+m(ix).*(1-m(ix))/(s(i)+1).*(2*(ix==iz)-1);
    M = M.*M2;
end
ix = bitand(bitshift(0:Ns-1,-N),1);
iy = bitand(0:Ns-1,Ns/2-1)+1;
m = [1-ml ml];
mc = mc';
mo = mo.*((1-mc(iy)).^(1-ix)).*(mc(iy).^ix).*m(ix+1);
iy = ones(Ns,1)*iy;
ix = ones(Ns,1)*ix;

iy2 = iy';
ix2 = ix';


M2 = mc(iy).*(1-mc(iy))./(sc(iy)+1).*((-1).^ix).*((-1).^ix2).*(iy==iy2)+((1-mc(iy)).^(1-ix)).*(mc(iy).^ix).*((1-mc(iy2)).^(1-ix2)).*(mc(iy2).^ix2);
M = M.*M2;

M2 = ml.*(1-ml)/(sl+1).*((-1).^ix).*((-1).^ix2)+m(ix+1).*m(ix2+1);
M = M.*M2;



ix = bitand(bitshift(0:Ns-1,parent-N),1);

mzxo = sum(mo(ix==1));
mznxo = sum(mo(ix==0));

ix = ones(Ns,1)*ix;
ix2 = ix';

varzx = max(sum(M((ix==1)&(ix2==1)))-mzxo^2,0);  %max put in due to numerical percision errors

varznx = max(sum(M((ix==0)&(ix2==0)))-mznxo^2,0);

cross = max(sum(M((ix==1)&(ix2==0)))-mzxo*mznxo,0);

mzx = mzxo/(mzxo+mznxo);

szx = 1/mzx/(1-mzx)/max(varzx/mzxo^2+varznx/mznxo^2-2*cross/mzxo/mznxo,0)-1; %max put in due to numerical percision errors

szx = max([szx 1/mzx 1/(1-mzx)]);

wout = [[mzx 1-mzx]-1/szx 2/szx];


