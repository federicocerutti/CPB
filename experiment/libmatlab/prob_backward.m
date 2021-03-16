function mout = prob_backward(mc,m_in,parent,ml)

if nargin<4,
    ml = 0.5;
end


N = length(m_in);
Ns = 2^(N+1);

mo = ones(1,Ns);
for i=[1:parent-1 parent+1:N],
    ix = bitand(bitshift(0:Ns-1,i-N),1);
    mo = mo.*((1-m_in(i)).^(1-ix)).*(m_in(i).^ix);
end
ix = bitand(bitshift(0:Ns-1,-N),1);
iy = bitand(0:Ns-1,Ns/2-1)+1;
m = [1-ml ml];
mc = mc';
mo = mo.*((1-mc(iy)).^(1-ix)).*(mc(iy).^ix).*m(ix+1);

ix = bitand(bitshift(0:Ns-1,parent-N),1);

mzxo = sum(mo(ix==1));
mznxo = sum(mo(ix==0));


mout = mzxo/(mzxo+mznxo);


