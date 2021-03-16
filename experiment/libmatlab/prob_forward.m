function m_out = prob_forward(mc,m_in)

N = length(m_in);
Ns = 2^N;

mo = ones(1,Ns);

for i=1:N,
    ix = bitand(bitshift(0:Ns-1,i-N),1);
    mo = mo.*((1-m_in(i)).^(1-ix)).*(m_in(i).^ix);
end

m_out = mc'*mo';






    
    