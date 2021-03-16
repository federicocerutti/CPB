function [wout,M,M2] = forwardprop(wcond,wmarg)

mc = wcond*[1;0;0.5];
sc = 2./wcond(:,3);

M = mc*mc'+diag(mc.*(1-mc)./(sc+1));
m_in = wmarg*[1;0;0.5];
N = length(m_in);
Ns = 2^N;
s = 2./wmarg(:,3);


mo = ones(1,Ns);

for i=1:N,
    ix = bitand(bitshift(0:Ns-1,i-N),1);
    mo = mo.*((1-m_in(i)).^(1-ix)).*(m_in(i).^ix);
    ix = ones(Ns,1)*(ix+1);
    m = [1-m_in(i); m_in(i)];
    iz = ix';
    M2 = m(iz).*m(ix)+m(ix).*(1-m(ix))/(s(i)+1).*(2*(ix==iz)-1);
    M = M.*M2;
end

m_out = mc'*mo';

v_out = sum(M(:));


s_out = (m_out-v_out)/max(v_out-m_out^2,0); %max operation put in for numerical percision

s_out = max([s_out 1/m_out 1/(1-m_out)]);

wout = [[m_out 1-m_out]-1/s_out 2/s_out]; 





    
    