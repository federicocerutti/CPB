function m_out = prob_fuse(m_in)

a = prod(m_in);
b = prod(1-m_in);
m_out = a/(a+b);