function wout = approx_likefusion_all(win)

epsilon = 2;

m = win*[1;0;0.5];
m1 = prod(m);
m2 = prod(1-m);
mout = m1/(m1+m2);

%s = 2./win(:,3);
u  = win(:,3);

%sout = 1/(mout*(1-mout)*sum(1./m./(1-m)./(s+1)))-1;

if (mout==1)||(mout==0),
    uout = 0;
else
    partial = mout*(1-mout)*sum(1./m./(1-m).*u./(u+2));
    if partial<=1/(1+epsilon),
        uout = 2*partial/(1-partial);
    else
        uout = 2/epsilon;
    end
end
uout = min([uout 2*mout 2*(1-mout)]);
wout = [[mout 1-mout]-uout/2 uout];
%sout = max([sout 1/mout 1/(1-mout)]);
% = [[mout 1-mout]-1/sout 2/sout];