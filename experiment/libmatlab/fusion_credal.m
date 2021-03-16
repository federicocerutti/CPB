function wout = fusion_credal(win)

pmin = win(:,1);
pmax = win*[1;0;1];

A = prod(pmin);
B = prod(1-pmin);

pfmin = A/(A+B);

A = prod(pmax);
B = prod(1-pmax);

pfmax = A/(A+B);

if isnan(pfmin)&&(pfmax==1),
    pfmin = 1;
elseif isnan(pfmax)&&(pfmin==0),
    pfmax = 0;
end

wout = [pfmin 1-pfmax pfmax-pfmin];
