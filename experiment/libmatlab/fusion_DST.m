function wout = fusion_DST(win)

wout = zeros(1,3);

wout(3) = prod(win(:,3));

wout(1) = prod(1-win(:,2))-wout(3);
wout(2) = prod(1-win(:,1))-wout(3);

wout = wout/sum(wout);