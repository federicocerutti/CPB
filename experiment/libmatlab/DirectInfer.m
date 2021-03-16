function node = DirectInfer(node,val)

Ns = size(val,1);
loc = ones(Ns,1);

N = length(node);
for i=1:N,
    if ~isempty(node(i).value),
        loc = loc.*(val(:,i)==node(i).value);
    end
end
val = val(loc==1,:);
Ns = size(val,1);
for i=1:N,
    if isempty(node(i).value),
        r = sum(val(:,i));
        node(i).we = [r Ns-r 2]/(Ns+2);
    end
end
        