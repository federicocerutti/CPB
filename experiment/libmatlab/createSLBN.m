function [node, val] = createSLBN(node,ntrain)

N = length(node);
queue = 1:N;
activated = zeros(N,1);
order = [];

while ~isempty(queue),
    for i=queue, 
        if sum(activated(node(i).parents))==length(node(i).parents),
            activated(i) = 1;
            order = [order i];
        end
    end
    queue = setdiff(queue,order);
end

val = zeros(ntrain,N);

for i = order,
    np = length(node(i).parents);
    if np==0,
        val(:,i) = rand(ntrain,1)<node(i).p;
    else
        val(:,i) = rand(ntrain,1)<node(i).p(val(:,node(i).parents)*(2.^(np-1:-1:0)')+1);
    end
end

for i=1:N,
    np = length(node(i).parents);
    if np==0,
        rcount = sum(val(:,i));
        scount = sum(val(:,i)==0);
        node(i).w = [rcount scount 2]/(rcount+scount+2);
    else
        idx = val(:,node(i).parents)*2.^(np-1:-1:0)'+1;
        for j=1:2^np,
            rcount = sum(val(idx==j,i));
            scount = sum(val(idx==j,i)==0);
            node(i).w(j,:) = [rcount scount 2]/(rcount+scount+2);
        end
    end
end


            
     