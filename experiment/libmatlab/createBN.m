function node = createBN(node)

N = length(node);
for i=1:N,
    nparents = length(node(i).parents);
    node(i).p = rand(2^nparents,1);
end