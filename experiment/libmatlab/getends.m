function loc = getends(node)

N = length(node);
loc = [];
for i=1:N,
    %if (length(node(i).children)==0)||(length(node(i).parents)==0),
    if length(node(i).children)+length(node(i).parents)==1,
        loc = [loc; i];
    end
end