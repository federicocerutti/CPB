function node = makenetwork(filename)

fp = fopen(filename,'r');


node(1).children = [];
node(1).parents = [];

index = fscanf(fp,'%d');
frewind(fp);

N = max(index);

node(N).children = [];

while ~feof(fp),
    index = fscanf(fp,'%d',2);
    node(index(1)).children = [node(index(1)).children index(2)];
    node(index(2)).parents = [node(index(2)).parents index(1)];
end

for i=1:N,
    node(i).parents = sort(node(i).parents);
    node(i).children = sort(node(i).children);
end

fclose(fp);
    