%{
Copyright (c) 2018 Federico Cerutti <CeruttiF@cardiff.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
%}

function res = BNfromPY(jsongraph, bn)
warning('off','all');

n = jsongraph;

node(1).children = [];
node(1).parents = [];

N = max(size(n));
node(N).children = [];

for i=1:N,
    node(i).children = n{i}.children;
    node(i).parents = n{i}.parents;
    node(i).children = sort(node(i).children);
    node(i).parents = sort(node(i).parents);
    node(i).p = n{i}.p';
    node(i).value = n{i}.value;
    if bn ~= 1
        if size(n{i}.w, 2) == 1
            node(i).w = n{i}.w;
        else
            node(i).w = n{i}.w;
        end
    end
end

N = length(node);
end_loc = getends(node);
Nends = length(end_loc);
interior_loc = setdiff((1:N)',end_loc);
Nint = length(interior_loc);


res = {};

if bn == 1,
    node = prob_infer(node);

    for j=1:Nint,
        res = [res, node(interior_loc(j)).pe];
    end
elseif bn == 0,
    node = inferSLBN(node);
    for j=1:Nint,
        size(node(interior_loc(j)).we);

        res = [res; node(interior_loc(j)).we];
    end
elseif bn == 2,
    node = inferCREDAL(node);
    for j=1:Nint,
        size(node(interior_loc(j)).we);

        res = [res; node(interior_loc(j)).we];
    end
elseif bn ==3,
    node = inferGBT(node);
    for j=1:Nint,
        size(node(interior_loc(j)).we);

        res = [res; node(interior_loc(j)).we];
    end
elseif bn ==4,
    node = createVBS(node);
    node = inferVBSfull(node);
    for j=1:Nint,
        size(node(interior_loc(j)).we);

        res = [res; node(interior_loc(j)).we];
    end
end
