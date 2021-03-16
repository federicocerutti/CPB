function node = createVBS(node)

global fusionfunction 

fusionfunction = 'fusion_VBS_DST_full';

N = length(node);

for i=1:N,
    node(i).wv = balloon_vbs_full(node(i).w);
end