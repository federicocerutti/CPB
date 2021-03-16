function node = inferGBT(node)

global backfunction fusionfunction forwardfunction

backfunction = 'backprop_gbt';
forwardfunction = 'forwardprop_gbt';
fusionfunction = 'fusion_DST';

node = inferMODULE_belief(node);
