function node = inferVBS(node)

global backfunction fusionfunction forwardfunction

backfunction = 'backprop_vbs_full';
forwardfunction = 'forwardprop_vbs_full';
fusionfunction = 'fusion_VBS_DST_full';

node = inferMODULEv(node);
