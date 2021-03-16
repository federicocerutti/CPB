function node = prob_infer(node)

N = length(node);

active = [];
for i=1:N,
    node(i).sendchildren = zeros(size(node(i).children));
    node(i).recchildren = zeros(size(node(i).children));
    node(i).sendparents = zeros(size(node(i).parents));
    node(i).recparents =  zeros(size(node(i).parents));
    node(i).childmsg = zeros(size(node(i).children));
    node(i).parentmsg = zeros(size(node(i).parents));
    node(i).forward = [];
    node(i).back = [];
    if (length(node(i).recparents)+length(node(i).recchildren)==1),
        active = [active i];
    end
end

while ~isempty(active),
    id = min(active);
    recparent = find(node(id).recparents==0);
    recchild = find(node(id).recchildren==0);
    if length(recparent)+length(recchild)==1,
        if isempty(recparent),
            sendid = recchild;
            schild = 1;
        else
            sendid = recparent;
            schild = 0;
        end
    else
        sendid = find(node(id).sendparents==0,1);
        schild = 0;
        if isempty(sendid),
            sendid = find(node(id).sendchildren==0,1);
            schild = 1;
        end
    end
    if schild==0,
        if isempty(node(id).back),
            if isempty(node(id).value),
                node(id).back = prob_fuse(node(id).childmsg);
            else
                node(id).back = node(id).value;
            end
        end
        
        node(id).sendparents(sendid) = 1;
        idx = node(id).parents(sendid);
        loc = find(node(idx).children==id);
        node(idx).childmsg(loc) = prob_backward(node(id).p,node(id).parentmsg,sendid,node(id).back);
        node(idx).recchildren(loc) = 1;
        leftsend = sum(node(idx).sendparents==0)+sum(node(idx).sendchildren==0);
        leftrec = sum(node(idx).recchildren==0)+sum(node(idx).recparents==0);
        if (leftrec<=1)&&(leftsend>0),
            active = union(active,idx);
        end     
    else 
        if isempty(node(id).forward),
            if isempty(node(id).value),
                node(id).forward = prob_forward(node(id).p,node(id).parentmsg);
            else
                node(id).forward = node(id).value;
            end
        end
        
        node(id).sendchildren(sendid) = 1;
        idx = node(id).children(sendid);
        loc = find(node(idx).parents==id);
        node(idx).parentmsg(loc) = prob_fuse([node(id).forward node(id).childmsg([1:sendid-1 sendid+1:end])]);
        node(idx).recparents(loc) = 1;
    end
    
    leftsend = sum(node(idx).sendparents==0)+sum(node(idx).sendchildren==0);
    leftrec = sum(node(idx).recchildren==0)+sum(node(idx).recparents==0);
    if (leftrec<=1)&&(leftsend>0),
        active = union(active,idx);
    end  
    
    leftsend = sum(node(id).sendparents==0)+sum(node(id).sendchildren==0);
    leftrec = sum(node(id).recparents==0)+sum(node(id).recchildren==0);
    if (leftrec==1)||(leftsend==0),
        loc = find(active==id);
        active(loc) = [];
    end
end
    
    
    
for id=1:N,
    if isempty(node(id).forward),
        node(id).forward = prob_forward(node(id).p,node(id).parentmsg);
    end
    if isempty(node(id).back),
        node(id).back = prob_fuse(node(id).childmsg);
    end
    node(id).pe = prob_fuse([node(id).forward node(id).back]);
end