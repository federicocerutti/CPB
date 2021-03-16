function Z = fusion_VBS_DST_full(X,Y)

if nargin<2,
    N = size(X,1);
    Z = X(1,:);
    for i=2:N,
        Z = fusion_VBS_DST(Z,X(i,:));
    end
else

    belief = zeros(size(X));

    Nv = round(log(log(length(X(:)))/log(2))/log(2));

    Ns = 2^Nv;

    X = X(:);
    Y = Y(:);
    loc1 = find(X~=0)';
    loc2 = find(Y~=0)';
    for i=loc1,
        ix1 = double(dec2bin(i,Ns))-48;
        for j=loc2,
            ix2 = double(dec2bin(j,Ns))-48;
            noconflict = sum(ix1&ix2);
            if noconflict,
                ix3 = ix1&ix2;
                k = bin2dec(char(ix3+48));
                belief(k) = belief(k)+X(i)*Y(j);
            end
        end
    end
    Z = belief/sum(belief(:));
end
            
        
    

