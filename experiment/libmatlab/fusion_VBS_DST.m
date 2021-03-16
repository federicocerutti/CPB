function Z = fusion_VBS_DST(X,Y)

if nargin<2,
    N = size(X,1);
    Z = X(1,:);
    for i=2:N,
        Z = fusion_VBS_DST(Z,X(i,:));
    end
else

    belief = zeros(size(X));

    Nv = round(log(length(X(:)))/log(3));


    Ne = 3^Nv;

    X = X(:);
    Y = Y(:);

    for i=0:Ne-1,
        ix1 = double(dec2base(i,3,Nv))-48;
        ix1 = ix1(Nv:-1:1);
        for j=0:Ne-1,
            ix2 = double(dec2base(j,3,Nv))-48;
            ix2 = ix2(Nv:-1:1);
            noconflict = prod((ix2==ix1)|(ix1==2)|(ix2==2));
            if noconflict,
                ix3 = min(ix1,ix2);
                ix3 = ix3(Nv:-1:1);
                k = base2dec(char(ix3+48),3)+1;
                belief(k) = belief(k)+X(i+1)*Y(j+1);
            end
        end
    end
    Z = belief./sum(belief(:));
end
            
        
    

