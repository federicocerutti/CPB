function [pa,pp] = perfbound_beta(p2cg,w2c)

sz = size(p2cg);
sz = [sz(1) 3 prod(sz(2:end))];
w2c = reshape(w2c,sz);


sz = numel(p2cg);

pe = reshape(w2c(:,1,:)+w2c(:,3,:)*.5,[sz 1]);
p2cg = p2cg(:);
alpha1 = reshape(2*w2c(:,1,:)./w2c(:,3,:),[sz 1])+1;
alpha2 = reshape(2*w2c(:,2,:)./w2c(:,3,:),[sz 1])+1;

loc = (~isinf(alpha1)&~isinf(alpha2));

alpha1 = alpha1(loc);
alpha2 = alpha2(loc);

mn = alpha1./(alpha1+alpha2);


er = p2cg-pe;

er = er(loc);

N = length(er);
%sz = N;


pp = 0:.01:.99;
pa = zeros(size(pp));
%hd = waitbar(0);
lb = zeros(N,1);
ub = zeros(N,1);
loc1 = find(((alpha1>1)&(alpha2>1))|((alpha1==1))&(alpha2==1));
loc2 = find((alpha1<1)&(alpha2<1));
loc3 = find((alpha1<=1)&(alpha2>1));
lb(loc3) = -mn(loc3);
loc4 = find((alpha1>1)&(alpha2<=1));
ub(loc4) = 1-mn(loc4);
N1 = length(loc1);
N2 = length(loc2);
N3 = length(loc3);
N4 = length(loc4);

%loc1 = 1:N;
%N1 = N; N2 = 0; N3 = 0; N4 = 0;


for i=1:length(pp),
    %waitbar(i/length(pp),hd);
    pp2 = (1-pp(i))/2*ones(N1,1);
    %mn = alpha1./(alpha1+alpha2);
    lb(loc1) = betaincinv(pp2,alpha1(loc1),alpha2(loc1))-mn(loc1);
    ub(loc1) = betaincinv(1-pp2,alpha1(loc1),alpha2(loc1))-mn(loc1);
    if N2>0,
        %pp2 = (0.5-pp(i)/2)*ones(N2,1);
        pp2 = pp(i)/2*ones(N2,1);
        lb(loc2) = betaincinv(pp2,alpha1(loc2),alpha2(loc2))-mn(loc2);
        ub(loc2) = betaincinv(1-pp2,alpha1(loc2),alpha2(loc2))-mn(loc2);
    end
    if N3>0,
        pp2 = pp(i)*ones(N3,1);
        ub(loc3) = betaincinv(pp2,alpha1(loc3),alpha2(loc3))-mn(loc3);
    end
    if N4>0,
        pp2 = (1-pp(i))*ones(N4,1);
        lb(loc4) = betaincinv(pp2,alpha1(loc4),alpha2(loc4))-mn(loc4);
    end
    
    in = (er<=ub).*(er>=lb);
    in(loc2) = 1-in(loc2);
    pa(i) = sum(in)/N;
end
%delete(hd);
