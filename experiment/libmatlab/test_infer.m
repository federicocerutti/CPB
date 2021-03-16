Ntrain = 10;
Nmonte = 2;
Nnetworks = 5;

load state.mat
rand('state',st);

Nruns = Nmonte*Nnetworks;

hd = waitbar(0);


node = makenetwork('net3.file');
%node = makenetwork('C:\Users\Lance\Documents\ARL\SubLogic\SubjectiveNetwork\Experiments\graph2.txt');

N = length(node);
end_loc = getends(node);
Nends = length(end_loc);
interior_loc = setdiff((1:N)',end_loc);
Nint = length(interior_loc);

pgt = zeros(Nruns*Nint,1);
we = zeros(Nruns*Nint,3);
wec = we;
weg = we;
wed = we;
wev = we;

tm = zeros(Nruns,1);
tc = tm;
tg = tm;
td = tm;
tv = tm;
ts = tm;

cnt = 0;
for i=1:Nruns,
    waitbar(i/Nruns,hd,sprintf('%d/%d',i,Nruns));
    if mod(i-1,Nmonte)==0,
        node = createBN(node);
    end
    loc = getends(node);
    [node,val] = createSLBN(node,Ntrain);
    node = createVBS(node);
    for j=1:Nends,
        node(end_loc(j)).value = rand(1)<0.5;
    end
    tic,
    node = prob_infer(node);
    tm(i) = toc;
    nodec = node;
    tic,
    nodec = inferCREDAL(nodec);
    tc(i) = toc;
    noded = node;
    tic,
    noded = DirectInfer(noded,val);
    td(i) = toc;
    nodeg = node;
    tic,
    nodeg = inferGBT(nodeg);
    tg(i) = toc;
    nodev = node;
    tic,
    nodev = inferVBSfull(nodev);
    tv(i) = toc;
    tic,
    node = inferSLBN(node);
    ts(i) = toc;
    for j=1:Nint,
        cnt = cnt+1;
        pgt(cnt) = node(interior_loc(j)).pe;
        we(cnt,:) = node(interior_loc(j)).we;
        wec(cnt,:) = nodec(interior_loc(j)).we;
        wed(cnt,:) = noded(interior_loc(j)).we;
        weg(cnt,:) = nodeg(interior_loc(j)).we;
        wev(cnt,:) = nodev(interior_loc(j)).we;
        %if isnan(wec(cnt,1)),
%         if sum(wec(cnt,:)<0)>0, %isnan(wec(cnt,1))||isnan(wec(cnt,2)),%sum(wec(cnt,:)<0)>0,
%             break,
%         end
    end
%     if sum(wec(cnt,:)<0)>0, %isnan(wec(cnt,1))||isnan(wec(cnt,2)),%sum(wec(cnt,:)<0)>0,
%         break
%     end
end
delete(hd);

[pa,pp] = perfbound_beta(pgt,we);
[pac,pp] = perfbound_beta(pgt,wec);
[pad,pp] = perfbound_beta(pgt,wed);
[pag,pp] = perfbound_beta(pgt,weg);
[pav,pp] = perfbound_beta(pgt,wev);

hd = plot(pp,pa,'b-',[pp 1],[pac 1],'r-',[pp 1],[pag 1],'g-',[pp 1],[pad 1],'c-',[pp 1],[pav 1],'m',[pp 1],[pp 1],'k:');

% plot(pp,pa,'b-',[pp 1],[pac 1],'r-',[pp 1],[pag 1],'g-',[pp 1],[pad 1],'c-',[pp 1],[pp 1],'k:');

set(hd,'linewidth',2);
set(gca,'linewidth',2,'fontsize',18,'fontweight','bold');
xlabel('Desired Confidence');
ylabel('Actual Confidence');
hl = legend('SLBN','Credal','GBT','Direct','VBS');
set(hl,'location','northwest');


pe = we*[1;0;0.5];
er = sqrt(mean((pe-pgt).^2));
sig = sqrt(mean(pe.*(1-pe).*we(:,3)./(2+we(:,3))));

pec = wec*[1;0;0.5];
erc = sqrt(mean((pec-pgt).^2));
sigc = sqrt(mean(pec.*(1-pec).*wec(:,3)./(2+wec(:,3))));

ped = wed*[1;0;0.5];
erd = sqrt(mean((ped-pgt).^2));
sigd = sqrt(mean(ped.*(1-ped).*wed(:,3)./(2+wed(:,3))));

peg = weg*[1;0;0.5];
erg = sqrt(mean((peg-pgt).^2));
sigg = sqrt(mean(peg.*(1-peg).*weg(:,3)./(2+weg(:,3))));

pev = wev*[1;0;0.5];
erv = sqrt(mean((pev-pgt).^2));
sigv = sqrt(mean(pev.*(1-pev).*wev(:,3)./(2+wev(:,3))));

disp(sprintf('SLBN \t Credal \t Direct \t GBT \t\t VBS'));
disp('Actual RMSE');
disp(sprintf('%5.3f \t %5.3f \t\t %5.3f \t\t %5.3f \t\t %5.3f',er,erc,erd,erg,erv));
disp('Predicted Error');
disp(sprintf('%5.3f \t %5.3f \t\t %5.3f \t\t %5.3f \t\t %5.3f',sig,sigc,sigd,sigg,sigv));
disp('Average Run Time (ms)');
disp(sprintf('%5.3f \t %5.3f \t\t %5.3f \t\t %5.3f \t %5.3f',mean(ts*1000),mean(tc*1000),mean(td*1000),mean(tg*1000),mean(tv*1000)));



    
    
    
