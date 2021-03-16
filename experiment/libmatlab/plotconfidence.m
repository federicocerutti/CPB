function plotconfidence(pa)

pp = 0:.01:.99;
hd = plot(pp,pa,'b-',[pp 1],[pp 1],'k:');


% plot(pp,pa,'b-',[pp 1],[pac 1],'r-',[pp 1],[pag 1],'g-',[pp 1],[pad 1],'c-',[pp 1],[pp 1],'k:');

set(hd,'linewidth',2);
set(gca,'linewidth',2,'fontsize',18,'fontweight','bold');
xlabel('Desired Confidence');
ylabel('Actual Confidence');
%hl = legend('SLBN','Credal','GBT','Direct','VBS');
%set(hl,'location','northwest');