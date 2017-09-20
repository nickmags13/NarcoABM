%%%%%%%%%%%%% Get Results %%%%%%%%%%%%%%%%%%%%

cd X:\model_results\SupplyChain_091817
fnames=dir;
fnamescell=struct2cell(fnames);
h=strncmp('supplychain_results_',fnamescell(1,:),20);
hind=find(h==1);

TMAX=180;
MRUNS=30;
ERUNS=5;
batchind=[reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1) ...
    reshape(repmat(1:MRUNS,1,ERUNS),MRUNS*ERUNS,1)];

sl_max=[160 240 320 400 480];
sl_min=ceil(sl_max/6);

NACTNODES=cell(MRUNS*ERUNS,1);
SLPEREVENT=cell(MRUNS*ERUNS,1);
dtoBDGT=cell(MRUNS*ERUNS,1,2);


for mr=1:length(hind)   % MRUNS*EXPTRUNS
    h=strcmp(sprintf('supplychain_results_091817_%d_%d.mat',...
        batchind(mr,1),batchind(mr,2)),fnamescell(1,:));
    filename=fnamescell{1,h};
    load(filename)
    
    NACTNODES(mr)=mat2cell(nactnodes,1,TMAX);
    SLPEREVENT(mr)=mat2cell(slperevent,1,TMAX);
    dtoBDGT(mr,1,1)=mat2cell(DTOBDGT(1,:),1,TMAX);
    dtoBDGT(mr,1,2)=mat2cell(DTOBDGT(2,:),1,TMAX);
    
    clear nactnodes slperevent
end

slevents=cell2mat(SLPEREVENT);
actrtes=cell2mat(NACTNODES);
dtobdgt_1=cell2mat(dtoBDGT(:,:,1));
dtobdgt_2=cell2mat(dtoBDGT(:,:,2));

meanrtes=zeros(length(hind),1);
medianrtes=zeros(length(hind),1);
minrtes=zeros(length(hind),1);
maxrtes=zeros(length(hind),1);
cvrtes=zeros(length(hind),1);
meansl=zeros(length(hind),1);
mediansl=zeros(length(hind),1);
minsl=zeros(length(hind),1);
maxsl=zeros(length(hind),1);
cvsl=zeros(length(hind),1);
maxdto=zeros(length(hind),2);

for mr=1:length(hind)
    istart=find(slevents(mr,:)>0,1,'first');
    meanrtes(mr)=mean(actrtes(mr,istart:TMAX));
    medianrtes(mr)=median(actrtes(mr,istart:TMAX));
    minrtes(mr)=min(actrtes(mr,istart:TMAX));
    maxrtes(mr)=max(actrtes(mr,istart:TMAX));
    cvrtes(mr)=var(actrtes(mr,istart:TMAX))/mean(actrtes(mr,istart:TMAX));
    meansl(mr)=mean(slevents(mr,istart:TMAX));
    mediansl(mr)=median(slevents(mr,istart:TMAX));
    minsl(mr)=min(slevents(mr,istart:TMAX));
    maxsl(mr)=max(slevents(mr,istart:TMAX));
    cvsl(mr)=var(slevents(mr,istart:TMAX))/mean(slevents(mr,istart:TMAX));
    maxdto(mr,1)=max(dtobdgt_1(mr,istart:TMAX));
    maxdto(mr,2)=max(dtobdgt_2(mr,istart:TMAX));
end

summroutes=cell(length(sl_max),5);
summsl=cell(length(sl_max),5);
for g=1:length(sl_max)
    ind=batchind(:,1) == g;
    
    [rtemu,rtesigma,rtemuci,~]=normfit(medianrtes(ind));
    rtequant=quantile(medianrtes(ind),[0.025 0.25 0.5 0.75 0.975]);
    summroutes{g,1}=rtemu;
    summroutes{g,2}=rtesigma;
    summroutes(g,3)=mat2cell(rtemuci,2,1);
    summroutes{g,4}=rtequant(3);
    summroutes(g,5)=mat2cell([rtequant(2) rtequant(4)]',2,1);
    
    [slmu,slsigma,slmuci,~]=normfit(mediansl(ind));
    slquant=quantile(mediansl(ind),[0.025 0.25 0.5 0.75 0.975]);
    summsl{g,1}=slmu;
    summsl{g,2}=slsigma;
    summsl(g,3)=mat2cell(slmuci,2,1);
    summsl{g,4}=slquant(3);
    summsl(g,5)=mat2cell([slquant(2) slquant(4)]',2,1);
end

% % h2_1=figure;
% set(h2_1,'Color','white')
% [hAx,hl1,hl2]=plotyy(1:TMAX,actrtes(61,:),1:TMAX,slevents(61,:));
% ylabel(hAx(1),'Number of Routes')
% ylabel(hAx(2),'Average S&L Volume (kg)')
% xlim(hAx(1),[1 TMAX])
% xlim(hAx(2),[1 TMAX])
% xlabel('Month')
% legend('Active Routes','S&L Volume','Orientation','horizontal','Location','southoutside')