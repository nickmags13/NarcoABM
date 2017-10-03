%%%%%%%%%%%%% Get Results %%%%%%%%%%%%%%%%%%%%

cd X:\model_results\SupplyChain_091817
fnames=dir;
fnamescell=struct2cell(fnames);
h=strncmp('supplychain_results_',fnamescell(1,:),20);
hind=find(h==1);

TSTART=1;
TMAX=180;
MRUNS=30;
ERUNS=5;
ndto=2;
batchind=[reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1) ...
    reshape(repmat(1:MRUNS,1,ERUNS),MRUNS*ERUNS,1)];

sl_max=[160 240 320 400 480];
sl_min=ceil(sl_max/6);

% NACTNODES=cell(ndto,MRUNS*ERUNS);
% SLPEREVENT=cell(ndto,MRUNS*ERUNS);
NACTNODES=zeros(MRUNS*ERUNS,TMAX,ndto);
SLPEREVENT=zeros(MRUNS*ERUNS,TMAX,ndto);
dtoBDGT=cell(MRUNS*ERUNS,1,ndto);
TMOV=cell(MRUNS*ERUNS,1);


for mr=1:length(hind)   % MRUNS*EXPTRUNS
    h=strcmp(sprintf('supplychain_results_091817_%d_%d.mat',...
        batchind(mr,1),batchind(mr,2)),fnamescell(1,:));
    filename=fnamescell{1,h};
    load(filename)
    
    %     NACTNODES(:,mr)=mat2cell(nactnodes,ndto,TMAX);
    %     SLPEREVENT(:,mr)=mat2cell(slperevent,ndto,TMAX);
    for idt=1:ndto
        NACTNODES(mr,:,idt)=nactnodes(idt,:);
        SLPEREVENT(mr,:,idt)=slperevent(idt,:);
    end
    dtoBDGT(mr,1,1)=mat2cell(DTOBDGT(1,:),1,TMAX);
    dtoBDGT(mr,1,2)=mat2cell(DTOBDGT(2,:),1,TMAX);
    TMOV(mr,1)=mat2cell(t_firstmov',1,length(t_firstmov));
    
    clear nactnodes slperevent
end


% slevents=cell2mat(SLPEREVENT);
% actrtes=cell2mat(NACTNODES);
slevents=SLPEREVENT;
actrtes=NACTNODES;
dtobdgt_1=cell2mat(dtoBDGT(:,:,1));
dtobdgt_2=cell2mat(dtoBDGT(:,:,2));
tmov=cell2mat(TMOV);

meanrtes=zeros(ndto,length(hind));
medianrtes=zeros(ndto,length(hind));
minrtes=zeros(ndto,length(hind));
maxrtes=zeros(ndto,length(hind));
cvrtes=zeros(ndto,length(hind));
meansl=zeros(ndto,length(hind));
mediansl=zeros(ndto,length(hind));
minsl=zeros(ndto,length(hind));
maxsl=zeros(ndto,length(hind));
cvsl=zeros(ndto,length(hind));
maxdto=zeros(length(hind),2);
meanfirstmov=zeros(length(hind),1);

for mr=1:length(hind)
    for idt=1:ndto
        istart=find(slevents(mr,:,idt)>0,1,'first');
        meanrtes(idt,mr)=mean(actrtes(mr,istart:TMAX,idt));
        medianrtes(idt,mr)=median(actrtes(mr,istart:TMAX,idt));
        minrtes(idt,mr)=min(actrtes(mr,istart:TMAX,idt));
        maxrtes(idt,mr)=max(actrtes(mr,istart:TMAX,idt));
        cvrtes(idt,mr)=var(actrtes(mr,istart:TMAX,idt))/mean(actrtes(mr,istart:TMAX,idt));
        meansl(idt,mr)=mean(slevents(mr,istart:TMAX,idt));
        mediansl(idt,mr)=median(slevents(mr,istart:TMAX,idt));
        minsl(idt,mr)=min(slevents(mr,istart:TMAX,idt));
        maxsl(idt,mr)=max(slevents(mr,istart:TMAX,idt));
        cvsl(idt,mr)=var(slevents(mr,istart:TMAX,idt))/mean(slevents(mr,istart:TMAX,idt));
    end
    maxdto(mr,1)=max(dtobdgt_1(mr,istart:TMAX));
    maxdto(mr,2)=max(dtobdgt_2(mr,istart:TMAX));
    meanfirstmov(mr)=mean(tmov(mr,:));
end

summroutes=cell(length(sl_max),5);
summsl=cell(length(sl_max),5);
for g=1:length(sl_max)
    ind=batchind(:,1) == g;
    
    [rtemu,rtesigma,rtemuci,~]=normfit(sum(medianrtes(:,ind),1));
    rtequant=quantile(sum(medianrtes(:,ind),1),[0.025 0.25 0.5 0.75 0.975]);
    summroutes{g,1}=rtemu;
    summroutes{g,2}=rtesigma;
    summroutes(g,3)=mat2cell(rtemuci,2,1);
    summroutes{g,4}=rtequant(3);
    summroutes(g,5)=mat2cell([rtequant(2) rtequant(4)]',2,1);
    
    [slmu,slsigma,slmuci,~]=normfit(sum(mediansl(:,ind),1));
    slquant=quantile(sum(mediansl(:,ind),1),[0.025 0.25 0.5 0.75 0.975]);
    summsl{g,1}=slmu;
    summsl{g,2}=slsigma;
    summsl(g,3)=mat2cell(slmuci,2,1);
    summsl{g,4}=slquant(3);
    summsl(g,5)=mat2cell([slquant(2) slquant(4)]',2,1);
end

% % h2_1=figure;
% set(h2_1,'Color','white')
% subplot(2,1,1)
% [hAx,hl1,hl2]=plotyy(1:TMAX,nactnodes(1,:),1:TMAX,slperevent(1,:));
% ylabel(hAx(1),'Number of Routes')
% ylabel(hAx(2),'Average S&L Volume (kg)')
% xlim(hAx(1),[1 TMAX])
% xlim(hAx(2),[1 TMAX])
% legend('Active Routes','S&L Volume','Orientation','vertical','Location','NorthWest')
% subplot(2,1,2)
% [hAx2,h21,h22]=plotyy(1:TMAX,nactnodes(2,:),1:TMAX,slperevent(2,:));
% ylabel(hAx2(1),'Number of Routes')
% ylabel(hAx2(2),'Average S&L Volume (kg)')
% xlim(hAx2(1),[1 TMAX])
% xlim(hAx2(2),[1 TMAX])
% xlabel('Month')

% [hAx,hl1,hl2]=plotyy(1:TMAX,sum(actrtes(61,:,:),3)),1:TMAX,sum(slevents(61,:,:),3));
% ylabel(hAx(1),'Number of Routes')
% ylabel(hAx(2),'Average S&L Volume (kg)')
% xlim(hAx(1),[1 TMAX])
% xlim(hAx(2),[1 TMAX])
% xlabel('Month')
% legend('Active Routes','S&L Volume','Orientation','horizontal','Location','southoutside')

h1=figure;
set(h1,'color','white','Visible','off')
boxplot(mediansl,batchind(:,1),'color','r','labels',...
    {'Low','Low Mod.','Mod.','High Mod.','High'},...
    'datalim',[0 6],'extrememode','compress')
xlabel('Interdiction Capacity')
ylabel('Median Annual Volume Seized')
saveas(h1,'Annual_SL_varcptcy','png')

h1_1=figure;
set(h1_1,'color','white','Visible','off')
boxplot(medianrtes,batchind(:,1),'color','b','labels',...
    {'Low','Low Mod.','Mod.','High Mod.','High'})
xlabel('Interdiction Capacity')
ylabel('Median Annual Number of Active Routes')
saveas(h1_1,'Annual_Routes_varcpcty','png')


h3=figure;
set(h3,'color','white','Visible','off')
geoshow(CAadm0,'FaceColor',[1 1 1])
hold on
for n=2:nnodes-1
    if t_firstmov(n) == 0
        continue
    else
    plot(nodelon(n),nodelat(n),'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',...
        [1-t_firstmov(n)/max(t_firstmov) t_firstmov(n)/max(t_firstmov) 0])
    end
end
title('First Time Step of Cocaine Movements')
xlabel('Longitude')
ylabel('Latitude')
saveas(h3,'Earliest_Movement','png')

%%% Locate Department statistics
cntrynames={'Guatemala','Panama','Panama','Costa Rica','Honduras','Panama','Nicaragua'};
deptnames={'Pet','Dari','Ember','Puntarenas','Grac','Col','Atlantico Norte'};
cntryind=zeros(length(CAattr1),length(deptnames));
deptind=zeros(length(CAattr1),length(deptnames));
hind=zeros(1,length(deptnames));
slctnodes=cell(1,length(deptnames));
deptflows=zeros(length(deptnames),TMAX);
deptflows_ts=zeros(length(deptnames),TMAX/12);
tsind=ceil((1:TMAX)./12);
for jj=1:length(cntrynames)
    for ii=1:length(CAattr1)
    cntryind(ii,jj)=strncmp(cntrynames(jj),CAattr1(ii).ADM0_NAME,length(cntrynames{jj}));
    deptind(ii,jj)=strncmp(deptnames(jj),CAattr1(ii).ADM1_NAME,length(deptnames{jj}));
    end
    hind(jj)=find(cntryind(:,jj) == 1 & deptind(:,jj) == 1);
    nodeids=NodeTable.ID(NodeTable.DeptCode == CAattr1(hind(jj)).ADM1_CODE);
    slctnodes(jj)=mat2cell(nodeids,length(nodeids),1);
    deptflows(jj,:)=sum(OUTFLOW(slctnodes{jj},:),1);
    for ts=1:TMAX/12
        deptflows_ts(jj,ts)=sum(deptflows(jj,tsind==ts));
    end
end





