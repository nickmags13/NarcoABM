%%%%%%%%%%%%% Get Results %%%%%%%%%%%%%%%%%%%%

cd C:\Users\nrmagliocca\'Box Sync'\'Data Drive'\model_results\SupplyChain_full_011618
fnames=dir;
fnamescell=struct2cell(fnames);
h=strncmp('supplychain_results_',fnamescell(1,:),20);
hind=find(h==1);

[CAadm0,CAattr0]=shaperead('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmericaData\GADM\g2015_2014_0\CAadm0.shp',...
            'UseGeoCoords',true);
[CAadm1,CAattr1]=shaperead('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmericaData\GADM\g2015_2014_1\CAadm1.shp',...
            'UseGeoCoords',true);

load(fnamescell{1,hind(1)})
NNODES=length(NodeTable.ID);

TSTART=1;
TMAX=180;
MRUNS=30;
ERUNS=1;
ndto=2;
batchind=[reshape(repmat(1:ERUNS,MRUNS,1),MRUNS*ERUNS,1) ...
    reshape(repmat(1:MRUNS,1,ERUNS),MRUNS*ERUNS,1)];

% sl_max=[80 160 240 320 400];
sl_max=400;
sl_min=ceil(sl_max/6);

% NACTNODES=cell(ndto,MRUNS*ERUNS);
% SLPEREVENT=cell(ndto,MRUNS*ERUNS);
NACTNODES=zeros(MRUNS*ERUNS,TMAX,ndto);
SLPEREVENT=zeros(MRUNS*ERUNS,TMAX,ndto);
VALPEREVENT=zeros(MRUNS*ERUNS,TMAX,ndto);
dtoBDGT=cell(MRUNS*ERUNS,1,ndto);
TMOV=cell(MRUNS*ERUNS,1);

SENDFLOW=zeros(NNODES,TMAX,MRUNS*ERUNS);
SLVOL=zeros(NNODES,TMAX,MRUNS*ERUNS);
INTRATE=zeros(NNODES,TMAX,MRUNS*ERUNS);
DEPTFLOWS=cell(7,length(hind));
CNTRYFLOWS=cell(7,length(hind));
    
for mr=1:length(hind)   % MRUNS*EXPTRUNS
    h=strcmp(sprintf('supplychain_results_011618_%d_%d.mat',...
        batchind(mr,1),batchind(mr,2)),fnamescell(1,:));
    filename=fnamescell{1,h};
    load(filename)
    
    % %%% Locate Department statistics
    cntrynames={'Guatemala','Panama','Panama','Costa Rica','Honduras','Panama','Nicaragua'};
    deptnames={'Pet','Dari','Ember','Puntarenas','Grac','Col','Atlantico Norte'};
    cntryind=zeros(length(CAattr1),length(deptnames));
    deptind=zeros(length(CAattr1),length(deptnames));
    hind=zeros(1,length(deptnames));
    chind=cell(1,length(deptnames));
    slctnodes=cell(1,length(deptnames));
    cslctnodes=cell(1,length(deptnames));
    deptflows=zeros(length(deptnames),TMAX);
    deptflows_ts=zeros(length(deptnames),TMAX/12);
    cntryflows=zeros(length(deptnames),TMAX);
    cntryflows_ts=zeros(length(deptnames),TMAX/12);
    deptslvol=zeros(length(deptnames),TMAX);
    deptslvol_ts=zeros(length(deptnames),TMAX/12);
    deptintrt=zeros(length(deptnames),TMAX);
    deptintrt_ts=zeros(length(deptnames),TMAX/12);
    tsind=ceil((1:TMAX)./12);
    for jj=1:length(cntrynames)
        for ii=1:length(CAattr1)
            cntryind(ii,jj)=strncmp(cntrynames(jj),CAattr1(ii).ADM0_NAME,length(cntrynames{jj}));
            deptind(ii,jj)=strncmp(deptnames(jj),CAattr1(ii).ADM1_NAME,length(deptnames{jj}));
        end
        hind(jj)=find(cntryind(:,jj) == 1 & deptind(:,jj) == 1);
        chind(jj)=mat2cell(find(cntryind(:,jj) == 1),length(find(cntryind(:,jj) == 1)),1);
        nodeids=NodeTable.ID(NodeTable.DeptCode == CAattr1(hind(jj)).ADM1_CODE);
        cnodeids=NodeTable.ID(ismember(NodeTable.DeptCode,cat(1,CAattr1(chind{jj}).ADM1_CODE)));
        slctnodes(jj)=mat2cell(nodeids,length(nodeids),1);
        cslctnodes(jj)=mat2cell(cnodeids,length(cnodeids),1);
        deptflows(jj,:)=sum(OUTFLOW(slctnodes{jj},:),1);    %single run code
        cntryflows(jj,:)=sum(OUTFLOW(cslctnodes{jj},:),1);
        test=reshape(sum(slsuccess,1),NNODES,TMAX); %single run code
        deptslvol(jj,:)=sum(test(slctnodes{jj},:),1);  %single run code
        %     deptflows(jj,:)=sum(mediandlvr(slctnodes{jj},:),1);
        %     deptslvol(jj,:)=sum(medianslvol(slctnodes{jj},:),1);
        deptintrt(jj,:)=deptslvol(jj,:)./(deptslvol(jj,:)+deptflows(jj,:));
        for ts=1:TMAX/12
            deptflows_ts(jj,ts)=sum(deptflows(jj,tsind==ts));
            cntryflows_ts(jj,ts)=sum(cntryflows(jj,tsind==ts));
            deptslvol_ts(jj,ts)=sum(deptslvol(jj,tsind==ts));
            deptintrt_ts(jj,:)=deptslvol_ts(jj,:)./(deptslvol_ts(jj,:)+deptflows_ts(jj,:));
        end
        DEPTFLOWS(jj,mr)=mat2cell(deptflows_ts(jj,:),1,ts);
        CNTRYFLOWS(jj,mr)=mat2cell(cntryflows_ts(jj,:),1,ts);
    end
    
    SENDFLOW(:,:,mr)=OUTFLOW;
    SLVOL(:,:,mr)=reshape(sum(slsuccess,1),NNODES,TMAX);
    INTRATE(:,:,mr)=SLVOL(:,:,mr)./(SLVOL(:,:,mr)+SENDFLOW(:,:,mr));
    %     NACTNODES(:,mr)=mat2cell(nactnodes,ndto,TMAX);
    %     SLPEREVENT(:,mr)=mat2cell(slperevent,ndto,TMAX);
    for idt=1:ndto
        NACTNODES(mr,:,idt)=nactnodes(idt,:);
        SLPEREVENT(mr,:,idt)=slperevent(idt,:);
        VALPEREVENT(mr,:,idt)=slval(idt,:);
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
totactrtes=sum(NACTNODES,3);
totslevents=sum(SLPEREVENT,3);
totvalevents=sum(VALPEREVENT./max(SLPEREVENT,1),3);
dtobdgt_1=cell2mat(dtoBDGT(:,:,1));
dtobdgt_2=cell2mat(dtoBDGT(:,:,2));
tmov=cell2mat(TMOV);

meanrtes=zeros(ndto,TMAX,ERUNS);
medianrtes=zeros(ndto,TMAX,ERUNS);
minrtes=zeros(ndto,TMAX,ERUNS);
maxrtes=zeros(ndto,TMAX,ERUNS);
cvrtes=zeros(ndto,TMAX,ERUNS);
meansl=zeros(ndto,TMAX,ERUNS);
mediansl=zeros(ndto,TMAX,ERUNS);
minsl=zeros(ndto,TMAX,ERUNS);
maxsl=zeros(ndto,TMAX,ERUNS);
cvsl=zeros(ndto,TMAX,ERUNS);
maxdto=zeros(length(hind),2);
meanfirstmov=zeros(ERUNS,NNODES);

medianslvol=zeros(NNODES,TMAX,ERUNS);
meanslvol=zeros(NNODES,TMAX,ERUNS);
mediandlvr=zeros(NNODES,TMAX,ERUNS);
meandlvr=zeros(NNODES,TMAX,ERUNS);
medianintrt=zeros(NNODES,TMAX,ERUNS);
meanintrt=zeros(NNODES,TMAX,ERUNS);
slvolQnt=cell(7,ERUNS);    % quantiles [0.025 0.05 0.25 0.5 0.75 0.95 0.975] for NNODES x TMAX
dlvrQnt=cell(7,ERUNS);
intrtQnt=cell(7,ERUNS);

for erun=1:ERUNS
    medianslvol(:,:,erun)=median(SLVOL(:,:,batchind(:,1)==erun),3); %volume seized
    meanslvol(:,:,erun)=mean(SLVOL(:,:,batchind(:,1)==erun),3); %volume seized
    mediandlvr(:,:,erun)=median(SENDFLOW(:,:,batchind(:,1)==erun),3); %volume delivered
    meandlvr(:,:,erun)=mean(SENDFLOW(:,:,batchind(:,1)==erun),3); %volume delivered
    medianintrt(:,:,erun)=median(INTRATE(:,:,batchind(:,1)==erun),3); %interdiction rate
    meanintrt(:,:,erun)=mean(INTRATE(:,:,batchind(:,1)==erun),3); %interdiction rate
    
    slvol_qnt=quantile(SLVOL(:,:,batchind(:,1)==erun),[0.025 0.05 0.25 0.5 0.75 0.95 0.975],3);
    dlvr_qnt=quantile(SENDFLOW(:,:,batchind(:,1)==erun),[0.025 0.05 0.25 0.5 0.75 0.95 0.975],3);
    intrt_qnt=quantile(INTRATE(:,:,batchind(:,1)==erun),[0.025 0.05 0.25 0.5 0.75 0.95 0.975],3);
    
    for k=1:7
        slvolQnt(k,erun)=mat2cell(slvol_qnt(:,:,k),NNODES,TMAX);
        dlvrQnt(k,erun)=mat2cell(slvol_qnt(:,:,k),NNODES,TMAX);
        intrtQnt(k,erun)=mat2cell(slvol_qnt(:,:,k),NNODES,TMAX);
    end

    for idt=1:ndto
        %         istart=find(slevents(mr,:,idt)>0,1,'first');
        %         meanrtes(mr,:,idt)=mean(actrtes(mr,istart:TMAX,idt));
        %         medianrtes(idt,mr)=median(actrtes(mr,istart:TMAX,idt));
        %         minrtes(idt,mr)=min(actrtes(mr,istart:TMAX,idt));
        %         maxrtes(idt,mr)=max(actrtes(mr,istart:TMAX,idt));
        %         cvrtes(idt,mr)=var(actrtes(mr,istart:TMAX,idt))/mean(actrtes(mr,istart:TMAX,idt));
        %         meansl(idt,mr)=mean(slevents(mr,istart:TMAX,idt));
        %         mediansl(idt,mr)=median(slevents(mr,istart:TMAX,idt));
        %         minsl(idt,mr)=min(slevents(mr,istart:TMAX,idt));
        %         maxsl(idt,mr)=max(slevents(mr,istart:TMAX,idt));
        %         cvsl(idt,mr)=var(slevents(mr,istart:TMAX,idt))/mean(slevents(mr,istart:TMAX,idt));

        meanrtes(idt,:,erun)=mean(actrtes(batchind(:,1)==erun,:,idt),1);
        medianrtes(idt,:,erun)=median(actrtes(batchind(:,1)==erun,:,idt));
        minrtes(idt,:,erun)=min(actrtes(batchind(:,1)==erun,:,idt));
        maxrtes(idt,:,erun)=max(actrtes(batchind(:,1)==erun,:,idt));
        cvrtes(idt,:,erun)=var(actrtes(batchind(:,1)==erun,:,idt))./mean(actrtes(batchind(:,1)==erun,:,idt));
        meansl(idt,:,erun)=mean(slevents(batchind(:,1)==erun,:,idt));
        mediansl(idt,:,erun)=median(slevents(batchind(:,1)==erun,:,idt));
        minsl(idt,:,erun)=min(slevents(batchind(:,1)==erun,:,idt));
        maxsl(idt,:,erun)=max(slevents(batchind(:,1)==erun,:,idt));
        cvsl(idt,:,erun)=var(slevents(batchind(:,1)==erun,:,idt))./mean(slevents(batchind(:,1)==erun,:,idt));
    end
    
    for j=1:length(tmov(1,:))
        isubmov=tmov(batchind(:,1)==erun,j) ~= 0;
        meanfirstmov(erun,j)=mean(tmov(isubmov,j),1);
    end
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

        
%%
deptrefvecs{1}=[5500 6300 2250 3450 1050 1375 625 625 1600];    %Peten, 2005-2013
deptrefvecs{2}=[29750 26812 37034 49043 41399 41535];   %Darien, 2009-2014
% deptrefvecs{3}=[4000 886 32782 26225.6 11400 25765 28684 34926];    %Puntarenas, 2007-2014
deptrefvecs{3}=[770 1275 3450 5000 23450 23779 92661 159368 114905 92472 ...
    73836 191770 209586 269160];    %Costa Rica, 2001-2014
deptrefvecs{4}=[3930 12800 4463 6906 63395 109875 219354 214087 107219 55268]; %gracias a dios, 2005-2014
% deptrefvecs{4}=[1550 12800 6075 7872 22500 67802 83100 66305 64756 29033];
deptrefvecs{5}=[560 6766 7045 17690 9850 3715 6100 25258];  %colon, 2006-2013
deptrefvecs{6}=[5575 775 32750 37650 48038 9400 9117 3970]; %atlantico norte, 2007-2014
cntryrefvecs{1}=[20776 16324 56392 23744 23002 39326 65296 59360 78652 146916 102396 76426 81620 152852]; %Guatemal, 2001-2014

% Peten Flows
h10=figure;
set(h10,'color','white')
labels={'05','06','07','08','09','10','11','12','13'};
subflows=cat(1,DEPTFLOWS{1,:});
p1=subflows(:,5:13);
p2=deptrefvecs{1};
line1=parallelcoords(p1,'labels',labels,'quantile',0.25);
line1(1).LineWidth=2;
line1(2).LineWidth=1.5;
line1(3).LineWidth=1.5;
hold on
line2=parallelcoords(p2);
hold off
line2.Color='red';
line2.LineWidth=2;
lgd=legend([line1(1) line2],'Model Output','CCDB Dept. Data','Location',...
    'NorthWest');
lgd.FontSize=14;
set(gca,'FontSize',12)
ylabel('Cocaine Shipment Volume(kg)')
xlabel('Year')
saveas(h10,'Flows_Peten_coords','png')

% Darien Flows
h10_1=figure;
set(h10_1,'color','white')
labels={'09','10','11','12','13','14'};
subflows=cat(1,DEPTFLOWS{2,:});
p1=subflows(:,9:14);
p2=deptrefvecs{2};
line1=parallelcoords(p1,'labels',labels,'quantile',0.25);
line1(1).LineWidth=2;
line1(2).LineWidth=1.5;
line1(3).LineWidth=1.5;
hold on
line2=parallelcoords(p2);
hold off
line2.Color='red';
line2.LineWidth=2;
lgd=legend([line1(1) line2],'Model Output','CCDB Dept. Data','Location',...
    'NorthWest');
lgd.FontSize=12;
set(gca,'FontSize',14)
ylabel('Cocaine Shipment Volume(kg)')
xlabel('Year')
saveas(h10_1,'Flows_Darien_coords','png')

%Costa Rica Flows
h10_2=figure;
set(h10_2,'color','white')
labels={'01','02','03','04','05','06','07','08','09','10','11','12','13','14'};
subflows=cat(1,CNTRYFLOWS{4,:});
p1=subflows(:,1:14);
p2=deptrefvecs{3};
line1=parallelcoords(p1,'labels',labels,'quantile',0.25);
line1(1).LineWidth=2;
line1(2).LineWidth=1.5;
line1(3).LineWidth=1.5;
hold on
line2=parallelcoords(p2);
hold off
line2.Color='red';
line2.LineWidth=2;
lgd=legend([line1(1) line2],'Model Output','CCDB Dept. Data','Location',...
    'NorthWest');
lgd.FontSize=12;
set(gca,'FontSize',14)
ylabel('Cocaine Shipment Volume(kg)')
xlabel('Year')
saveas(h10_2,'Flows_CostaRica_coords','png')

% Gracias Flows
h10_3=figure;
set(h10_3,'color','white')
labels={'05','06','07','08','09','10','11','12','13','14'};
subflows=cat(1,DEPTFLOWS{5,:});
p1=subflows(:,5:14);
p2=deptrefvecs{4};
line1=parallelcoords(p1,'labels',labels,'quantile',0.25);
line1(1).LineWidth=2;
line1(2).LineWidth=1.5;
line1(3).LineWidth=1.5;
hold on
line2=parallelcoords(p2);
hold off
line2.Color='red';
line2.LineWidth=2;
lgd=legend([line1(1) line2],'Model Output','CCDB Dept. Data','Location',...
    'NorthWest');
lgd.FontSize=12;
set(gca,'FontSize',14)
ylabel('Cocaine Shipment Volume(kg)')
xlabel('Year')
saveas(h10_3,'Flows_Gracias_coords','png')

% Norte Flows
h10_4=figure;
set(h10_4,'color','white')
labels={'07','08','09','10','11','12','13','14'};
subflows=cat(1,DEPTFLOWS{7,:});
p1=subflows(:,7:14);
p2=deptrefvecs{6};
line1=parallelcoords(p1,'labels',labels,'quantile',0.25);
line1(1).LineWidth=2;
line1(2).LineWidth=1.5;
line1(3).LineWidth=1.5;
hold on
line2=parallelcoords(p2);
hold off
line2.Color='red';
line2.LineWidth=2;
lgd=legend([line1(1) line2],'Model Output','CCDB Dept. Data','Location',...
    'NorthWest');
lgd.FontSize=12;
set(gca,'FontSize',14)
ylabel('Cocaine Shipment Volume(kg)')
xlabel('Year')
saveas(h10_4,'Flows_Norte_coords','png')

% %%% Single run
% % % h2_1=figure;
% % set(h2_1,'Color','white')
subplot(2,1,1)
[hAx,hl1,hl2]=plotyy(1:TMAX,nactnodes(1,:),1:TMAX,slperevent(1,:));
ylabel(hAx(1),'Number of Routes')
ylabel(hAx(2),'Average S&L Volume (kg)')
% [hAx,hl1,hl2]=plotyy(1:TMAX,nactnodes(1,:),1:TMAX,slval(1,:)./1000000);
% ylabel(hAx(1),'Number of Routes')
% ylabel(hAx(2),'Average S&L Value ($Mil)')
xlim(hAx(1),[1 TMAX])
xlim(hAx(2),[1 TMAX])
legend('Active Routes','S&L Volume','Orientation','vertical','Location','NorthWest')
subplot(2,1,2)
[hAx2,h21,h22]=plotyy(1:TMAX,nactnodes(2,:),1:TMAX,slperevent(2,:));
ylabel(hAx2(1),'Number of Routes')
ylabel(hAx2(2),'Average S&L Volume (kg)')
% [hAx2,h21,h22]=plotyy(1:TMAX,nactnodes(2,:),1:TMAX,slval(2,:)./1000000);
% ylabel(hAx2(1),'Number of Routes')
% ylabel(hAx2(2),'Average S&L Value ($Mil)')
xlim(hAx2(1),[1 TMAX])
xlim(hAx2(2),[1 TMAX])
xlabel('Month')
% 
%%% Value per shipment
subplot(2,1,1)
% [hAx,hl1,hl2]=plotyy(1:TMAX,nactnodes(1,:),1:TMAX,slperevent(1,:));
% ylabel(hAx(1),'Number of Routes')
% ylabel(hAx(2),'Average S&L Volume (kg)')
valkg=slval./slperevent;
valkg(isnan(valkg))=0;
[hAx,hl1,hl2]=plotyy(1:TMAX,nactnodes(1,:),1:TMAX,valkg(1,:));
ylabel(hAx(1),'Number of Routes')
ylabel(hAx(2),'Average S&L Value ($/kilo)')
xlim(hAx(1),[1 TMAX])
xlim(hAx(2),[1 TMAX])
legend('Active Routes','S&L VAlue','Orientation','vertical','Location','NorthWest')
subplot(2,1,2)
% [hAx2,h21,h22]=plotyy(1:TMAX,nactnodes(2,:),1:TMAX,slperevent(2,:));
% ylabel(hAx2(1),'Number of Routes')
% ylabel(hAx2(2),'Average S&L Volume (kg)')
[hAx2,h21,h22]=plotyy(1:TMAX,nactnodes(2,:),1:TMAX,valkg(2,:));
ylabel(hAx2(1),'Number of Routes')
ylabel(hAx2(2),'Average S&L Value ($/kilo)')
xlim(hAx2(1),[1 TMAX])
xlim(hAx2(2),[1 TMAX])
xlabel('Month')
% 
% h3=figure;
% set(h3,'color','white','Visible','off')
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
%red is early, green is late
title('First Time Step of Cocaine Movements')
xlabel('Longitude')
ylabel('Latitude')
% saveas(h3,'Earliest_Movement','png')

%%%% Single run from results file
NNODES=size(activeroute,1);
h3=figure;
set(h3,'color','white','Visible','off')
geoshow(CAadm0,'FaceColor',[1 1 1])
hold on
for n=2:NNODES-1
    if t_firstmov(n) == 0
        continue
    else
    plot(NodeTable.Lon(n),NodeTable.Lat(n),'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',...
        [1-t_firstmov(n)/max(t_firstmov) t_firstmov(n)/max(t_firstmov) 0])
    end
end
title('Average First Time Step of Cocaine Movements')
xlabel('Longitude')
ylabel('Latitude')

% 
% %%% Locate Department statistics
cntrynames={'Guatemala','Panama','Panama','Costa Rica','Honduras','Panama','Nicaragua'};
deptnames={'Pet','Dari','Ember','Puntarenas','Grac','Col','Atlantico Norte'};
cntryind=zeros(length(CAattr1),length(deptnames));
deptind=zeros(length(CAattr1),length(deptnames));
hind=zeros(1,length(deptnames));
chind=cell(1,length(deptnames));
slctnodes=cell(1,length(deptnames));
cslctnodes=cell(1,length(deptnames));
deptflows=zeros(length(deptnames),TMAX);
deptflows_ts=zeros(length(deptnames),TMAX/12);
cntryflows=zeros(length(deptnames),TMAX);
cntryflows_ts=zeros(length(deptnames),TMAX/12);
deptslvol=zeros(length(deptnames),TMAX);
deptslvol_ts=zeros(length(deptnames),TMAX/12);
deptintrt=zeros(length(deptnames),TMAX);
deptintrt_ts=zeros(length(deptnames),TMAX/12);
tsind=ceil((1:TMAX)./12);
for jj=1:length(cntrynames)
    for ii=1:length(CAattr1)
    cntryind(ii,jj)=strncmp(cntrynames(jj),CAattr1(ii).ADM0_NAME,length(cntrynames{jj}));
    deptind(ii,jj)=strncmp(deptnames(jj),CAattr1(ii).ADM1_NAME,length(deptnames{jj}));
    end
    hind(jj)=find(cntryind(:,jj) == 1 & deptind(:,jj) == 1);
    chind(jj)=mat2cell(find(cntryind(:,jj) == 1),length(find(cntryind(:,jj) == 1)),1);
    nodeids=NodeTable.ID(NodeTable.DeptCode == CAattr1(hind(jj)).ADM1_CODE);
    cnodeids=NodeTable.ID(ismember(NodeTable.DeptCode,cat(1,CAattr1(chind{jj}).ADM1_CODE)));
    slctnodes(jj)=mat2cell(nodeids,length(nodeids),1);
    cslctnodes(jj)=mat2cell(cnodeids,length(cnodeids),1);
    deptflows(jj,:)=sum(OUTFLOW(slctnodes{jj},:),1);    %single run code
    cntryflows(jj,:)=sum(OUTFLOW(cslctnodes{jj},:),1);
    test=reshape(sum(slsuccess,1),nnodes,TMAX); %single run code
    deptslvol(jj,:)=sum(test(slctnodes{jj},:),1);  %single run code
%     deptflows(jj,:)=sum(mediandlvr(slctnodes{jj},:),1);
%     deptslvol(jj,:)=sum(medianslvol(slctnodes{jj},:),1);
    deptintrt(jj,:)=deptslvol(jj,:)./(deptslvol(jj,:)+deptflows(jj,:));
    for ts=1:TMAX/12
        deptflows_ts(jj,ts)=sum(deptflows(jj,tsind==ts));
        cntryflows_ts(jj,ts)=sum(cntryflows(jj,tsind==ts));
        deptslvol_ts(jj,ts)=sum(deptslvol(jj,tsind==ts));
        deptintrt_ts(jj,:)=deptslvol_ts(jj,:)./(deptslvol_ts(jj,:)+deptflows_ts(jj,:));
    end
end

h4=figure;
% set(h4,'Color','white','Visible','off')
set(h4,'Color','white')
plot(2005:2013,deptrefvecs{1},'-b','LineWidth',3)
hold on
plot(2005:2013,deptflows_ts(1,5:13),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4,'Flows_Peten','png')

h4_1=figure;
% set(h4_1,'Color','white','Visible','off')
set(h4_1,'Color','white')
plot(2009:2014,deptrefvecs{2},'-b','LineWidth',3)
hold on
plot(2009:2014,deptflows_ts(2,9:14),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4_1,'Flows_Darien','png')

h4_2=figure;
% set(h4_2,'Color','white','Visible','off')
set(h4_2,'Color','white')
plot(2001:2014,deptrefvecs{3},'-b','LineWidth',3)
hold on
plot(2001:2014,cntryflows_ts(4,1:14),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4_2,'Flows_Punta','png')

h4_3=figure;
% set(h4_3,'Color','white','Visible','off')
set(h4_3,'Color','white')
plot(2005:2014,deptrefvecs{4},'-b','LineWidth',3)
hold on
plot(2005:2014,deptflows_ts(5,5:14),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4_3,'Flows_Gracias','png')

h4_4=figure;
set(h4_4,'Color','white','Visible','off')
plot(2006:2013,colon(1,:),'-b','LineWidth',3)
hold on
plot(2006:2013,deptflows_ts(6,6:13),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4_4,'Flows_Colon','png')

h4_5=figure;
% set(h4_5,'Color','white','Visible','off')
set(h4_5,'Color','white')
plot(2007:2014,deptflows_ts(7,7:14),'-r','LineWidth',3)
hold on
plot(2007:2014,deptrefvecs{6},'-b','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('Model Output','CCDB Data','Location','NorthWest')
saveas(h4_5,'Flows_Norte','png')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Multiple runs
for erun=1:ERUNS
    h2_1=figure;
    set(h2_1,'Color','white','Visible','off')
    subplot(2,1,1)
    [hAx,hl1,hl2]=plotyy(1:TMAX,medianrtes(1,:,erun),1:TMAX,mediansl(1,:,erun));
    ylabel(hAx(1),'Number of Routes')
    ylabel(hAx(2),'Average S&L Volume (kg)')
    % [hAx,hl1,hl2]=plotyy(1:TMAX,nactnodes(1,:),1:TMAX,slval(1,:)./1000000);
    % ylabel(hAx(1),'Number of Routes')
    % ylabel(hAx(2),'Average S&L Value ($Mil)')
    xlim(hAx(1),[1 TMAX])
    xlim(hAx(2),[1 TMAX])
    legend('Active Routes','S&L Volume','Orientation','vertical','Location','NorthWest')
    subplot(2,1,2)
    [hAx2,h21,h22]=plotyy(1:TMAX,medianrtes(2,:,erun),1:TMAX,mediansl(2,:,erun));
    ylabel(hAx2(1),'Number of Routes')
    ylabel(hAx2(2),'Average S&L Volume (kg)')
    % [hAx2,h21,h22]=plotyy(1:TMAX,nactnodes(2,:),1:TMAX,slval(2,:)./1000000);
    % ylabel(hAx2(1),'Number of Routes')
    % ylabel(hAx2(2),'Average S&L Value ($Mil)')
    xlim(hAx2(1),[1 TMAX])
    xlim(hAx2(2),[1 TMAX])
    xlabel('Month')
    saveas(h2_1,sprintf('MedianVolume_%d',erun),'png')
end

mediansl_ts=reshape(sum(mediansl,1),ERUNS,TMAX);
medianrtes_ts=reshape(sum(medianrtes,1),ERUNS,TMAX);
h1=figure;
set(h1,'color','white','Visible','off')
boxplot(reshape(median(totslevents,2),MRUNS,ERUNS),1:ERUNS,'color','r','labels',...
    {'Low','Low Mod.','Mod.','High Mod.','High'},'extrememode','compress')
xlabel('Interdiction Capacity')
ylabel('Median Volume Seized (kg)')
saveas(h1,'Annual_SL_varcptcy','png')

valuedata=zeros(MRUNS*ERUNS,1);
for m=1:length(totvalevents(:,1))
    valuedata(m)=median(totvalevents(m,totvalevents(m,:)~=0));
%     valuedata(m)=mean(totvalevents(m,totvalevents(m,:)~=0));
end
h1_2=figure;
set(h1_2,'color','white','Visible','off')
boxplot(reshape(valuedata,MRUNS,ERUNS),1:ERUNS,'color','k','labels',...
    {'Low','Low Mod.','Mod.','High Mod.','High'},'extrememode','compress')
xlabel('Interdiction Capacity')
ylabel('Median Value Seized ($/kg)')
saveas(h1_2,'Annual_VAL_varcptcy','png')


h1_1=figure;
set(h1_1,'color','white','Visible','off')
boxplot(reshape(median(totactrtes,2),MRUNS,ERUNS),1:ERUNS,'color','b','labels',...
    {'Low','Low Mod.','Mod.','High Mod.','High'})
xlabel('Interdiction Capacity')
ylabel('Median Annual Number of Active Routes')
saveas(h1_1,'Annual_Routes_varcpcty','png')


for er=1:ERUNS
h3=figure;
set(h3,'color','white','Visible','off')
geoshow(CAadm0,'FaceColor',[1 1 1])
hold on
for n=2:NNODES-1
    if isnan(meanfirstmov(er,n)) == 1
        continue
    else
    plot(NodeTable.Lon(n),NodeTable.Lat(n),'o','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',...
        [1-meanfirstmov(er,n)/max(meanfirstmov(er,:)) meanfirstmov(er,n)/max(meanfirstmov(er,:)) 0])
    end
end
title('Average First Time Step of Cocaine Movements')
xlabel('Longitude')
ylabel('Latitude')
saveas(h3,sprintf('Earliest_Movement_%d',er),'png')
clf
end

%%% Compare trajectories
idmat=repmat(1:NNODES,MRUNS*ERUNS,1)';
tmovseq=zeros(NNODES,MRUNS,ERUNS);
corrmat=zeros(MRUNS,MRUNS,ERUNS);
for erun=1:ERUNS
    subdata=tmov(batchind(:,1) == erun,:);
    subid=idmat(:,batchind(:,1) == erun);
    for c=1:MRUNS
        datavec=[subdata(c,:)' subid(:,c)];
        vecsorted=sortrows(datavec,1);
        tmovseq(:,c,erun)=vecsorted(:,2);
    end
    [RHO,PVAL]=corr(tmovseq(:,:,erun));
    corrmat(:,:,erun)=RHO;
end

[idx,C]=kmeans(RHO,4);

%%% Locate Department statistics
cntrynames={'Guatemala','Panama','Panama','Costa Rica','Honduras','Panama','Nicaragua'};
deptnames={'Pet','Dari','Ember','Puntarenas','Grac','Col','Atlantico Norte'};
cntryind=zeros(length(CAattr1),length(deptnames));
deptind=zeros(length(CAattr1),length(deptnames));
hind=zeros(1,length(deptnames));
chind=cell(1,length(deptnames));
cslctnodes=cell(1,length(deptnames));
cntryflows=zeros(length(deptnames),TMAX);
cntryflows_ts=zeros(length(deptnames),TMAX/12);
slctnodes=cell(1,length(deptnames));
deptflows=zeros(length(deptnames),TMAX);
deptflows_ts=zeros(length(deptnames),TMAX/12);
deptslvol=zeros(length(deptnames),TMAX);
deptslvol_ts=zeros(length(deptnames),TMAX/12);
deptintrt=zeros(length(deptnames),TMAX);
deptintrt_ts=zeros(length(deptnames),TMAX/12);
tsind=ceil((1:TMAX)./12);
for jj=1:length(cntrynames)
    for ii=1:length(CAattr1)
    cntryind(ii,jj)=strncmp(cntrynames(jj),CAattr1(ii).ADM0_NAME,length(cntrynames{jj}));
    deptind(ii,jj)=strncmp(deptnames(jj),CAattr1(ii).ADM1_NAME,length(deptnames{jj}));
    end
    hind(jj)=find(cntryind(:,jj) == 1 & deptind(:,jj) == 1);
    chind(jj)=mat2cell(find(cntryind(:,jj) == 1),length(find(cntryind(:,jj) == 1)),1);
    nodeids=NodeTable.ID(NodeTable.DeptCode == CAattr1(hind(jj)).ADM1_CODE);
    cnodeids=NodeTable.ID(ismember(NodeTable.DeptCode,cat(1,CAattr1(chind{jj}).ADM1_CODE)));
    slctnodes(jj)=mat2cell(nodeids,length(nodeids),1);
    cslctnodes(jj)=mat2cell(cnodeids,length(cnodeids),1);
%     deptflows(jj,:)=sum(mediandlvr(slctnodes{jj},:),1);
%     deptslvol(jj,:)=sum(medianslvol(slctnodes{jj},:),1);
    deptflows(jj,:)=sum(meandlvr(slctnodes{jj},:),1);
    deptslvol(jj,:)=sum(meanslvol(slctnodes{jj},:),1);
    cntryflows(jj,:)=sum(meandlvr(cslctnodes{jj},:),1);
    deptintrt(jj,:)=deptslvol(jj,:)./(deptslvol(jj,:)+deptflows(jj,:));
    for ts=1:TMAX/12
        deptflows_ts(jj,ts)=sum(deptflows(jj,tsind==ts));
        cntryflows_ts(jj,ts)=sum(cntryflows(jj,tsind==ts));
        deptslvol_ts(jj,ts)=sum(deptslvol(jj,tsind==ts));
        deptintrt_ts(jj,:)=deptslvol_ts(jj,:)./(deptslvol_ts(jj,:)+deptflows_ts(jj,:));
    end
end
%%% CCDB Dept flows data
% peten=[3300 3000 2070 1140 1675 1727 625 1112.5 1600; 2005:2013];
peten=[6100 1600 5500 6300 2250 3450 1050 1375 625 625 1600; 2003:2013];
darien=[3000 1000 1700 29750 26812 37034 49043 41399 41535; 2006:2014];
% puntarenas=[4000 886 32782 26225.6 11400 25765 28684 34926; 2007:2014];
puntarenas=[770 1275 3450 5000 23450 23779 92661 159368 114905 92472 ...
    73836 191770 209586 269160];
colon=[560 6766 7045 17690 9850 3715 6100 25258; 2006:2013];
% gracias=[1550 12800 6075 7872 22500 67802 83100 66305 64756 29033; 2005:2014];
gracias=[3930 12800 4463 6906 63395 109875 219354 214087 107219 55268; 2005:2014];
norte=[1673 1223 4891 NaN 5575 775 32750 37650 48038 9400 9117 3970; 2003:2014];

h4=figure;
set(h4,'Color','white','Visible','off')
plot(2003:2013,peten(1,:),'-b','LineWidth',3)
hold on
plot(2003:2013,deptflows_ts(1,3:13),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4,'Flows_Peten','png')

h4_1=figure;
set(h4_1,'Color','white','Visible','off')
plot(2006:2014,darien(1,:),'-b','LineWidth',3)
hold on
plot(2006:2014,deptflows_ts(2,6:14),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4_1,'Flows_Darien','png')

h4_2=figure;
set(h4_2,'Color','white','Visible','off')
plot(2001:2014,puntarenas(1,:),'-b','LineWidth',3)
hold on
plot(2001:2014,deptflows_ts(4,1:14),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4_2,'Flows_Punta','png')

h4_3=figure;
set(h4_3,'Color','white','Visible','off')
plot(2005:2014,gracias(1,:),'-b','LineWidth',3)
hold on
plot(2005:2014,deptflows_ts(5,5:14),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4_3,'Flows_Gracias','png')

h4_4=figure;
set(h4_4,'Color','white','Visible','off')
plot(2006:2013,colon(1,:),'-b','LineWidth',3)
hold on
plot(2006:2013,deptflows_ts(6,6:13),'-r','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('CCDB Data','Model Output','Location','NorthWest')
saveas(h4_4,'Flows_Colon','png')

h4_5=figure;
set(h4_5,'Color','white','Visible','off')
plot(2003:2014,deptflows_ts(7,3:14),'-r','LineWidth',3)
hold on
plot(2003:2014,norte(1,:),'-b','LineWidth',3)
ylabel('Trafficking Volume (kg)')
xlabel('Year')
legend('Model Output','CCDB Data','Location','NorthWest')
saveas(h4_5,'Flows_Norte','png')

% ccdb_dept=table(peten(1,:)',peten(2,:)',darien(1,:)',darien(2,:)',...
%     puntarenas(1,:)',puntarenas(2,:)',colon(1,:)',colon(2,:)',...
%     gracias(1,:)',gracias(2,:)',norte(1,:)',norte(2,:)','VariableNames',...
%     {'Peten' 'Year' 'Darien' 'Year' 'Puntarenas' 'Year' 'Colon' 'Year' ...
%     'Gracias' 'Year' 'Norte' 'Year'});




