%%%%%%%%%%%%%% Make network map %%%%%%%%%%%%%%%%%%%%%%
cd C:\Users\nrmagliocca\'Box Sync'\'Data Drive'\model_results\SupplyChain_full_021618
load supplychain_results_021618_1_6.mat

[CAadm0,CAattr0]=shaperead('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmericaData\GADM\g2015_2014_0\CAadm0.shp',...
            'UseGeoCoords',true);
[CAadm1,CAattr1]=shaperead('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmericaData\GADM\g2015_2014_1\CAadm1.shp',...
            'UseGeoCoords',true);
        
TMAX=180;
NNODES=size(activeroute,1);
nnodes=NNODES;
% %%% Locate Department statistics
cntrynames={'Guatemala','Panama','Panama','Costa Rica','Honduras','Panama','Nicaragua','Panama','Nicaragua','Guatemala','Honduras'};
    deptnames={'Pet','Dari','Ember','Puntarenas','Grac','Col','Atlantico Norte','Veraguas','Atlantico Sur','Izabal','Col'};
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
        
%%%% Single run from results file

h3=figure;
set(h3,'color','white')
geoshow(CAadm0,'FaceColor',[1 1 1])
hold on
geoshow(CAadm1(hind(1)),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
geoshow(CAadm1(hind(2)),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
geoshow(CAadm1(hind(3)),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
geoshow(CAadm0(6),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
geoshow(CAadm1(hind(5)),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
geoshow(CAadm1(hind(6)),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
geoshow(CAadm1(hind(7)),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
geoshow(CAadm1(hind(9)),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
geoshow(CAadm1(hind(11)),'FaceColor','blue','FaceAlpha',0.3,'LineWidth',2)
plot(NodeTable.Lon,NodeTable.Lat,'o','MarkerSize',8,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.5 0.5 0.5])
tstep=18;
iflow=find(OUTFLOW(:,tstep)>0);
for n=1:length(iflow)
    plot(NodeTable.Lon(iflow(n)),NodeTable.Lat(iflow(n)),'o','MarkerSize',8,...
        'MarkerEdgeColor','k','MarkerFaceColor','red')
    rtes=activeroute{iflow(n),tstep};
    for i=1:length(rtes)
        plot([NodeTable.Lon(iflow(n)) NodeTable.Lon(rtes(i))],...
            [NodeTable.Lat(iflow(n)) NodeTable.Lat(rtes(i))],':k')
    end
end

title('Average First Time Step of Cocaine Movements')
xlabel('Longitude')
ylabel('Latitude')