%%%%%%%%%% Build empirical S&L layers %%%%%%%%%%%%%%%%%%%%
function [empSLPROB,slctnodes] = build_SLemp(nnodes,TMAX,CAattr1,NodeTable,ADJ,ccdb)

% [CAadm0,CAattr0]=shaperead('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmericaData\GADM\g2015_2014_0\CAadm0.shp',...
%             'UseGeoCoords',true);
% [CAadm1,CAattr1]=shaperead('C:\Users\nrmagliocca\Box Sync\Data Drive\CentralAmericaData\GADM\g2015_2014_1\CAadm1.shp',...
%             'UseGeoCoords',true);
empSLPROB=zeros(nnodes,nnodes,TMAX);
        
cntrynames={'Guatemala','Guatemala','Guatemala','Panama','Honduras',...
    'Honduras','Honduras','Nicaragua','Nicaragua'};
deptnames={'Pet','Izabal','Alta Verapaz','Ember','Grac','Col','Olancho',...
    'Atlantico Norte','Atlantico Sur'};
cntryind=zeros(length(CAattr1),length(deptnames));
deptind=zeros(length(CAattr1),length(deptnames));
hind=zeros(1,length(deptnames));
slctnodes=cell(1,length(deptnames));

for jj=1:length(cntrynames)
    for ii=1:length(CAattr1)
    cntryind(ii,jj)=strncmp(cntrynames(jj),CAattr1(ii).ADM0_NAME,length(cntrynames{jj}));
    deptind(ii,jj)=strncmp(deptnames(jj),CAattr1(ii).ADM1_NAME,length(deptnames{jj}));
    end
    hind(jj)=find(cntryind(:,jj) == 1 & deptind(:,jj) == 1);
    nodeids=NodeTable.ID(NodeTable.DeptCode == CAattr1(hind(jj)).ADM1_CODE);
    slctnodes(jj)=mat2cell(nodeids,length(nodeids),1);
end

for T=1:length(ccdb(1,:))   %starting at 2001, rather than 2000 (TSTART)
    isl=(ccdb(:,T) > 0);
    trgtnodes=cat(1,slctnodes{isl});
    empSLPROB(trgtnodes,:,(T-1)*12+1)=ADJ(trgtnodes,:);
    empSLPROB(:,trgtnodes,(T-1)*12+1)=ADJ(:,trgtnodes);
    empSLPROB(:,:,(T-1)*12+2:(T-1)*12+12)=ones(1,1,11).*...
        empSLPROB(:,:,(T-1)*12+1);
end