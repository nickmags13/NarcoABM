%%%%%%%%%%%%%%%%%    Build trafficking network    %%%%%%%%%%%%%%%%%%%%%%%%%
function [NodeTable,EdgeTable]=build_network(ca_adm0,Rcagrid_cntry,dptgrid,...
    Rdptgrid,LANDSUIT,dptcodes,dptorder,savedState,stock_0)

% Set-up producer and end supply nodes
strow=size(ca_adm0,1);
stcol=size(ca_adm0,2);
edrow=100;
edcol=100;
[startlat,startlon]=pix2latlon(Rcagrid_cntry,strow,stcol);
pstart=geopoint(startlat,startlon,'NodeName',{'Start Node'});
[endlat,endlon]=pix2latlon(Rcagrid_cntry,edrow,edcol);
pend=geopoint(endlat,endlon,'NodeName',{'End Node'});

nodeid=1;
noderow=strow;
nodecol=stcol;
nodelat=[pstart.Latitude];
nodelon=[pstart.Longitude];
nodecode=1;
nodestck=stock_0;
% nodestck=0;
% nodecptl=startvalue*stock_0;
nodecptl=0;
% nodename={'startnode'};
nodetcov=0;     % node tree cover
nodepopsuit=0;
nodedcsuit=0;
nodedbsuit=0;
nodeslpsuit=0;
nodemktsuit=0;
nodelusuit=0;
nodelsuit=0;
nodedto=0;
% Create empty table
snode=ones([],1);
tnode=[];
weights=ones([],1);
flows=ones([],1);
%                 cpcty=2*stock_0*ones(length(randnode),1); %currently all the same capacity, but could introduce heterogeneity
cpcty=[];
EdgeTable=table([snode tnode],weights,flows,cpcty,'VariableNames',...
    {'EndNodes' 'Weight' 'Flows' 'Capacity'});

dtoassign=[1; 2; 1; 2];

% Allocate nodes based on suitability
nodepct=0.00005; %percentage of high suitability cells that contain possible nodes
nodequant=quantile(LANDSUIT(~isnan(LANDSUIT)),[0.025 0.50 0.66 0.75 0.99]);
%         inodepick=find(LANDSUIT > nodequant(4));
inodepick=find(LANDSUIT > 0.8);
% inodepick=find(LANDSUIT > 0.5);     %Use with RAT suitability

avgnodealloc=ceil((length(inodepick)*nodepct)/length(dptcodes));
pctdptsuit=zeros(length(dptorder(:,1)),1);
for dc=1:length(pctdptsuit)
    dptsuit=LANDSUIT(dptgrid == dptorder(dc,1));
    pctdptsuit(dc)=length(find(dptsuit > 0.8))/length(dptsuit);
%     pctdptsuit(dc)=length(find(dptsuit > 0.5))/length(dptsuit); %Use with RAT suitability
end
allocnodes=round(1.75*pctdptsuit/mean(pctdptsuit));
%         allocnodes=round(pctdptsuit/mean(pctdptsuit));
for i=1:length(dptorder)
    if pctdptsuit(i)==0
        continue
    else
        %place nodes based on LANDSUIT
        idptmnt=find(dptgrid == dptorder(i,1));
        ipotnode=find(ismember(idptmnt,inodepick)==1);
        % allocate nodes per department based on suitability within department
        randnode=idptmnt(ipotnode(randperm(length(ipotnode),...
            allocnodes(i))));
        
        %                 randnode=idptmnt(ipotnode(randperm(length(ipotnode),...
        %                     min(4,length(ipotnode)))));
        
        [nrow,ncol]=ind2sub(size(dptgrid),randnode);
        [nlat,nlon]=pix2latlon(Rdptgrid,nrow,ncol);
        nodeid=[nodeid length(nodeid)+(1:length(randnode))];
        noderow=[noderow; nrow];
        nodecol=[nodecol; ncol];
        nodelat=[nodelat; nlat];
        nodelon=[nodelon; nlon];
        nodecode=[nodecode; dptorder(i,1)*ones(length(randnode),1)];
        nodestck=[nodestck; zeros(length(randnode),1)];
        nodecptl=[nodecptl; zeros(length(randnode),1)];
        nodetcov=[nodetcov; zeros(length(randnode),1)];
        nodepopsuit=[nodepopsuit; zeros(length(randnode),1)];
        nodedcsuit=[nodedcsuit; zeros(length(randnode),1)];
        nodedbsuit=[nodedbsuit; zeros(length(randnode),1)];
        nodeslpsuit=[nodeslpsuit; zeros(length(randnode),1)];
        nodemktsuit=[nodemktsuit; zeros(length(randnode),1)];
        nodelusuit=[nodelusuit; zeros(length(randnode),1)];
        nodelsuit=[nodelsuit; LANDSUIT(randnode)];
        nodedto=[nodedto; zeros(length(randnode),1)];
        %                 nodedto=[nodedto; dtoassign(1:min(length(dtoassign),length(ipotnode)))];
        if i == 1
            snode=ones(length(randnode),1);
            tnode=(1+(1:length(randnode)))';
            weights=ones(length(randnode),1);
            flows=ones(length(randnode),1);
            %                 cpcty=2*stock_0*ones(length(randnode),1); %currently all the same capacity, but could introduce heterogeneity
            cpcty=1000000*ones(length(randnode),1);
            %                     EdgeTable=table([snode tnode],weights,flows,cpcty,'VariableNames',...
            %                         {'EndNodes' 'Weight' 'Flows' 'Capacity'});
            EdgeTable=table([EdgeTable.EndNodes; snode tnode],[EdgeTable.Weight; ...
                weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
                'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
        end
        if i == length(dptcodes)
            ineicode=ismember(nodecode,unique(dptgrid(ca_adm0 == 23 | ca_adm0 ==94)));
            snode=nodeid(ineicode)';
            tnode=(length(nodeid)+1)*ones(length(snode),1);
            nodeid=[nodeid max(nodeid)+1];  %add end node
            noderow=[noderow; edrow];
            nodecol=[nodecol; edcol];
            nodelat=[nodelat; pend.Latitude];
            nodelon=[nodelon; pend.Longitude];
            nodecode=[nodecode; 2];
            nodestck=[nodestck; 0];
            nodecptl=[nodecptl; 0];
            nodetcov=[nodetcov; 0];
            nodepopsuit=[nodepopsuit; 0];
            nodedcsuit=[nodedcsuit; 0];
            nodedbsuit=[nodedbsuit; 0];
            nodeslpsuit=[nodeslpsuit; 0];
            nodemktsuit=[nodemktsuit; 0];
            nodelusuit=[nodelusuit; 0];
            nodelsuit=[nodelsuit; 0];
            nodedto=[nodedto; 0];
            weights=ones(length(snode),1);
            flows=ones(length(snode),1);
            %                 cpcty=2*stock_0*ones(length(snode),1);
            cpcty=1000000*ones(length(snode),1);
            EdgeTable=table([EdgeTable.EndNodes; snode tnode],[EdgeTable.Weight; ...
                weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
                'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
        end
    end
end
NodeTable=table(nodeid',noderow,nodecol,nodelat,nodelon,nodecode,nodestck,...
    nodecptl,nodetcov,nodepopsuit,nodedcsuit,nodedbsuit,nodeslpsuit,...
    nodemktsuit,nodelusuit,nodelsuit,nodedto,'VariableNames',{'ID','Row','Col','Lat',...
    'Lon','DeptCode','Stock','Capital','TreeCover','PopSuit',...
    'DistCoastSuit','DistBorderSuit','SlopeSuit','MktAccSuit','LandUseSuit',...
    'LandSuit','DTO'});

rng(savedState);
hitrngstate=rand(height(NodeTable),1);

for k=1:height(NodeTable)-1
    if k == 1
%         newedges=2:height(NodeTable);
        newedges=2:height(NodeTable)-1; %eliminate direct producer to consumer edge
        nodechk=ismember(newedges,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2)); %check for redundant edges
        newedges=newedges(~nodechk);
        weights=ones(length(newedges),1);
        flows=ones(length(newedges),1);
        %                 cpcty=2*stock_0*ones(length(newedges),1);
        cpcty=1000000*ones(length(newedges),1);
        EdgeTable=table([EdgeTable.EndNodes; k*ones(length(newedges),1) newedges'],...
            [EdgeTable.Weight; weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
            'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
    else
        nodeset=k+1:height(NodeTable);
        nnewedges=ceil(0.1*length(nodeset)*rand(1));    %generate new edges
        newedges=nodeset(randperm(length(nodeset),min(nnewedges,length(nodeset))));
        nodechk=ismember(newedges,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2)); %check for redundant edges
        newedges=newedges(~nodechk);
        weights=ones(length(newedges),1);
        flows=ones(length(newedges),1);
        %                 cpcty=2*stock_0*ones(length(newedges),1);
        cpcty=1000000*ones(length(newedges),1);
        EdgeTable=table([EdgeTable.EndNodes; k*ones(length(newedges),1) newedges'],...
            [EdgeTable.Weight; weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
            'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
    end
end

% Make sure all nodes connect to end node
iendnode=NodeTable.ID(NodeTable.DeptCode == 2);
% newedges=1:height(NodeTable)-1;
newedges=2:height(NodeTable)-1; %eliminate direct producer to consumer edge
nodechk=ismember(newedges,EdgeTable.EndNodes(EdgeTable.EndNodes(:,2)==iendnode,1)); %check for redundant edges
newedges=newedges(~nodechk);
weights=ones(length(newedges),1);
flows=ones(length(newedges),1);
%         cpcty=2*stock_0*ones(length(newedges),1);
cpcty=1000000*ones(length(newedges),1);
EdgeTable=table([EdgeTable.EndNodes; newedges' iendnode*ones(length(newedges),1)],...
    [EdgeTable.Weight; weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
    'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});

% Calculate interception probability based on number of edges
p_intcpt=zeros(height(NodeTable),1);
for q=2:height(NodeTable)
%     n_out=length(find(EdgeTable.EndNodes(:,1) == q));
    n_in=length(find(EdgeTable.EndNodes(:,2) == q));
    p_intcpt(q)=1/n_in;
end
NodeTable.pintcpt=p_intcpt;