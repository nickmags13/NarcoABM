%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   NarcoLogic ABM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic         % start run timer
cd C:\Users\nmagliocca\Documents\Matlab_code\NarcoLogic
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@ Procedures @@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
TSTART=1;
TMAX=180;   % 15 years at monthly time steps

rng default

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@ Environment @@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

% Load Central America shapefiles and rasters
[CAadm0,CAattr0]=shaperead('X:\CentralAmericaData\GADM\g2015_2014_0\CAadm0.shp',...
    'UseGeoCoords',true);  %polygons
% calat=cat(1,CAmap(:).Lat);
% calon=CAmap.Lon;
% cabox=CAmap.BoundingBox;
caadmid0=cat(1,CAattr0.ADM0_CODE);
maxlat=18.49656;
minlat=5.49908990000006;
maxlon=-77.163943085;
minlon=-92.231320038;

[CAadm1,CAattr1]=shaperead('X:\CentralAmericaData\GADM\g2015_2014_1\CAadm1.shp',...
    'UseGeoCoords',true);  %polygons
% calat=cat(1,CAmap(:).Lat);
% calon=CAmap.Lon;
% cabox=CAmap.BoundingBox;
caadmid1=cat(1,CAattr1.ADM1_CODE);

[CAcntr,CAcntrattr]=shaperead('X:\CentralAmericaData\CentralAmerica\Vector\CAcentroids.shp','UseGeoCoords',...
    true);
cntrlat=cat(1,CAcntr.Lat);
cntrlon=cat(1,CAcntr.Lon);
CApts=geopoint(CAcntr);

%%% Raster layers %%%
% Central America adminstrative boundaries level 2, objectid
[cagrid,Rcagrid]=geotiffread('X:\CentralAmericaData\CentralAmerica\cagd_adm2_500.tif');
cagrid(cagrid == 2147483647)=0;     %remove No Data value
ca_adm2=double(cagrid);

% Central America adminstrative boundaries level 0, objectid
[cagrid_cntry,Rcagrid_cntry]=geotiffread('X:\CentralAmericaData\CentralAmerica\cagd_adm0_500.tif');
cagrid_cntry(cagrid_cntry == 255)=0; %remove No Data value
ca_adm0=double(cagrid_cntry);
cellsize=Rcagrid.CellExtentInLatitude;

cntrycodes=unique(cagrid_cntry);    %subset landscape by country to place nodes
% Belize(23),Costa
% Rica(55),Panama(173),Guatemala(94),Honduras(101),Nicuragua(161),El
% Salvador(70)
cntrycodes=cntrycodes(cntrycodes ~= 0);
orderccodes=[173 55 161 94 70 101 23]; %order countries in desired order of nodes

% Tree cover
[tcov,Rtcov]=geotiffread('X:\CentralAmericaData\CentralAmerica\trcv_2000_500.tif');
treecov=tcov;
treecov(treecov==255 | cagrid_cntry == 0)=-1;
treecov=double(treecov);
itreecov=find(treecov > 0); %identify high forest cover areas
treecovpct=quantile(treecov(itreecov),[0.025 0.25 0.50 0.75 0.975]);
itreehigh=find(treecov >= 97);
avgtcov=zeros(length(cntrycodes),1);    %average tree cover per county, weighting for generating trade nodes
for cc=1:length(orderccodes)
    avgtcov(cc)=mean(treecov(ca_adm0 == orderccodes(cc)));
end

LANDSUIT=zeros(size(treecov));  % land suitability based on forest loss (and narco variable) predictors
LANDSUIT(itreecov)=treecov(itreecov)./max(treecov(itreecov));   % for now, just based on tree cover (high tree cover = high suitability)

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@ Agent Attributes @@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%%% Interdiction Agent %%%
slprob_0=0.001;     % baseline probability of seisure and loss event
delta_sl=0.05;      % reinforcement learning rate for S&L vents (i.e., weight on new information)

%%% Network Agent %%%
stock_0=100000;     %initial cocaine stock at producer node
startvalue=385; %producer price, $385/kg
deltavalue=8;   %added value for distance traveled $8/kilo/km
nodeloss=0;     % amount of cocaine that is normally lost (i.e., non-interdiction) at each node
ctrans_ground=3.5;  % transportation costs (kg/km) over-ground
ctrans_air=2;     % transportation costs (kg/km) via plane or boat
% Set-up producer and end supply nodes
strow=2700;
stcol=3600;
edrow=100;
edcol=100;
[startlat,startlon]=pix2latlon(Rcagrid,strow,stcol);
pstart=geopoint(startlat,startlon,'NodeName',{'Start Node'});
[endlat,endlon]=pix2latlon(Rcagrid,edrow,edcol);
pend=geopoint(endlat,endlon,'NodeName',{'End Node'});

%%% Nodes agents %%%
% perceived risk model
alpharisk=2;  %baseline
timewght_0=0.91;     %time discounting for subjective risk perception (Gallagher, 2014), range[0,1.05]
betarisk=alpharisk/slprob_0-alpharisk;

% cntrycpcty=[0.1 0.1 0.1 0.1 0.1 0.1 0.1];   %country-specific, per node trafficking capacity


%%%%%%%%%%%%%  Set-up Trafficking Network  %%%%%%%%%%%%%%%%%%
% G=digraph;  
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
for i=1:length(cntrycodes)
    icntry=find(ca_adm0 == orderccodes(i));
    ipotnode=find(ismember(icntry,itreehigh)==1);   %place nodes based on treecover
    randnode=icntry(ipotnode(randperm(length(ipotnode),...
        round(10*avgtcov./median(avgtcov)))));
    [nrow,ncol]=ind2sub(size(ca_adm0),randnode);
    [nlat,nlon]=pix2latlon(Rcagrid,nrow,ncol);
    nodeid=[nodeid length(nodeid)+(1:length(randnode))];
    noderow=[noderow; nrow];
    nodecol=[nodecol; ncol];
    nodelat=[nodelat; nlat];
    nodelon=[nodelon; nlon];
    nodecode=[nodecode; orderccodes(i)*ones(length(randnode),1)];
    nodestck=[nodestck; zeros(length(randnode),1)];
    nodecptl=[nodecptl; zeros(length(randnode),1)];
    if i == 1
        snode=ones(length(randnode),1);
        tnode=(1+(1:length(randnode)))';
        weights=ones(length(randnode),1);
        flows=ones(length(randnode),1);
        cpcty=0.1*stock_0*ones(length(randnode),1); %currently all the same capacity, but could introduce heterogeneity
        EdgeTable=table([snode tnode],weights,flows,cpcty,'VariableNames',...
            {'EndNodes' 'Weight' 'Flows' 'Capacity'});
    end
    if i == length(cntrycodes)
        inei=(nodecode == 23 | nodecode == 94);
        snode=nodeid(inei)';
%         snode=(length(nodeid)-(length(randnode)-1):length(nodeid))';
        tnode=(length(nodeid)+1)*ones(length(snode),1);
        nodeid=[nodeid max(nodeid)+1];  %add end node
        noderow=[noderow; edrow];
        nodecol=[nodecol; edcol];
        nodelat=[nodelat; pend.Latitude];
        nodelon=[nodelon; pend.Longitude];
        nodecode=[nodecode; 2];
        nodestck=[nodestck; 0];
        nodecptl=[nodecptl; 0];
        weights=ones(length(snode),1);
        flows=ones(length(snode),1);
        cpcty=0.1*stock_0*ones(length(snode),1);
        EdgeTable=table([EdgeTable.EndNodes; snode tnode],[EdgeTable.Weight; ...
            weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
            'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
    end 
end
NodeTable=table(nodeid',noderow,nodecol,nodelat,nodelon,nodecode,nodestck,nodecptl,...
    'VariableNames',{'ID','Row','Col','Lat','Lon','CountryCode','Stock','Capital'});
nnodes=height(NodeTable);
ADJ=zeros(nnodes);      % adjacency matrix for trafficking network
DIST=zeros(nnodes);     % geographic distance associated with edges
ADDVAL=zeros(nnodes);    % added value per edge in trafficking network
WGHT=zeros(nnodes);     % dynamic weighting of edges
FLOW=zeros(nnodes,nnodes,TMAX);     % dynamic flows of cocaine between nodes
SLRISK=zeros(nnodes);     % dynamic perceived risk of seisure and loss per edge by node agent
INTRISK=zeros(nnodes,TMAX); % dynamic perceived risk of interdiction at each node
CPCTY=zeros(nnodes);     % maximum flow possible between nodes
CTRANS=zeros(nnodes);   % transportation costs between nodes
rentcap=0.3*ones(nnodes,1);     % proportion of value of shipments 'captured' by nodes

slevent=zeros(nnodes,nnodes,TMAX);  % occurrence of S&L event
slsuccess=zeros(nnodes,nnodes,TMAX);    % S&L events in which cocaine was seized
SLPROB=zeros(nnodes,nnodes,TMAX);   % dynamic probability of S&L event per edge
intrdevent=zeros(nnodes,TMAX);
INTRDPROB=zeros(nnodes,TMAX);

STOCK=zeros(nnodes,TMAX);       %dynamic cocaine stock at each node
PRICE=zeros(nnodes,TMAX);       % $/kilo at each node
INFLOW=zeros(nnodes,TMAX);       % dynamic stock of cocaine coming into at each node
OUTFLOW=zeros(nnodes,TMAX);       % dynamic stock of cocaine leaving from each node
TOTCPTL=zeros(nnodes,TMAX);     % total value of cocaine at each node
ICPTL=zeros(nnodes,TMAX);       % dynamic illicit capital accumulated at each node
LCPTL=zeros(nnodes,TMAX);       % dynamic legitimate capital accumulated at each node
LEAK=zeros(nnodes,TMAX);        % dynamic amount of cocaine leaked at each node

for k=1:nnodes-1
    if k == 1
        newedges=2:nnodes;
        nodechk=ismember(newedges,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2)); %check for redundant edges
        newedges=newedges(~nodechk);
        weights=ones(length(newedges),1);
        flows=ones(length(newedges),1);
        cpcty=0.1*stock_0*ones(length(newedges),1);
        EdgeTable=table([EdgeTable.EndNodes; k*ones(length(newedges),1) newedges'],...
            [EdgeTable.Weight; weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
            'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
        
        % Create adjacency matrix (without graph toolbox)
        ADJ(k,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2))=1;
    else
        nodeset=k+1:nnodes;
        nnewedges=ceil(0.1*length(nodeset)*rand(1));    %generate new edges
        newedges=nodeset(randperm(length(nodeset),min(nnewedges,length(nodeset))));
        nodechk=ismember(newedges,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2)); %check for redundant edges
        newedges=newedges(~nodechk);
        weights=ones(length(newedges),1);
        flows=ones(length(newedges),1);
        cpcty=0.1*stock_0*ones(length(newedges),1);
        EdgeTable=table([EdgeTable.EndNodes; k*ones(length(newedges),1) newedges'],...
            [EdgeTable.Weight; weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
            'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
        
        % Create adjacency matrix (without graph toolbox)
        ADJ(k,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2))=1;
    end
end
% Make sure all nodes connect to end node
iendnode=NodeTable.ID(NodeTable.CountryCode == 2);
newedges=1:nnodes-1;
nodechk=ismember(newedges,EdgeTable.EndNodes(EdgeTable.EndNodes(:,2)==iendnode,1)); %check for redundant edges
newedges=newedges(~nodechk);
weights=ones(length(newedges),1);
flows=ones(length(newedges),1);
cpcty=0.05*stock_0*ones(length(newedges),1);
EdgeTable=table([EdgeTable.EndNodes; newedges' iendnode*ones(length(newedges),1)],...
    [EdgeTable.Weight; weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
    'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});

% Create adjacency matrix (without graph toolbox)
ADJ(EdgeTable.EndNodes(EdgeTable.EndNodes(:,2)==iendnode,1),iendnode)=1;

for j=1:nnodes
    %Create weight and capacity matrices
    WGHT(j,ADJ(j,:)==1)=EdgeTable.Weight(ADJ(j,:)==1);
    CPCTY(j,ADJ(j,:)==1)=EdgeTable.Capacity(ADJ(j,:)==1);
     
    % Create distance (in km) matrix
    latlon2=[NodeTable.Lat(ADJ(j,:)==1) NodeTable.Lon(ADJ(j,:)==1)];
    latlon1=repmat([NodeTable.Lat(j) NodeTable.Lon(j)],length(latlon2(:,1)),1);
    [d1km,d2km]=lldistkm(latlon1,latlon2);
    DIST(j,ADJ(j,:)==1)=d1km;
    
    % Create transport cost matrix (USD kg km-1)
    idist_ground=(DIST(j,:)  >0 & DIST(j,:) <= 500);
    idist_air=(DIST(j,:) > 500);
    CTRANS(j,idist_ground)=ctrans_ground.*DIST(j,idist_ground);
    CTRANS(j,idist_air)=ctrans_air.*DIST(j,idist_air);
    
    %Create added value matrix (USD) and price per node
    ADDVAL(j,ADJ(j,:)==1)=deltavalue.*DIST(j,ADJ(j,:)==1);
    if j == 1
        PRICE(j,TSTART)=startvalue;
%         PRICE(ADJ(j,:)==1,TSTART)=startvalue+ADDVAL(j,ADJ(j,:)==1);
    else
        isender=EdgeTable.EndNodes(EdgeTable.EndNodes(:,2) == j,1);
%         isender=find(ADJ(:,j) == 1);
        PRICE(j,TSTART)=mean(PRICE(isender,TSTART)+ADDVAL(isender,j));
    end
end

%%% Initialize Interdiction agent
SLPROB(:,:,TSTART)=(DIST./max(max(DIST)));   % dynamic probability of seisure and loss at edges
SLPROB(:,:,TSTART+1)=SLPROB(:,:,TSTART);
INTRDPROB(:,TSTART+1)=slprob_0*ones(nnodes,1); % dynamic probability of interdiction at nodes

%%% Initialize Node agents
STOCK(:,TSTART)=NodeTable.Stock(:);
TOTCPTL(:,TSTART)=NodeTable.Capital(:);
PRICE(:,TSTART+1)=PRICE(:,TSTART);

%%% Set-up node and network risk perceptions
SLRISK(:,:)=(DIST./max(max(DIST)));
INTRISK(:,TSTART:TSTART+1)=slprob_0.*ones(nnodes,2);

% subjective risk perception with time distortion
twght=timewght_0*ones(nnodes,1);    % time weighting for dynamic, subjective perceived risk of interdiction event
    
%%% Set-up trafficking netowrk benefit-cost logic  %%%%%%%%%%%%
ltcoeff=ones(nnodes,1);

%%% Define Node Investment Choice Sets

%%% Set-up figure for trafficking movie
MOV=zeros(nnodes,nnodes,TMAX);

%%
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@ Dynamics @@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

for t=TSTART+1:TMAX
    %%% S&L and interdiction events
    slevent(:,:,t)=(SLPROB(:,:,t) > rand(size(ADJ)));
    intrdevent(:,t)=(INTRDPROB(:,t) > rand(nnodes,1));
    MOV(:,1,t)=NodeTable.Stock(:);
    
    for n=1:nnodes-1 %exclude end node
       
      %%%%%  Route cocaine shipmments %%%%%
      STOCK(n,t)=STOCK(n,t-1)+STOCK(n,t);
      TOTCPTL(n,t)=TOTCPTL(n,t-1)+TOTCPTL(n,t);
%       ICPTL(n,t)=ICPTL(n,t-1)+ICPTL(n,t);
       if STOCK(n,t) > 0
         if n > 1
             LEAK(n,t)=nodeloss*STOCK(n,t); %drugs 'leaked' at each node
             STOCK(n,t)=STOCK(n,t)-LEAK(n,t);
         end
         inei=find(ADJ(n,:)==1);
         %%% Procedure for selecting routes based on expected profit %%%
         c_trans=CTRANS(n,inei);
         p_sl=SLRISK(n,inei);
         p_int=INTRISK(n,inei); %currently not used
         y_node=ADDVAL(n,inei);
         q_node=min(floor(STOCK(n,t)./length(inei)),CPCTY(n,inei));
         lccf=ltcoeff(n);
         totstock=STOCK(n,t);
         totcpcty=CPCTY(n,inei);
         
         [neipick,neivalue]=calc_neival(c_trans,p_sl,y_node,q_node,lccf,...
             totstock,totcpcty);
         inei=inei(neipick);
         
%          WGHT(n,inei)=SLRISK(n,inei).*ADDVAL(n,inei) ... something like
%          that
%          if length(inei) >  10
%              inei=randperm(length(inei),10);
%          end

         %%% !!! Put checks in to make sure buying node has enough capital
         FLOW(n,inei,t)=min(WGHT(n,inei).*floor(STOCK(n,t)/length(inei)),CPCTY(n,inei));
         % Check for S%L event
         if isempty(find(ismember(find(slevent(n,:,t)),inei),1)) == 0
             isl=(slevent(n,inei,t)==1);
             slsuccess(n,inei(isl),t)=FLOW(n,inei(isl),t);
             OUTFLOW(n,t)=sum(FLOW(n,inei,t));
             STOCK(n,t)=STOCK(n,t)-OUTFLOW(n,t);
             FLOW(n,inei(isl),t)=0;     % remove from trafficking route due to S&L event
             STOCK(inei,t)=STOCK(inei,t)+FLOW(n,inei,t)';
             TOTCPTL(n,t)=TOTCPTL(n,t)+sum(FLOW(n,inei,t).*ADDVAL(n,inei));
             TOTCPTL(inei,t)=TOTCPTL(inei,t-1)-(FLOW(n,inei,t)'.*PRICE(inei,t));
             ICPTL(n,t)=rentcap(n)*sum(FLOW(n,inei).*ADDVAL(n,inei));
             
         else
             OUTFLOW(n,t)=sum(FLOW(n,inei,t));
             STOCK(n,t)=STOCK(n,t)-OUTFLOW(n,t);
             STOCK(inei,t)=STOCK(inei,t)+FLOW(n,inei,t)';
             TOTCPTL(n,t)=TOTCPTL(n,t)+sum(FLOW(n,inei,t).*ADDVAL(n,inei));
             TOTCPTL(inei,t)=TOTCPTL(inei,t-1)-(FLOW(n,inei,t)'.*PRICE(inei,t));
             ICPTL(n,t)=rentcap(n)*sum(FLOW(n,inei).*ADDVAL(n,inei));
         end
         %!!!!!!!!!!!
%          ICPTL(n,t)=ICPTL(n,t)-OUTFLOW(n,t)*VALUE(  %account for value retained at node

         NodeTable.Stock(:)=STOCK(:,t);
         NodeTable.Capital(:)=TOTCPTL(:,t);
      end

      %%%% Update perceived risk in response to S&L and Interdiction events
      timeweight=twght(n);
      % identify neighbors in network (without network toolbox)
      bcknei=EdgeTable.EndNodes(EdgeTable.EndNodes(:,2) == n,1)';
      fwdnei=inei;
      sloccur=reshape(slevent(n,[bcknei fwdnei],TSTART+1:t),length(TSTART+1:t),...
          length([bcknei fwdnei]));
      intrdoccur=intrdevent([bcknei fwdnei],TSTART+1:t);
      [sl_risk,intrd_risk,slevnt,intrdevnt,tmevnt]=calc_intrisk(sloccur,...
          intrdoccur,t,TSTART,alpharisk,betarisk,timeweight);
      SLRISK(n,[bcknei fwdnei])=sl_risk;
      INTRISK(n,t+1)=mean(intrd_risk);  %node-specific risk is the average of neighbor risks
      
      %%% Updating interdiction event probability
      SLPROB(:,:,t+1)=max((1-delta_sl).*SLPROB(:,:,t)+delta_sl.*...
          (slsuccess(:,:,t)./max(max(slsuccess(:,:,t)))),SLPROB(:,:,TSTART));
      INTRDPROB(:,t+1)=INTRDPROB(:,t);
      
      
      %%% Make trafficking movie
      MOV(:,n,t)=STOCK(:,t);      % Capture stock data after each node iteration
    end
    STOCK(1,t+1)=stock_0;    %additional production to enter network next time step
    STOCK(nnodes,t+1)=0;    %remove stock at end node for next time step
    NodeTable.Stock(1)=stock_0;
    NodeTable.Stock(nnodes)=0;
end
toc     % stop run timer
% %%
% % %%% Visualization %%%
% 
% %%% Trafficking movie
% writerObj = VideoWriter('trafficking_risk_v2.mp4','MPEG-4');
% writerObj.FrameRate=2;
% open(writerObj);
% 
% h1=figure;
% set(h1,'Color','white','Visible','off')
% 
% for tt=TSTART+1:t
%     geoshow(CAadm0,'FaceColor',[1 1 1])
%     hold on
% %     plot(nodelon(1),nodelat(1),'r.','MarkerSize',ceil(MOV(1,1,tt)/1000))
%     plot(nodelon(1),nodelat(1),'r.','MarkerSize',ceil(stock_0/1000))
%     fedge=find(FLOW(1,:,tt) > 0);
%     for g=1:length(fedge)
%         plot([nodelon(1); nodelon(fedge(g))],[nodelat(1); nodelat(fedge(g))],'-k')
%     end
%     plot(nodelon(2:nnodes),nodelat(2:nnodes),'b.','MarkerSize',3)
%     xlabel('Longitude')
%     ylabel('Latitude')
%     title(sprintf('Timestep(month) = %d',tt-1))
%     frame = getframe(h1);
%     writeVideo(writerObj,frame);
%     % set(gca,'nextplot','replacechildren');
%     % set(gcf,'Renderer','zbuffer');
%     % ax=gca;
%     % movfilename='testmov.gif';
%     % cmap=get(h1,'ColorMap');
%     % ax.NextPlot='replaceChildren';
%     % MOV(nnodes-1) = struct('cdata',[],'colormap',[]);
%     for mm=1:nnodes-1
%         if mm == 1
%             clf
%             geoshow(CAadm0,'FaceColor',[1 1 1])
%             hold on
%             for nn=1:nnodes
%                 if MOV(nn,mm,tt) > 0
%                     plot(nodelon(nn),nodelat(nn),'r.','MarkerSize',ceil(MOV(nn,mm,tt)./1000))
%                 else
%                     plot(nodelon(nn),nodelat(nn),'b.','MarkerSize',3)
%                 end
%             end
%             continue
%         end
%         if isempty(find(MOV(mm,mm-1,tt) > 0,1)) == 1
%             continue
%         else
%             clf
%             geoshow(CAadm0,'FaceColor',[1 1 1])
%             hold on
%             for nn=1:nnodes
%                 if MOV(nn,mm,tt) > 0
%                     plot(nodelon(nn),nodelat(nn),'r.','MarkerSize',ceil(MOV(nn,mm,tt)./1000))
%                 else
%                     plot(nodelon(nn),nodelat(nn),'b.','MarkerSize',3)
%                 end
%             end
%             fedge=find(FLOW(mm,:,tt) > 0);
%             if isempty(find(fedge,1)) == 1
%                 plot(nodelon(mm),nodelat(mm),'kx','MarkerSize',ceil(MOV(mm,mm-1,tt)./1000))
%             else
%                 for g=1:length(fedge)
%                     plot([nodelon(mm); nodelon(fedge(g))],[nodelat(mm); nodelat(fedge(g))],'-k')
%                 end
%             end
%         end
%         xlabel('Longitude')
%         ylabel('Latitude')
%         title(sprintf('Timestep(month) = %d',tt-1))
%         frame = getframe(h1);
%         writeVideo(writerObj,frame);
%     end
%     clf
% end
% close(writerObj);
% 
% % im = frame2im(frame);
% % [A,cmap] = rgb2ind(im,256);
% % if n == 1;
% %     imwrite(A,cmap,movfilename,'gif','LoopCount',Inf,'DelayTime',0.5);
% % else
% %     imwrite(A,cmap,movfilename,'gif','WriteMode','append','DelayTime',0.5);
% % end
% % 
% % figure
% % axes('Position',ax.Position)
% % movie(MOV,1,1)
% % geoshow(CAmap)
% % hold on
% % geoshow(startlat,startlon,'DisplayType','point','MarkerEdgeColor','g','Marker', '+')
% % geoshow(endlat,endlon,'DisplayType','point','MarkerEdgeColor','r','Marker', '+')
% % geoshow(NodeTable.Lat,NodeTable.Lon,'DisplayType','point','MarkerEdgeColor','b','Marker', '.')
% % %quick trafficking network viz
% % figure
% % plot(nodelon,nodelat,'k.','MarkerSize',1)
% % hold on
% % for i=1:length(nodelat)
% %     %     fedge=find(EdgeTable.EndNodes(:,1)==i);
% %     fedge=find(ADJ(i,:)==1);
% %     for g=1:length(fedge)
% %         plot([nodelon(i); nodelon(fedge(g))],[nodelat(i); nodelat(fedge(g))],'-k')
% %     end
% %     plot(nodelon(i),nodelat(i),'r.','MarkerSize',ceil(PRICE(i,1)./1000))
% % end
% % % plot(nodelon,nodelat,'r.','MarkerSize',ceil(PRICE(:,1)./1000))
% % % % Trafficking volume per node
% % figure
% % plot(nodelon,nodelat,'k.','MarkerSize',1)
% % hold on
% % for i=1:length(nodelat)
% %     if STOCK(i,t) > 0
% %         plot(nodelon(i),nodelat(i),'r.','MarkerSize',ceil(STOCK(i,t)./1000))
% %     else
% %         plot(nodelon(i),nodelat(i),'k.','MarkerSize',1)
% %     end
% % end
% %%% Command to use
% % digraph, maxflow, nearest