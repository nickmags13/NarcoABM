%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   NarcoLogic ABM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic         % start run timer
cd C:\Users\nmagliocca\Documents\Matlab_code\NarcoLogic
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@ Procedures @@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
TSTART=1;
% TMAX=180;   % 15 years at monthly time steps
TMAX=120;
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
% simplify geometry
latin=extractfield(CAadm0,'Lat')';
lonin=extractfield(CAadm0,'Lon')';
% [CAadm0_latrdc,CAadm0_lonrdc]=reducem(latin,lonin);
    
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

% Spatial narco vars by administrative departments
[dptvars,dptvarsattr]=shaperead('X:\CentralAmericaData\GADM\CA_ALLt_UTM\CA_ALLt_narcovars.shp',...
    'UseGeoCoords',true);
dptcode=cat(1,dptvarsattr.ADM1_CODE);
intlbrdrdmmy=cat(1,dptvarsattr.MAX_1);
coastdmmy=cat(1,dptvarsattr.COASTDMMY);
% meanlat=cat(1,dptvarsattr.MEAN_1);

%%% Palm oil mills
% startRow=1;
% endRow=Inf;
% filename='X:\CentralAmericaData\Model_inputs\pomills_latlon.txt';
% [pomilllat,pomilllon] = import_pomills(filename, startRow, endRow);
% pomilllat=pomilllat(2:length(pomilllat));
% pomilllon=pomilllon(2:length(pomilllon));
load pomilldist.mat

%%% Raster layers %%%
% Central America adminstrative boundaries level 2, objectid
[dptgrid,Rdptgrid]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\dptgrid_clp.tif');
dptnodataval=-9999;
dptgrid=double(dptgrid);
% dptgrid(dptgrid == dptnodataval)=NaN;     %remove No Data value
dptcodes=unique(dptgrid);
dptcodes=dptcodes(dptcodes ~= dptnodataval);

% Central America adminstrative boundaries level 0, objectid
[cagrid_cntry,Rcagrid_cntry]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\ca_cntry_clp.tif');
cntrynodataval=255;
ca_adm0=double(cagrid_cntry);
ca_adm0(cagrid_cntry == cntrynodataval)=NaN; %remove No Data value

cellsize=Rcagrid_cntry.CellExtentInLatitude;
%%%% reconciled to 'ca_slope_250.tif'

cntrycodes=unique(cagrid_cntry);    %subset landscape by country to place nodes
% Belize(23),Costa
% Rica(55),Panama(173),Guatemala(94),Honduras(101),Nicuragua(161),El
% Salvador(70)
cntrycodes=cntrycodes(cntrycodes ~= 0 & cntrycodes ~= cntrynodataval);
dbrdr_suit=zeros(size(dptgrid));
for iadm=1:length(dptcodes)
    admind=(dptgrid == dptcodes(iadm));
    dbrdr_suit(admind)=intlbrdrdmmy(iadm);
end
%order countries in desired order of nodes based on macroeconomics
% ordernodes=[173 55 161 94 70 101 23]; 
% ordernodes=sortrows([dptcode meanlat],2);
% ordernodes(ismember(ordernodes,sortpanlon(:,1)),1)=sortpanlon(:,1);
% meanlat=zeros(length(dptcodes),1);
% meanlon=zeros(length(dptcodes),1);
maxlat=zeros(length(dptcodes),1);
maxlon=zeros(length(dptcodes),1);
for j=1:length(dptcodes)
    [nrow,ncol]=ind2sub(size(dptgrid),find(dptgrid == dptcodes(j)));
    [nlat,nlon]=pix2latlon(Rdptgrid,nrow,ncol);
%     meanlat(j)=mean(nlat);
%     meanlon(j)=mean(nlon);
    maxlat(j)=max(nlat);
    maxlon(j)=max(nlon);
end
% nodevec=sqrt(meanlat.^2+meanlon.^2);
nodevec=sqrt(0.9*maxlat.^2+0.1*maxlon.^2);
nodemat=[dptcodes nodevec];
nodeorder=sortrows(nodemat,2);

% Tree cover
[tcov,Rtcov]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\treecov_clp.tif');
treecov=tcov;
tcovnodataval=255;
treecov=double(treecov);
treecov(cagrid_cntry==cntrynodataval)=NaN;
itreecov=find(treecov > 0); %identify high forest cover areas
treecovpct=quantile(treecov(itreecov),[0.025 0.25 0.50 0.75 0.975]);
itreepick=find(treecov >= treecovpct(2) & treecov < treecovpct(3));
avgtcov=zeros(length(cntrycodes),1);    %average tree cover per county, weighting for generating trade nodes
for cc=1:length(nodeorder(:,1))
    avgtcov(cc)=mean(treecov(ca_adm0 == nodeorder(cc,1)));
end

% Distance to coast and country borders
[dcoast,Rdcoast]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\dcoast_clp.tif');
dcoastnodataval=-9999;
dcoast(cagrid_cntry==cntrynodataval)=NaN;
dcoast(dcoast == dcoastnodataval)=NaN;
dcoast_suit=1-dcoast./max(max(dcoast));

% Population density as a proxy for remoteness
[popden,Rpopden]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\popden_clp.tif');
popnodataval=-1;
popden=double(popden);
popden(cagrid_cntry==cntrynodataval)=NaN;
popden(popden == popnodataval)=NaN;
popq=quantile(reshape(popden,size(popden,1)*size(popden,2),1),...
    [0.25 0.5 0.75 0.95]);
pop_suit=zeros(size(popden));
pop_suit(popden > popq(4))=0;
pop_suit(popden <= popq(4))=1-popden(popden <= popq(4))./popq(4);

% Topography
[slope,Rslope]=geotiffread('X:\CentralAmericaData\Model_inputs\ca_slope_250m.tif');
slope(cagrid_cntry==cntrynodataval)=NaN;
slopeclass=[8 16 30 31; 0 25 50 100]';   % GAEZ (see Magliocca et al., 2013, PLOS ONE)
slp_suit=zeros(size(slope));
slp_suit(slope < slopeclass(1,1))=1-slopeclass(1,2)/100;
slp_suit(slope >= slopeclass(1,1) & slope < slopeclass(2,1))=1-slopeclass(2,2)/100;
slp_suit(slope >= slopeclass(2,1) & slope < slopeclass(3,1))=1-slopeclass(3,2)/100;
slp_suit(slope >= slopeclass(4,1))=1-slopeclass(4,2)/100;

% Market Access
[mktacc,Rmktacc]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\mktacc_clp.tif');
manodataval=-9999;
mktacc=double(mktacc);
mktacc(cagrid_cntry==cntrynodataval)=NaN;
mktacc(mktacc == manodataval)=NaN;
mktacc_suit=zeros(size(mktacc));
% submasuit=mktacc./median(mktacc(~isnan(mktacc)));

% Maize Yield
[mazyld,Rmazyld]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\mazyld_clp.tif');
maznodataval=-9999;
mazyld=double(mazyld);
mazyld(cagrid_cntry==cntrynodataval)=NaN;
mazyld(mazyld == maznodataval)=NaN;

% Oil Palm Yield
[plmyld,Rplmyld]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\plmyld_clp.tif');
plmnodataval=-9999;
plmyld=double(plmyld);
plmyld(cagrid_cntry==cntrynodataval)=NaN;
plmyld(plmyld == plmnodataval)=NaN;

% Cattle density
[ctlden,Rctlden]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\ctlden_clp.tif');
plmnodataval=-9999;
ctlden=double(ctlden);
ctlden(cagrid_cntry==cntrynodataval)=NaN;
ctlden(ctlden == plmnodataval)=NaN;

% Protected Areas
[protarea,Rprotarea]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\protarea_clp.tif');
protnodataval=255;
protarea=double(protarea);
protarea(cagrid_cntry==cntrynodataval)=NaN;
protarea(protarea == protnodataval)=NaN;
protsuit=1-isnan(protarea);

% Land Prices
[landprice,Rlandprice]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\landprice_clp.tif');
landprice=double(landprice);
ldnprcnodataval=min(min(landprice));
landprice(cagrid_cntry==cntrynodataval)=NaN;
landprice(landprice < 0)=NaN;
landprice=landprice/28.36; %conversion to USD

% Initial land use
% Land cover classes: 
% 1. built-up
% 2. cropland — row crop agriculture
% 3. shrubs
% 4. trees
% 5. pastureland/grassland — grazing land or natural grassland
% 6. bare
% 7. plantation — citrus, vineyard, coffee, etc.
% 8. water
% 9. plantation tree — eucalyptus, pine, etc.
[luint,Rluint]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\luint_clp.tif');
lunodataval=0;
luint=double(luint);
luint(cagrid_cntry==cntrynodataval)=NaN;
luint(luint == lunodataval)=NaN;
lu_suit=zeros(size(luint));
lu_suit(luint == 1 | luint == 6 | luint == 7 | luint == 8 | luint ==9)=0;
lu_suit(luint == 2)=0.5;
lu_suit(luint == 3 | luint == 4 | luint == 5)=1;

%%% Land-base investment value
invst_suit=zeros(size(luint));
invst_suit(luint == 1 | luint == 6 | luint == 7 | luint == 8 | luint ==9)=0;
invst_suit(luint == 2)=0.5;
invst_suit(luint == 3 | luint == 5)=0.75;
invst_suit(luint == 4)=1;

%%% Weight each landscape attribute
tcwght=0;       % tree cover
brdwght=0;      % distance to country border
dcstwght=0;     % distance to coast
mktwght=1;      % market access
popwght=1;      % population density - proxy for remoteness
slpwght=1;      % slope-constrained land suitability
luwght=1;       % suitability based on initial land use
invstwght=1;    % investment potential of initial land use
protwght=1;         % protected area

wghts=[tcwght brdwght dcstwght mktwght popwght slpwght luwght invstwght protwght]./...
    sum([tcwght brdwght dcstwght mktwght popwght slpwght luwght invstwght protwght]);

%%% Null Model
LANDSUIT=wghts(1).*treecov./100+wghts(2).*dbrdr_suit+wghts(3).*dcoast_suit+...
    wghts(4).*mktacc_suit+wghts(5).*pop_suit+wghts(6).*slp_suit+wghts(6).*...
    lu_suit+wghts(7).*invst_suit+wghts(8)*(1-protsuit);  % land suitability based on biophysical and narco variable predictors

% %%% Narco model
% LANDSUIT=wghts(1).*treecov./100+wghts(2).*dbrdr_suit+wghts(3).*dcoast_suit+...
%     wghts(4).*(1-mktacc_suit)+wghts(5).*pop_suit+wghts(6).*slp_suit+wghts(6).*...
%     lu_suit+wghts(7).*invst_suit+wghts(8)*(1-protsuit);  % land suitability based on biophysical and narco variable predictors

%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@ Agent Attributes @@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%%% Interdiction Agent %%%
slprob_0=0.02;     % baseline probability of seisure and loss event
slcpcty_0=30;         % assumed number of S&L events that can be carried out per time step
delta_sl=0.5;      % reinforcement learning rate for S&L vents (i.e., weight on new information)
losstol=0.9;        % tolerance threshold for loss due to S&L, triggers route fragmentation
%%% Network Agent %%%
stock_0=100;     %initial cocaine stock at producer node
startvalue=385; %producer price, $385/kg
deltavalue=8;   %added value for distance traveled $8/kilo/km
nodeloss=0;     % amount of cocaine that is normally lost (i.e., non-interdiction) at each node
ctrans_inland=3.5;  % transportation costs (kg/km) over-ground
ctrans_coast=2;     % transportation costs (kg/km) via plane or boat
delta_rt=0.5;       % reinforcement learning rate for network agent 
                    %(i.e., weight on new information for successful routes)
                    % note: faster learning rate than for interdiction
                    % agent

% Set-up producer and end supply nodes
strow=size(ca_adm0,1);
stcol=size(ca_adm0,2);
edrow=100;
edcol=100;
[startlat,startlon]=pix2latlon(Rcagrid_cntry,strow,stcol);
pstart=geopoint(startlat,startlon,'NodeName',{'Start Node'});
[endlat,endlon]=pix2latlon(Rcagrid_cntry,edrow,edcol);
pend=geopoint(endlat,endlon,'NodeName',{'End Node'});

%%% Nodes agents %%%
% perceived risk model
alpharisk=0.001;  %baseline
timewght_0=0.91;     %time discounting for subjective risk perception (Gallagher, 2014), range[0,1.05]
% betarisk=alpharisk/slprob_0-alpharisk;
betarisk=0.5;

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
nodetcov=0;     % node tree cover
nodepopsuit=0;
nodedcsuit=0;
nodedbsuit=0;
nodeslpsuit=0;
nodemktsuit=0;
nodelusuit=0;
nodelsuit=0;
noderoad=0;
nodequant=quantile(LANDSUIT(~isnan(LANDSUIT)),[0.025 0.50 0.66 0.75 0.99]);
inodepick=find(LANDSUIT > nodequant(4));
% for i=1:length(cntrycodes)
for i=1:length(dptcodes)
%     icntry=find(ca_adm0 == nodeorder(i));
%     ipotnode=find(ismember(icntry,inodepick)==1);   %place nodes based on LANDSUIT
    idptmnt=find(dptgrid == nodeorder(i,1));
    ipotnode=find(ismember(idptmnt,inodepick)==1);
    % Select number of nodes per admin boundary based on drug intensity
    % index
    allocnodes=1;
%     randnode=icntry(ipotnode(randperm(length(ipotnode),...
%         round(10*avgtcov(i)./median(avgtcov)))));
    if isempty(find(ipotnode,1)) == 1
        subinodepick=find(LANDSUIT > nodequant(2));
        ipotnode=find(ismember(idptmnt,subinodepick)==1);
    end
    randnode=idptmnt(ipotnode(randperm(max(length(ipotnode),1),...
        allocnodes)));
%     [nrow,ncol]=ind2sub(size(ca_adm0),randnode);
%     [nlat,nlon]=pix2latlon(Rcagrid_cntry,nrow,ncol);
    [nrow,ncol]=ind2sub(size(dptgrid),randnode);
    [nlat,nlon]=pix2latlon(Rdptgrid,nrow,ncol);
    nodeid=[nodeid length(nodeid)+(1:length(randnode))];
    noderow=[noderow; nrow];
    nodecol=[nodecol; ncol];
    nodelat=[nodelat; nlat];
    nodelon=[nodelon; nlon];
    nodecode=[nodecode; nodeorder(i,1)*ones(length(randnode),1)];
    nodestck=[nodestck; zeros(length(randnode),1)];
    nodecptl=[nodecptl; zeros(length(randnode),1)];
    nodetcov=[nodetcov; treecov(randnode)];
    nodepopsuit=[nodepopsuit; pop_suit(randnode)];
    nodedcsuit=[nodedcsuit; dcoast_suit(randnode)];
    nodedbsuit=[nodedbsuit; dbrdr_suit(randnode)];
    nodeslpsuit=[nodeslpsuit; slp_suit(randnode)];
    nodemktsuit=[nodemktsuit; mktacc_suit(randnode)];
    nodelusuit=[nodelusuit; lu_suit(randnode)];
    nodelsuit=[nodelsuit; LANDSUIT(randnode)];
    noderoad=[noderoad; 0];
    if i == 1
        snode=ones(length(randnode),1);
        tnode=(1+(1:length(randnode)))';
        weights=ones(length(randnode),1);
        flows=ones(length(randnode),1);
        cpcty=2*stock_0*ones(length(randnode),1); %currently all the same capacity, but could introduce heterogeneity
        EdgeTable=table([snode tnode],weights,flows,cpcty,'VariableNames',...
            {'EndNodes' 'Weight' 'Flows' 'Capacity'});
    end
    if i == length(dptcodes)
%     if i == length(cntrycodes)
%         inei=(nodecode == 23 | nodecode == 94);
        inei=ismember(nodecode,unique(dptgrid(ca_adm0 == 23 | ca_adm0 ==94)));
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
        nodetcov=[nodetcov; 0];
        nodepopsuit=[nodepopsuit; 0];
        nodedcsuit=[nodedcsuit; 0];
        nodedbsuit=[nodedbsuit; 0];
        nodeslpsuit=[nodeslpsuit; 0];
        nodemktsuit=[nodemktsuit; 0];
        nodelusuit=[nodelusuit; 0];
        nodelsuit=[nodelsuit; 0];
        noderoad=[noderoad; 0];
        weights=ones(length(snode),1);
        flows=ones(length(snode),1);
        cpcty=2*stock_0*ones(length(snode),1);
        EdgeTable=table([EdgeTable.EndNodes; snode tnode],[EdgeTable.Weight; ...
            weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
            'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
    end 
end
NodeTable=table(nodeid',noderow,nodecol,nodelat,nodelon,nodecode,nodestck,...
    nodecptl,nodetcov,nodepopsuit,nodedcsuit,nodedbsuit,nodeslpsuit,...
    nodemktsuit,nodelusuit,nodelsuit,noderoad,'VariableNames',{'ID','Row','Col','Lat',...
    'Lon','DeptCode','Stock','Capital','TreeCover','PopSuit',...
    'DistCoastSuit','DistBorderSuit','SlopeSuit','MktAccSuit','LandUseSuit',...
    'LandSuit','RoadFlag'});
nnodes=height(NodeTable);
ADJ=zeros(nnodes);      % adjacency matrix for trafficking network
DIST=zeros(nnodes);     % geographic distance associated with edges
ADDVAL=zeros(nnodes);    % added value per edge in trafficking network
WGHT=zeros(nnodes);     % dynamic weighting of edges
FLOW=zeros(nnodes,nnodes,TMAX);     % dynamic flows of cocaine between nodes
SLRISK=slprob_0*ones(nnodes);     % dynamic perceived risk of seisure and loss per edge by node agent
INTRISK=zeros(nnodes,TMAX); % dynamic perceived risk of interdiction at each node
CPCTY=zeros(nnodes);     % maximum flow possible between nodes
CTRANS=zeros(nnodes);   % transportation costs between nodes
RMTFAC=zeros(nnodes);   % landscape factor (remoteness) influencing S&L risk
COASTFAC=zeros(nnodes);   % landscape factor (distance to coast) influencing S&L risk
LATFAC=zeros(nnodes);   % decreased likelihood of S&L moving north to reflect greater DTO investment
rentcap=0.4*ones(nnodes,1);     % proportion of value of shipments 'captured' by nodes

NEIHOOD=zeros(size(ca_adm0,1),size(ca_adm0,2),nnodes);

routepref=zeros(nnodes,nnodes,TMAX);   % weighting by network agent of successful routes
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
activeroute=cell(nnodes,TMAX);  % track active routes
avgslrisk=cell(nnodes,TMAX);    % average S&L risk at each node given active routes
totslrisk=zeros(1,TMAX);        % network-wide average S&L risk
slcpcty=zeros(1,TMAX);
for k=1:nnodes-1
    if k == 1
        newedges=2:nnodes;
        nodechk=ismember(newedges,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2)); %check for redundant edges
        newedges=newedges(~nodechk);
        weights=ones(length(newedges),1);
        flows=ones(length(newedges),1);
        cpcty=2*stock_0*ones(length(newedges),1);
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
        cpcty=2*stock_0*ones(length(newedges),1);
        EdgeTable=table([EdgeTable.EndNodes; k*ones(length(newedges),1) newedges'],...
            [EdgeTable.Weight; weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
            'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
        
        % Create adjacency matrix (without graph toolbox)
        ADJ(k,EdgeTable.EndNodes(EdgeTable.EndNodes(:,1)==k,2))=1;
    end
end
% Make sure all nodes connect to end node
iendnode=NodeTable.ID(NodeTable.DeptCode == 2);
newedges=1:nnodes-1;
nodechk=ismember(newedges,EdgeTable.EndNodes(EdgeTable.EndNodes(:,2)==iendnode,1)); %check for redundant edges
newedges=newedges(~nodechk);
weights=ones(length(newedges),1);
flows=ones(length(newedges),1);
cpcty=2*stock_0*ones(length(newedges),1);
EdgeTable=table([EdgeTable.EndNodes; newedges' iendnode*ones(length(newedges),1)],...
    [EdgeTable.Weight; weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
    'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});

%%% Node Attributes
% forest cover as proxy for remoteness; the higher the forest cover, the
% more remote and lower the S&L risk. Start and end node unchanged.
% remotefac=[0; 1-NodeTable.PopSuit(NodeTable.PopSuit~=0); 0];
remotefac=[0; 1-NodeTable.PopSuit(2:nnodes-1); 0];
% proximity to the coast also increases risk of S&L event
% Find node distance to coast
% lats_in=NodeTable.Lat;
% lons_in=NodeTable.Lon;
% [dists_min,lats_closest,lons_closest]=dist_from_coast(lats_in,...
%     lons_in);
icoastdist=sub2ind(size(dcoast),NodeTable.Row,NodeTable.Col);
coastdist=dcoast(icoastdist);  %convert to km
NodeTable.CoastDist=coastdist;
% coastfac=2-NodeTable.CoastDist(:)./max(NodeTable.CoastDist(:));
coastfac=[0; NodeTable.CoastDist(2:nnodes-1)./max(NodeTable.CoastDist(:)); 0];
nwvec=sqrt(0.9.*NodeTable.Lat(2:nnodes-1).^2+0.1.*NodeTable.Lon(2:nnodes-1).^2);
latfac=[0; 1-nwvec./max(nwvec); 0];
% Create adjacency matrix (without graph toolbox)
ADJ(EdgeTable.EndNodes(EdgeTable.EndNodes(:,2)==iendnode,1),iendnode)=1;
iedge=find(ADJ == 1);
for j=1:nnodes
    %Create weight and capacity matrices
    WGHT(j,ADJ(j,:)==1)=EdgeTable.Weight(ADJ(j,:)==1);
    CPCTY(j,ADJ(j,:)==1)=EdgeTable.Capacity(ADJ(j,:)==1);
     
    % Create distance (in km) matrix
    latlon2=[NodeTable.Lat(ADJ(j,:)==1) NodeTable.Lon(ADJ(j,:)==1)];
    latlon1=repmat([NodeTable.Lat(j) NodeTable.Lon(j)],length(latlon2(:,1)),1);
    [d1km,d2km]=lldistkm(latlon1,latlon2);
    DIST(j,ADJ(j,:)==1)=d1km;
    
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
    RMTFAC(j,ADJ(j,:)==1)=remotefac(ADJ(j,:)==1);
    COASTFAC(j,ADJ(j,:)==1)=coastfac(ADJ(j,:)==1);
    LATFAC(j,ADJ(j,:)==1)=latfac(ADJ(j,:)==1);
    
    % Transportation costs
    idist_ground=(DIST(j,:)  >0 & DIST(j,:) <= 500);
    idist_air=(DIST(j,:) > 500);
    if NodeTable.CoastDist(j) < 20
        ireceiver=EdgeTable.EndNodes(EdgeTable.EndNodes(:,1) == j,2);
        idist_coast=(NodeTable.CoastDist(ireceiver) < 20);
        idist_inland=(NodeTable.CoastDist(ireceiver) >= 20);
%         CTRANS(j,ireceiver(idist_coast))=ctrans_coast.*...
%             COASTFAC(j,ireceiver(idist_coast)).*DIST(j,ireceiver(idist_coast));
%         CTRANS(j,ireceiver(idist_inland))=ctrans_inland.*...
%             RMTFAC(j,ireceiver(idist_inland)).*DIST(j,ireceiver(idist_inland));
        CTRANS(j,ireceiver(idist_coast))=ctrans_coast.*...
            DIST(j,ireceiver(idist_coast));
        CTRANS(j,ireceiver(idist_inland))=ctrans_inland.*...
            DIST(j,ireceiver(idist_inland));
    else
        ireceiver=EdgeTable.EndNodes(EdgeTable.EndNodes(:,1) == j,2);
%         CTRANS(j,ireceiver)=ctrans_inland.*RMTFAC(j,ireceiver).*...
%             DIST(j,ireceiver);
        CTRANS(j,ireceiver)=ctrans_inland.*DIST(j,ireceiver);
    end

    % Create gridded distance from node for calculating land use
    % neighborhood
    if j > 1 || j < nnodes
        hdir=repmat([(NodeTable.Col(j)-1):-1:1 0 1:(size(ca_adm0,2)-...
            NodeTable.Col(j))],size(ca_adm0,1),1);
%         hdir(luint == 0 | luint == 8)=NaN;
        vdir=repmat([((NodeTable.Row(j)-1):-1:1)'; 0; (1:(size(ca_adm0,1)-...
            NodeTable.Row(j)))'],1,size(ca_adm0,2));
%         vdir(luint == 0 | luint == 8)=NaN;
        cmpdir=zeros(size(hdir,1),size(hdir,2),2);
        cmpdir(:,:,1)=hdir;
        cmpdir(:,:,2)=vdir;
        NEIHOOD(:,:,j)=mean(cmpdir,3);
    end
end


%%% Initialize Interdiction agent
%%% Initialize Interdiction agent
facmat=LATFAC;
facmat(:,:,2)=COASTFAC;
facmat(:,:,3)=RMTFAC;
% SLPROB(:,:,TSTART)=max(min(DIST./max(max(DIST)),1),0);
% SLPROB(:,:,TSTART)=max(min(COASTFAC.*RMTFAC+(DIST./max(max(DIST))),1),0);
SLPROB(:,:,TSTART)=max(min(max(facmat,[],3)+DIST./max(max(DIST)),1),0);   % dynamic probability of seisure and loss at edges
slmin=SLPROB(:,:,TSTART);
SLPROB(:,:,TSTART+1)=SLPROB(:,:,TSTART);
INTRDPROB(:,TSTART+1)=slprob_0*ones(nnodes,1); % dynamic probability of interdiction at nodes

%%% Initialize Node agents
STOCK(:,TSTART)=NodeTable.Stock(:);
TOTCPTL(:,TSTART)=NodeTable.Capital(:);
PRICE(:,TSTART+1)=PRICE(:,TSTART);
slcpcty(TSTART+1)=slcpcty_0;
%%% Set-up node and network risk perceptions
% SLRISK(:,:)=(DIST./max(max(DIST)));
% SLRISK(:,:)=SLPROB(:,:,TSTART);
% INTRISK(:,TSTART:TSTART+1)=slprob_0.*ones(nnodes,2);

% subjective risk perception with time distortion
twght=timewght_0*ones(nnodes,1);    % time weighting for dynamic, subjective perceived risk of interdiction event
    
%%% Set-up trafficking network benefit-cost logic  %%%%%%%%%%%%
ltcoeff=ones(nnodes,1);
margprofit=ADDVAL-CTRANS;
routepref(1,:,TSTART+1)=(margprofit(1,:)==max(margprofit(1,:)));
totslrisk(TSTART+1)=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   Node Land Use Choices   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Land use choices are made annually, unlike trafficking choices
%%% (monthly)
%%% Land uses based on Aide data
% 1. built-up
% 2. cropland — row crop agriculture
% 4. trees
% 5. pastureland/grassland — grazing land or natural grassland
% (speculative)
% 7. plantation — citrus, vineyard, coffee, etc.
% 8. water
% 9. Livestock - additional to Aide data
% 10. Road building - additional to Aide data
Nuse=length(unique(luint(~isnan(luint))))+2; % 
iluse=[unique(luint(~isnan(luint))); 9; 10];
luprice=zeros(Nuse,TMAX);   % price time series for each land use
LU=zeros(size(luint,1),size(luint,2),(TMAX/12)+1);
LU(:,:,TSTART)=luint;

discount=0.05;
celldist=0.25;  % 250 meter resolution, in km
cell2ha=6.25;   % 250 meter resolution equals 6.25 ha
bu2kg=27.22;    %bushel of wheat = 27.22kg
lbs2kg=0.45359237;  %1lbs = 0.453 ... kg
mi2km=1.60934;  % miles to km, 1 mile equals ... km
crop2meat=1/7;   %Conversion factor for beef, Trostle (2010)
avgcowwght=round(638*lbs2kg); %ISU extension

%%%% Per cell costs %%%%
tree2past=cell2ha*456;      % land clearing from forest
x2lvstck=cell2ha*(218000+1037*200)/400;    %assuming 1 head/2ha; captial cost, herd costs (all PV)
lvstckcost=cell2ha*(711/2);    %assuming 1 head/2ha; annual ownership and operating costs (PV)
x2plnt=cell2ha*1960-tree2past;    %assuming 8000 ha farm
plntcost=cell2ha*(36704/25);  %operating costs over assumed lifetime of 25 years
x2crop=cell2ha*743;     % fixed costs and variable costs, including cash rent equivalent, ISU extension and outreach
cropcost=cell2ha*(100+373);  % fixed costs (less land) and variable costs, ISU extension and outreach
roadcost=(plntcost*8000)/cell2ha;   %assuming capital to build plantation is related, cost to build infrastructure
% Land use costs include fixed, capital costs and variable operation and
% labor costs
lucost=[...
    0 0                  0 0         0                  0 0 0;
    0 cropcost           0 0         x2plnt             0 x2lvstck 0;
    0 tree2past+x2crop   0 tree2past tree2past+x2plnt   0 tree2past+x2lvstck 0;
    0 x2crop             0 0         x2plnt             0 x2lvstck 0;
    0 0                  0 0         plntcost           0 0 0;
    0 0                  0 0         0                  0 0 0;
    0 0                  0 0         0                  0 lvstckcost 0;
    0 0                  0 0         0                  0 0 roadcost];

%%%% Land conversion thresholds (minimum capital needed for assumed size)
cropminsize=round(100/cell2ha);
cropthresh=x2crop*cropminsize;
plntminsize=ceil(8000/cell2ha);
palmthresh=x2plnt*plntminsize;
cattleminsize=ceil(1000/cell2ha);
cattlethresh=x2lvstck*cattleminsize;
minsize=zeros(Nuse,1);
minsize(2)=cropminsize;
minsize(5)=plntminsize;
minsize(7)=cattleminsize;
minsize(10)=plntminsize;    %need to have cleared enough land for plantation to build road

%%%% Transport costs %%%%
lvstcktrans=celldist*mi2km*10; % per kg/cell transport costs, ISU extension
palmtrans=celldist*mi2km*60; %assuming this includes offsite processing costs
maizetrans=celldist*mi2km*2;
transdist=dcoast.*(1-mktacc);

%%% Code used to generate distance to palm oil mill points
% palmtransdist=zeros(size(ca_adm0,1),size(ca_adm0,2),length(pomilllat));
% for k=1:length(pomilllat)
%     [porow,pocol]=latlon2pix(Rslope,pomilllat(k),pomilllon(k));
%     porow=round(porow);
%     pocol=round(pocol);
%     hdir=repmat([(pocol-1):-1:1 0 1:(size(ca_adm0,2)-...
%         pocol)],size(ca_adm0,1),1);
%     %         hdir(luint == 0 | luint == 8)=NaN;
%     vdir=repmat([((porow-1):-1:1)'; 0; (1:(size(ca_adm0,1)-...
%         porow))'],1,size(ca_adm0,2));
%     %         vdir(luint == 0 | luint == 8)=NaN;
%     cmpdir=zeros(size(hdir,1),size(hdir,2),2);
%     cmpdir(:,:,1)=hdir;
%     cmpdir(:,:,2)=vdir;
%     palmtransdist(:,:,k)=mean(cmpdir,3);
%     % latlon2=[pomilllat pomilllon];
%     % latlon1=repmat([NodeTable.Lat(k) NodeTable.Lon(k)],length(latlon2(:,1)),1);
%     % [d1km,d2km]=lldistkm(latlon1,latlon2);
%     % palmtransdist(k)=min(d1km);
% end
% pomilltransdist=min(palmtransdist,[],3);

%%%% Commodity Prices %%%%
stckrate=diff([280/554 2312/3500; 19 34],1,2);  %change in revenue with stocking rate
plvstck=ones(size(ca_adm0));    %$/kg, FAOSTAT
pmaize=ones(size(ca_adm0));    %$/kg, FAOSTAT
ppalm=ones(size(ca_adm0));    %$/kg, FAOSTAT
plvstck(ca_adm0 == 23)=0;     % Belize(23)
pmaize(ca_adm0 == 23)=0.1;
ppalm(ca_adm0 == 23)=0;
plvstck(ca_adm0 == 55)=2.16;    % Costa Rica(55),
pmaize(ca_adm0 == 55)=0.14;
ppalm(ca_adm0 == 55)=0.38;
plvstck(ca_adm0 == 70)=8;   % El Salvador(70)
pmaize(ca_adm0 == 70)=0.16;
ppalm(ca_adm0 == 70)=0.64;
plvstck(ca_adm0 == 94)=0;   %Guatemala(94)
pmaize(ca_adm0 == 94)=0.7;
ppalm(ca_adm0 == 94)=0.38;
plvstck(ca_adm0 == 101)=1.71;   %Honduras(101)
pmaize(ca_adm0 == 101)=0.21;
ppalm(ca_adm0 == 101)=0.42;
plvstck(ca_adm0 == 161)=1.79;   %Nicuragua(161),
pmaize(ca_adm0 == 161)=0.15;
ppalm(ca_adm0 == 161)=0.72;
plvstck(ca_adm0 == 173)=0;  %Panama(173)
pmaize(ca_adm0 == 173)=0;
ppalm(ca_adm0 == 173)=0.47;


%%%% Land Use Productivity %%%%
YIELDS=zeros(size(ca_adm0,1),size(ca_adm0,2),Nuse);
LUPROD=zeros(size(ca_adm0,1),size(ca_adm0,2),Nuse);
PROD=cell(nnodes,TMAX/12);       %[cellid,existing lu,scalefac,profit,yield,cost]
EXPTPROD=cell(nnodes,Nuse,TMAX);   % Expected productivity

%%% Derive regression relationship to fill-out missing areas in official
%%% data - correlate slope suitability and yields
% mazmdl=fitlm(mktacc(~isnan(mazyld)),...
%     mazyld(~isnan(mazyld)),'linear','RobustOpts','on');
% ctlmdl=fitlm(slope(~isnan(ctlden)),...
%     ctlden(~isnan(ctlden)),'linear','RobustOpts','on');
load yieldmdls.mat

maizeyld=mazyld;
maizeyld(isnan(mazyld))=predict(mazmdl,mktacc(isnan(mazyld)));
maizeyld(luint == 8)=NaN;
palmyld=plmyld;
for cc=1:length(nodeorder(:,1))
    palmyld(ca_adm0 == nodeorder(cc,1))=nanmean(plmyld(ca_adm0 == ...
        nodeorder(cc,1)));
end
palmyld(luint == 8)=NaN;
cattleyld=ctlden;
cattleyld(isnan(ctlden))=predict(ctlmdl,slope(isnan(ctlden)));
cattleyld(luint ==8)=NaN;

YIELDS(:,:,2)=maizeyld*1000;  %convert to kg per cell
YIELDS(:,:,5)=palmyld*1000;  %convert to kg per cell
YIELDS(:,:,7)=avgcowwght*cattleyld;   %multiply this by head of cattle, assumes 1 head/2 ha

LUPROD(:,:,2)=pmaize.*YIELDS(:,:,2)-maizetrans*transdist;
LUPROD(:,:,5)=ppalm.*YIELDS(:,:,5)-palmtrans*celldist*pomilltransdist;
LUPROD(:,:,7)=plvstck.*YIELDS(:,:,7)-lvstcktrans*transdist;

OWN=zeros(size(LANDSUIT));  % node agent land ownership
IOWN=cell(nnodes,TMAX/12);     % dynamic list of owned parcels
for n=2:nnodes-1
    inodeown=sub2ind(size(LANDSUIT),NodeTable.Row(n),NodeTable.Col(n));
    IOWN{n,TSTART}=inodeown;
    OWN(inodeown)=n;
end

%%% Set-up figure for trafficking movie
MOV=zeros(nnodes,nnodes,TMAX);

%%
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@ Dynamics @@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

for t=TSTART+1:TMAX
    %%%%%% S&L and interdiction events %%%%%%
%     %%% Fully random S&L events
%     rndslevents=ones(size(ADJ));
%     rndslevents(iedge(randperm(length(iedge),slcpcty)))=rand(slcpcty,1);
%     slevent(:,:,t)=(SLPROB(:,:,t) > rndslevents);
    %%% Random number of events, selection of highest probability nodes (in addition to p=1) 
%     slcpcty(t)=slcpcty(t-1);
%     rndslevents=ceil(slcpcty(t)*rand(1));
    subslevent=slevent(:,:,t);
    subslprob=reshape(SLPROB(:,:,t),nnodes*nnodes,1);
%     subslprob=subslprob(subslprob~=1);
    [subslprobsort,isubslprobsort]=sort(subslprob,'descend');
    islevent=isubslprobsort(1:length(find(subslprob==1))+slcpcty(t));
%     islevent=find(ismember(SLPROB(:,:,t),subslprobsort(1:rndslevents))==1);
    subslevent(islevent)=1;
    slevent(:,:,t)=subslevent;
%     slevent(:,:,t)=(SLPROB(:,:,t) == 1);
    intrdevent(:,t)=(INTRDPROB(:,t) > rand(nnodes,1));
    MOV(:,1,t)=NodeTable.Stock(:);
    
    PRICE(:,t)=PRICE(:,t-1);    %price constant for now
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
%          inei=find(ADJ(n,:)==1);
         inei=find(ADJ(n,:) == 1 & routepref(n,:,t) > 0);
         %%% Procedure for selecting routes based on expected profit %%%
         c_trans=CTRANS(n,inei);
         p_sl=SLRISK(n,inei);
         p_int=INTRISK(n,inei); %currently not used
         y_node=ADDVAL(n,inei);
         q_node=min(floor(STOCK(n,t)./length(inei)),CPCTY(n,inei));
         lccf=ltcoeff(n);
         totstock=STOCK(n,t);
         totcpcty=CPCTY(n,inei);
         tslrisk=totslrisk(t);
         rtpref=routepref(n,inei,t);

%          [neipick,neivalue]=calc_neival(c_trans,p_sl,y_node,q_node,lccf,...
%              totstock,totcpcty,tslrisk);
         [neipick,neivalue]=calc_neival(c_trans,p_sl,y_node,q_node,lccf,...
             rtpref,tslrisk);
         inei=inei(neipick);
         activeroute(n,t)=mat2cell(inei',length(inei),1);
         
         % weight according to salience value fuction
         if isempty(find(neivalue > 0,1)) == 1
             WGHT(n,inei)=1;
         else
%              WGHT(n,inei)=1+(abs(routepref(inei,t).*neivalue)./...
%                  sum(abs(routepref(inei,t).*neivalue))-1/length(inei));
             WGHT(n,inei)=1+(abs(neivalue)./sum(abs(neivalue))-1/length(inei));
         end
         
         %%% !!! Put checks in to make sure buying node has enough capital
%          FLOW(n,inei,t)=min(floor(WGHT(n,inei).*(STOCK(n,t)/length(inei))),CPCTY(n,inei));
         FLOW(n,inei,t)=min(WGHT(n,inei).*(STOCK(n,t)/length(inei)),CPCTY(n,inei));
         % Check for S%L event
         if isempty(find(ismember(find(slevent(n,:,t)),inei),1)) == 0
             isl=(slevent(n,inei,t)==1);
             slsuccess(n,inei(isl),t)=FLOW(n,inei(isl),t);
             if slsuccess(n,inei(isl),t) == 0
                 slevent(n,inei(isl),t)=0;
             end
             OUTFLOW(n,t)=sum(FLOW(n,inei,t));
             STOCK(n,t)=STOCK(n,t)-OUTFLOW(n,t);
             FLOW(n,inei(isl),t)=0;     % remove from trafficking route due to S&L event
             STOCK(inei,t)=STOCK(inei,t)+FLOW(n,inei,t)';
             %              TOTCPTL(n,t)=TOTCPTL(n,t)+sum(FLOW(n,inei,t).*ADDVAL(n,inei));
             TOTCPTL(n,t)=TOTCPTL(n,t)+sum(FLOW(n,inei,t)'.*PRICE(inei,t));
             TOTCPTL(inei,t)=TOTCPTL(inei,t-1)-(FLOW(n,inei,t)'.*PRICE(inei,t));
             %              ICPTL(n,t)=rentcap(n)*sum(FLOW(n,inei).*ADDVAL(n,inei));
             ICPTL(n,t)=rentcap(n)*TOTCPTL(n,t);
         else
             OUTFLOW(n,t)=sum(FLOW(n,inei,t));
             STOCK(n,t)=STOCK(n,t)-OUTFLOW(n,t);
             STOCK(inei,t)=STOCK(inei,t)+FLOW(n,inei,t)';
%              TOTCPTL(n,t)=TOTCPTL(n,t)+sum(FLOW(n,inei,t).*ADDVAL(n,inei));
             TOTCPTL(n,t)=TOTCPTL(n,t)+sum(FLOW(n,inei,t)'.*PRICE(inei,t));
             TOTCPTL(inei,t)=TOTCPTL(inei,t-1)-(FLOW(n,inei,t)'.*PRICE(inei,t));
%              ICPTL(n,t)=rentcap(n)*sum(FLOW(n,inei).*ADDVAL(n,inei));
             ICPTL(n,t)=rentcap(n)*TOTCPTL(n,t);
         end
         %%%% Update perceived risk in response to S&L and Interdiction events
         timeweight=twght(n);
         % identify neighbors in network (without network toolbox)
         %       bcknei=EdgeTable.EndNodes(EdgeTable.EndNodes(:,2) == n,1)';
         fwdnei=inei;
         if t == TSTART+1
             %           sloccur=slevent(n,[bcknei fwdnei],TSTART+1:t);
             sloccur=slevent(n,fwdnei,TSTART+1:t);
         elseif t > TSTART+1 && length(fwdnei) == 1
             %           sloccur=squeeze(slevent(n,[bcknei fwdnei],TSTART+1:t))';
             sloccur=squeeze(slevent(n,fwdnei,TSTART+1:t));
         else
             sloccur=squeeze(slevent(n,fwdnei,TSTART+1:t))';
         end
         %       intrdoccur=intrdevent([bcknei fwdnei],TSTART+1:t);
         intrdoccur=intrdevent(fwdnei,TSTART+1:t);
         [sl_risk,intrd_risk,slevnt,intrdevnt,tmevnt]=calc_intrisk(sloccur,...
             intrdoccur,t,TSTART,alpharisk,betarisk,timeweight);
         %       SLRISK(n,[bcknei fwdnei])=sl_risk;
         SLRISK(n,fwdnei)=sl_risk;
         
         if isempty(find(sl_risk,1)) == 0
             avgslrisk(n,t)=mat2cell(SLRISK(n,activeroute{n,t}),1,...
                 length(activeroute{n,t}));
         end
         INTRISK(n,t+1)=mean(intrd_risk);  %node-specific risk is the average of neighbor risks
         
         NodeTable.Stock(:)=STOCK(:,t);
         NodeTable.Capital(:)=TOTCPTL(:,t);
       end
       
       %%% Make trafficking movie
       MOV(:,n,t)=STOCK(:,t);      % Capture stock data after each node iteration
    end
    totslrisk(t+1)=mean(cat(2,avgslrisk{:,t}));
    %%% Updating interdiction event probability
    subslsuc=slsuccess(:,:,t);
    subslprob=SLPROB(:,:,t);
    islcheck=(slevent(:,:,t) == 1);
%     SLPROB(:,:,t+1)=max((1-delta_sl).*SLPROB(:,:,t)+delta_sl.*...
%         (slsuccess(:,:,t)./max(max(slsuccess(:,:,t)))),SLPROB(:,:,TSTART));
    subslprob(islcheck)=max((1-delta_sl).*subslprob(islcheck)+delta_sl.*...
        (subslsuc(islcheck) > 0),slmin(islcheck));
    SLPROB(:,:,t+1)=subslprob;
    %%% Eventually, this should be tied to perception of negative
    %%% consequences of trafficking (e.g., violence, lost profits from
    %%% licit markets, etc.)
    slcpcty(t+1)=max(slcpcty(t)+ceil(delta_sl*(sum(sum(slsuccess(:,:,t)))-...
        sum(sum(slsuccess(:,:,t-1))))),slcpcty_0);
    
    INTRDPROB(:,t+1)=INTRDPROB(:,t);
    
    % Reinforcement learning for successful routes
    iactivenode=find(OUTFLOW(2:nnodes-1,t) > 0)+1;
    avgflow=STOCK(iendnode,t)/length(iactivenode);
%     subroutepref=zeros(nnodes,1);
%     for nn=2:nnodes-1
%         if isempty(activeroute{nn,t}) == 1
% %             routepref(nn,t+1)=(1-delta_rt).*routepref(nn,t-1);
%             subroutepref(nn)=(1-delta_rt).*routepref(nn,t-1);
%         else
%             rtwght=mean(FLOW(nn,activeroute{nn,t},t)./avgflow);
% %             routepref(nn,t+1)=(1-delta_rt).*routepref(nn,t)+delta_rt.*rtwght;
%             subroutepref(nn)=(1-delta_rt).*routepref(nn,t)+delta_rt.*rtwght;
%         end
%     end
%     routepref(:,t+1)=subroutepref./max(subroutepref);
% %     routepref(iendnode,t+1)=1;
% %     routepref(iendnode,t+1)=1.1*max(routepref(:,t+1));

    subroutepref=routepref(:,:,t);
    activenodes=unique(cat(1,activeroute{:,t}));
    subflow=FLOW(:,:,t);
    %%% calculate losses from S&L events
%     volume-based - does not matter where in supply chain
%     supplyfit=STOCK(iendnode,t)/stock_0;
%     losstolval=losstol*stock_0;

    % value-based - price varies with location in supply chain
    ipossl=find(slsuccess(:,:,t)>0);
    [nrow,ncol]=ind2sub(size(slsuccess(:,:,t)),ipossl);
%     supplyfit=PRICE(nnodes,t)*(STOCK(iendnode,t)/stock_0);
    supplyfit=stock_0*PRICE(nnodes,t)-sum(subslsuc(ipossl).*PRICE(ncol,t));  %value-based loss calc
    losstolval=losstol*stock_0*PRICE(nnodes,t); %value-based loss threshold
    
    %call top-down route optimization
    newroutepref=optimizeroute(nnodes,subflow,supplyfit,activenodes,...
        subroutepref,EdgeTable,SLRISK,ADDVAL,CTRANS,losstolval);
    routepref(:,:,t+1)=newroutepref;

    PRICE(:,t+1)=PRICE(:,t);
    STOCK(1,t+1)=stock_0;    %additional production to enter network next time step
    STOCK(nnodes,t+1)=0;    %remove stock at end node for next time step
    NodeTable.Stock(1)=stock_0;
    NodeTable.Stock(nnodes)=0;
    
    if rem(t,12) == 0       % land use decision-making loop
        lt=t/12;    %land use time indexing
        subcropprod=LUPROD(:,:,2);
        subplntprod=LUPROD(:,:,5);
        subctlprod=LUPROD(:,:,7);
        sublu=LU(:,:,lt);
        isublu=sublu;
        for m=1:Nuse
            isublu(sublu == iluse(m))=m;
        end
        for n=2:nnodes-1
            if TOTCPTL(n,t) < 0
                LU(:,:,lt+1)=sublu;
                if lt > 1
                    IOWN(n,lt)=IOWN(n,lt-1);
                end
                continue
            end
%             if lt>1
%                 keyboard
%             end
            % Update licit revenues from land use
            neihood=NEIHOOD(:,:,n);     % find closet, contiguous cells
            neihood(isnan(mktacc) | sublu==1 | sublu== 8)=10000;            
            if lt > 1 && isempty(PROD{n,lt-1}) == 0 
                % need to recalculate with operating costs
                subprod=PROD{n,lt-1};   %[cellid,existing lu,scalefac,profit,yield,cost]
                ulu=unique(subprod(:,2));
                for iu=1:length(ulu)
                    ilu=(subprod(:,2)==ulu(iu));
                    econfac=subprod(ilu,3).*1/(1+(length(subprod(ilu,1))-...
                        minsize(ulu(iu) == iluse)));
                    cellprod=subprod(ilu,5)-econfac.*...
                        lucost(isublu(subprod(ilu,1)),ulu(iu) == iluse);
                    LCPTL(n,t)=LCPTL(n,t)+sum(cellprod);
                    TOTCPTL(n,t)=TOTCPTL(n,t)-sum(lucost(isublu(subprod(ilu,1)),ulu(iu) == iluse));
                    subprod(ilu,4)=cellprod;
                    subprod(ilu,6)=lucost(isublu(subprod(ilu,1)),ulu(iu) == iluse);
                    neihood(subprod(:,1))=10000;
                end
                PROD(n,lt)=mat2cell(subprod,size(subprod,1),size(subprod,2));
            end
            neiind=[(1:(size(sublu,1)*size(sublu,2)))' ...
                reshape(neihood,size(sublu,1)*size(sublu,2),1)];
            inanland=find(isnan(sublu) == 1 | OWN == n);
            neiind=neiind(~ismember(neiind(:,1),inanland),:);
%             proxcells=sort(reshape(neihood,size(neihood,1)*size(neihood,2),...
%                 1),'ascend');
            proxcells=sortrows(neiind,2);
            
            %%%% Establishing new land uses %%%%
            potprod=[];
            if LUPROD(NodeTable.Row(n),NodeTable.Col(n),2) > 0 && ...
                    TOTCPTL(n,t) >= cropthresh && ...
                    NodeTable.MktAccSuit(n) >= 0.05 && ...
                    length(find(sublu(IOWN{n,lt-1}) == 7)) >= minsize(10)
                % row crop land use
%                 icropcells=find(neihood <= proxcells(cropminsize) & ...
%                     subcropprod(~isnan(subcropprod)));
                icropcells=neiind(neiind(:,2) <= proxcells(minsize(2),2),1);
                cropscale=1-max((sum(cropcost*ones(minsize(2),1))-...
                    sum(sum(subcropprod(icropcells(1:minsize(2))))./...
                    (1+discount).^(1:25)))/sum(cropcost*ones(minsize(2),1)),0);
                cropecon=cropscale*1/(1+(length(icropcells)-minsize(2)));
                pland=landprice(icropcells);
                pland(ismember(icropcells,IOWN{n,lt-1}))=0;
                cellprod=subcropprod(icropcells)-pland-...
                    cropecon*lucost(isublu(icropcells),2);
                potcropprod=[icropcells iluse(2)*ones(length(icropcells),1) ...
                    sublu(icropcells) cropscale*ones(length(icropcells),1) ...
                    cellprod subcropprod(icropcells) lucost(isublu(icropcells),2)];
                icrop=find(cumsum(potcropprod(:,7)) >= TOTCPTL(n,t),1,'first');
                if isempty(find(icrop,1)) == 1
                    icrop=length(potcropprod(:,1));
                end
                potprod=[potprod; potcropprod(1:icrop,:)];
            end
            
            if LUPROD(NodeTable.Row(n),NodeTable.Col(n),5) > 0 && ...
                    TOTCPTL(n,t) >= palmthresh && ...
                    NodeTable.MktAccSuit(n) >= 0.05 && ...
                    length(find(sublu(IOWN{n,lt-1}) == 7)) >= minsize(10)
                % plantation (palm oil) land use
%                 iplntcells=find(neihood <= proxcells(minsize(5)) & ...
%                     subplntprod(~isnan(subplntprod)));
                iplntcells=neiind(neiind(:,2) <= proxcells(minsize(5),2),1);
                plntscale=1-max((sum(plntcost*ones(minsize(5),1))-...
                    sum(sum(subplntprod(iplntcells(1:minsize(5))))./...
                    (1+discount).^(1:25)))/sum(plntcost*ones(minsize(5),1)),0);
                plntecon=plntscale*1/(1+(length(iplntcells)-minsize(5)));
                pland=landprice(iplntcells);
                pland(ismember(iplntcells,IOWN{n,lt-1}))=0;
                cellprod=subplntprod(iplntcells)-pland-...
                    plntecon*lucost(isublu(iplntcells),5);
                potplntprod=[iplntcells iluse(5)*ones(length(iplntcells),1) ...
                    sublu(iplntcells) plntscale*ones(length(iplntcells),1) ...
                    cellprod subplntprod(iplntcells) lucost(isublu(iplntcells),5)];
                iplnt=find(cumsum(potplntprod(:,7)) >= TOTCPTL(n,t),1,'first');
                if isempty(find(iplnt,1)) == 1
                    iplnt=length(potplntprod(:,1));
                end
                potprod=[potprod; potplntprod(1:iplnt,:)];
            end
            
            if LUPROD(NodeTable.Row(n),NodeTable.Col(n),7) > 0 && ...
                    TOTCPTL(n,t) >= cattlethresh
                %cattle ranching
%                 ictlcells=find(neihood <= proxcells(minsize(7)) & ...
%                     subctlprod(~isnan(subctlprod)));
                ictlcells=neiind(neiind(:,2) <= proxcells(minsize(7),2),1);
                ctlscale=1-max((sum(lvstckcost*ones(minsize(7),1))-...
                    sum(sum(subctlprod(ictlcells(1:minsize(7))))./...
                    (1+discount).^(1:25)))/sum(lvstckcost*ones(minsize(7),1)),0);
                ctlecon=ctlscale*1/(1+(length(ictlcells)-minsize(7)));
                pland=landprice(ictlcells);
                pland(ismember(ictlcells,IOWN{n,lt-1}))=0;
                cellprod=subctlprod(ictlcells)-pland-...
                    ctlecon*lucost(isublu(ictlcells),7);
                potctlprod=[ictlcells iluse(7)*ones(length(ictlcells),1) ...
                    sublu(ictlcells) ctlscale*ones(length(ictlcells),1) ...
                    cellprod subctlprod(ictlcells) lucost(isublu(ictlcells),7)];
                ictl=find(cumsum(potctlprod(:,7)) >= TOTCPTL(n,t),1,'first');
                if isempty(find(ictl,1)) == 1
                    ictl=length(potctlprod(:,1));
                end
                potprod=[potprod; potctlprod(1:ictl,:)];
            end
            
            %%% Add road building decision
            if NodeTable.RoadFlag(n) == 0 && TOTCPTL(n,t) >= roadcost && ...
                    length(find(sublu(IOWN{n,lt-1}) == 7)) >= minsize(10)
                display('Check extent of cleared land')
                NodeTable.RoadFlag(n)=1;
                NodeTable.MktAccSuit(n)=max(0.05,NodeTable.MktAccSuit(n));
                keyboard
            end
            
            if isempty(find(potprod,1)) == 1
                LU(:,:,lt+1)=sublu;
                if lt > 1
                    IOWN(n,lt)=IOWN(n,lt-1);
                end
                continue
            end
            [maxprod,~]=max([sum(sum(potprod(potprod(:,2)==iluse(2),5))./(1+discount).^(1:25)) ...
                sum(sum(potprod(potprod(:,2)==iluse(5),5))./(1+discount).^(1:25)) ...
                sum(sum(potprod(potprod(:,2)==iluse(7),5))./(1+discount).^(1:25))],...
                [],1);
            luset=iluse([2 5 7]);
            maxlu=luset(maxprod == max(maxprod));
            if length(maxlu) > 1
                maxlu=iluse(7);
            end
            imaxlu=(potprod(:,2) == maxlu);
            
            subprod=PROD{n,lt};
            if isempty(subprod) == 1 || isempty(find(subprod(:,2) == maxlu,1)) == 1
                TOTCPTL(n,t)=TOTCPTL(n,t)-sum(potprod(imaxlu,7));
                PROD(n,lt)=mat2cell([subprod; potprod(imaxlu,[1 2 4 5 6 7])],...
                    length(subprod)+size(potprod(imaxlu,[1 2 4 5 6 7]),1),...
                    size(potprod(imaxlu,[1 2 4 5 6 7]),2));   %[cellid,new lu,scalefac,profit,yield,cost]
                sublu(potprod(imaxlu,1))=potprod(imaxlu,2);
                iadd2own=~ismember(potprod(imaxlu,1),IOWN{n,lt-1});
                IOWN(n,lt)=mat2cell([IOWN{n,lt-1}; potprod(imaxlu(iadd2own),...
                    1)],length(IOWN{n,lt-1})+size(potprod(imaxlu(iadd2own),...
                    [1 2 4 5 6 7]),1));
            elseif isempty(find(subprod(:,2) == maxlu,1)) == 0
                %pick highest profit and max cost cells
                ichoose=(potprod(imaxlu,5) >=0);
                TOTCPTL(n,t)=TOTCPTL(n,t)-sum(potprod(imaxlu(ichoose),7));
                PROD(n,lt)=mat2cell([subprod; potprod(imaxlu(ichoose),[1 2 4 5 6 7])],...
                    length(subprod(:,1))+size(potprod(imaxlu(ichoose),[1 2 4 5 6 7]),1),...
                    size(potprod(imaxlu(ichoose),[1 2 4 5 6 7]),2));   %[cellid,new lu,scalefac,profit,yield,cost]
                sublu(potprod(imaxlu(ichoose),1))=potprod(imaxlu(ichoose),2);
                iadd2own=~ismember(potprod(imaxlu(ichoose),1),IOWN{n,lt-1});
                IOWN(n,lt)=mat2cell([IOWN{n,lt-1}; ...
                    potprod(imaxlu(ichoose(iadd2own)),1)],length(IOWN{n,lt-1})+...
                    size(potprod(imaxlu(ichoose(iadd2own)),[1 2 4 5 6 7]),1));
            end
            LU(:,:,lt+1)=sublu;
            OWN(IOWN{n,lt})=n;
            %%%% Rules needed
            % check if agent is already engaged in land use, update revenue
            % once size thresh passed, how long to wait to convert?
            % calculate profitability (costs) via conversion costs given
            % current land use
            % rank profitability of all cells within a given neighborhood
            % adjust illicit capital (use TOTCPTL) after conversion, yearly
            
            
        end
    end
end
cd X:\model_results\NarcoLogic_null_070417
save('narcologic_results_null_070417','EdgeTable','NodeTable','LU','MOV','FLOW',...
    'TOTCPTL','ICPTL','LCPTL','slsuccess','PROD','LUPROD','activeroute','STOCK','-v7.3')

toc     % stop run timer
%%
%%% Create LUC map for viz
difflumap=zeros(size(ca_adm0));
iforest=find(LU(:,:,1)==4);
difflumap(iforest)=4;
LUCMAP=zeros(size(ca_adm0,1),size(ca_adm0,2),16);
LUCMAP(:,:,1)=difflumap;
for lt=2:(TMAX/12)+1
    sublumap=LU(:,:,lt);
    subprvmap=LU(:,:,1);
    ichange=((sublumap-subprvmap) ~= 0);
    LUCMAP(:,:,lt)=difflumap;
    sublucmap=LUCMAP(:,:,lt);
    
%     ichange=(sublumap(iforest) ~= 4);
    sublucmap(ichange)=sublumap(ichange);
    LUCMAP(:,:,lt)=sublucmap;
end
% % %%% Visualization %%%
clrmap=[1 1 1;  %built-up, nodata
    1 1 1;  %crop
    1 1 1;  %placeholder
    0.2 0.5 0.2];    %forest

clrmap2=[1 1 1;  %built-up, nodata
    0 1 0;  %crop
    1 1 1;  %placeholder
    0.2 0.5 0.2;    %forest
    0.7 1 0;    %pasture
    1 1 1;  %placeholder
    1 0 1;  %plantation
    0 0 1]; %water

% % %%% Trafficking movie
% writerObj = VideoWriter('trafficking_lucombo_full_v1.mp4','MPEG-4');
% writerObj.FrameRate=10;
% open(writerObj);
% 
% h1=figure;
% set(h1,'Color','white','Visible','off')
% geoshow(zeros(size(LUCMAP(:,:,1))),Rcagrid_cntry,'CData',LUCMAP(:,:,1),'DisplayType','surface')
% colormap(clrmap)
% hold on
% for tt=144:180
%     if rem(tt,12) == 0
%         lt=tt/12;
%         sublumap=LU(:,:,lt);
%         subprvmap=LU(:,:,1);
%         ichange=((sublumap-subprvmap) ~= 0);
%         LUCMAP(:,:,lt)=difflumap;
%         sublucmap=LUCMAP(:,:,lt);
%         
%         %     ichange=(sublumap(iforest) ~= 4);
%         sublucmap(ichange)=sublumap(ichange);
%         LUCMAP(:,:,lt)=sublucmap;
%         geoshow(zeros(size(LUCMAP(:,:,tt/12))),Rcagrid_cntry,'CData',...
%             LUCMAP(:,:,tt/12),'DisplayType','surface')
%         colormap(clrmap2)
%     end
%     %     geoshow(CAadm0,'FaceColor',[1 1 1])
%     geoshow(CAadm0,'FaceColor','none')
%     hold on
%     %     plot(nodelon(1),nodelat(1),'r.','MarkerSize',ceil(MOV(1,1,tt)/1000))
%     plot(nodelon(1),nodelat(1),'r.','MarkerSize',stock_0)
% % %     fedge=find(FLOW(1,:,tt) > 0);
%     
%     fedge=activeroute{1,tt};
%     %
%     islevent=find(slevent(1,fedge,tt) == 1);
%     inoslevent=find(slevent(1,fedge,tt) == 0);
%     if isempty(islevent) == 1
%         for g=1:length(fedge)
%             plot([nodelon(1); nodelon(fedge(g))],[nodelat(1); nodelat(fedge(g))],'-k')
%         end
%     else
%         for g=1:length(fedge(islevent))
%             plot([nodelon(1); nodelon(fedge(islevent(g)))],[nodelat(1); ...
%                 nodelat(fedge(islevent(g)))],'-k')
%             plot(nodelon(fedge(islevent(g))),nodelat(fedge(islevent(g))),...
%               'kx','MarkerSize',ceil(slsuccess(1,fedge(islevent(g)),tt)))
%         end
%         for k=1:length(fedge(inoslevent))
%             plot([nodelon(1); nodelon(fedge(inoslevent(k)))],[nodelat(1); ...
%                 nodelat(fedge(inoslevent(k)))],'-k')
%         end
%     end
% %     for g=1:length(fedge)
% %         plot([nodelon(1); nodelon(fedge(g))],[nodelat(1); nodelat(fedge(g))],'-k')
% %     end
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
%             if rem(tt,12) == 0
%                 geoshow(zeros(size(LUCMAP(:,:,tt/12))),Rcagrid_cntry,'CData',...
%                     LUCMAP(:,:,tt/12),'DisplayType','surface')
%                 colormap(clrmap2)
%             else
%                 geoshow(zeros(size(LUCMAP(:,:,1))),Rcagrid_cntry,'CData',LUCMAP(:,:,1),'DisplayType','surface')
%                 colormap(clrmap)
%             end
%             hold on
%             geoshow(CAadm0,'FaceColor','none')
% %             hold on
%             for nn=1:nnodes
%                 if MOV(nn,mm,tt) > 0
%                     plot(nodelon(nn),nodelat(nn),'r.','MarkerSize',ceil(MOV(nn,mm,tt)))
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
%             if rem(tt,12) == 0
%                 geoshow(zeros(size(LUCMAP(:,:,tt/12))),Rcagrid_cntry,'CData',...
%                     LUCMAP(:,:,tt/12),'DisplayType','surface')
%                 colormap(clrmap2)
%             else
%                 geoshow(zeros(size(LUCMAP(:,:,1))),Rcagrid_cntry,'CData',LUCMAP(:,:,1),'DisplayType','surface')
%                 colormap(clrmap)  
%             end
%             hold on
%             geoshow(CAadm0,'FaceColor','none')
% %             hold on
%             for nn=1:nnodes
%                 if MOV(nn,mm,tt) > 0
%                     plot(nodelon(nn),nodelat(nn),'r.','MarkerSize',ceil(MOV(nn,mm,tt)))
%                 else
%                     plot(nodelon(nn),nodelat(nn),'b.','MarkerSize',3)
%                 end
%             end
% %             fedge=find(FLOW(mm,:,tt) > 0);
%             fedge=activeroute{mm,tt};
% %             if isempty(find(fedge,1)) == 1
% %                 plot(nodelon(mm),nodelat(mm),'kx','MarkerSize',ceil(MOV(mm,mm-1,tt)./1000))
% %             else
% %                 for g=1:length(fedge)
% %                     plot([nodelon(mm); nodelon(fedge(g))],[nodelat(mm); nodelat(fedge(g))],'-k')
% %                 end
% %             end
%             islevent=find(slevent(mm,fedge,tt) == 1);
%             inoslevent=find(slevent(mm,fedge,tt) == 0);
%             if isempty(islevent) == 1
%                 for g=1:length(fedge)
%                     plot([nodelon(mm); nodelon(fedge(g))],[nodelat(mm); nodelat(fedge(g))],'-k')
%                 end
%             else
%                 for g=1:length(fedge(islevent))
%                 plot([nodelon(mm); nodelon(fedge(islevent(g)))],[nodelat(mm); ...
%                     nodelat(fedge(islevent(g)))],'-k')
%                 plot(nodelon(fedge(islevent(g))),nodelat(fedge(islevent(g))),...
%                     'kx','MarkerSize',ceil(slsuccess(mm,fedge(islevent(g)),tt)))
%                 end
%                 for k=1:length(fedge(inoslevent))
%                     plot([nodelon(mm); nodelon(fedge(inoslevent(k)))],[nodelat(mm); nodelat(fedge(inoslevent(k)))],'-k')
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
% % Plot some node attribute
% figure
% geoshow(CAadm0,'FaceColor',[1 1 1])
% hold on
% for i=1:length(nodelat)
%         plot(nodelon(i),nodelat(i),'r.','MarkerSize',round(NodeTable.CoastDist(i)/10))
% end
%%% Plot time series of flows and S&L
h1_1=figure;
set(h1_1,'Color','white')
% plot(1:TMAX,STOCK(1,1:TMAX),'-b')
plot(1:TMAX,STOCK(nnodes,1:TMAX),'--k')
hold on
sltot=sum(sum(slsuccess(:,:,1:TMAX),1));
plot(1:TMAX,reshape(sltot,1,TMAX),'-r')
xlim([0 TMAX])
ylabel('Cocaine Volume kg/month')
xlabel('Month')
legend('Consumer','S&L','Orientation','horizontal','Location','southoutside')
saveas(h1_1,'Flows_vs_SL_null.png')
% 
% %%%%%% Diagnostics %%%%%%%
% slrecord=sum(slsuccess(:,:,1:t),3);
% [slr,slc]=ind2sub(size(slrecord),find(slrecord > 0));

% %%%% Land use maps
% % Initial land use
% clrmap=[0 0 0;  %built-up, nodata
%     0 1 0;  %crop
%     0 0 0;  %placeholder
%     0.2 0.5 0.2;    %forest
%     0.7 1 0;    %pasture
%     0 0 0;  %placeholder
%     1 0 1;  %plantation
%     0 0 1]; %water

% %%% change map
% 
h1_2=figure;
set(h1_2,'Color','white')
% difflumap=zeros(size(ca_adm0));
difflumap=zeros(size(LU(:,:,1)));
iforest=find(LU(:,:,1)==4);
difflumap(iforest)=4;
% LUCMAP=zeros(size(ca_adm0,1),size(ca_adm0,2),16);
LUCMAP=zeros(size(LU(:,:,1),1),size(LU(:,:,1),2),16);
LUCMAP(:,:,1)=difflumap;
for lt=2:(size(activeroute,2)/12)+1
    sublumap=LU(:,:,lt);
    subprvmap=LU(:,:,1);
    ichange=((sublumap-subprvmap) ~= 0);
    LUCMAP(:,:,lt)=difflumap;
    sublucmap=LUCMAP(:,:,lt);
    
%     ichange=(sublumap(iforest) ~= 4);
    sublucmap(ichange)=sublumap(ichange);
    LUCMAP(:,:,lt)=sublucmap;
end
imagesc(LUCMAP(:,:,lt));
clrmap=[1 1 1;  %built-up, nodata
    0 1 0;  %crop
    1 1 1;  %placeholder
    0.2 0.5 0.2;    %forest
    0.7 1 0;    %pasture
    1 1 1;  %placeholder
    1 0 1;  %plantation
    0 0 1]; %water
colormap(clrmap);
set(gca,'Visible','off')
saveas(h1_2,'LUmap_null.png')

%%% Plot time series of active nodes
nactnodes=zeros(1,TMAX);
nslsuccess=zeros(1,TMAX);
slperevent=zeros(1,TMAX);
for z=1:TMAX
    nactnodes(z)=length(cat(1,activeroute{:,z}));
    subslsccss=slsuccess(:,:,z);
    nslsuccess(z)=length(find(subslsccss > 0));
    if isempty(find(subslsccss > 0,1)) == 0
        slperevent(z)=mean(subslsccss(subslsccss > 0)./...
            length(find(subslsccss > 0)));
    else
        slperevent(z)=0;
    end
end 
h2_1=figure;
set(h2_1,'Color','white')
[hAx,hl1,hl2]=plotyy(1:TMAX,nactnodes,1:TMAX,slperevent);
% plot(1:TMAX,nactnodes,'-b')
ylabel(hAx(1),'Number of Routes')
% hold on
% yyaxis right
% plot(1:TMAX,slpervent,'--r')
ylabel(hAx(2),'Average S&L Volume (kg)')
xlim([0 TMAX])
xlabel('Month')
legend('Active Routes','S&L Volume','Orientation','horizontal','Location','southoutside')
% saveas(h2_1,'Nodes_vs_SL_null.png')

% %%% Command to use
% % digraph, maxflow, nearest