%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%   NarcoLogic ABM   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tic         % start run timer
% cd C:\Users\nmagliocca\Documents\Matlab_code\NarcoLogic
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@ Procedures @@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
TSTART=1;
TMAX=180;   % 15 years at monthly time steps

% rng default
rng(mrun)

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


% [POmills_pts,POattr0]=shaperead('X:\CentralAmericaData\Model_inputs\ca_pomill_pts.SHP','UseGeoCoords',...
%     true);


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
% cntryorder=[173 55 161 101 70 94];
cntrycodes=cntrycodes(cntrycodes ~= 0 & cntrycodes ~= cntrynodataval);
% pandptcodes=unique(dptgrid(ca_adm0 == 173));
% pandptcodes=pandptcodes(pandptcodes ~= -9999);
% meanpanlon=zeros(length(pandptcodes),1);
% for j=1:length(pandptcodes)
%     [nrow,ncol]=ind2sub(size(dptgrid),find(dptgrid == pandptcodes(j)));
%     [nlat,nlon]=pix2latlon(Rdptgrid,nrow,ncol);
%     meanpanlon(j)=mean(nlon);
% end
% sortpanlon=sortrows([pandptcodes meanpanlon],-2);

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
dptvec=sqrt(0.9*maxlat.^2+0.1*maxlon.^2);
dptmat=[dptcodes dptvec];
dptorder=sortrows(dptmat,2);

% Tree cover
[tcov,Rtcov]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\treecov_clp.tif');
treecov=tcov;
tcovnodataval=255;
treecov=double(treecov);
treecov(cagrid_cntry==cntrynodataval)=NaN;
itreecov=find(treecov > 0); %identify high forest cover areas
treecovpct=quantile(treecov(itreecov),[0.025 0.25 0.50 0.75 0.975]);
itreepick=find(treecov >= treecovpct(2) & treecov < treecovpct(3));
avgtcov=zeros(length(dptcodes),1);    %average tree cover per county, weighting for generating trade nodes
for cc=1:length(dptcodes)
    avgtcov(cc)=mean(treecov(ca_adm0 == dptcodes(cc)));
end

% Distance to coast and country borders
[dcoast,Rdcoast]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\dcoast_clp.tif');
dcoastnodataval=-9999;
dcoast(cagrid_cntry==cntrynodataval)=NaN;
dcoast(dcoast == dcoastnodataval)=NaN;
dcoast_suit=1-dcoast./max(max(dcoast));

% [dbrdr,Rdbrdr]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\dbrdr_clp.tif');
% brdrnodataval=-9999;
% dbrdr(cagrid_cntry==cntrynodataval)=NaN;
% dbrdr(dbrdr == brdrnodataval)=NaN;
% dbrdr_suit=1-dbrdr./max(max(dbrdr));
    
% Population density as a proxy for remoteness
[popden,Rpopden]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\popden_clp.tif');
popnodataval=-1;
popden=double(popden);
popden(cagrid_cntry==cntrynodataval)=NaN;
popden(popden == popnodataval)=NaN;
popq=quantile(reshape(popden,size(popden,1)*size(popden,2),1),...
    [0.25 0.5 0.75 0.95]);
pop_suit=zeros(size(popden));
pop_suit(popden > popq(3))=0;
pop_suit(popden <= popq(3))=1-popden(popden <= popq(3))./popq(3);

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
mktacc_suit=mktacc;
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
luint(cagrid_cntry==cntrynodataval)=NaN;
luint(luint == lunodataval)=NaN;
lu_suit=zeros(size(luint));
lu_suit(luint == 1 | luint == 6 | luint == 7 | luint == 8 | luint ==9)=0;
lu_suit(luint == 2)=0.5;
lu_suit(luint == 3 | luint == 4 | luint == 5)=1;

% Protected Areas
[protarea,Rprotarea]=geotiffread('X:\CentralAmericaData\Model_inputs\clipped\protarea_clp.tif');
protnodataval=255;
protarea=double(protarea);
protarea(cagrid_cntry==cntrynodataval)=NaN;
protarea(protarea == protnodataval)=NaN;
protsuit=1-isnan(protarea);

%%% Land-base investment value
invst_suit=zeros(size(luint));
invst_suit(luint == 1 | luint == 6 | luint == 7 | luint == 8 | luint ==9)=0;
invst_suit(luint == 2)=0.5;
invst_suit(luint == 3 | luint == 5)=0.75;
invst_suit(luint == 4)=1;

%%% Weight each landscape attribute
tcwght=1;       % tree cover
brdwght=1;      % distance to country border
dcstwght=0;     % distance to coast
mktwght=1;      % market access
popwght=1;      % population density - proxy for remoteness
slpwght=1;      % slope-constrained land suitability
luwght=0;       % suitability based on initial land use
invstwght=0;    % investment potential of initial land use
protwght=1;

% LANDSUIT=tcwght.*treecov./100+brdwght.*dbrdr_suit+dcstwght.*dcoast_suit+...
%     mktwght.*(1-mktacc_suit)+popwght.*pop_suit+slpwght.*slp_suit+luwght.*...
%     lu_suit+invstwght.*invst_suit+protwght.*(1-protsuit);  % land suitability based on biophysical and narco variable predictors

wghts=[tcwght brdwght dcstwght mktwght popwght slpwght luwght invstwght protwght]./...
    sum([tcwght brdwght dcstwght mktwght popwght slpwght luwght invstwght protwght]);

% %%% Null Model
% LANDSUIT=wghts(1).*treecov./100+wghts(2).*dbrdr_suit+wghts(3).*dcoast_suit+...
%     wghts(4).*mktacc_suit+wghts(5).*pop_suit+wghts(6).*slp_suit+wghts(6).*...
%     lu_suit+wghts(7).*invst_suit+wghts(8)*(1-protsuit);  % land suitability based on biophysical and narco variable predictors

% %%% Full model
LANDSUIT=wghts(1).*treecov./100+wghts(2).*dbrdr_suit+wghts(3).*dcoast_suit+...
    wghts(4).*(1-mktacc_suit)+wghts(5).*pop_suit+wghts(6).*slp_suit+wghts(6).*...
    lu_suit+wghts(7).*invst_suit+wghts(8)*protsuit;  % land suitability based on biophysical and narco variable predictors


%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%@@@@@@@@@@ Agent Attributes @@@@@@@@@@@@
%@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
%%% Interdiction Agent %%%
slprob_0=0.02;     % baseline probability of seisure and loss event
delta_sl=0.5;      % reinforcement learning rate for S&L vents (i.e., weight on new information)

%%% Network Agent %%%
ndto=2;         %initial number of DTOs
dtocutflag=zeros(ndto,1);
DTOBDGT=zeros(ndto,TMAX);
losstol=0.9;        % tolerance threshold for loss due to S&L, triggers route fragmentation
stock_0=200;     %initial cocaine stock at producer node
startvalue=385; %producer price, $385/kg
deltavalue=8;   %added value for distance traveled $8/kilo/km
nodeloss=0;     % amount of cocaine that is normally lost (i.e., non-interdiction) at each node
ctrans_inland=3.5;  % transportation costs (kg/km) over-ground
ctrans_coast=1.5;     % transportation costs (kg/km) via plane or boat
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
alpharisk=2;  %baseline is 0.001
timewght_0=0.91;     %time discounting for subjective risk perception (Gallagher, 2014), range[0,1.05]
% betarisk=alpharisk/slprob_0-alpharisk;
betarisk=0.5;
nodepct=0.0001; %percentage of high suitability cells that contain possible nodes
% cntrycpcty=[0.1 0.1 0.1 0.1 0.1 0.1 0.1];   %country-specific, per node trafficking capacity
rentcap=0.4;     % proportion of value of shipments 'captured' by nodes
bribepct=0.3;       % Annual proportion of gross profits from drug trafficking that go towards securing node territory
bribethresh=12;      % Maximum number of months a node can o without bribes to maintain control
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
nodedto=0;
nodequant=quantile(LANDSUIT(~isnan(LANDSUIT)),[0.025 0.50 0.66 0.75 0.99]);
inodepick=find(LANDSUIT > nodequant(4));
avgnodealloc=ceil((length(inodepick)*nodepct)/length(dptcodes));
pctdptsuit=zeros(length(dptorder(:,1)),1);
for dc=1:length(pctdptsuit)
    dptsuit=LANDSUIT(dptgrid == dptorder(dc,1));
    pctdptsuit(dc)=length(find(dptsuit > nodequant(4)))/length(dptsuit);
end
% for i=1:length(cntrycodes)
for i=1:length(dptorder)
%     icntry=find(ca_adm0 == nodeorder(i));
%     ipotnode=find(ismember(icntry,inodepick)==1);   %place nodes based on LANDSUIT
    idptmnt=find(dptgrid == dptorder(i,1));
    ipotnode=find(ismember(idptmnt,inodepick)==1);
    % allocate nodes per department based on suitability within department
    allocnodes=round(avgnodealloc*pctdptsuit(i)./median(pctdptsuit(pctdptsuit~=0)));
    if isempty(find(ipotnode,1)) == 1
        subinodepick=find(LANDSUIT(idptmnt) > nodequant(2));
        if isempty(find(subinodepick,1)) == 1
            continue
        end
%         ipotnode=find(ismember(idptmnt,subinodepick)==1);
        ipotnode=subinodepick;
    end
    
    randnode=idptmnt(ipotnode(randperm(max(length(ipotnode),...
        length(ipotnode)),allocnodes)));

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
    nodetcov=[nodetcov; treecov(randnode)];
    nodepopsuit=[nodepopsuit; pop_suit(randnode)];
    nodedcsuit=[nodedcsuit; dcoast_suit(randnode)];
    nodedbsuit=[nodedbsuit; dbrdr_suit(randnode)];
    nodeslpsuit=[nodeslpsuit; slp_suit(randnode)];
    nodemktsuit=[nodemktsuit; mktacc_suit(randnode)];
    nodelusuit=[nodelusuit; lu_suit(randnode)];
    nodelsuit=[nodelsuit; LANDSUIT(randnode)];
    nodedto=[nodedto; zeros(length(randnode),1)];
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
        cpcty=2*stock_0*ones(length(snode),1);
        EdgeTable=table([EdgeTable.EndNodes; snode tnode],[EdgeTable.Weight; ...
            weights],[EdgeTable.Flows; flows],[EdgeTable.Capacity; cpcty],...
            'VariableNames',{'EndNodes' 'Weight' 'Flows' 'Capacity'});
    end 
end
NodeTable=table(nodeid',noderow,nodecol,nodelat,nodelon,nodecode,nodestck,...
    nodecptl,nodetcov,nodepopsuit,nodedcsuit,nodedbsuit,nodeslpsuit,...
    nodemktsuit,nodelusuit,nodelsuit,nodedto,'VariableNames',{'ID','Row','Col','Lat',...
    'Lon','DeptCode','Stock','Capital','TreeCover','PopSuit',...
    'DistCoastSuit','DistBorderSuit','SlopeSuit','MktAccSuit','LandUseSuit',...
    'LandSuit','DTO'});
nnodes=height(NodeTable);
ADJ=zeros(nnodes);      % adjacency matrix for trafficking network
TRRTY=zeros(nnodes);    % control of nodes by each DTO
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

NEIHOOD=cell(nnodes,2);

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
BRIBE=zeros(nnodes,TMAX);       % annual bribe payments made at each node to maintain control
MARGIN=zeros(nnodes,TMAX);      % gross profit per node after purchasing, trafficking, and selling
RENTCAP=zeros(nnodes,TMAX);     % portion of MARGIN retained at node as profit
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
% remotefac=[0; 1-NodeTable.TreeCover(2:nnodes-1)./100; 0];                             
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
subneihood=zeros(size(LANDSUIT));
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
        idist_ground=(DIST(j,ireceiver)  >0 & DIST(j,ireceiver) <= 500);
        idist_air=(DIST(j,ireceiver) > 500);
%         CTRANS(j,ireceiver)=ctrans_inland.*RMTFAC(j,ireceiver).*...
%             DIST(j,ireceiver);
        CTRANS(j,ireceiver(idist_ground))=ctrans_inland.*DIST(j,ireceiver(idist_ground));
        CTRANS(j,ireceiver(idist_air))=ctrans_coast.*DIST(j,ireceiver(idist_air));
    end
    
%     % Create gridded distance from node for calculating land use
%     % neighborhood
%     if j > 1 && j < nnodes
%         idptcode=find(dptgrid == NodeTable.DeptCode(j));
%         [dptrows,dptcols]=ind2sub(size(dptgrid),idptcode);
%         hdir=repmat([(NodeTable.Col(j)-1):-1:min(dptcols) 0 1:(max(dptcols)-...
%             NodeTable.Col(j))],max(dptrows)-min(dptrows)+1,1);
% %         hdir(luint == 0 | luint == 8)=NaN;
%         vdir=repmat([((NodeTable.Row(j)-1):-1:min(dptrows))'; 0; (1:(max(dptrows)-...
%             NodeTable.Row(j)))'],1,max(dptcols)-min(dptcols)+1);
% %         vdir(luint == 0 | luint == 8)=NaN;
%         cmpdir=zeros(size(hdir,1),size(hdir,2),2);
%         cmpdir(:,:,1)=hdir;
%         cmpdir(:,:,2)=vdir;
%         dirholder=max(hdir,vdir);
%         subneihood(min(dptrows):max(dptrows),min(dptcols):max(dptcols))=dirholder;
%         NEIHOOD(j,1)=mat2cell(idptcode,length(idptcode),1);
%         NEIHOOD(j,2)=mat2cell(subneihood(idptcode),length(idptcode),1);
%     end
end

%%% Assign nodes to initials DTOs
for nn=2:nnodes-1
    westdir=NodeTable.Col(nn)-find(isnan(dcoast(NodeTable.Row(nn),1:...
        NodeTable.Col(nn)-1))==1,1,'last');
    eastdir=find(isnan(dcoast(NodeTable.Row(nn),NodeTable.Col(nn)+1:...
        size(LANDSUIT,2)))==1,1,'first');
    northdir=NodeTable.Row(nn)-find(isnan(dcoast(1:NodeTable.Row(nn)-1,...
        NodeTable.Col(nn)))==1,1,'last');
    southdir=find(isnan(dcoast(NodeTable.Row(nn)+1:size(LANDSUIT,1),...
        NodeTable.Col(nn)))==1,1,'first');
    [mindist,imindist]=min([westdir eastdir northdir southdir]);
    if imindist == 1 || imindist == 4
        NodeTable.DTO(nn)=1;
    else
        NodeTable.DTO(nn)=2;
    end
end

%%% Initialize Interdiction agent
facmat=LATFAC;
facmat(:,:,2)=COASTFAC;
facmat(:,:,3)=RMTFAC;
% SLPROB(:,:,TSTART)=max(min((COASTFAC./max(max(COASTFAC))).*...
%     (RMTFAC./max(max(RMTFAC))).*(DIST./max(max(DIST))),1),0);
SLPROB(:,:,TSTART)=max(min(max(facmat,[],3)+DIST./max(max(DIST)),1),0);   % dynamic probability of seisure and loss at edges
slmin=SLPROB(:,:,TSTART);
SLPROB(:,:,TSTART+1)=SLPROB(:,:,TSTART);
INTRDPROB(:,TSTART+1)=slprob_0*ones(nnodes,1); % dynamic probability of interdiction at nodes

%%% Initialize Node agents
STOCK(:,TSTART)=NodeTable.Stock(:);
TOTCPTL(:,TSTART)=NodeTable.Capital(:);
PRICE(:,TSTART+1)=PRICE(:,TSTART);
% slcpcty_0=ceil(length(find(SLPROB(:,:,TSTART) == 1))/6);    % assumed number of S&L events that can be carried out per time step
% slcpcty_max=ceil(length(find(SLPROB(:,:,TSTART) == 1)));
slcpcty_0=sl_min(erun);
slcpcty_max=sl_max(erun);
slcpcty(TSTART+1)=slcpcty_0;

% subjective risk perception with time distortion
twght=timewght_0*ones(nnodes,1);    % time weighting for dynamic, subjective perceived risk of interdiction event
    
%%% Set-up trafficking netowrk benefit-cost logic  %%%%%%%%%%%%
ltcoeff=ones(nnodes,1);
% routepref(:,:,TSTART+1)=ADJ;
margprofit=ADDVAL-CTRANS;
for nd=1:ndto
    idto=(NodeTable.DTO == nd);
    routepref(1,idto,TSTART+1)=(margprofit(1,idto)==max(margprofit(1,idto)));
end
routepref(:,nnodes,TSTART+1)=1;
totslrisk(TSTART+1)=1;

OWN=zeros(size(LANDSUIT));  % node agent land ownership
IOWN=cell(nnodes,TMAX);     % dynamic list of owned parcels


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
    %%% Specified number of events, selection of highest probability nodes (in addition to p=1) 
%     newslevents=ceil(slcpcty(t)*rand(1));
    subslevent=slevent(:,:,t);
    subslprob=reshape(SLPROB(:,:,t),nnodes*nnodes,1);
    irevisit=find(slsuccess(:,:,t-1) > 0);
    irevisit=irevisit(randperm(length(irevisit),min(slcpcty(t),...
        length(irevisit))));
    %     subslprob=subslprob(subslprob~=1);
    if isempty(find(irevisit,1)) == 1 && ...
            sum(sum(sum(slsuccess(:,:,TSTART+1:t-1))),3) > 0
        slquant=quantile(subslprob(subslprob~=0),[0.5 0.75 0.9]);
        sleligible=find(subslprob > slquant(1));
    else
        slquant=quantile(subslprob(subslprob~=0),[0.5 0.75 0.9]);
        sleligible=find(subslprob > slquant(2));
    end
    sleligible=sleligible(~ismember(sleligible,irevisit));
    islevent=[irevisit; sleligible(randperm(length(sleligible),...
        min(max(slcpcty(t)-length(irevisit),0),length(sleligible))))];
    
%     [subslprobsort,isubslprobsort]=sort(subslprob,'descend');
%     islevent=isubslprobsort(1:length(find(subslprob==1))+slcpcty(t));
%     if length(find(subslprob==1)) >= slcpcty(t)
%         islevent=isubslprobsort(randperm(length(find(subslprob==1)),...
%             slcpcty(t)));
%     else
%         islevent=isubslprobsort(1:slcpcty(t));
%     end
    subslevent(islevent)=1;
    slevent(:,:,t)=subslevent;
%     slevent(:,:,t)=(SLPROB(:,:,t) == 1);
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
%             inei=zeros([],1);
            if n == 1
%                 if t >= 166
%                     keyboard
%                 end
                inei=find(ADJ(n,:) == 1 & routepref(n,:,t) > 0);
                for nd=1:length(unique(NodeTable.DTO(2:nnodes-1)))
                    if isempty(find(NodeTable.DTO(inei) == nd,1)) == 1
                        idtombr=(NodeTable.DTO == nd);
                        subinei=find(ADJ(n,idtombr) == 1 & routepref(n,idtombr,t) > 0);
                        if isempty(find(subinei,1)) == 1
                            subinei=find(ADJ(n,idtombr) == 1 & routepref(n,idtombr,t) == ...
                                max(routepref(n,idtombr,t)));
                        end
                        inei=[inei subinei];
                    end
                end
            else
%                 inei=find(ADJ(n,:) == 1 & routepref(n,:,t) > 0 & ...
%                     ((NodeTable.DTO' == NodeTable.DTO(n)) == 1 | ...
%                     (NodeTable.DTO' == 0) == 1));
                inei=find(ADJ(n,:) == 1 & routepref(n,:,t) > 0);
                inei=inei(ismember(inei,[find(NodeTable.DTO == ...
                    NodeTable.DTO(n)); nnodes]));
                if isempty(find(inei,1)) == 1
                    %                     inei=find(ADJ(n,:) == 1 & routepref(n,:,t) == ...
                    %                         max(routepref(n,NodeTable.DTO' == NodeTable.DTO(n),t)));
                    inei=find(ADJ(n,:) == 1 & routepref(n,:,t) == ...
                        max(routepref(n,:,t)));
                    inei=inei(ismember(inei,[find(NodeTable.DTO == ...
                        NodeTable.DTO(n)); nnodes]));
                end
            end
            %%% Procedure for selecting routes based on expected profit %%%
            c_trans=CTRANS(n,inei);
            p_sl=SLRISK(n,inei);
            %          p_int=INTRISK(n,inei); %currently not used
            y_node=ADDVAL(n,inei);
%             q_node=min(floor(STOCK(n,t)./length(inei)),CPCTY(n,inei));
            q_node=min(STOCK(n,t)./length(inei),CPCTY(n,inei));
            lccf=ltcoeff(n);
            totstock=STOCK(n,t);
            totcpcty=CPCTY(n,inei);
            tslrisk=totslrisk(t);
            rtpref=routepref(n,inei,t);
            dtonei=NodeTable.DTO(inei);
            cutflag=dtocutflag(unique(dtonei(dtonei~=0)));
            
            % ********** Need to include check on each DTO's performance so
            % that after some number of time steps (3?) without receiving
            % shipment at end node from a DTO's network, that DTO can be
            % cut out **************************************************
            
            %          [neipick,neivalue]=calc_neival(c_trans,p_sl,y_node,q_node,lccf,...
            %              totstock,totcpcty,tslrisk);
            [neipick,neivalue]=calc_neival(c_trans,p_sl,y_node,q_node,lccf,...
                rtpref,tslrisk,dtonei,cutflag);
            inei=inei(neipick);
            activeroute(n,t)=mat2cell(inei',length(inei),1);
            
            neiset=unique(NodeTable.DTO(inei));
            if length(neiset(neiset~=0)) > 1 && n ~= 1
                display('Check DTO membership of neighbors')
                keyboard
            end
            
            % weight according to salience value fuction
%             if isempty(find(neivalue > 0,1)) == 1
%                 WGHT(n,inei)=1;
%             else
            if isempty(find(neivalue > 0,1)) == 0
                WGHT(n,inei)=max(neivalue,0)./sum(max(neivalue,0)); 
            end
%             if t >= 115 && n == 1
%                 keyboard
%             end
            %%% !!! Put checks in to make sure buying node has enough capital
            FLOW(n,inei,t)=min(WGHT(n,inei).*STOCK(n,t),CPCTY(n,inei));
            %          FLOW(n,inei,t)=min(WGHT(n,inei).*(STOCK(n,t)/length(inei)),CPCTY(n,inei));
            
            % Check for S%L event
            if isempty(find(ismember(find(slevent(n,:,t)),inei),1)) == 0
                isl=(slevent(n,inei,t)==1);
                slsuccess(n,inei(isl),t)=FLOW(n,inei(isl),t);
                if slsuccess(n,inei(isl),t) == 0
                    slevent(n,inei(isl),t)=0;
                end
                OUTFLOW(n,t)=sum(FLOW(n,inei,t));
                STOCK(n,t)=STOCK(n,t)-OUTFLOW(n,t);
                nodecosts=sum(FLOW(n,inei,t).*CTRANS(n,inei));
%                 TOTCPTL(n,t)=TOTCPTL(n,t)-sum(FLOW(n,inei,t).*CTRANS(n,inei));
                FLOW(n,inei(isl),t)=0;     % remove from trafficking route due to S&L event
                STOCK(inei,t)=STOCK(inei,t)+FLOW(n,inei,t)';
%                 TOTCPTL(n,t)=TOTCPTL(n,t)+sum(FLOW(n,inei,t).*ADDVAL(n,inei));
                noderevenue=sum(FLOW(n,inei,t).*PRICE(inei,t)');
%                 TOTCPTL(n,t)=TOTCPTL(n,t)+sum(FLOW(n,inei,t).*PRICE(inei,t)');
                TOTCPTL(inei,t)=TOTCPTL(inei,t)-(FLOW(n,inei,t)'.*PRICE(inei,t));
                ICPTL(n,t)=rentcap*sum(FLOW(n,inei).*ADDVAL(n,inei));
                MARGIN(n,t)=noderevenue-nodecosts+min(TOTCPTL(n,t),0);
                if n > 1
                    BRIBE(n,t)=max(bribepct*MARGIN(n,t),0);
                    if MARGIN(n,t) > 0
                        RENTCAP(n,t)=rentcap*(MARGIN(n,t)-BRIBE(n,t));
                    else
                        RENTCAP(n,t)=MARGIN(n,t);
                    end
                    TOTCPTL(n,t)=max(TOTCPTL(n,t),0)+RENTCAP(n,t);  % losses on top of debt capture by MARGIN
                    DTOBDGT(NodeTable.DTO(n),t)=DTOBDGT(NodeTable.DTO(n),t)+...
                        (1-rentcap)*max(MARGIN(n,t)-BRIBE(n,t),0);
                else
                    RENTCAP(n,t)=MARGIN(n,t);
                    TOTCPTL(n,t)=TOTCPTL(n,t)+RENTCAP(n,t);
                end
            else
                OUTFLOW(n,t)=sum(FLOW(n,inei,t));
                STOCK(n,t)=STOCK(n,t)-OUTFLOW(n,t);
                STOCK(inei,t)=STOCK(inei,t)+FLOW(n,inei,t)';
                nodecosts=sum(FLOW(n,inei,t).*CTRANS(n,inei));
                noderevenue=sum(FLOW(n,inei,t).*PRICE(inei,t)');
                TOTCPTL(inei,t)=TOTCPTL(inei,t)-(FLOW(n,inei,t)'.*PRICE(inei,t));
                ICPTL(n,t)=rentcap*sum(FLOW(n,inei).*ADDVAL(n,inei));
                MARGIN(n,t)=noderevenue-nodecosts+min(TOTCPTL(n,t),0);
                if n > 1
                    BRIBE(n,t)=max(bribepct*MARGIN(n,t),0);
                    if MARGIN(n,t) > 0
                        RENTCAP(n,t)=rentcap*(MARGIN(n,t)-BRIBE(n,t));
                    else
                        RENTCAP(n,t)=MARGIN(n,t);
                    end
                    TOTCPTL(n,t)=max(TOTCPTL(n,t),0)+RENTCAP(n,t);
                    DTOBDGT(NodeTable.DTO(n),t)=DTOBDGT(NodeTable.DTO(n),t)+...
                        (1-rentcap)*max(MARGIN(n,t)-BRIBE(n,t),0);
                else
                    RENTCAP(n,t)=MARGIN(n,t);
                    TOTCPTL(n,t)=TOTCPTL(n,t)+RENTCAP(n,t);
                end
            end
            %%%% Update perceived risk in response to S&L and Interdiction events
            timeweight=twght(n);
            % identify neighbors in network (without network toolbox)
            %       bcknei=EdgeTable.EndNodes(EdgeTable.EndNodes(:,2) == n,1)';
            fwdnei=inei;
            if t == TSTART+1
                sloccur=slevent(n,fwdnei,TSTART+1:t);
            elseif t > TSTART+1 && length(fwdnei) == 1
                sloccur=squeeze(slevent(n,fwdnei,max(TSTART+1,t-12):t));
            else
                sloccur=squeeze(slevent(n,fwdnei,max(TSTART+1,t-12):t))';
            end
            intrdoccur=intrdevent(fwdnei,max(TSTART+1,t-12):t);
            [sl_risk,intrd_risk,slevnt,intrdevnt,tmevnt]=calc_intrisk(sloccur,...
                intrdoccur,t,TSTART,alpharisk,betarisk,timeweight);
            %       SLRISK(n,[bcknei fwdnei])=sl_risk;
            SLRISK(n,fwdnei)=sl_risk;
            
            if isempty(find(sl_risk,1)) == 0
                avgslrisk(n,t)=mat2cell(SLRISK(n,activeroute{n,t}),1,...
                    length(activeroute{n,t}));
            end
            %       INTRISK(n,t+1)=mean(intrd_risk);  %node-specific risk is the average of neighbor risks
            
            %!!!!!!!!!!!
            %          ICPTL(n,t)=ICPTL(n,t)-OUTFLOW(n,t)*VALUE(  %account for value retained at node
            
            NodeTable.Stock(:)=STOCK(:,t);
            NodeTable.Capital(:)=TOTCPTL(:,t);
        end
        
        %%%%%%%%% Loss of node control %%%%%%%%%%%
%         if n ~= 1 && n ~= nnodes && isempty(find(OUTFLOW(n,1:t),1)) == 0 && ...
%                 t > find(OUTFLOW(n,1:t),1,'first')+bribethresh && ...
%                 length(find(BRIBE(n,find(BRIBE(n,:),1,'last'):t) == 0)) > bribethresh
%             NodeTable.DTO(n)=0;
%             WGHT(n,:)=0;
%             routepref(n,:,t+1)=0;
%         end
    
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
%     slcpcty(t+1)=min(max(ceil((1-delta_sl)*slcpcty(t)+delta_sl*...
%         (sum(sum(slsuccess(:,:,t)))-sum(sum(slsuccess(:,:,t-1))))),...
%         slcpcty_0),slcpcty_max);
    slcpcty(t+1)=min(max(ceil((1-delta_sl)*slcpcty(t)+delta_sl*...
        (sum(sum(slsuccess(:,:,t)))-sum(sum(mean(slsuccess(:,:,TSTART+1:t-1),3))))),...
        slcpcty_0),slcpcty_max);
%     slcpcty(t+1)=max(slcpcty(t)+ceil(delta_sl*(sum(sum(slsuccess(:,:,t)))-...
%         sum(sum(mean(slsuccess(:,:,TSTART+1:t-1),3))))),slcpcty_0);
%     INTRDPROB(:,t+1)=INTRDPROB(:,t);
    
    % Reinforcement learning for successful routes
    iactivenode=find(OUTFLOW(2:nnodes-1,t) > 0)+1;
    avgflow=STOCK(iendnode,t)/length(iactivenode);

    activenodes=unique(cat(1,activeroute{:,t}));
    actedge=activeroute(:,t);
    
    %%%%%%%%% Route Optimization %%%%%%%%%%%
    for dt=1:ndto
        idto=find(NodeTable.DTO == dt);
%         DTOBDGT(dt,t)=sum(TOTCPTL(idto,t)); %total DTO funds for expansion/viability
%         if t > 3 && isempty(find(DTOBDGT(dt,t-3:t) > 0,1)) == 1
%             dtocutflag(dt)=1;
%             display('dto cut')
%             keyboard
%         end
        dtorefvec=[1; idto; nnodes];
        subnnodes=length(idto);
        subroutepref=routepref(dtorefvec,dtorefvec,t);
        subactivenodes=activenodes(ismember(activenodes,idto));
        subactedges=actedge(dtorefvec);
        subflow=FLOW(dtorefvec,dtorefvec,t);
        dtoslsuc=slsuccess(dtorefvec,dtorefvec,t);
        
        sendedge=ismember(EdgeTable.EndNodes(:,1),dtorefvec);
%         recedge=ismember(EdgeTable.EndNodes(:,2),dtorefvec);
%         dtoEdgeTable=EdgeTable(sendedge == 1 & recedge == 1,:);
        dtoEdgeTable=EdgeTable(sendedge,:);
        dtoEdgeTable=dtoEdgeTable(ismember(dtoEdgeTable.EndNodes(:,2),dtorefvec),:);
        dtoSLRISK=SLRISK(dtorefvec,dtorefvec);
        dtoADDVAL=ADDVAL(dtorefvec,dtorefvec);
        dtoCTRANS=CTRANS(dtorefvec,dtorefvec);
        %%% calculate losses from S&L events
        % volume-based - does not matter where in supply chain
        %     supplyfit=STOCK(iendnode,t)/stock_0;
        %     losstolval=losstol*stock_0;
        
        % value-based - price varies with location in supply chain
        ipossl=find(dtoslsuc > 0);
        [nrow,ncol]=ind2sub(size(dtoslsuc),ipossl);
        %     losstolval=losstol*stock_0*PRICE(nnodes,t); %value-based loss threshold
        
        supplyfit=sum(dtoslsuc(ipossl).*PRICE(dtorefvec(ncol),t));  %value-based loss calc
        losstolval=(1-losstol)*sum(subflow(1,:)+dtoslsuc(1,:))*PRICE(nnodes,t);
        
        if isempty(find(supplyfit ~= 0,1)) == 1 && isempty(find(losstolval ~= 0,1)) == 1
            supplyfit=1;
        end
        %call top-down route optimization
        %     newroutepref=optimizeroute(nnodes,subflow,supplyfit,activenodes,...
        %         subroutepref,EdgeTable,SLRISK,ADDVAL,CTRANS,losstolval);
        newroutepref=optimizeroute_multidto(dtorefvec,subflow,supplyfit,subactivenodes,...
            subroutepref,dtoEdgeTable,dtoSLRISK,dtoADDVAL,dtoCTRANS,losstolval);
        
        routepref(dtorefvec,dtorefvec,t+1)=newroutepref;
        if isempty(find(routepref(1,dtorefvec(1:length(dtorefvec)-1),t+1),1)) ==1
            display('check route optimization for each dto')
%             keyboard
        end
    end
%     routepref(:,:,t+1)=routepref(:,:,t);
    PRICE(:,t+1)=PRICE(:,t);
    STOCK(1,t+1)=stock_0;    %additional production to enter network next time step
    STOCK(nnodes,t+1)=0;    %remove stock at end node for next time step
    NodeTable.Stock(1)=stock_0;
    NodeTable.Stock(nnodes)=0;
end
% toc     % stop run timer
%%
% % %%% Visualization %%%
% 
% %%% Trafficking movie
% writerObj = VideoWriter('trafficking_risk_dtos_short_v1.mp4','MPEG-4');
% writerObj.FrameRate=3;
% open(writerObj);
% 
% h1=figure;
% set(h1,'Color','white','Visible','off')
% 
% for tt=TSTART+1:36
%     geoshow(CAadm0,'FaceColor',[1 1 1])
%     hold on
%     %     plot(nodelon(1),nodelat(1),'r.','MarkerSize',ceil(MOV(1,1,tt)/1000))
%     plot(nodelon(1),nodelat(1),'m.','MarkerSize',stock_0)
% %     fedge=find(FLOW(1,:,tt) > 0);
%     fedge=activeroute{1,tt};
%     %
%     islevent=find(slevent(1,fedge,tt) == 1);
%     inoslevent=find(slevent(1,fedge,tt) == 0);
%     if isempty(islevent) == 1
%         for g=1:length(fedge)
%             if NodeTable.DTO(fedge(g)) == 1
%                 plot([nodelon(1); nodelon(fedge(g))],[nodelat(1); nodelat(fedge(g))],'-r')
%             elseif NodeTable.DTO(fedge(g)) == 2
%                 plot([nodelon(1); nodelon(fedge(g))],[nodelat(1); nodelat(fedge(g))],'-b')
%             else
%                 plot([nodelon(1); nodelon(fedge(g))],[nodelat(1); nodelat(fedge(g))],'--k')
%             end
%         end
%     else
%         for g=1:length(fedge(islevent))
%             if NodeTable.DTO(fedge(islevent(g))) == 1
%                 plot([nodelon(1); nodelon(fedge(islevent(g)))],[nodelat(1); ...
%                     nodelat(fedge(islevent(g)))],'-r')
%                 plot(nodelon(fedge(islevent(g))),nodelat(fedge(islevent(g))),...
%                     'rx','MarkerSize',ceil(slsuccess(1,fedge(islevent(g)),tt)./10))
%             elseif NodeTable.DTO(fedge(islevent(g))) == 2
%                 plot([nodelon(1); nodelon(fedge(islevent(g)))],[nodelat(1); ...
%                     nodelat(fedge(islevent(g)))],'-b')
%                 plot(nodelon(fedge(islevent(g))),nodelat(fedge(islevent(g))),...
%                     'bx','MarkerSize',ceil(slsuccess(1,fedge(islevent(g)),tt)./10))
%             else
%                 plot([nodelon(1); nodelon(fedge(islevent(g)))],[nodelat(1); ...
%                     nodelat(fedge(islevent(g)))],'--k')
%                 plot(nodelon(fedge(islevent(g))),nodelat(fedge(islevent(g))),...
%                     'kx','MarkerSize',ceil(slsuccess(1,fedge(islevent(g)),tt)./10))
%             end
%         end
%         for k=1:length(fedge(inoslevent))
%             if NodeTable.DTO(fedge(inoslevent(k))) == 1
%                 plot([nodelon(1); nodelon(fedge(inoslevent(k)))],[nodelat(1); ...]
%                     nodelat(fedge(inoslevent(k)))],'-r')
%             elseif NodeTable.DTO(fedge(inoslevent(k))) == 2
%                 plot([nodelon(1); nodelon(fedge(inoslevent(k)))],[nodelat(1); ...]
%                     nodelat(fedge(inoslevent(k)))],'-b')
%             else
%                 plot([nodelon(1); nodelon(fedge(inoslevent(k)))],[nodelat(1); ...]
%                     nodelat(fedge(inoslevent(k)))],'--k')
%             end
%         end
%     end
%     plot(nodelon(2:nnodes),nodelat(2:nnodes),'k.','MarkerSize',3)
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
%                     plot(nodelon(nn),nodelat(nn),'m.','MarkerSize',ceil(MOV(nn,mm,tt)))
%                 else
%                     plot(nodelon(nn),nodelat(nn),'k.','MarkerSize',3)
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
%                     plot(nodelon(nn),nodelat(nn),'m.','MarkerSize',ceil(MOV(nn,mm,tt)))
%                 else
%                     plot(nodelon(nn),nodelat(nn),'k.','MarkerSize',3)
%                 end
%             end
%             fedge=activeroute{mm,tt};
%             islevent=find(slevent(mm,fedge,tt) == 1);
%             inoslevent=find(slevent(mm,fedge,tt) == 0);
%             if isempty(islevent) == 1
%                 for g=1:length(fedge)
%                     if NodeTable.DTO(mm) == 1
%                         plot([nodelon(mm); nodelon(fedge(g))],[nodelat(mm); nodelat(fedge(g))],'-r')
%                     elseif NodeTable.DTO(mm) == 2
%                         plot([nodelon(mm); nodelon(fedge(g))],[nodelat(mm); nodelat(fedge(g))],'-b')
%                     else
%                         plot([nodelon(mm); nodelon(fedge(g))],[nodelat(mm); nodelat(fedge(g))],'--k')
%                     end
%                 end
%             else
%                 for g=1:length(fedge(islevent))
%                     if NodeTable.DTO(mm) == 1
%                         plot([nodelon(mm); nodelon(fedge(islevent(g)))],[nodelat(mm); ...
%                             nodelat(fedge(islevent(g)))],'-r')
%                         plot(nodelon(fedge(islevent(g))),nodelat(fedge(islevent(g))),...
%                             'rx','MarkerSize',ceil(slsuccess(mm,fedge(islevent(g)),tt)))
%                     elseif NodeTable.DTO(mm) == 2
%                         plot([nodelon(mm); nodelon(fedge(islevent(g)))],[nodelat(mm); ...
%                             nodelat(fedge(islevent(g)))],'-b')
%                         plot(nodelon(fedge(islevent(g))),nodelat(fedge(islevent(g))),...
%                             'bx','MarkerSize',ceil(slsuccess(mm,fedge(islevent(g)),tt)))
%                     else
%                         plot([nodelon(mm); nodelon(fedge(islevent(g)))],[nodelat(mm); ...
%                             nodelat(fedge(islevent(g)))],'--k')
%                         plot(nodelon(fedge(islevent(g))),nodelat(fedge(islevent(g))),...
%                             'kx','MarkerSize',ceil(slsuccess(mm,fedge(islevent(g)),tt)))
%                     end
%                 end
%                 for k=1:length(fedge(inoslevent))
%                     if NodeTable.DTO(mm) == 1
%                         plot([nodelon(mm); nodelon(fedge(inoslevent(k)))],...
%                             [nodelat(mm); nodelat(fedge(inoslevent(k)))],'-r')
%                     elseif NodeTable.DTO(mm) == 2
%                         plot([nodelon(mm); nodelon(fedge(inoslevent(k)))],...
%                             [nodelat(mm); nodelat(fedge(inoslevent(k)))],'-b')
%                     else
%                         plot([nodelon(mm); nodelon(fedge(inoslevent(k)))],...
%                             [nodelat(mm); nodelat(fedge(inoslevent(k)))],'--k')
%                     end
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

% %%%%%% Plot some node attribute
% % Map of node locations and DTO membership
% geoshow(CAadm0,'FaceColor',[1 1 1])
% hold on
% for nn=1:nnodes
%     if NodeTable.DTO(nn) == 1
%         plot(nodelon(nn),nodelat(nn),'.r','MarkerSize',3)
%     elseif NodeTable.DTO(nn) == 2
%         plot(nodelon(nn),nodelat(nn),'.b','MarkerSize',3)
%     else
%         plot(nodelon(nn),nodelat(nn),'.k','MarkerSize',3)
%     end
% end

% figure
% geoshow(CAadm0,'FaceColor',[1 1 1])
% hold on
% for i=1:length(nodelat)
%         plot(nodelon(i),nodelat(i),'r.','MarkerSize',round(NodeTable.CoastDist(i)/10))
% end
%%% Plot time series of flows and S&L
% plot(1:t,STOCK(1,1:t),'-b')
% h2=figure;
% plot(1:t,STOCK(nnodes,1:t),'--k')
% hold on
% sltot=sum(sum(slsuccess(:,:,1:t),1));
% plot(1:t,reshape(sltot,1,t),'-r')
% xlim([0 t])
% ylabel('Cocaine Volume kg/month')
% xlabel('Month')
% legend('Consumer','S&L','Orientation','horizontal','Location','southoutside')
% 
% %%%%%% Diagnostics %%%%%%%
% slrecord=sum(slsuccess(:,:,1:t),3);
% [slr,slc]=ind2sub(size(slrecord),find(slrecord > 0));
% 
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
% h2_1=figure;
% set(h2_1,'Color','white')
% [hAx,hl1,hl2]=plotyy(1:TMAX,nactnodes,1:TMAX,slperevent);
% % plot(1:TMAX,nactnodes,'-b')
% ylabel(hAx(1),'Number of Routes')
% % hold on
% % yyaxis right
% % plot(1:TMAX,slpervent,'--r')
% ylabel(hAx(2),'Average S&L Volume (kg)')
% xlim([0 TMAX])
% xlabel('Month')
% legend('Active Routes','S&L Volume','Orientation','horizontal','Location','southoutside')

sltot=sum(sum(slevent(:,:,1:TMAX),1));
% plot(1:TMAX,reshape(sltot,1,TMAX),'--r')
% xlim([0 TMAX])
% ylabel('Number of Routes & SL Events')
% xlabel('Month')
% legend('Active Routes','S&L Events','Orientation','horizontal','Location','southoutside')
% saveas(h2_1,'Nodes_vs_SL_null.png')

% %%% Command to use
% % digraph, maxflow, nearest