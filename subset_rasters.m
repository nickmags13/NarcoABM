%%%%%%%%% Subset geotiffs for alignment %%%%%%%%%%%%%%%%%%%%

[cagrid_cntry,Rcagrid_cntry]=geotiffread('X:\CentralAmericaData\Model_inputs\ca_cntry_250.tif');
cntrynodataval=255;

[treecov,Rtreecov]=geotiffread('X:\CentralAmericaData\Model_inputs\treecov_250.tif');
treenodataval=255;

% Distance to coast and country borders
[dcoast,Rdcoast]=geotiffread('X:\CentralAmericaData\Model_inputs\coastdist_250.tif');
dcoastnodataval=-9999;

[dbrdr,Rdbrdr]=geotiffread('X:\CentralAmericaData\Model_inputs\borderdist_250.tif');
brdrnodataval=-9999;

% Population density as a proxy for remoteness
[popden,Rpopden]=geotiffread('X:\CentralAmericaData\Model_inputs\lspop2001_250.tif');
popnodataval=-1;

% Topography
[slope,Rslope]=geotiffread('X:\CentralAmericaData\Model_inputs\ca_slope_250m.tif');

% Market Access
[mktacc,Rmktacc]=geotiffread('X:\CentralAmericaData\Model_inputs\mkt_acc_250.tif');
manodataval=-9999;

% Initial land use
[luint,Rluint]=geotiffread('X:\CentralAmericaData\Model_inputs\luca_2001.tif');
lunodataval=-1;

% Yields
[mazyld,Rmazyld]=geotiffread('X:\CentralAmericaData\Model_inputs\ca_mazyld_250.tif');
maznodataval=-9999;

[plmyld,Rplmyld]=geotiffread('X:\CentralAmericaData\Model_inputs\ca_plmyld_250.tif');
plmnodataval=-9999;

[ctlden,Rctlden]=geotiffread('X:\CentralAmericaData\Model_inputs\ca_ctlden_250.tif');
ctlnodataval=-9999;

%%% Reconcile to slope because it is the smallest raster
firstrowlat=Rslope.LatitudeLimits(2);
lastrowlat=Rslope.LatitudeLimits(1);
firstcollon=Rslope.LongitudeLimits(1);
lastcollon=Rslope.LongitudeLimits(2);
cellext=Rslope.CellExtentInLatitude;

% Distance to coast
topdiff=round((Rdcoast.LatitudeLimits(2)-firstrowlat)/cellext);
bottomdiff=round((lastrowlat-Rdcoast.LatitudeLimits(1))/cellext);
leftdiff=round((Rdcoast.LongitudeLimits(1)-firstcollon)/cellext);
rghtdiff=round((lastcollon-Rdcoast.LongitudeLimits(2))/cellext);
dcoastdim=size(dcoast);
dcoast=dcoast(1+topdiff:dcoastdim(1)-bottomdiff,...
    1+leftdiff:dcoastdim(2)+rghtdiff);
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\dcoast_clp.tif',dcoast,Rslope);

% Distance to border
topdiff=round((Rdbrdr.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rdbrdr.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rdbrdr.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rdbrdr.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
dbrdrdim=size(dbrdr);
dbrdr=dbrdr(1+max(topdiff,0):dbrdrdim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):dbrdrdim(2)+rghtdiff);
dbrdr=[zeros(length(dbrdr(:,1)),lftadd) dbrdr];
dbrdr=[zeros(abs(topadd),length(dbrdr(1,:))); dbrdr];
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\dbrdr_clp.tif',dbrdr,Rslope);

% Country grid
topdiff=round((Rcagrid_cntry.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rcagrid_cntry.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rcagrid_cntry.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rcagrid_cntry.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
cagrid_cntrydim=size(cagrid_cntry);
cagrid_cntry=cagrid_cntry(1+max(topdiff,0):cagrid_cntrydim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):cagrid_cntrydim(2)+rghtdiff);
cagrid_cntry=[cntrynodataval*ones(length(cagrid_cntry(:,1)),lftadd) cagrid_cntry];
cagrid_cntry=[cntrynodataval*ones(abs(topadd),length(cagrid_cntry(1,:))); cagrid_cntry];
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\ca_cntry_clp.tif',cagrid_cntry,Rslope);

% Market Access
topdiff=round((Rmktacc.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rmktacc.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rmktacc.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rmktacc.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
mktaccdim=size(mktacc);
mktacc=mktacc(1+max(topdiff,0):mktaccdim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):mktaccdim(2)+rghtdiff);
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\mktacc_clp.tif',mktacc,Rslope);

% Initial land use
topdiff=round((Rluint.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rluint.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rluint.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rluint.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
luintdim=size(luint);
luint=luint(1+max(topdiff,0):luintdim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):luintdim(2)+rghtdiff);
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\luint_clp.tif',luint,Rslope);

% population density
topdiff=round((Rpopden.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rpopden.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rpopden.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rpopden.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
popdendim=size(popden);
popden=popden(1+max(topdiff,0):popdendim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):popdendim(2)+rghtdiff);
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\popden_clp.tif',popden,Rslope);

% Tree cover
topdiff=round((Rtreecov.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rtreecov.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rtreecov.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rtreecov.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
treecovdim=size(treecov);
treecov=treecov(1+max(topdiff,0):treecovdim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):treecovdim(2)+rghtdiff);
treecov=[treenodataval*ones(length(treecov(:,1)),lftadd) treecov];
treecov=[treenodataval*ones(abs(topadd),length(treecov(1,:))); treecov];
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\treecov_clp.tif',treecov,Rslope);

% Maize yield
topdiff=round((Rmazyld.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rmazyld.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rmazyld.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rmazyld.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
mazylddim=size(mazyld);
mazyld=mazyld(1+max(topdiff,0):mazylddim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):mazylddim(2)+rghtdiff);
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\mazyld_clp.tif',mazyld,Rslope);

% Oil Palm yield
topdiff=round((Rplmyld.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rplmyld.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rplmyld.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rplmyld.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
plmylddim=size(plmyld);
plmyld=plmyld(1+max(topdiff,0):plmylddim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):plmylddim(2)+rghtdiff);
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\plmyld_clp.tif',plmyld,Rslope);

% Cattle density
topdiff=round((Rctlden.LatitudeLimits(2)-firstrowlat)/cellext);
topadd=min(topdiff,0);
bottomdiff=round((lastrowlat-Rctlden.LatitudeLimits(1))/cellext);
btmadd=min(bottomdiff,0);
leftdiff=round((Rctlden.LongitudeLimits(1)-firstcollon)/cellext);
lftadd=max(leftdiff,0);
rghtdiff=round((lastcollon-Rctlden.LongitudeLimits(2))/cellext);
rghtadd=max(rghtdiff,0);
ctldendim=size(ctlden);
ctlden=ctlden(1+max(topdiff,0):ctldendim(1)-max(bottomdiff,0),...
    1-min(leftdiff,0):ctldendim(2)+rghtdiff);
geotiffwrite('X:\CentralAmericaData\Model_inputs\clipped\ctlden_clp.tif',ctlden,Rslope);
