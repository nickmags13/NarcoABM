%%%%%%%%%%%%%%   Process distance matricies   %%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
cd X:\CentralAmericaData\Model_inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Load template geotiff distance to coast   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script finds the minimum linear distance to the coastline in each
% cardinal direction, and takes the global minimum to create a single
% distance to coast raster
[Z,R]=geotiffread('ca_cntry_250.tif');
Z=double(Z);
% %%%% Make a test raster
subZ=Z(1350:1650,3700:4000);
imagesc(subZ)
subR = R;
subR.RasterSize=size(subZ);
subR.LatitudeLimits=[R.LatitudeLimits(2)-1800*R.CellExtentInLatitude ...
    R.LatitudeLimits(2)-1400*R.CellExtentInLatitude];
subR.LongitudeLimits=[R.LongitudeLimits(1)+3600*R.CellExtentInLongitude ...
    R.LongitudeLimits(1)+4000*R.CellExtentInLongitude];
% subR.RasterExtentInLatitude=subR.LatitudeLimits(2)-subR.LatitudeLimits(1);
% subR.RasterExtentInLongiitude=subR.LongitudeLimits(1)-subR.LongitudeLimits(1);
subR.ColumnsStartFrom=R.ColumnsStartFrom;
subR.RowsStartFrom=R.RowsStartFrom;
subR.CellExtentInLatitude=R.CellExtentInLatitude;
subR.CellExtentInLongitude=R.CellExtentInLongitude;
% testcoastdist=zeros(size(subZ));
NoDataVal=255;
comparestack=repmat(zeros(size(subZ)),1,1,4);
celldist=0.250; %250 meter resolution dataset
% matrix to minimize
maxval=max(size(subZ));
indmat=(maxval+1)*ones(size(subZ));
% edge check
% can I modify this to deal with irregular coastlines?
w_edgechk=repmat((subZ(:,1) == NoDataVal),1,length(subZ(1,:)));
e_edgechk=repmat((subZ(:,length(subZ(1,:))) == NoDataVal),1,length(subZ(1,:)));
n_edgechk=repmat((subZ(1,:) == NoDataVal),length(subZ(:,1)),1);
s_edgechk=repmat((subZ(length(subZ(:,1)),:) == NoDataVal),length(subZ(:,1)),1);
% calc distance from edges

% deal with irregular coastlines
testvec=west2east(25,1:50);
diffvec=diff(testvec,1,2);
idiff=find(diffvec == 0);
testvec(1+idiff)=testvec(idiff+1)-testvec(idiff);


testvec=north2south(1:50,25);
diffvec=diff(testvec,1,1);
idiff=find(diffvec == 0);
testvec(1+idiff)=testvec(idiff+1)-testvec(idiff);

logiclyr=(subZ ~= NoDataVal);
west2east=cumsum(logiclyr,2);
wdiff=diff(west2east,1,2);
iwdiff=find(wdiff == 0);
[rowdiff,coldiff]=ind2sub(size(wdiff),iwdiff);
ishft=sub2ind(size(west2east),rowdiff,coldiff+1);
west2east(ishft)=west2east(ishft)-west2east(iwdiff);
comparestack(:,:,1)=west2east;

east2west=cumsum(logiclyr(:,length(logiclyr(1,:)):-1:1),2);
ediff=diff(east2west,1,2);
iediff=find(ediff == 0);
[rowdiff,coldiff]=ind2sub(size(ediff),iediff);
ishft=sub2ind(size(east2west),rowdiff,coldiff+1);
east2west(ishft)=east2west(iediff)-east2west(ishft);
east2west=fliplr(east2west);
comparestack(:,:,2)=east2west;

north2south=flipud(rot90(cumsum(logiclyr,1),1));
ndiff=diff(north2south,1,2);
indiff=find(ndiff == 0);
[rowdiff,coldiff]=ind2sub(size(ndiff),indiff);
ishft=sub2ind(size(north2south),rowdiff,coldiff+1);
north2south(ishft)=north2south(ishft)-north2south(indiff);
north2south=rot90(flipud(north2south),-1);
comparestack(:,:,3)=north2south;

south2north=rot90(flipud(cumsum(logiclyr(length(logiclyr(:,1)):-1:1,:),1)),-1);
sdiff=diff(south2north,1,2);
isdiff=find(sdiff == 0);
[rowdiff,coldiff]=ind2sub(size(sdiff),isdiff);
ishft=sub2ind(size(south2north),rowdiff,coldiff+1);
south2north(ishft)=south2north(ishft)-south2north(isdiff);
south2north=rot90(south2north,1);
comparestack(:,:,4)=south2north;
% filter based on edge check and select minimum distance for each cell
w2e=indmat;
w2e(w_edgechk == 1)=west2east(w_edgechk == 1);
comparestack(:,:,1)=w2e;
e2w=indmat;
e2w(e_edgechk == 1)=east2west(e_edgechk == 1);
comparestack(:,:,2)=e2w;
n2s=indmat;
n2s(n_edgechk == 1)=north2south(n_edgechk == 1);
comparestack(:,:,3)=n2s;
s2n=indmat;
s2n(s_edgechk == 1)=south2north(s_edgechk == 1);
comparestack(:,:,4)=s2n;
testcoastdist=celldist*min(comparestack,[],3);

NoDataVal=255;
comparestack=repmat(zeros(size(Z)),1,1,4);
celldist=0.250; %250 meter resolution dataset
% matrix to minimize
maxval=max(size(Z));
indmat=(maxval+1)*ones(size(Z));
% edge check
w_edgechk=repmat((Z(:,1) == NoDataVal),1,length(Z(1,:)));
e_edgechk=repmat((Z(:,length(Z(1,:))) == NoDataVal),1,length(Z(1,:)));
n_edgechk=repmat((Z(1,:) == NoDataVal),length(Z(:,1)),1);
s_edgechk=repmat((Z(length(Z(:,1)),:) == NoDataVal),length(Z(:,1)),1);
% calc distance from edges
logiclyr=(Z ~= NoDataVal);
west2east=cumsum(logiclyr,2);
east2west=fliplr(cumsum(logiclyr(:,length(logiclyr(1,:)):-1:1),2));
north2south=cumsum(logiclyr,1);
south2north=flipud(cumsum(logiclyr(length(logiclyr(:,1)):-1:1,:),1));

% filter based on edge check and select minimum distance for each cell
w2e=indmat;
w2e(w_edgechk == 1)=west2east(w_edgechk == 1);
comparestack(:,:,1)=w2e;
e2w=indmat;
e2w(e_edgechk == 1)=east2west(e_edgechk == 1);
comparestack(:,:,2)=e2w;
n2s=indmat;
n2s(n_edgechk == 1)=north2south(n_edgechk == 1);
comparestack(:,:,3)=n2s;
s2n=indmat;
s2n(s_edgechk == 1)=south2north(s_edgechk == 1);
comparestack(:,:,4)=s2n;
coastdist=celldist*min(comparestack,[],3);

% iland=find(Z ~= NoDataVal);
% for i=1:length(iland)
%     [row,col]=ind2sub(size(Z),iland(i));
%     lats_in=R.LatitudeLimits(2)-row*R.CellExtentInLatitude;
%     lons_in=R.LongitudeLimits(1)+col*R.CellExtentInLongitude;
%     [dists_min,lats_closest,lons_closest]=dist_from_coast(lats_in,...
%         lons_in);
%     coastdist(iland(i))=dists_min/1000;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%   Create distance to national border   %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cntryid=unique(Z);
cntryid=cntryid(cntryid ~= NoDataVal);
borderdist=zeros(size(Z));
ndmat=NoDataVal*ones(size(Z));
celldist=0.250; %250 meter resolution dataset
% matrix to minimize
cntrymaxval=max(size(Z));
cntryindmat=(maxval+1)*ones(size(Z));
for c=1:length(cntryid)
    cntrymat=ndmat;
    cntrymat(Z == cntryid(c))=cntryid(c);
    
    cntrystck=repmat(zeros(size(Z)),1,1,4);
    % edge check
    w_edgechk=repmat((cntrymat(:,1) == NoDataVal),1,length(cntrymat(1,:)));
    e_edgechk=repmat((cntrymat(:,length(cntrymat(1,:))) == ...
        NoDataVal),1,length(cntrymat(1,:)));
    n_edgechk=repmat((cntrymat(1,:) == NoDataVal),length(cntrymat(:,1)),1);
    s_edgechk=repmat((cntrymat(length(cntrymat(:,1)),:) == ...
        NoDataVal),length(cntrymat(:,1)),1);
    % calc distance from edges
    logiclyr=(cntrymat ~= NoDataVal);
    west2east=cumsum(logiclyr,2);
    east2west=fliplr(cumsum(logiclyr(:,length(logiclyr(1,:)):-1:1),2));
    north2south=cumsum(logiclyr,1);
    south2north=flipud(cumsum(logiclyr(length(logiclyr(:,1)):-1:1,:),1));
    
    % filter based on edge check and select minimum distance for each cell
    w2e=cntryindmat;
    w2e(w_edgechk == 1)=west2east(w_edgechk == 1);
    cntrystck(:,:,1)=w2e;
    e2w=cntryindmat;
    e2w(e_edgechk == 1)=east2west(e_edgechk == 1);
    cntrystck(:,:,2)=e2w;
    n2s=cntryindmat;
    n2s(n_edgechk == 1)=north2south(n_edgechk == 1);
    cntrystck(:,:,3)=n2s;
    s2n=cntryindmat;
    s2n(s_edgechk == 1)=south2north(s_edgechk == 1);
    cntrystck(:,:,4)=s2n;
    mindist=celldist*min(cntrystck,[],3);
    borderdist(Z == cntryid(c))=celldist*mindist(Z == cntryid(c));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  Write geotiffs   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd X:\CentralAmericaData\Model_inputs
geotiffwrite('coastdist_250.tif',coastdist,R);
geotiffwrite('borderdist_250.tif',borderdist,R);
toc