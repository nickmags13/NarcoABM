%%%%%%%%%%%%%%% Aster Global DEM Processing %%%%%%%%%%%%%%%%%%%%
cd X:\CentralAmericaData\ASTER\CAextent

fnames=dir;
fnamescell=struct2cell(fnames);
% (2) change the prefix of the results file names
h=strncmp('AST',fnamescell(1,:),3);
hind=find(h==1);

for j=1:length(hind)
    filename=fnamescell{1,hind(j)};
    [Z,R]=geotiffread(filename);
    dims=size(Z);
    maxy=dims(1);
    maxx=dims(2);
    ipeak=find(Z > 4240);% max elevation in CA is 4220 m +/- 20m vertical accuracy
    for k=1:length(ipeak)
        [row,col]=ind2sub(size(Z),ipeak(k));
        ltslp=Z(row,max(col-1,1))+mean(diff(Z(row,max(col-3:col-1,1))));
        rtslp=Z(row,min(col+1,maxx))-mean(diff(Z(row,min(col+1:col+3,maxx))));
        tpslp=Z(max(row-1,1),col)+mean(diff(Z(max(row-3:row-1,1),col)));
        btslp=Z(min(row+1,maxy),col)-mean(diff(Z(min(row+1:row+3,maxy),col)));
        Z(ipeak(k))=min(mean([ltslp rtslp btslp tpslp]),4240);
    end
%     %%% Fill-in no data (-9999) pixels with interpolation
%     % Would need to find full extent of no data patches, widen subsample to
%     % encompass, then set query points within to run interpolation. Edges
%     % would be a challenge.
%     igap=find(Z == -9999);
%     minz=min(Z(Z ~= -9999));
%     for g=1:length(igap)
%         [row,col]=ind2sub(size(Z),igap(g));
%         ltslp=Z(row,max(col-1,1))+mean(diff(Z(row,max(col-3:col-1,1))));
%         rtslp=Z(row,min(col+1,maxx))-mean(diff(Z(row,min(col+1:col+3,maxx))));
%         tpslp=Z(max(row-1,1),col)+mean(diff(Z(max(row-3:row-1,1),col)));
%         btslp=Z(min(row+1,maxy),col)-mean(diff(Z(min(row+1:row+3,maxy),col)));
%         Z(igap(k))=max(min(mean([ltslp rtslp btslp tpslp]),4240),minz);
%     end
    % write new dem as geotiff
    cd X:\CentralAmericaData\ASTER\CAextent\ca_dem_processed_tiles
    
    geotiffwrite(filename,Z,R)
    cd X:\CentralAmericaData\ASTER\CAextent
end

%%%% Calculate slope for CA dem
cd X:\CentralAmericaData\Model_inputs

[Z,R]=geotiffread('ca_dem_250m.tif');
Zclean=Z;
Zclean(Z == -9999)=NaN;
[aspect,slope,gradN,gradE]=gradientm(Z,R);

%%% clean up edge effects
slope(gradE < -19)=0;
slope(gradE > 19)=0;
slope(gradN > 19)=0;
slope(gradN < -19)=0;

geotiffwrite('ca_slope_250m',slope,R);
% %%% slice array to process separately
% [nw_aspect,nw_slope,nw_gradN,nw_gradE] = gradientm(...
%     Zclean(1:ceil(length(Zclean(:,1))/2),1:ceil(length(Zclean(1,:))/2)),R);