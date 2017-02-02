%%%% Calculate distance to the coast %%%%
gridsize=size(ca_adm0);
ROWS=gridsize(1);
COLS=gridsize(2);
coastdist=zeros(size(ca_adm0));
save coastdistdata coastdist
icell=find(ca_adm0 > 0);
poolobj=parpool(12);
addAttachedFiles(poolobj,{'coastdistsave.m','load_cdist.m'});
parfor i=1:length(icell)
   index=icell(i);
   cdistfname='C:\Users\nmagliocca\Documents\Matlab_code\NarcoLogic\coastdistdata.mat';
%    load_cdist(cdistfname);
   [irow,icol]=ind2sub(size(ca_adm0),index); 
   northvec=ca_adm0(max(irow-1,1):-1:1,icol);
   southvec=ca_adm0(min(irow+1,ROWS):1:ROWS,icol);
   westvec=ca_adm0(irow,max(icol-1,1):-1:1);
   eastvec=ca_adm0(irow,min(icol+1,COLS):1:COLS);
   northdist=find(northvec == 0,1,'first');
   southdist=find(southvec == 0,1,'first');
   eastdist=find(eastvec == 0,1,'first');
   westdist=find(westvec == 0,1,'first');
   distvec=[northdist southdist eastdist westdist];
   savefname=('coastdistdata.mat');
   subcdist=cdist.coastdist;
   coastdistsave(distvec,subcdist,index,savefname)
%    coastdist(icell(i))=min(distvec);
end
% save coastdistdata coastdist