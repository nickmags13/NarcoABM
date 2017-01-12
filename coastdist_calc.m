%%%% Calculate distance to the coast %%%%
gridsize=size(ca_adm0);
ROWS=gridsize(1);
COLS=gridsize(2);
coastdist=zeros(size(ca_adm0));
icell=find(ca_adm0 > 0);
for i=1:length(icell)
   [irow,icol]=ind2sub(size(ca_adm0),icell(i)); 
   northvec=ca_adm0(max(irow-1,1):-1:1,icol);
   southvec=ca_adm0(min(irow+1,ROWS):1:ROWS,icol);
   westvec=ca_adm0(irow,max(icol-1,1):-1:1);
   eastvec=ca_adm0(irow,min(icol+1,COLS):1:COLS);
   northdist=find(northvec == 0,1,'first');
   southdist=find(southvec == 0,1,'first');
   eastdist=find(eastvec == 0,1,'first');
   westdist=find(westvec == 0,1,'first');
   distvec=[northdist southdist eastdist westdist];
   coastdist(icell)=min(distvec);
end

save coastdist coastdist