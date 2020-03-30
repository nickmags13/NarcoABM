%%%%%%%%%%   Interdiction events from optimization model   %%%%%%%%%%%%%%%%
function [intrdct_events]=optimize_interdiction(t,ADJ)

cd C:\Users\nrmagliocca\'Google Drive'\NSF_EAGER_Models

readflag=0;
fnames=dir;
fnamescell=struct2cell(fnames);
t1=double(t >= 100);
if t>= 100
    t2=floor((t-100)/10);
else
    t2=floor(t/10);
end
t3=mod(t,10);
trgtfile=sprintf('%d%d%d_INT_3_MCI.txt',t1,t2,t3);
while readflag == 0
    
   h=strcmp(trgtfile,fnamescell(1,:));
   hind=find(h == 1);
   readflag=~isempty(find(hind,1));
end
Tintevent=readmatrix(trgtfile);
intrdct_events=zeros(size(ADJ));
for j=1:length(Tintevent)
    iupstream=(ADJ(:,Tintevent(j))==1);
    intrdct_events(iupstream,Tintevent(j))=1;
end
cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM