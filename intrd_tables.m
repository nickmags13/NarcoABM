%%%%%% Interdiction Initialization %%%%%%%%%%
function [Tflow,Tintrd]=intrd_tables(FLOW,SLPROB,EdgeTable,t)

varTypes={'double','double','double'};
Tflow=table('Size',[height(EdgeTable) 3],'VariableTypes',varTypes,'VariableNames',{'End_Node','Start_Node','IntitFlow'});
startFLOW=FLOW(:,:,t);
% iflow=find(startFLOW > 0);
% [snd,rec]=ind2sub(size(startFLOW),iflow);
% Tflow(1:length(rec),1)=array2table(rec);
% Tflow(:,2)=array2table(snd);
% Tflow(:,3)=array2table(startFLOW(iflow));

Tintrd=table('Size',[height(EdgeTable) 3],'VariableTypes',varTypes,'VariableNames',{'End_Node','Start_Node','IntitProb'});
if t == 1
    startSLPROB=SLPROB(:,:,1);
else
    startSLPROB=SLPROB(:,:,t-1);
end
sumprob=sum(sum(startSLPROB));
for g=1:height(EdgeTable)
    edge=table2array(EdgeTable(g,1));
    Tflow(g,1)={edge(2)};
    Tflow(g,2)={edge(1)};
    Tflow(g,3)={startFLOW(edge(1),edge(2))};
    
    Tintrd(g,1)={edge(2)};
    Tintrd(g,2)={edge(1)};
    Tintrd(g,3)={startSLPROB(edge(1),edge(2))/sumprob};
%     iedge=find(startSLPROB > 0);
%     tindex=height(Tintrd(:,1))+1:height(Tintrd(:,1))+length(iedge);
%     Tintrd(tindex,1)=array2table(iedge');
%     Tintrd(tindex,2)=array2table(ones(length(iedge),1)*i);
%     Tintrd(tindex,3)=array2table(SLPROB(i,iedge)');
end
t1=double(t >= 100);
if t>= 100
    t2=floor((t-100)/10);
else
    t2=floor(t/10);
end
t3=mod(t,10);
cd C:\Users\nrmagliocca\'Google Drive'\NSF_EAGER_Models
writetable(Tflow,sprintf('%d%d%d_SL_3_MCI.txt',t1,t2,t3))
writetable(Tintrd,sprintf('%d%d%d_SLPROB_3_MCI.txt',t1,t2,t3)) %apply to next time step's interdiction actions
cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM