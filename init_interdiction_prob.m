%%%%%% Create table of interdiction probabilities %%%%%%%%%%

cd \\asfs.asnet.ua-net.ua.edu\users$\home\nrmagliocca\'My Documents'\MATLAB\NarcoLogic\NarcoABM

load intprob_080719.mat

sz=[0 3];
varTypes={'double','double','double'};
T=table('Size',sz,'VariableTypes',varTypes,'VariableNames',{'End_Node','Start_Node','IntProb'});
for i=1:size(INTPROB,1)-1
    iedge=find(INTPROB(i,:) > 0);
    tindex=height(T(:,1))+1:height(T(:,1))+length(iedge);
    T(tindex,1)=array2table(iedge');
    T(tindex,2)=array2table(ones(length(iedge),1)*i);
    T(tindex,3)=array2table(INTPROB(i,iedge)');
end

writetable(T,'init_interdiction_prob.csv')