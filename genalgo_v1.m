%%%%%%%%%%%   NarcoLogic Genetic Algorithm Parameter Searh   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
poolobj=parpool(10);

PARMS=7;
POP=30;
GEN=10;
parmset=zeros(POP,PARMS,GEN);
% parmset = [nfwageparm farmcostparm nfarmcostparm croppriceparm]
fitness=zeros(POP,2,GEN);   %[fitness_score id]
outcomeset=zeros(POP,6,GEN);
evalset=zeros(POP,6,GEN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng default
% Parameter set
parmset(:,1,1)=round(100+(140-100)*rand(POP,1)); %maximum capacity
parmset(:,2,1)=0.4+(0.6-0.4)*rand(POP,1);   %baserisk
parmset(:,3,1)=0.5+(0.8-0.5)*rand(POP,1);   %sl_learn
parmset(:,4,1)=0.2+(0.5-0.2)*rand(POP,1);   %rt_learn
parmset(:,5,1)=0.035+(0.045-0.035)*rand(POP,1); %losslim
parmset(:,6,1)=0.25+(0.35-0.25)*rand(POP,1);   %targetseize
parmset(:,7,1)=round(8+(10-8)*rand(POP,1));  %expandmax

pset=(1:PARMS);

genset=zeros([],1);
finalset=zeros([],PARMS);
finalevalset=zeros([],6);
finaloutset=zeros([],6);

% Indicator variables
tmovcorr=zeros(POP,GEN);
tmovpval=zeros(POP,GEN);
deptcorr=zeros(POP,6,GEN);
deptpval=zeros(POP,6,GEN);

% addAttachedFiles(poolobj,{'calc_neival.m','optimizeroute_multidto.m',...
%     'calc_intrisk.m','load_expmntl_parms_ga.m','Master_supplychain_genalgo.m',...
%     'run_master_file.m'});
parmfname=sprintf('%sparmsfile_%d','C:\Users\nrmagliocca\Box Sync\Data Drive\model_results\SupplyChain_013018_ga\',1);
g_id=1;
save(parmfname,'g_id','parmset');

for g=1:GEN
    [sub_tmovcorr,sub_tmovpval,sub_deptcorr,sub_deptpval]=...
        run_master_file(g,parmfname,poolobj,POP);
    
    
    tmovcorr(:,g)=sub_tmovcorr;
    tmovpval(:,g)=sub_tmovpval;
    deptcorr(:,:,g)=sub_deptcorr;
    deptpval(:,:,g)=sub_deptpval;
    
    %Fitness evaluation
%     evalset(:,:,g)=([tmovpval(:,g) deptpval(:,:,g)] < 0.1); %p-value based
%     evalset(:,:,g)=[tmovcorr(:,g) deptcorr(:,:,g)];     %correlation coefficient based
    evalset(:,:,g)=[tmovcorr(:,g) deptcorr(:,2:6,g)];
%     outcomeset(:,:,g)=[tmovpval(:,g) deptpval(:,:,g)];
%     outcomeset(:,:,g)=[tmovcorr(:,g) deptcorr(:,:,g)];
    outcomeset(:,:,g)=[tmovcorr(:,g) deptcorr(:,2:6,g)];
    
    fitness(:,:,g)=[sum(evalset(:,:,g),2) (1:POP)'];
    sortfit=sortrows(fitness(:,:,g),-1);
    fitcut=round(POP*0.1667);
%     fitcut=length(find(fitness(:,1,g) >= max(reshape(fitness(:,1,:),POP*GEN,1))));

    fitset=(1:fitcut);
    bestset=parmset(sortfit(1:fitcut,2),:,g);

    % crossover
    crsspt=ceil((PARMS-1)*rand(fitcut,1));
    pairs=round(1+(fitcut-1)*rand(fitcut,2));
    
    %mutant
    mutes=ceil((PARMS-2)*rand(10,1));

    for i=1:5
        
        if pairs(i,1) == pairs(i,2)
            altset=find(fitset ~= pairs(i,2));
            randpick=ceil(length(altset)*rand(1));
            pairs(i,2)=fitset(altset(randpick));
        end
        parmset(i,:,g+1)=[bestset(pairs(i,1),1:crsspt(i)) ...
            bestset(pairs(i,2),crsspt(i)+1:PARMS)];
    end
        
%     for i=1:length(mutes)
    for i=1:5
        mutpts=ceil(PARMS*rand(1,mutes(i)));
        mutset=round(1+(length(bestset(:,1))-1)*rand(1,mutes(i)));
        mutant=bestset(ceil(length(bestset(:,1))*rand(1)),:);
        for im=1:mutes(i)
            mutant(mutpts(im))=bestset(mutset(im),mutpts(im));
        end
        parmset(i+fitcut,:,g+1)=mutant;
    end

    % New, random generation
    parmset(11:POP,1,g+1)=round(100+(140-100)*rand(20,1));
    parmset(11:POP,2,g+1)=0.4+(0.6-0.4)*rand(20,1);
    parmset(11:POP,3,g+1)=0.5+(0.8-0.5)*rand(20,1);
    parmset(11:POP,4,g+1)=0.2+(0.5-0.2)*rand(20,1);
    parmset(11:POP,5,g+1)=0.035+(0.045-0.035)*rand(20,1);
    parmset(11:POP,6,g+1)=0.25+(0.35-0.25)*rand(20,1);
    parmset(11:POP,7,g+1)=round(8+(10-8)*rand(20,1));
    
    parmset(:,:,g+1)=[max(min(parmset(:,1,g+1),140),100) ...
        max(min(parmset(:,2,g+1),0.6),0.4) ...
        max(min(parmset(:,3,g+1),0.8),0.5) ...
        max(min(parmset(:,4,g+1),0.5),0.2) ...
        max(min(parmset(:,5,g+1),0.05),0.03) ...
        max(min(parmset(:,6,g+1),0.35),0.25) ...
        max(min(parmset(:,7,g+1),10),8)];
    
%     parmfname=sprintf('parmsfile_%d',g);
    parmfname=sprintf('%sparmsfile_%d','C:\Users\nrmagliocca\Box Sync\Data Drive\model_results\SupplyChain_013018_ga\',g+1);
    g_id=g+1;
    save(parmfname,'g_id','parmset');
    
%     ikeep=find(fitness(:,1,g) == max(fitness(:,1,g)));
    ikeep=find(fitness(:,1,g) >= 0.95*max(reshape(fitness(:,1,:),POP*GEN,1)));
    genset(length(genset)+1:length(genset)+length(ikeep),1)=g;
    finalset(size(finalset,1)+1:size(finalset,1)+length(ikeep),:)=parmset(ikeep,:,g);
    finalevalset(size(finalevalset,1)+1:size(finalevalset,1)+length(ikeep),:)=evalset(ikeep,:,g);
    finaloutset(size(finaloutset,1)+1:size(finaloutset,1)+length(ikeep),:)=outcomeset(ikeep,:,g);
end
fname=sprintf('%sga_results_013018','C:\Users\nrmagliocca\Box Sync\Data Drive\model_results\SupplyChain_013018_ga\');
save(fname,'finalset','finalevalset','finaloutset','genset','deptcorr','deptpval','evalset','outcomeset')
toc
delete(poolobj)