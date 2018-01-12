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
outcomeset=zeros(POP,7,GEN);
evalset=zeros(POP,7,GEN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng default
% Parameter set
parmset(:,1,1)=round(100+(200-100)*rand(POP,1)); %maximum capacity
parmset(:,2,1)=0.3+(0.75-0.3)*rand(POP,1);   %baserisk
parmset(:,3,1)=0.25+(0.55-0.25)*rand(POP,1);   %sl_learn
parmset(:,4,1)=0.1+(0.75-0.1)*rand(POP,1);   %rt_learn
parmset(:,5,1)=0.01+(0.05-0.01)*rand(POP,1); %losslim
parmset(:,6,1)=0.1+(0.25-0.1)*rand(POP,1);   %targetseize
parmset(:,7,1)=round(8+(16-1)*rand(POP,1));  %expandmax

pset=(1:PARMS);

genset=zeros([],1);
finalset=zeros([],PARMS);
finalevalset=zeros([],7);
finaloutset=zeros([],7);

% Indicator variables
tmovcorr=zeros(POP,GEN);
tmovpval=zeros(POP,GEN);
deptcorr=zeros(POP,6,GEN);
deptpval=zeros(POP,6,GEN);

% addAttachedFiles(poolobj,{'calc_neival.m','optimizeroute_multidto.m',...
%     'calc_intrisk.m','load_expmntl_parms_ga.m','Master_supplychain_genalgo.m',...
%     'run_master_file.m'});
parmfname=sprintf('%sparmsfile_%d','C:\Users\nrmagliocca\Box Sync\Data Drive\model_results\SupplyChain_011218_ga\',1);
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
    evalset(:,:,g)=([tmovpval(:,g) deptpval(:,:,g)] < 0.1); %p-value based
    %     evalset(:,:,g)=[tmovcorr(:,g) deptcorr(:,:,g)];     %correlation coefficient based
    outcomeset(:,:,g)=[tmovpval(:,g) deptpval(:,:,g)];
    %     outcomeset(:,:,g)=[tmovcorr(:,g) deptcorr(:,:,g)];
    
    fitness(:,:,g)=[sum(evalset(:,:,g),2) (1:POP)'];
    sortfit=sortrows(fitness(:,:,g),-1);
    fitcut=round(POP*0.1667);

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
    parmset(11:POP,1,g+1)=round(100+(200-100)*rand(20,1));
    parmset(11:POP,2,g+1)=0.3+(0.75-0.3)*rand(20,1);
    parmset(11:POP,3,g+1)=0.25+(0.55-0.25)*rand(20,1);
    parmset(11:POP,4,g+1)=0.1+(0.75-0.1)*rand(20,1);
    parmset(11:POP,5,g+1)=0.01+(0.05-0.01)*rand(20,1);
    parmset(11:POP,6,g+1)=0.1+(0.25-0.1)*rand(20,1);
    parmset(11:POP,7,g+1)=round(8+(16-8)*rand(20,1));
    
    parmset(:,:,g+1)=[max(min(parmset(:,1,g+1),200),100) ...
        max(min(parmset(:,2,g+1),0.75),0.3) ...
        max(min(parmset(:,3,g+1),0.55),0.25) ...
        max(min(parmset(:,4,g+1),0.75),0.1) ...
        max(min(parmset(:,5,g+1),0.05),0.01) ...
        max(min(parmset(:,6,g+1),0.25),0.1) ...
        max(min(parmset(:,7,g+1),16),8)];
    
%     parmfname=sprintf('parmsfile_%d',g);
    parmfname=sprintf('%s_parmsfile_%d','C:\Users\nrmagliocca\Box Sync\Data Drive\model_results\SupplyChain_011218_ga',g+1);
    g_id=g+1;
    save(parmfname,'g_id','parmset');
    
    ikeep=find(fitness(:,1,g) == max(fitness(:,1,g)));
    genset(length(genset)+1:length(genset)+length(ikeep),1)=g;
    finalset(size(finalset,1)+1:size(finalset,1)+length(ikeep),:)=parmset(ikeep,:,g);
    finalevalset(size(finalevalset,1)+1:size(finalevalset,1)+length(ikeep),:)=evalset(ikeep,:,g);
    finaloutset(size(finaloutset,1)+1:size(finaloutset,1)+length(ikeep),:)=outcomeset(ikeep,:,g);
end
toc
delete(poolobj)