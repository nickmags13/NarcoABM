function [sub_tmovcorr,sub_tmovpval,sub_deptcorr,sub_deptpval]=...
    run_master_file(g,parmfname,poolobj,POP)
load savedrngstate.mat
addAttachedFiles(poolobj,{'calc_neival.m','optimizeroute_multidto.m',...
    'calc_intrisk.m','load_expmntl_parms_ga.m','fitness_calc.m',...
    'parsave_illicit_supplychain_ga.m'});
Master_supplychain_genalgo
% 
% savefname=sprintf('%sga_results_%d_%d','C:\Users\nrmagliocca\Box Sync\Data Drive\model_results\SupplyChain_011118_ga\',g);
% save(savefname,'t_firstmovcorr','t_firstmovpval','deptmovcorr','deptmovpval')

end