function [sl_max,sl_min,baserisk,riskmltplr,startstock,sl_learn,rt_learn,...
    losslim,prodgrow,targetseize,intcpctymodel,profitmodel,endstock,...
    growthmdl,timewght,locthink,expandmax,empSLflag]=load_expmntl_parms_ga(parmfname,p)
    
load(parmfname);

empSLflag=0;    %determines is empirical (1) or artificial (0) S&L schedule used

sl_max=parmset(p,1,g_id);       %baseline; maximum interdiction capacity
sl_min=ceil(sl_max/6);          %baseline; minimum interdiction capacity
% sl_min=ceil(sl_max/2);

baserisk=parmset(p,2,g_id);     %baseline; threshold for risk premium, Caulkins et al. (1993)
riskmltplr=2;     %baseline; risk multiplier for risk premium, Caulkins et al. (1993)
% riskcomp=12000*ones(1,ERUNS);   %baseline; risk compensation factor ofr overland transport, Caulkins et al. (1993)

startstock=200;       %baseline; stock at production node
endstock=150000;

sl_learn=parmset(p,3,g_id);     %baseline; rate of interdiction learning
rt_learn=parmset(p,4,g_id);     %baseline; rate of network agent learning

losslim=parmset(p,5,g_id);      %baseline; loss tolerance

growthmdl=2;       %production function; 1 = linear, 2 = logistic
prodgrow=0.75;   %baseline; rate of volume increase into CA per year, based on mean value for valume max in Gracias a Dios from CCDB


timewght=1;     %time discounting for subjective risk perception (Gallagher, 2014), range[0,1.05]
locthink=0.5;     %'Local thinker' coefficient for salience function (Bordalo et al., 2012)

targetseize=parmset(p,6,g_id);  %baseline; target portion of total flow to seize
intcpctymodel=1;    %decreaseing(1) or increasing(2) capacity response to missing target seizures 

profitmodel=1;    %profit maximization model for node selection: 1 is standard, 2 is cumulative

expandmax=parmset(p,7,g_id);     %number of new nodes established per month

