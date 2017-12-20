function [sl_max,sl_min,baserisk,riskmltplr,startstock,sl_learn,rt_learn,...
    losslim,prodgrow,targetseize,intcpctymodel,profitmodel,endstock,...
    growthmdl,timewght,locthink,expandmax,empSLflag]=load_expmntl_parms(ERUNS)
    
empSLflag=zeros(1,ERUNS);    %determines is empirical (1) or artificial (0) S&L schedule used

sl_max=400*ones(1,ERUNS);       %baseline; maximum interdiction capacity
% sl_max=[160 240 320 400 480];
sl_min=ceil(sl_max/6);          %baseline; minimum interdiction capacity
% sl_min=ceil(sl_max/2);

baserisk=0.5*ones(1,ERUNS);     %baseline; threshold for risk premium, Caulkins et al. (1993)
riskmltplr=2*ones(1,ERUNS);     %baseline; risk multiplier for risk premium, Caulkins et al. (1993)
% riskcomp=12000*ones(1,ERUNS);   %baseline; risk compensation factor ofr overland transport, Caulkins et al. (1993)

startstock=200*ones(1,ERUNS);       %baseline; stock at production node
endstock=150000*ones(1,ERUNS);

sl_learn=0.5*ones(1,ERUNS);     %baseline; rate of interdiction learning
rt_learn=0.5*ones(1,ERUNS);     %baseline; rate of network agent learning

losslim=0.1*ones(1,ERUNS);      %baseline; loss tolerance

% growthmdl=1*ones(1,ERUNS);       %production function; 1 = linear, 2 = logistic
% prodgrow=1000*ones(1,ERUNS);   %baseline; rate of volume increase into CA per year, based on mean value for valume max in Gracias a Dios from CCDB
growthmdl=2*ones(1,ERUNS);       
prodgrow=0.75*ones(1,ERUNS);

timewght=ones(1,ERUNS);     %time discounting for subjective risk perception (Gallagher, 2014), range[0,1.05]
locthink=0.5*ones(1,ERUNS);     %'Local thinker' coefficient for salience function (Bordalo et al., 2012)

targetseize=0.3*ones(1,ERUNS);  %baseline; target portion of total flow to seize
intcpctymodel=ones(1,ERUNS);    %decreaseing(1) or increasing(2) capacity response to missing target seizures 

profitmodel=ones(1,ERUNS);    %profit maximization model for node selection: 1 is standard, 2 is cumulative

expandmax=10*ones(1,ERUNS);     %number of new nodes established per month

