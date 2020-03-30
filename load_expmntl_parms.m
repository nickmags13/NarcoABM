function [sl_max,sl_min,baserisk,riskmltplr,startstock,sl_learn,rt_learn,...
    losslim,prodgrow,targetseize,intcpctymodel,profitmodel,endstock,...
    growthmdl,timewght,locthink,expandmax,empSLflag,optSLflag,suitflag]=load_expmntl_parms(ERUNS)
    
empSLflag=zeros(1,ERUNS);    %determines is empirical (1) or artificial (0) S&L schedule used
% optSLflag=zeros(1,ERUNS);   %use interdiction events from optimization model (1)
optSLflag=ones(1,ERUNS);
suitflag=zeros(1,ERUNS);    % use RAT suitability (1) or build from covariates (0)
% suitflag=ones(1,ERUNS);
sl_max=125*ones(1,ERUNS);       %baseline; maximum interdiction capacity
% sl_max=[75 100 125 150 175];
% sl_max=114*ones(1,ERUNS);
sl_min=ceil(sl_max/6);          %baseline; minimum interdiction capacity
% sl_min=ceil(sl_max/2);

% % baserisk=0.4592*ones(1,ERUNS);     %baseline; threshold for risk premium, Caulkins et al. (1993)
% baserisk=[0.33 0.38 0.43 0.48 0.53];
baserisk=0.43*ones(1,ERUNS);
riskmltplr=2*ones(1,ERUNS);     %baseline; risk multiplier for risk premium, Caulkins et al. (1993)
% riskcomp=12000*ones(1,ERUNS);   %baseline; risk compensation factor ofr overland transport, Caulkins et al. (1993)

startstock=200*ones(1,ERUNS);       %baseline; stock at production node (MT)
endstock=111500*ones(1,ERUNS);
% endstock=71250*ones(1,ERUNS);

sl_learn=0.6*ones(1,ERUNS);     %baseline; rate of interdiction learning
% sl_learn=[0.4 0.5 0.6 0.7 0.8];
% sl_learn=0.7477*ones(1,ERUNS);
rt_learn=0.3558*ones(1,ERUNS);     %basline; rate of network agent learning
% rt_learn=[0.25 0.3 0.35 0.4 0.45];
% rt_learn=0.2408*ones(1,ERUNS);

losslim=0.0455*ones(1,ERUNS);      %baseline; loss tolerance
% losslim=0.0399*ones(1,ERUNS);
% losslim=[0.01 0.025 0.05 0.1 0.2];

% growthmdl=1*ones(1,ERUNS);       %production function; 1 = linear, 2 = logistic
% prodgrow=1000*ones(1,ERUNS);   %baseline; rate of volume increase into CA per year, based on mean value for valume max in Gracias a Dios from CCDB
growthmdl=2*ones(1,ERUNS);       
prodgrow=0.75*ones(1,ERUNS);

timewght=ones(1,ERUNS);     %time discounting for subjective risk perception (Gallagher, 2014), range[0,1.05]
locthink=0.5*ones(1,ERUNS);     %'Local thinker' coefficient for salience function (Bordalo et al., 2012)

targetseize=0.3417*ones(1,ERUNS);  %baseline; target portion of total flow to seize
% targetseize=[0.25 0.3 0.35 0.4 0.45];
% targetseize=0.2669*ones(1,ERUNS);
intcpctymodel=ones(1,ERUNS);    %decreaseing(1) or increasing(2) capacity response to missing target seizures 

profitmodel=ones(1,ERUNS);    %profit maximization model for node selection: 1 is standard, 2 is cumulative

expandmax=8*ones(1,ERUNS);     %number of new nodes established per month
% expandmax=9*ones(1,ERUNS);
% expandmax=[6 7 8 9 10];

