function [sl_max,sl_min,baserisk,riskmltplr,startstock,sl_learn,rt_learn,...
    losslim,prodgrow,targetseize,intcpctymodel,profitmodel,endstock,...
    growthmdl,timewght,locthink,expandmax,empSLflag]=load_expmntl_parms(ERUNS)
    
empSLflag=zeros(1,ERUNS);    %determines is empirical (1) or artificial (0) S&L schedule used

sl_max=113*ones(1,ERUNS);       %baseline; maximum interdiction capacity
% sl_max=111*ones(1,ERUNS);
sl_min=ceil(sl_max/6);          %baseline; minimum interdiction capacity
% sl_min=ceil(sl_max/2);

baserisk=0.4017*ones(1,ERUNS);     %baseline; threshold for risk premium, Caulkins et al. (1993)
% baserisk=0.5114*ones(1,ERUNS);
riskmltplr=2*ones(1,ERUNS);     %baseline; risk multiplier for risk premium, Caulkins et al. (1993)
% riskcomp=12000*ones(1,ERUNS);   %baseline; risk compensation factor ofr overland transport, Caulkins et al. (1993)

startstock=200*ones(1,ERUNS);       %baseline; stock at production node
endstock=150000*ones(1,ERUNS);

sl_learn=0.5450*ones(1,ERUNS);     %baseline; rate of interdiction learning
% sl_learn=0.7405*ones(1,ERUNS);
rt_learn=0.6551*ones(1,ERUNS);     %baseline; rate of network agent learning
% rt_learn=0.5154*ones(1,ERUNS);

losslim=0.0473*ones(1,ERUNS);      %baseline; loss tolerance
% losslim=0.1282*ones(1,ERUNS);
% losslim=[0.01 0.05 0.1 0.2 0.3];

% growthmdl=1*ones(1,ERUNS);       %production function; 1 = linear, 2 = logistic
% prodgrow=1000*ones(1,ERUNS);   %baseline; rate of volume increase into CA per year, based on mean value for valume max in Gracias a Dios from CCDB
growthmdl=2*ones(1,ERUNS);       
prodgrow=0.75*ones(1,ERUNS);

timewght=ones(1,ERUNS);     %time discounting for subjective risk perception (Gallagher, 2014), range[0,1.05]
locthink=0.5*ones(1,ERUNS);     %'Local thinker' coefficient for salience function (Bordalo et al., 2012)

targetseize=0.2463*ones(1,ERUNS);  %baseline; target portion of total flow to seize
% targetseize=0.2413*ones(1,ERUNS);
intcpctymodel=ones(1,ERUNS);    %decreaseing(1) or increasing(2) capacity response to missing target seizures 

profitmodel=ones(1,ERUNS);    %profit maximization model for node selection: 1 is standard, 2 is cumulative

expandmax=8*ones(1,ERUNS);     %number of new nodes established per month
% expandmax=14*ones(1,ERUNS);
% expandmax=[1 5 10 15 20];

