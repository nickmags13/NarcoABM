function [neipick,neivalue]=calc_neival(c_trans,p_sl,y_node,q_node,lccf,...
    totstock,totcpcty)

pay_noevent=zeros(length(c_trans),1);
pay_event=zeros(length(c_trans),1);
xpay_noevent=zeros(length(c_trans),1);
ypay_noevent=zeros(length(c_trans),1);
xpay_event=zeros(length(c_trans),1);
ypay_event=zeros(length(c_trans),1);
value_noevent=zeros(length(c_trans),1);
value_event=zeros(length(c_trans),1);
ival_noevent=zeros(length(c_trans),1);
ival_event=zeros(length(c_trans),1);
dwght_noevent=zeros(length(c_trans),1);
dwght_event=zeros(length(c_trans),1);
salwght_noevent=zeros(length(c_trans),1);
salwght_event=zeros(length(c_trans),1);
valuex=zeros(length(c_trans),1);
valuey=zeros(length(c_trans),1);
iset=1:length(c_trans);
for i=1:length(c_trans)
    pay_noevent(i)=y_node(i)*q_node(i)-c_trans(i)*q_node(i);  % payoff with no S&L event
    pay_event(i)=y_node(i)*q_node(i)-c_trans(i)*q_node(i)-y_node(i)*q_node(i);    % payoff with S&L event
    
    xpay_noevent(i)=pay_noevent(i); % payoff for route A with no S&L event
    xpay_event(i)=pay_event(i);     % payoff for route A with S&L event
    ypay_noevent(i)=mean(pay_noevent(iset ~= i));   % payoff for all other routes with no S&L event
    ypay_event(i)=mean(pay_event(iset ~= i));       % payoff for all other routes with S&L event
    value_noevent(i)=abs(ypay_noevent(i)-xpay_noevent(i))/...
        (abs(ypay_noevent(i))+abs(xpay_noevent(i)));
    value_event(i)=abs(ypay_event(i)-xpay_event(i))/...
        (abs(ypay_event(i))+abs(xpay_event(i)));
    
    [~,ipntlval]=sort([value_noevent(i) value_event(i)],2,'descend');
    ival_noevent(i)=ipntlval(1);
    ival_event(i)=ipntlval(2);
    dwght_noevent(i)=(lccf.^ival_noevent(i))./((lccf.^ival_noevent(i))*(1-p_sl(i))+...
        (lccf.^ival_event(i))*p_sl(i));
    dwght_event(i)=(lccf.^ival_event(i))./((lccf.^ival_noevent(i))*(1-p_sl(i))+...
        (lccf.^ival_event(i))*p_sl(i));
    salwght_noevent(i)=(1-p_sl(i))*dwght_noevent(i);
    salwght_event(i)=p_sl(i)*dwght_event(i);
    valuey(i)=salwght_noevent(i)*ypay_noevent(i)+salwght_event(i)*ypay_event(i);
    valuex(i)=salwght_noevent(i)*xpay_noevent(i)+salwght_event(i)*xpay_event(i);
end
[ineivalue,ineipick]=max([valuex valuey],[],2);
iroute=find(ineipick == 1);
rankroute=sortrows([ineivalue(iroute) totcpcty(iroute)' q_node(iroute)' iroute],-1);  %rank trafficking routes by salient payoff
% rankroute=sortrows([ineivalue(iroute) p_sl(iroute)' q_node(iroute)' iroute],2);  %rank trafficking routes by risk level
icut=find(cumsum(rankroute(:,2)) <= totstock);  % select route based on total capcity
neipick=rankroute(icut,4);
neivalue=rankroute(icut,1);