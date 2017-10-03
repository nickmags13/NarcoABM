function [neipick,neivalue,valuex]=calc_neival(c_trans,p_sl,y_node,q_node,lccf,...
    rtpref,tslrisk,dtonei,profmdl,cutflag)

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
    xpay_event(i)=pay_event(i); % payoff for route A with S&L event
end
for i=1:length(c_trans)
    ypay_noevent(i)=mean(pay_noevent(iset ~= i));   % payoff for all other routes with no S&L event
    ypay_event(i)=mean(pay_event(iset ~= i));       % payoff for all other routes with S&L event
    value_noevent(i)=abs(ypay_noevent(i)-xpay_noevent(i))/...
        (abs(ypay_noevent(i))+abs(xpay_noevent(i))+1);
    value_event(i)=abs(ypay_event(i)-xpay_event(i))/...
        (abs(ypay_event(i))+abs(xpay_event(i))+1);
    
    [~,ipntlval]=sort([value_noevent(i) value_event(i)],2,'descend');
    ival_noevent(i)=ipntlval(1);
    ival_event(i)=ipntlval(2);
    dwght_noevent(i)=(lccf.^ival_noevent(i))./((lccf.^ival_noevent(i))*(1-p_sl(i))+...
        (lccf.^ival_event(i))*p_sl(i));
    dwght_event(i)=(lccf.^ival_event(i))./((lccf.^ival_noevent(i))*(1-p_sl(i))+...
        (lccf.^ival_event(i))*p_sl(i));
    salwght_noevent(i)=(1-p_sl(i))*dwght_noevent(i);
    salwght_event(i)=p_sl(i)*dwght_event(i);
%     valuey(i)=rtpref(i)*salwght_noevent(i)*ypay_noevent(i)+...
%         (1-rtpref(i))*salwght_event(i)*ypay_event(i);
%     valuex(i)=rtpref(i)*salwght_noevent(i)*xpay_noevent(i)+...
%         (1-rtpref(i))*salwght_event(i)*xpay_event(i);
    valuey(i)=salwght_noevent(i)*ypay_noevent(i)+...
        salwght_event(i)*ypay_event(i);
    valuex(i)=salwght_noevent(i)*xpay_noevent(i)+...
        salwght_event(i)*xpay_event(i);
end
%%% Select only the most salient route options
%[ineivalue,ineipick]=max([valuex valuey],[],2);
% iroute=find(ineipick == 1);
% % % Select all route options that are within a specified loss tolerance
% % iroute=find((1+losstol).*valuex > valuey);
% if isempty(find(iroute,1)) == 1
%     [~,iroute]=max(ineivalue,[],1);
% end
% rankroute=sortrows([ineivalue(iroute) min(totstock/length(iroute),...
%     totcpcty(iroute))' q_node(iroute)' iroute],-1);  %rank trafficking routes by salient payoff
% % rankroute=sortrows([ineivalue(iroute) p_sl(iroute)' q_node(iroute)' iroute],2);  %rank trafficking routes by risk level
% icut=find(cumsum(rankroute(:,2)) <= totstock);  % select route based on total capcity

%%% Selection based on maximize profits while less than average S&L risk
rankroute=sortrows([rtpref'.*valuex p_sl' q_node' iset' dtonei],-1);  %rank trafficking routes by salient payoff
dtos=unique(dtonei(dtonei~=0));
if length(dtos) > 1
    icut=[];
    for j=1:length(dtos)
        idto=find(rankroute(:,5) == dtos(j));
        if profmdl == 1
            if isempty(find(valuex(dtonei == dtos(j)) > 0,1)) == 1
                [~,subicut]=min(rankroute(idto,2),[],1);
            elseif isempty(find(rankroute(idto,1) > 0,1)) == 1
                subicut=find(rankroute(idto,1) >= 0);
            else
                subicut=find(rankroute(idto,1) > 0);
            end
            if cutflag(dtos(j)) == 1
                subicut=[];
            end
            icut=[icut; idto(subicut)];
        elseif profmdl == 2
            if isempty(find(valuex(dtonei == dtos(j)) > 0,1)) == 1
                [~,subicut]=min(rankroute(idto,2),[],1);
            elseif isempty(find(cumsum(rankroute(idto,1)) > 0,1)) == 1
                subicut=find(cumsum(rankroute(idto,1)) >= 0);
            else
                subicut=find(cumsum(rankroute(idto,1)) > 0);
            end
            if cutflag(dtos(j)) == 1
                subicut=[];
            end
            icut=[icut; idto(subicut)];
        end
    end
    if rankroute(rankroute(:,5) == 0,1) > min(rankroute(icut,1))
        icut=[icut; find(rankroute(:,5) == 0)];
    end
else
    if profmdl == 1
        if isempty(find(valuex > 0,1)) == 1
            [~,icut]=min(rankroute(:,2),[],1);
        elseif isempty(find(rankroute(:,1) > 0,1)) == 1
            icut=find(rankroute(:,1) >= 0);
        else
            icut=find(rankroute(:,1) > 0);
        end
    elseif profmdl == 2
        if isempty(find(valuex > 0,1)) == 1
            [~,icut]=min(rankroute(:,2),[],1);
        elseif isempty(find(cumsum(rankroute(:,1)) > 0,1)) == 1
            icut=find(cumsum(rankroute(:,1)) >= 0);
        else
            icut=find(cumsum(rankroute(:,1)) > 0);
        end
    end
end
neipick=rankroute(icut,4);
neivalue=rankroute(icut,1);