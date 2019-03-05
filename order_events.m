function [events,flag] = order_events(events,loc_mode)

if nargin == 1
    loc_mode = 'walk';
end
lhs     = events.lhs;
rhs     = events.rhs;
lto     = events.lto;
rto     = events.rto;
switch loc_mode
    case 'walk'
        event_list  = [lhs;rto;rhs;lto];
        lhs_code    = ones(length(lhs),1)*1;
        rto_code    = ones(length(rto),1)*2;
        rhs_code    = ones(length(rhs),1)*3;
        lto_code    = ones(length(lto),1)*4;
        code        = [lhs_code;rto_code;rhs_code;lto_code];
        
        % make one data matrix, sort,
        data            = [code event_list];
        [data(:,2),ind] = sort(data(:,2));
        data(:,1)       = data(ind,1);
        
        % throw away everything before first left hs, and after last lto
        ind_begin   = find(data(:,1)==1,1,'first');
        ind_end     = find(data(:,1)==4,1,'last');
        data        = data(ind_begin(1):ind_end(end),:);
        
        % and rewrite the lhs lto etc :)
        events.lhs  = data(data(:,1)==1,2);
        events.rto  = data(data(:,1)==2,2);
        events.rhs  = data(data(:,1)==3,2);
        events.lto  = data(data(:,1)==4,2);
        
        % check if same amounf of each event, and if order is correct.
        flag = 0;
        if ~(length(events.lto)==length(events.rto)&&length(events.rhs)==length(events.rto)&&length(events.lhs)==length(events.rto))
            flag = 1;
        elseif any((events.rto-events.lhs)<0)||any((events.rhs-events.rto)<0)||any((events.lto-events.rhs)<0)||any((events.lhs(2:end)-events.lto(1:end-1))<0)
            flag = 2;
        end
    case 'run'
        event_list=[lhs;lto;rhs;rto];
        lhs_code=ones(length(lhs),1)*1;
        lto_code=ones(length(lto),1)*2;
        rhs_code=ones(length(rhs),1)*3;
        rto_code=ones(length(rto),1)*4;
        code=[lhs_code;lto_code;rhs_code;rto_code];
        % make one data matrix, sort,
        data=[code event_list];
        [data(:,2),ind]=sort(data(:,2));
        data(:,1)=data(ind,1);
        
        % throw away everything before first left hs, and after last lto
        ind_begin=find(data(:,1)==1,1,'first');
        ind_end=find(data(:,1)==4,1,'last');
        data=data(ind_begin:ind_end,:);
        % and rewrite the lhs lto etc :)
        events.lhs=data(data(:,1)==1,2);
        events.lto=data(data(:,1)==2,2);
        events.rhs=data(data(:,1)==3,2);
        events.rto=data(data(:,1)==4,2);
        flag = 0;
        if ~(length(events.lto)==length(events.rto)&&length(events.rhs)==length(events.rto)&&length(events.lhs)==length(events.rto))
            flag = 1;
            %         elseif any((events.rto-events.lhs)<0)||any((events.rhs-events.rto)<0)||any((events.lto-events.rhs)<0)||any((events.lhs(2:end)-events.lto(1:end-1))<0)
            %             flag = 2;
        end
end
