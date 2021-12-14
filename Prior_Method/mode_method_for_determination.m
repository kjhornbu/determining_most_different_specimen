    %Determine the "most change" specimen ID -- mode method of most
    %different for each contrast we are not doing that now. 
%     if isempty(o) %if no regions signficant in the given contrast just move on
%     elseif o==1 %if only one signficant region in the contrast;
%         Runno_low_mode_idx(n,:)=squeeze(Runno_idx_low(n,1,:));
%         Runno_high_mode_idx(n,:)=squeeze(Runno_idx_high(n,1,:));  
%     else
%         [~,~,low_idx]=unique(squeeze(Runno_idx_low(n,1:o,:)),'rows');
%         Runno_low_mode_idx(n,:)=Runno_idx_low(n,find(low_idx==mode(low_idx),1),:); %per contrast basis of mode finding -- sort of not a good at find the actual most specimen
%         
%         [~,~,high_idx]=unique(squeeze(Runno_idx_high(n,1:o,:)),'rows');
%         Runno_high_mode_idx(n,:)=Runno_idx_high(n,find(high_idx==mode(high_idx),1),:); %per contrast basis of mode finding -- sort of not a good at find the actual most specimen
%     end