function [gc_red_surr,gc_syn_surr]=CTE_surr(x,p,d, drivers_red,drivers_syn,nsurr, souce_index, target_index,condition_index,phase) 
%i -> j
gc_red_surr = zeros(nsurr,size(x,2));
gc_syn_surr = zeros(nsurr,size(x,2));
for dr=condition_index% 3:kmax+2
    xx=x(:,drivers_red(1:dr-1));
    for k=1:nsurr
        ys=surr_iaafft(x(:,drivers_red(dr)));
%         disp(std(ys))  % 看一下是否是常数
        xx(:,dr)=ys;
        xxx=xx;
        gc=CTE_iteration(xxx,p,d,souce_index,target_index,phase);
        Gs=mean(gc);
        gc_red_surr(k,dr)=Gs;

    end
end

for dr=condition_index %3:kmax+2
    xx=x(:,drivers_syn(1:dr-1));
    for k=1:nsurr
        ys=surr_iaafft(x(:,drivers_syn(dr)));
        xx(:,dr)=ys;
        xxx=xx;    
        gc=CTE_iteration(xxx,p,d,souce_index,target_index,phase);
        Gs=mean(gc);
        gc_syn_surr(k,dr)=Gs;
    end
end
