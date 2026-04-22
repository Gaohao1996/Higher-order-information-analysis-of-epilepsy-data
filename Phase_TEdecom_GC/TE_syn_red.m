function [drivers_red,drivers_syn,g_red,g_syn]=TE_syn_red(x,i,j,p,d,phase)
%transfer entropy i -> j 
%p=order of the model
%d=delay between the variables

xxx=size(x);
drivers_red=[i j];
ind=setdiff(1:xxx(2),[i j]); %the varying set of condition variables
m=length(ind);
%% TE Methods
TE= CTE_iteration(x(:,[i j]),p,d,i,j,phase); 
g=mean(TE);
g_red(1)=g;g_syn(1)=g;
icont=1;

while m > 0
    icont=icont+1;
    for h=1:m
        dt=[drivers_red ind(h)];
        TE= CTE_iteration(x(:,dt),p,d,i,j,phase); 
        TE_mean=mean(TE);
        deltas(h)=TE_mean;
    end
    [g_red(icont), ii]=min(deltas);
    drivers_red=[drivers_red ind(ii)];
    ind=setdiff(ind,ind(ii));m=length(ind);
    clear deltas;
end
%%find the subset  ZM (drivers_syn) that maximize the conditional TE: TE(x→y|ZM)
drivers_syn=[i j];
ind=setdiff(1:xxx(2),[i j]);
m=length(ind);
icont=1;
while m > 0
    icont=icont+1;
    for h=1:m
        dt=[drivers_syn ind(h)];
        TE=CTE_iteration(x(:,dt),p,d,i,j,phase); 
        TE_mean=mean(TE);
        deltas(h)=TE_mean;
    end
    [g_syn(icont), ii]=max(deltas);
    drivers_syn=[drivers_syn ind(ii)];
    ind=setdiff(ind,ind(ii));m=length(ind);
    clear deltas;
end
