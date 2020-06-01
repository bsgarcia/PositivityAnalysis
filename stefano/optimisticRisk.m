function [act, rew, P] = optimisticRisk(params,trials,contingencies)
% [a r] = stationary(params,trials,states,contingencies)
% Q-learning silulator, for two-armed bandit task, with reciprocal
% probabilites


%contingencies(2,:)=1-contingencies(1,:);


beta  =params(1);
alpha1=params(2);
alpha2=params(3);


Q  = zeros(1,2);

revpoint=floor(trials/2);


t=0;


    
    for i=1:trials;
        
        t=t+1;
        
      
        P(t)=1/(1+exp((Q(1,1)-Q(1,2))/(  beta )));
        
        act(t)=1+(P(t)>rand);
        
        
        if act(t)==2
%             if t==1
%               rew(t)=0
%             else
            rew(t)=((contingencies>rand)-0.5)*2;
           % end
        else
            rew(t)=0;
          %rew(t)=((rand)-.5)*2;
        end
        
        PE = rew(t) - Q(1,act(t));
        
        Q(1,act(t))   = Q(1,act(t)) +  alpha1 * PE * (PE>0) +  alpha2 * PE * (PE<0) ;
        
        % normalyzed q-learning
        % Q(1,3-act(t)) = Q(1,3-act(t)) -   Unc(1) * alpha * PE;
        
        %Unc(1) = Unc(1) + gamma * (abs(PE) - Unc(1));
        
        
        
        
        
        
        
    end
    


end

