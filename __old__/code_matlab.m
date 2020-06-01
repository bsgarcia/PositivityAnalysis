%% Model Free Analysis of the Impact of Past Reward and their Interactions on the Current Choice 


%% Load and Compute Behavioural Variables
% _______________________________________

    
    learning_data = load('kentaro.csv');
    
    nsub = 0;
    
    subjects = unique(learning_data(:,2))';
    
    for s = subjects
        
        nsub = nsub +1;
        
        cho{nsub} = learning_data(learning_data(:,2)  == s &...  subjects
                                  learning_data(:,9)  == 1 &...  valence     (positive)
                                  learning_data(:,10) == 0,...   information (factual only)
                                  11); % choice
                              
        out{nsub} = learning_data(learning_data(:,2)  == s &...  subjects
                                  learning_data(:,9)  == 1 &...  valence     (positive)
                                  learning_data(:,10) == 0,...   information (factual only)
                                  12); % outcome 
                              
        con{nsub} = learning_data(learning_data(:,2)  == s &...  subjects
                                  learning_data(:,9)  == 1 &...  valence     (positive)
                                  learning_data(:,10) == 0,...   information (factual only)
                                  5);  % state (condition)
                              
        con{nsub} = (con{nsub} == max(con{nsub}))+1; % Recoding the condition 1,2 instead of 1,5
                              
    end


%% Compute Variable of Interest for the logistic Regression
% _________________________________________________________


% We need a variable that includes all choice triplet for which the two first ones are identical.



% Re-arranging Choices & Rewards by Condition
%--------------------------------------------

for s = 1 : numel(subjects)
    
    for c = unique(con{s})'
        
        choCon{s}{c}    = cho{s}(con{s}==c);
        outCon{s}{c}    = out{s}(con{s}==c);
        
        
    end
    
end

% Computing the Relevant Triplets
%--------------------------------

    
for s = 1 : numel(subjects)
    
    nTriplet = 0;
    
    for c = unique(con{s})'

        for t = 1 : length(choCon{s}{c}) - 2 % number of trial minus two represent the first trial of the last possible triplet
            
            if choCon{s}{c}(t) == choCon{s}{c}(t+1)
                
                nTriplet = nTriplet +1;
                
                % Factual Outcomes (raw)
                %-----------------
                
                r1{s}(nTriplet,1)   = outCon{s}{c}(t+1);    % R1 is the reward one trial  ago relatively to the trial of interest t+2
                r2{s}(nTriplet,1)   = outCon{s}{c}(t);      % R2 is the reward two trials ago relatively to the trial of interest t+2
                
                
                
                % Choices
                %--------
                
                c1{s}(nTriplet,1)   = choCon{s}{c}(t+1);    % Control variable
                c2{s}(nTriplet,1)   = choCon{s}{c}(t);      % Control variable
                
                if choCon{s}{c}(t+2) == choCon{s}{c}(t+1)
                    
                    stay{s}(nTriplet,1) = 1;
                    
                else
                    
                    stay{s}(nTriplet,1) = 0;
                    
                end
                
            end
            
        end
        
    end

    
end



%% Big Matrix all Subjects
%-------------------------


for sub = 1 : numel(subjects)
    
    if sub == 1
        matrix = [         r1{sub} r2{sub} stay{sub} repmat(sub,size(r1{sub}))];
    else
        matrix = [matrix ; r1{sub} r2{sub} stay{sub} repmat(sub,size(r1{sub}))];
    end
    
    
end


%cd('/Users/Germain/Dropbox/Post-doc/Confirmation Bias in the Brain/R statistics/GLM')

%csvwrite(strcat('modelFree_asymetry_all_',expName,bin,'.csv'),matrix)

%cd('/Users/Germain/Dropbox/Post-doc/Confirmation Bias in the Brain/Matlab Workspace')

















