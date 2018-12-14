


clear
%make sure matlib is in the path
addpath(genpath('Apples v2'));

%first create group names. You can add multiple allStudyGroups but for this
%particular script I will stick with one. 
allStudyGroups = {{'Apples_%g', [1:6]}}; %expects the raw data to have file format Apples_YC0*.mat
groupName = {'AD'};
% Move data into an array (one cell for on and off), each row = 1 subject
%loop through allStudyGroups
for subs = 1:length(allStudyGroups)
    alldata=[];
    %loop through subjects
    for i = 1:length(allStudyGroups{subs}{2})
        data = load(sprintf(allStudyGroups{subs}{1},i));
        
        d=data.result.data; %retrieve trial data
        for t=126:135          %for the last 10 trials collecting force data
            AUC1{subs}(i,t) = nansum(d(t).data); % extract AUC data before removing data and 2 fields
            
            if d(t).Yestrial ==1 % if accepted
                if ~isnan(d(t).reward) % if data recorded
                    mF{subs}(i,t) = max(d(t).data(1:500));%500 is trial length
                else mF{subs}(i,t)=nan;
                end
            else mF{subs}(i,t)=nan;
            end
            
        end
        
        d=rmfield(d, 'data'); % don't process the actual squeeze data
        
        if isempty(alldata) % first subject: new structure
            alldata = d;
        else % subsequent subjects - make sure fields match
            d=ensureStructsAssignable(d, alldata);
            alldata=[alldata;d]; % and then add a new row
        end
    end
    d=transpIndex(alldata);
    alld{subs}=d;
    
    end
    
   
    
   D=alld{subs};
   
    save('AGT_DATA_AD','D');
    save('mF','mF')
    save('AUC1','AUC1')


%%
if ~exist('D','var') % if above has not been run earlier load data that is already in folder
    load AGT_DATA_Cam13
end

alld=D;
filter = allStudyGroups{1}{2}'; % select appropriate patients
subj = length(filter);

load AUC1;
AUC1{1} = AUC1{1}(filter,:);

%% This section loads lars questionnaire data and classifies patients as apathetic according to some threshold

%- can be adpated to load other questionnaire data for later
%correlations

% load lars_self_PD;load lars_CG_PD;load aes_PD;
% lars_self = lars_self(filter,:);lars_CG = lars_CG(filter,:);aes=aes(filter,:);
% subj = length(allStudyGroups{1}{2});
% lars=[]; %use this as main vector for apathy
% apVec=[]; % categorise as apathetic (1) or not
% cut=-17;
% for i=1:subj
%     if 0 % if you want to use CG lars, otherwise use self.
%         if ~isnan(lars_CG(i,5))
%              lars(i)=lars_CG(i,5);
%         else lars(i)=lars_self(i,5);
%         end
%     else lars(i)=lars_self(i,5);
%     end
%     if 1  % **************** CHOSEN MEASURE ******************
%         if lars(i)>-22
%             apVec(i)=1;
%         else apVec(i)=0;
%         end
%     elseif 0 % if you want to use aes
%         if aes(i) > 36
%             apVec(i)=1;
%         else apVec(i)=0;
%         end
%     else % if you want to take combo approach
%         if lars_self(i,5)>cut && lars_CG(i,5)>cut
%             apVec(i)=1;
%         elseif lars_self(i,5)>cut && isnan(lars_CG(i,5))
%             apVec(i)=1;
%         elseif lars_self(i,5)<(cut+1) && lars_CG(i,5)<(cut+1)
%             apVec(i)=0;
%         elseif lars_self(i,5)<(cut+1) && isnan(lars_CG(i,5))
%             apVec(i)=0;
%         elseif aes(i) > 36
%             apVec(i)=1;
%         else apVec(i)=0;
%         end
%     end
% end
% lars=lars';
% apVec=apVec';

%% Create MVC matrix
for i = 1:length(allStudyGroups{1}{2}) % for each sub
    mvc{1}(i,1)=max(D.MVC(i,:,1)); %extract the max MVC for the subject (from 10 trials at the end)
end
% and now reshape data into arrays

for i=1:size(D.R,1) % for each subject in each state
    es =  [0.16 0.32 0.48 0.64 0.80]; %effort levels
    ss =  [1 4 7 10 13];%stake
    
    
    y = d.Yestrial(i,:)'; % did they accept?
    y = y(1:125); %MV we only need the first 125 trials as last 10 are not decisions
    y(isnan(y))=0; % compensate for someone's incompetence
    dt{subs} = d.endChoice-d.startChoice; %MV decision time
    dt{subs} = dt{subs}(:,1:125); %dt{subs}(:,37:end) %exclude practice session
    maxF = d.maximumForce(:,126:135); %MV we only need last 10 trials
    AUC_R = AUC1{subs}(:,126:135);%right % AUC
    % added by YS on 23/8/2018
    yesLoc = d.yeslocation(i,1:125)';
   % H = d.hand(:,126:135); %pretty irrelevant as we only have 1 hand


    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Normalize maxF & AUC by each subject's MVC (per hand)
    %MV we can only do this for last 10 trials
    for k=1:10  %(for each trial)
        if ~isnan(maxF(i,k)) % if they had to squeeze
            
                maxFNorm(i,k) = maxF(i,k)./d.MVC(i,k,1);
                AUCNorm(i,k) = AUC_R(i,k)./d.MVC(i,k,1);
          
        elseif isnan(maxF(i,k))
            maxFNorm(i,k) = nan;
            AUCNorm(i,k) = nan;
        end
    end
    for effort = 1:5 % get the frequency "map"
        for reward = 1:5
            filter = (d.stake (i,1:125) == ss(reward)) ...
                & (d.effort(i,1:125) == es(effort));
            filter2 = (d.stake (i,126:135) == ss(reward)) ...
                & (d.effort(i,126:135) == es(effort));
            
            freqmap(effort,reward,i) = nanmean( y(filter) );
            if sum(filter)==5
                choicemap(effort, reward, i, :) = y(filter);
                decisiont(effort, reward, i, :) = dt{subs}(i, filter);
                yesMap(effort,reward,i,:) = yesLoc(filter);
            else
                % there are not 5 trials for this condition for this subject
                warning('bad data sub %g eff %g sta %g has %g', i, effort, reward, sum(filter));
                if sum(filter)==4  % if there were only 4 matching trials,
                    choicemap(effort, reward, i, :) = [nan; y(filter)]; % 4: add a nan for block 1
                    decisiont(effort, reward, i, :) = [nan; dt{subs}(i, filter)'];
                    
                    
                    
                    
                end
            end
        end
    end
    end
    % And now remove all trials where decision time was less than certain
    % value (e.g. 400ms) - replace them with nan
    for effort = 1:5
        for reward = 1:5
            for block = 1:5
                if decisiont(effort,reward,i,block) <= 0.4
                    choicemap(effort,reward,i,block) = nan;
                    maxF(effort,reward,i,block) = nan;
                    decisiont(effort,reward,i,block) = nan
                    yesMap(effort,reward,i,block) = nan;
                end
            end
        end
    end

grpD.freqmap{subs}=freqmap; % Note this has not had erroneous trials removed %%
grpD.choicemap{subs}=choicemap; % ALL_CHOICEMAP { subs } ( EFFORT, REWARD, SUBJECT, BLOCK )
grpD.decisiont{subs}=decisiont;
grpD.maxF{subs}=maxF;
maxFData{subs}=maxFNorm;
grpD.yesMap{subs}=yesMap;

choices = nanmean(grpD.choicemap{1},4);

D.Yestrial=D.Yestrial(:,1:125);
D.Yestrial(isnan(D.Yestrial))=0;
D.stake = D.stake(:,1:125);
D.effort = D.effort(:,1:125);

save('AGT_grpD_AD')
