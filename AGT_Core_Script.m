%% This will be the first version of the analysis script for the new AGT_CP 
%written by Youssuf Saleh September 24th 2018. 

%
% There are 4 phases to the AGT experiment
% 1  calibration - calculate their max grip strength (MVC) initial squeeze
%    and then 2 attempts at the yellow line (=110% then 105% MVC)
% 2.a. Physical Familiarisation - 2 rounds of squeezing at each force level to
% familiarise oneself with effort required
% 2.b. Decision Practice - 5 trials of just making decisions to get a feel
% for the task, especially to get used to the fact that yes/no randomly
% change in orientation after each trial. This can be confusing. 
% 3. Decisions - Self-paced (10s time out) choices of effort for reward, 
%    no squeezing, 5 effort x 5 reward levels, 25 trials x 5 blocks.
% 4. Execute 10 trials selected randomly (but actually fixed) from the
%    choices in part 3. Forced squeezing required for trials that were
%    accepted.

%%% SCRIPT structure: 
%1. Questionnaire correlations 
%2. 

%git test

clear
close all
load AGT_grpD % This should be generated by the last line of the script ...
% that is in the AGT2grpD script. THe purpose of that script is to reshape
% all the data into matrices that can be manipulated, plotted and then
% finally analysed using a Generalised mixed effects model. 

subj = size(D.R,1); %How many subjects?

% Create a binary vector which classifies all patients as apathetic('1')...
% or non-apathetic('0'). There are several ways to classify patients as
% apathetic. One way to do this is to use the LARS alone, 
% and another is to use a combined grouping as I have done below. Here I
% have used either LARS > -22 OR AES score > 37. If using the LARS only, we
% have previously used -22 as a cut off rather than -16 which is used
% clinically as it seems patients demonstrate apathetic behaviour earlier. 
 
% I have created a test vector ...
%called larsTest so the script can work when we dont have the LARS
%questionnaire data and we want to test the script anyways.
%Otherwise use a variable in place of larsTest(I use larsT for
%total scores and larsSub for my subscores).

load larsT aesVec larsSub % these are variable names I have created to ...
% represent the subjects' total Lars(larsT), AES (aesVec) and lars
% Subscores. Do feel free to ise your own variable names but effectively
% they are matlab vectors which I have saved. We will be hopefully adding
% more questionnaire data to this (AMI,GDS, etc with the new core
% protocol). 
larsT=larsTest

apVec=[]';
for i = 1:subj
  if AMI_T(i) > 1.91
    %if larsTest(i)>-22 %|| aesVec(i) > 37
    apVec(i)=1;
  else apVec(i)=0;
  end
end
apVec = apVec';

% draw a plot of aesVec and LARS correlations
if 0
    figure()
    h = scatterRegress(AMI_T,aesVec,'*','MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',5);
    h = ax;
    set(gca,'LineWidth',2)
    xlim([-36 11])
    ylabel('Apathy evaluation scale')
    xlabel('Lille apathy rating scale')
    title('Correlation between LARS and AES in 53 SVD patients')
    set(gca,'fontSize',16,'fontWeight','bold');hold on;
    plot(-22*ones(46,1),15:60,':k')
end

apVec = apVec';
%apVec(9)=1; %aesVec for this subject 40
if 1 % if generating MRI inputs therefore need to exclude subj 2 and 19
    mriApVec = apVec;
    
end

%set yo choice data so that it can be used for analyses. 

choices = nanmean(grpD.choicemap{1},4); %all subjects average choices with accidental squeezes removed (nan)
dt = D.endChoice - D.startChoice; % decision time info
decisionTime=dt; % exclude accidental squeezes below from this matrix


%Lets not look at offers Accepted 
yesTrial = D.Yestrial; % manipulate the new variable yesTrial

% this is to remove all trials where a decision is made accidentally. 
% this is defined as trials which are performed in less than 0.4s. 
for i=1:subj
    for t=1:125
        if dt(i,t)<0.4  % if squeezed accidentally
            yesTrial(i,t)=nan;
            decisionTime(i,t)=nan;
        end
    end
end

% in these matrices, every '1' is an accepted offer '0 is rejected and nan
% refers to an accidental squeeze. 
reward = D.stakeIx(:,1:125);

for i=1:subj
    grpD.accept(i) = length(find(yesTrial(i,:)==1));
    grpD.reject(i) = length(find(yesTrial(i,:)==0));
    grpD.mistake(i)= length(find(isnan(yesTrial(i,:))));
    grpD.trialsCorr(i)= length(find(~isnan(yesTrial(i,:)))); % how many trials after removing mistakes
    grpD.fail(i) = length(find(yesTrial(i,:)==1 & reward(i,:)==0));
    grpD.failCorr(i) = length(find(yesTrial(i,:)==1 & reward(i,:)==0))/grpD.trialsCorr(i);
    grpD.failHighEff(i) = length(find(yesTrial(i,:)==1 & reward(i,:)==0 & D.effort(i,:)>0.6))/grpD.trialsCorr(i); % this metric erroneous as divides by ALL trials
    grpD.failHighEff2(i) = length(find(yesTrial(i,:)==1 & reward(i,:)==0 & D.effort(i,:)>0.6));
end

accept=grpD.accept';

% And plot this
figure()
bar(1,mean(accept(apVec==0)),0.5);hold on;bar(2,mean(accept(apVec==1)),0.5);
errorbar([mean(accept(apVec==0)) mean(accept(apVec==1))],[std(accept(apVec==0))/sqrt(length(find(apVec==0))) ...
    std(accept(apVec==1))/sqrt(length(find(apVec==1)))],'r.','LineWidth',2);
ylabel('Raw number of offers accepted')
xlabel('apathy status (no/yes)')
title('apathetic patients accepted fewer total offers')
set(gca,'fontSize',14,'fontWeight','bold')
% and also look at cumulative offers accepted

% and plot ~rate of acceptance across experiment for 2 groups
sm=15;%what is the smoothing kernal
figure()
errorBarPlot(smoothn(yesTrial(apVec==0,:),sm),'b','LineWidth',1.5);hold on
errorBarPlot(smoothn(yesTrial(apVec==1,:),sm),'r','LineWidth',1.5);
xlim([0 179])
legend('No Apathy','Apathy')
xlabel('Trial number')
ylabel('Rate of Acceptance')
title('Smoothed (36) response rate for 2 groups across entire experiment')


if 0 % draw some other graphs of correlations.
    figure()
    scatterRegress(AMI_T,accept);
    xlabel('Increasing apathy severity --> (aesVec)')
    ylabel('Number of offers accepted');
    title('aesVec & Offers accepted');
    hold off
    figure()
    scatterRegress(larsSub(:,4),accept);hold on;
    xlabel('Increasing severity on LARS Action Initiation subscale')
    ylabel('Number of offers accepted');
    title('Lars-Action Initiation & Offers accepted');
    figure()
    scatterRegress(larsT,accept);hold on;
    xlabel('Increasing apathy severity --> LARS ')
    ylabel('Number of offers accepted');
    title('Lars-self rated & Offers accepted');
   
end


%% *************** RESULTS - Offers accepted - RAW **********************

% ************* Basic results - number of apples gathered, offers accepted
% etc******************************
close all
numT=125; % use this to change acceptances to a proportion, if want raw just make it 1.
aaa=accept/125;

figure();hold on;

bar(1,mean(aaa(apVec==0,1)),'b');hold on;errorbar(1,mean(aaa(apVec==0,1)),std(aaa(apVec==0,1))./sqrt(length(find(apVec==0))),'r.','LineWidth',3)
bar(2,mean(aaa(apVec==1,1)),'r');hold on;errorbar(2,mean(aaa(apVec==1,1)),std(aaa(apVec==1,1))./sqrt(length(find(apVec==1))),'b.','LineWidth',3)
if 0
  
  plot(1*ones(length(accept(apVec==0,1)),1),accept(apVec==0,1),'k.')
  plot(2*ones(length(accept(apVec==1,1)),1),accept(apVec==1,1),'k.')
end
xlim([0.5 2.5])

ax=gca;
set(ax,'XTick',[1 2],'fontWeight','bold','fontSize',20,'ylim',[.5 1],'YTick',[.5:.1:1],'YTickLabel',{'0.5','0.6','0.7','0.80','0.9','1.0'},'XTicklabel',{'No Apathy','Apathy'});
title('Mean Offers Accepted')

ylabel('Proportion offers accepted (%)')
hold off


% ********** Cumulative Acceptance ********************
cumAccept=(zeros(subj,1));
for i = 1:subj
    for t=1:125
        if yesTrial(i,t)==1
            cumAccept(i,t+1) = cumAccept(i,t) + 1;
        else cumAccept(i,t+1) = cumAccept(i,t);
        end
    end
end

figure();hold on
if 0 % If want to include HC individual cum accept plots - actually this makes things messier so don't use
% for i=1:subj_HC
%     plot(smoothn(cumAccept_HC(i,:)),'m:','LineWidth',3)
% end
end
for i = 1:subj
    if apVec(i)==0
        plot(smoothn(cumAccept(i,:)),'b','LineWidth',3)
    else
        plot(smoothn(cumAccept(i,:)),'r','LineWidth',3)
    end
end
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[0:20:125],'YTick',[0:20:125]);

xlabel('Trial Number');ylabel('Cumulative trials accepted');


% ****************   And plot response rate *****************

sm=36;%what is the smoothing kernal
figure()
%errorBarPlot(smoothn(2,yesTrial_HC,sm),'m','LineWidth',4);
hold on
errorBarPlot(smoothn(yesTrial(apVec==0,:),sm),'b','LineWidth',3);
errorBarPlot(smoothn(yesTrial(apVec==1,:),sm),'r','LineWidth',3);
xlim([18 162]);ylim([0.4 1.2])
legend('No Apathy','Apathy')
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'YTick',0.5:0.1:1);
xlabel('Trial number')
ylabel('Acceptance rate (smoothed)')
title('Smoothed (36) response rates')



% ****************  And plot correlation between lars_AI and acceptance *****************
figure();hold on
scatter(larsSub(:,3),accept/125,'MarkerEdgeColor',[0.5 0.5 0.5],'LineWidth',5)
ax=gca;
lsline(ax);
h = lsline(ax);
set(h,'LineWidth',3,'Color',[0.7 0.7 0.7])

ylabel('Proportion offers accepted')
xlabel('Action initiation subscale of LARS')
set(ax,'fontWeight','bold','FontSize',20,'YTick',0.2:0.2:1)
xlim([-4.5 1.5]);ylim([0.2 1.1])
axis square

if 0
figure()
errorBarPlot(cumAccept(apVec==0,:),'LineWidth',1);
hold on; errorBarPlot(cumAccept(apVec==1,:),'LineWidth',1);
xlim([0 200]); ylim ([0 200]);
title('cumulative offers accepted across experiment')
legend('No Apathy n=34','Apathy n=19')
end

%% *********** RESULTS - EFFECTS OF REWARD AND EFFORT *************
% Plot raw choice proportions in expanded form
close all
figure()
for i=1:5
    errorBarPlot(squeeze(choices(i,:,apVec==0))','--','LineWidth',5);hold on
end
set(gca,'ColorOrderIndex',1);
for i=1:5
    errorBarPlot(squeeze(choices(i,:,apVec==1))',':','LineWidth',5);hold on
end

%%
% generate the RM-ANOVA dataset - use the ARCSINE transformed dataset for
% thesis
temp=[]; % variable to store the choice proportions
tempArc=[];
for i=1:53 % each subject
    temp(i,:)=reshape(choices(:,:,i)',36,1);
    tempArc(i,:)=asin(temp(i,:));
end

%% **************** 2D plots ***********************
% using just errorbar function to avoid difficulties with errorBarPlot
close all

figure()
dat = squeeze(mean(choices,2))';
%dat_hc=squeeze(mean(choices_HC,2))';
%errorbar(1:6,mean(dat_hc),std(dat_hc)./sqrt(19),'--','Color',[0.7 0.7 0.7],'LineWidth',3); hold on;
errorbar(1:5,mean(dat(apVec==0,:)),std(dat(apVec==0,:))./sqrt(length(find(apVec==0))),'b','LineWidth',3); 
hold on
errorbar(1:5,mean(dat(apVec==1,:)),std(dat(apVec==1,:))./sqrt(length(find(apVec==1))),'r','LineWidth',3)
title('proportion of offers accepted as effort level increases')
axis square
ylim([0 1.1]);xlim([0 7])
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1:1:6],'XTickLabel',{'10','24','38','52','66','80'})
xlabel('Effort level (% MVC)')
hold off
legend({'SVD No Apathy','SVD Apathy'});
figure()
dat = squeeze(mean(choices,1))';
%dat_hc=squeeze(mean(choices_HC,1))';
%errorbar(1:6,mean(dat_hc),std(dat_hc)./sqrt(19),'--','Color',[0.7 0.7 0.7],'LineWidth',3); hold on;
errorbar(1:5,mean(dat(apVec==0,:)),std(dat(apVec==0,:))./sqrt(length(find(apVec==0))),'b','LineWidth',3); hold on;
errorbar(1:5,mean(dat(apVec==1,:)),std(dat(apVec==1,:))./sqrt(length(find(apVec==1))),'r','LineWidth',3)
%title('proportion of offers accepted as reward level increases')
legend('SVD No Apathy','SVD Apathy');
axis square
ylim([0 1.1]);xlim([0 7])
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1:1:6],'XTickLabel',{'1','3','6','9','12','15'})
xlabel('Reward level')
ylabel('Proportion of offers accepted')

%% Apathy - No Apathy plots (raw difference)
%  3D difference plot (2D plot not amazing...
close all
choiceDif=(mean(choices(:,:,apVec==0),3)-mean(choices(:,:,apVec==1),3));
h=surf(choiceDif);shading('interp');hold on;colormap('jet');%colorbar('Ticks',0:.05:.2)
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1:1:5],'YTickLabel',{'10','24','38','52','66','80'},'YTick',[1:1:6],'XTickLabel',{'1','3','6','9','12','15'},'ZTickLabel',{'','0','0.1','0.2','0.3','0.4','0.5'})
title('3D plot NoAp vs. AP')
ylabel('Effort (%MVC)')
xlabel('Reward')
zlabel('Diff. Proportion accepted')
hold on;
base=zeros(5,5);
hh=surf(base);
hh.FaceColor=[0.5 0.5 0.5];hh.FaceAlpha=1;
view(25,30)
if 1 % if want to add on grid lines
    for i=1:5
        plot3(1:5,(i)*ones(5,1),choiceDif(i,1:5),'k:','LineWidth',2)
        plot3((i)*ones(5,1),1:5,choiceDif(1:5,i),'k:','LineWidth',2)
    end
end

%% Decision Time
% Use decisionTime and decisionTime_HC - have values <0.4 'nan'
close all
for i=1:length(find(apVec==0))
    vecnoA(i) = 1+((rand-0.5)/5);
end
for i=1:length(find(apVec==1))
    vecA(i) = 2+((rand-0.5)/5);
end

for i=1:19
    vecHC(i) = 1+((rand-0.5)/5);
end
figure() % First just plot mean dt for the 3 groups

%bar(1,mean(nanmean(decisionTime_HC,2)),'FaceColor',[0.7 0.7 0.7]);hold on;
bar(1,mean(nanmean(decisionTime(apVec==0,:),2)),'b');hold on
bar(2,mean(nanmean(decisionTime(apVec==1,:),2)),'r');
if 1
 %   plot(vecHC,nanmean(decisionTime_HC,2),'.','MarkerSize',20,'MarkerEdgeColor',[0.5 0.5 0.5]);hold on
    plot(vecnoA,nanmean(decisionTime(apVec==0,:),2),'.','MarkerSize',20,'MarkerEdgeColor',[0.5 0.5 0.5])
    plot(vecA,nanmean(decisionTime(apVec==1,:),2),'.','MarkerSize',20,'MarkerEdgeColor',[0.5 0.5 0.5])
end
%errorbar(1,mean(nanmean(decisionTime_HC,2)),std(nanmean(decisionTime_HC,2))./sqrt(19),'k','LineWidth',3);
errorbar(1,mean(nanmean(decisionTime(apVec==0,:),2)),std(nanmean(decisionTime(apVec==0,:),2))./sqrt(length(find(apVec==0))),'k','LineWidth',3);
errorbar(2,mean(nanmean(decisionTime(apVec==1,:),2)),std(nanmean(decisionTime(apVec==1,:),2))./sqrt(length(find(apVec==1))),'k','LineWidth',3)
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',1:2,'XTickLabel',{'SVD-noAp','SVDAp'})
ylabel('Decision time (s)')
xlim([0.5 2.5])
title('Average Decision time is Longer in apathy')
% Now Yes versus No trials  ****************
dtYes=[]; % 1st column = yes, 2nd column = no
%dtYes_HC=[];
filter=yesTrial;filter(isnan(filter))=0;filter=logical(filter);
filter2=yesTrial;filter2(filter2==1)=nan;filter2(~isnan(filter2))=1;filter2(isnan(filter2))=0;filter2=logical(filter2);
for i=1:subj % each subject
    dtYes(i,1)=mean(decisionTime(i,filter(i,:)));
    dtYes(i,2)=nanmean(decisionTime(i,filter2(i,:)));
end

ap=apVec;

figure()
%bar(1,nanmean(dtYes_HC(:,1)),'FaceColor',[0.7 0.7 0.7]);hold on;errorbar(1,nanmean(dtYes_HC(:,1)),nanstd(dtYes_HC(:,1))./sqrt(19),'k','LineWidth',3);
%bar(2,nanmean(dtYes_HC(:,2)),'FaceColor',[0.3 0.3 0.3]);hold on;errorbar(2,nanmean(dtYes_HC(:,2)),nanstd(dtYes_HC(:,2))./sqrt(19),'k','LineWidth',3)
h1 = bar(1,nanmean(dtYes(ap==0,1)),'FaceColor',[0.7 0.7 0.7]);hold on;errorbar(1,nanmean(dtYes(ap==0,1)),nanstd(dtYes(ap==0,1))./sqrt(34),'k','LineWidth',3); 
h2 = bar(2,nanmean(dtYes(ap==0,2)),'FaceColor',[0.3 0.3 0.3]);hold on;errorbar(2,nanmean(dtYes(ap==0,2)),nanstd(dtYes(ap==0,2))./sqrt(34),'k','LineWidth',3); 
h3 = bar(4,nanmean(dtYes(ap==1,1)),'FaceColor',[0.7 0.7 0.7]);hold on;errorbar(4,nanmean(dtYes(ap==1,1)),nanstd(dtYes(ap==1,1))./sqrt(19),'k','LineWidth',3);
h4 = bar(5,nanmean(dtYes(ap==1,2)),'FaceColor',[0.3 0.3 0.3]);hold on;errorbar(5,nanmean(dtYes(ap==1,2)),nanstd(dtYes(ap==1,2))./sqrt(19),'k','LineWidth',3);
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1.5 4.5],'XTickLabel',{})
ylabel('Decision time (s)');
xlim([0 6]);
ylim([0 3]);
xticklabels({'NoAP','Ap'});
legend([h1 h2],{'Yes','No'});
title('Decision time vs Trial Type')

% And finally just plot the difference - YES - NO ie a within subject
% comparison rather than average time plots
figure()
temp=dtYes(:,2)-dtYes(:,1);%tempHC=dtYes_HC(:,2)-dtYes_HC(:,1);
%bar(1,nanmean(tempHC),'m');hold on;errorbar(1,nanmean(tempHC),nanstd(tempHC)./sqrt(19),'k','LineWidth',3);
bar(1,nanmean(temp(apVec==0)),'b');hold on;errorbar(1,nanmean(temp(apVec==0)),nanstd(temp(apVec==0))./sqrt(34),'k','LineWidth',3);
bar(2,nanmean(temp(apVec==1)),'r');hold on;errorbar(2,nanmean(temp(apVec==1)),nanstd(temp(apVec==1))./sqrt(19),'k','LineWidth',3);
if 1
    %plot(vecHC,nanmean(tempHC,2),'.','MarkerSize',12,'MarkerEdgeColor',[0.5 0.5 0.5])
    plot(vecnoA,nanmean(temp(apVec==0,:),2),'.','MarkerSize',12,'MarkerEdgeColor',[0.5 0.5 0.5])
    plot(vecA,nanmean(temp(apVec==1,:),2),'.','MarkerSize',12,'MarkerEdgeColor',[0.5 0.5 0.5])
    if 1
        ylim([-1 1.5]);plot(3,-.93,'.','MarkerSize',12,'MarkerEdgeColor',[0.5 0.5 0.5]);
    end
end
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1 2],'XTickLabel',{'noAp','Ap'})
ylabel('Difference in decision time NO-YES (s)');


%% NOW look at how DT relates to value of offer - to try and show all groups sensitive to this.
close all
tempDec=nanmean(grpD.decisiont{1},4);
easyDec=[];
hardDec=[];
for i=1:subj % each subject
    temp=[];
    temp2=[];
    c=1; %counter
    cc=1;
    for e=1:5
        for r=1:5
            if choices(e,r,i) <=0.25 || choices(e,r,i) >=0.75
                temp(c,1)=tempDec(e,r,i);
                c=c+1;
            else
                temp2(cc,1)=tempDec(e,r,i);
                cc=cc+1;
            end
        end
    end
    easyDec(i,1)=nanmean(temp);
    easyDec(i,2)=nanmean(temp2);
end

for i=1:subj % each subject
    temp=[];
    temp2=[];
    c=1; %counter
    cc=1;
    for e=1:5
        for r=1:5
            if choices(e,r,i) <=0.25 || choices(e,r,i) >=0.75
                temp(c,1)=tempDec(e,r,i);
                c=c+1;
            else
                temp2(cc,1)=tempDec(e,r,i);
                cc=cc+1;
            end
        end
    end
   % easyDecHC(i,1)=nanmean(temp);
    %easyDecHC(i,2)=nanmean(temp2);
end
if 1 % if want to exclude al values if one NaN
    easyDec(12,1)=nan;easyDec(16,1)=nan;easyDec(17,1)=nan;
   % easyDecHC(12,1)=nan;easyDecHC(16,1)=nan;easyDecHC(17,1)=nan;
end
figure()
%errorBarPlot(easyDecHC,'Color',[0.7 0.7 0.7],'LineWidth',3);hold on
errorBarPlot(easyDec(apVec==0,:),'b','LineWidth',3);
hold on 
errorBarPlot(easyDec(apVec==1,:),'r','LineWidth',3);
legend('SVD No Apathy','SVD Apathy')
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1 2],'XTickLabel',{'Easy','Hard'})
ylabel('Decision time (s)');
ylim([1 4])
xlim([0.5 2.5])

%% ******************** Block Effects   *************************
% grpD.choiceMap has already had accidental choices removed (nan)
close all
figure();hold on
for i=1:5 % for each block
    errorBarPlot(squeeze(nanmean(grpD.choicemap{1}(:,:,apVec==0,i),1))','LineWidth',3);
end
legend('Block 1','Block 2','Block 3','Block 4','Block 5');
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1:6])
ylim([0.2 1.1]);xlim([0 7]);
figure();hold on
for i=1:5 % for each block
    errorBarPlot(squeeze(nanmean(grpD.choicemap{1}(:,:,apVec==1,i),1))','LineWidth',3);
end
legend('Block 1','Block 2','Block 3','Block 4','Block 5');
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1:6])
ylim([0.2 1.1]);xlim([0 7]);
legend('Block 1','Block 2','Block 3','Block 4','Block 5');
ax=gca;
set(ax,'fontWeight','bold','fontSize',20,'XTick',[1:6])
ylim([0.2 1.1]);xlim([0 7]);


%% Set up for statistics - GLME models
%% Hierachical linear mixed effects model - use fitglme -
% move data into subData - [subject choice rew eff ap lars DA]
... stake, effort and Yestrial within the D array should already have had the practise block removed (as long as above parts run)
... stake, effort and Yestrial within the D array should already have had the practise block removed (as long as above parts run)
subData=[];
subDataT=[];
linear=0;

    for i=1:subj % each subject
        choicesVec = D.Yestrial(i,:)';
        choicesVec(isnan(choicesVec))=0; %change nans to 0
        reward = D.stake(i,:)';
        effort = D.effort(i,:)';
        decVec = dt(i,1:125)'; %pick individual subjects' dec times
        decVec(isnan(decVec))=0;
        a=1;
        for j = 1:length(decVec) %for each of the 125 trials
            if decVec(j) < 0.4 | choicesVec(j)==2
                removal(a) = j; %create vector of trials to remove
                a=a+1;
            end
        end
        if exist('removal')
            choicesVec(removal)=[];
            reward(removal)=[];
            effort(removal)=[];
            decVec(removal)=[];
        end
        clear removal
        if ~linear
            subData = [subData;i*ones(length(choicesVec),1) choicesVec reward effort apVec(i)*ones(length(choicesVec),1)];
            subDataT = [subDataT;i*ones(length(decVec),1) decVec reward effort apVec(i)*ones(length(decVec),1)];
        elseif linear
            unEff=unique(effort);
            unRew=unique(reward);
            for e=1:6
                for r=1:6
                    filter = effort==unEff(e) & reward==unRew(r);
                    %tempD=mean(choicesVec(filter));
                    tempD=choices{k}(e,r,i);
                    subData = [subData;i tempD unRew(r) unEff(e) apVec(i) k-1];
                    tempDT=choices{k}(e,r,i);
                    subDataT = [subDataT;i tempD unRew(r) unEff(e) apVec(i) k-1];
                    
                end
            end
            
        end
    end

%subject = categorical(subData(:,1));
subject = categorical(subData(:,1));
choice  = subData(:,2);
rew     = subData(:,3);       
eff     = subData(:,4);
ap      = subData(:,5);
dTime   = subDataT(:,2);
%lars    = subData(:,6);

if 0
    % ****** IF WANT Quadratic effort *******
    eff = eff.^2;
end
if 1
    rew=zscore(rew);
    eff=zscore(eff);
    ap = zscore(ap);
    %lars = zscore(lars);
end
if 0
    rew = categorical(rew);
    eff = categorical(eff);
    ap = categorical(ap);
    DA = categorical(DA);
    
    
end
if 0 % adjust each reward level based on group level random effect of reward (debatable...)
    load randRew % the estimates for reward effects to apply to the basic rewards
    unRew=unique(rew)
    
    for r = 1:6 % each reward level
        rew(rew==unRew(r))=unRew(r)+randRew(r);
        
    end
end
if linear
    if 1 % arcsine transform the data -  T = asin(sqrt(p))
        choice = asin(sqrt(choice));
    end
end
Design = table(choice,rew,eff,ap,subject);
Design1 = table(dTime,rew,eff,ap,subject);


clear aic glme_fit bic
linear =0; % which model type to run
models = {    
    

'choice ~ rew*eff*ap + (1|subject)'
'choice ~ rew*eff*ap + (1+rew| subject)' 
'choice ~ rew*eff*ap + (1+rew| subject)' 
'choice ~ rew*eff + rew*ap + eff*ap + (1|subject)'

% to generate reward and effort sensitivities. 
'choice ~ rew+eff + (rew+eff|subject)'

    };

if linear
    for i=1:length(models)
        sprintf('starting model %g',i)
        glme_fit{i}=fitglme(Design,models{i},'Distribution','normal');
        aic(i)=glme_fit{i}.ModelCriterion.AIC;
        bic(i)=glme_fit{i}.ModelCriterion.BIC;
    end
else
    for i=1:length(models)
        sprintf('starting model %g',i)
        glme_fit{i}=fitglme(Design,models{i},'Distribution','binomial','fitmethod','Laplace');
        aic(i)=glme_fit{i}.ModelCriterion.AIC;
        bic(i)=glme_fit{i}.ModelCriterion.BIC;
                

    end
end
aic=aic-min(aic);
bic=bic-min(bic);

if 0
    save('glme_fit_lin_cat','models','glme_fit')
end
if 0 % save outputs
save('glme_fitPD_Z','models','glme_fit')
end


% run the GLME with the two main outputs being B and BNames which
% correspond to the random effects parameters and their corresponding names
% in table form respectively. 
% Then run the GLME with the two main outputs being B and BNames which
% correspond to the random effects parameters and their corresponding names
% in table form respectively. 
model = { ...
    'choice ~ 1 + rew + eff + (rew+eff|subject)'
    };

% Then run the GLME with the two main outputs being B and BNames which
% correspond to the random effects parameters and their corresponding names
% in table form respectively. 
i = 1;
[B,BNames] = randomEffects(fitglme(Design,model{i},'Distribution','binomial','fitmethod','Laplace'));

% The output of this is a 159*1 matrix and a table whose dimensions is
% 159*3. The output format contains three random effects parameters per
% subject (53*3 = 159) and arranged so that each subject has 3 values
% (intercept,reward,effort) before moving onto the next subject. 
% I will now create three vectors to represent our variables of interest. 

%Intrinsic Motivation or int represents the intercept variation per
%subject. index into all intercept values.  
int = B(1:3:end);
% reward sensitivity encoded as rewSen and indexes into all reward
% parameters
rewSen = B(2:3:end);
%effort sensitivity encoded as effSen 
effSen = B(3:3:end);


M = fitglme(Design,model{i},'Distribution','binomial','fitmethod','Laplace');
fe=M.fixedEffects;
close all

subplot(2,2,1)

bar(1,fe(1)+(mean(int(apVec==0))),0.5);hold on;bar(2,fe(1)+(mean(int(apVec==1))),0.5);
errorbar(fe(1)+[(mean(int(apVec==0)))  (mean(int(apVec==1)))],...
    [std(int(apVec==0))/sqrt(length(int(apVec==0))) ...
    std(int(apVec==1))/sqrt(length(int(apVec==1)))],'k.','LineWidth',2);
ylabel('Parameter estimate')
xlabel('apathy status (no/yes)')
xticks([1 2]);
xticklabels({'noAp','Ap'})
title('Intrinsic Motivation')
set(gca,'fontSize',14,'fontWeight','bold')
 
 % reward Sensitivity parameter. 
 subplot(2,2,2)
 bar(1,fe(2)+(mean(rewSen(apVec==0))),0.5);hold on; bar(2,fe(2)+(mean(rewSen(apVec==1))),0.5);
 errorbar(fe(2)+[(mean(rewSen(apVec==0))) (mean(rewSen(apVec==1)))],...
     [std(rewSen(apVec==0))/sqrt(length(rewSen(apVec==0))) ...
     std(rewSen(apVec==1))/sqrt(length(rewSen(apVec==1)))],'k.','LineWidth',2);
 ylabel('Parameter estimate')
 xlabel('apathy status (no/yes)')
 xticks([1 2]);
 xticklabels({'noAp','Ap'})
 title('reward sensitivity')
 set(gca,'fontSize',14,'fontWeight','bold')
 
 
 % effort Sensitivity parameter. 
 subplot(2,2,3)
 bar(1,fe(3)+(mean(effSen(apVec==0))),0.5);hold on; bar(2,fe(3)+(mean(effSen(apVec==1))),0.5);
 errorbar(fe(3)+[(mean(effSen(apVec==0))) (mean(effSen(apVec==1)))],[std(effSen(apVec==0))/sqrt(length(effSen(apVec==0))) ...
     std(effSen(apVec==1))/sqrt(length(effSen(apVec==1)))],'k.','LineWidth',2);
 ylabel('Parameter estimate')
 xlabel('apathy status (no/yes)')
 xticks([1 2]);
 xticklabels({'noAp','Ap'})
 title('effort sensitivity')
 set(gca,'fontSize',14,'fontWeight','bold')

hold off 
hold off 
hold off

%% Can we now compute subjective value for every subject for every trial? 
% this is computed by V(r,e) = int - effSen*Eff^2 + rewSen*rew
clear sub v i;
% fixed effects
zrew = (unique(rew))';
zeff = (unique(eff))';
fe=M.fixedEffects; % these are constant across subjects
for sub = 1:subj
    for r = 1:5
        for e = 1:5
           
            v(sub,r,e) = fe(1) + int(sub) ...
                + (fe(2) + rewSen(sub))*zrew(r)' ...
                + (fe(3) + effSen(sub))*(zeff(e));
            
        end
    end
end
pp=1./(1 + exp(-v));
subplot(1,2,1)
imagesc(measn(pp),2)
colorbar
subplot(1,2,2)
%errorBarPlot(pp)
plot(pp')
%% now can we fit the chosen model to the data and compare modelled to raw?
clf
PLOT=3;
for i=1:subj
    subplot(4,4,i)
    if PLOT==1
        plot(sq(pp(i,:,:)));
        hold on
        set(gca,'colororderindex',1);
        plot(choices(:,:,i)','o:');
        hold off
    elseif PLOT==2
        imagesc(sq(pp(i,:,:)));
    elseif PLOT==3
        imagesc(sq(choices(:,:,i))');
    end
end

colormap(jet(256))
colorbar






