clear all

% "The role of pathogen mediated insect superabundance in the east-African emergence of a plant virus" Fig. 2 Fig. S1.1 Fig. S2.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMMENTS                                             % R. Donnelly 2021
%%% Script structured as follows
%%% Run three scenarios sequentially and plot each scenario as figure panel
%%% Within each scenario run odesolver for model in aggreg2_wCuttings.m
%%% run initially to whitefly steady-state
%%% run subsequently from steady-state to generate epidemic waves
%%%
%%% Subsequent runs have a stopping condition (midScapeEvents.m)
%%%  in order to produce comparable landscape snapshots for each scenario  
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% LIFE HISTORY PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
muu=1/50;        % vector mortality
a=1*(100*muu/2); % vector reproduction
K=40;            % carrying capacity
H=100;           % number plants in field
Pacq=0.032;      % rate acquire
Pinoc=0.032;     % rate inoculate
thet=1;          % dispersal rate
sig=0;           % assuming life long retention period
plE=1/360;       % plant mortality (exposed and susceptible)
plD=plE;         % addditional plant mortality e.g. roguing (diseased)
muuN=2*muu;      % vector mortality for nymphs

pMigrate=0.5;    % probability dispersal leaves focal field (to wider landscape)
pDisappear=0.5;  % probability dispersal that leaves focal field is lost from wider landscape

icbnMain=1/30;   % disease latency (30d exposed -> infected)
dNmain=1/25;     % vector development (25d nymph -> adult)

lh_pams=[a K thet H muu Pacq Pinoc plD sig pMigrate pDisappear muuN plE dNmain icbnMain];

numFieldsIn=120; % number of fields in landscape ring
tMaxSS=1000;     % initial steady-state run duration for solver
tMax=10000;      % epidemic run duration for solver (note, upr limit as there is stop condition)

% several additional parameters for infection through cuttings
propCleanSeed=0; % proportion of cuttings coming from completely virus-free source
cutFromWithin=1; % proportion of cuttings coming from within focal field
range=numFieldsIn-1; % number of fields for averaging cuttings infection pressure (if cuttings originate outside focal field)
is_discrim=1; % if growers are discriminating then cuttings might be exposed but not infected; otherwise exposed and infected
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE 1, to toggle between Fig 2 main text setup and SI2 cases:
% if propCleanSeed=1 then main text Fig 2 configuration regardless of cutFromWithin and range
% if propCleanSeed=0 and cutFromWithin=1; then SI2 case 1 and case 3
% if propCleanSeed=0 and cutFromWithin=0; and range=numFieldsIn-1;  then SI2 case 2 and case 4

% NOTE 2, to toggle between SI2 cases:
% 'aggreg2_wCuttings.m'
% is_discrim=1;     <--- discriminate cuttings (case 1 and 2 SI2)
% is_discrim=0;     <--- indiscriminate cuttings (case 3 and 4 SI2)
% above settings don't do anything when propCleanSeed=1 (i.e. as for main text results)

%%% SCENARIO 1 %%%
% additional scenario specific parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
eps1=20;                                                                   %%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
epsSP2=1;                                                                  %%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
epslscape=1;                                                               %%%% MODIFYING PATHOGEN %%%%%%%%%%%%%
sp2On=0;   % note SP2 switch must be on for invasive wf scenario... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MODIFYING PATHOGEN %%%%%%%%%%%%%



%%%%%% RUN first TO insect STEADY-STATE (no pathogen!) %%%%%%%%%%%%%%%%%%%%
invIndex=1; % field location of invader
initVals=[zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) K*(1-(muu/a))*ones(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn)];
[tnewSS,ynewSS] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,range,propCleanSeed,cutFromWithin,is_discrim,[lh_pams eps1 epsSP2 epslscape]),[0 tMaxSS],initVals);
disp('Scenario 1 ss run complete');
nymphsSS=ynewSS(:,0*numFieldsIn+(1:numFieldsIn));
insectsSS=ynewSS(:,3*numFieldsIn+(1:numFieldsIn));
nymphsSS2=ynewSS(:,11*numFieldsIn+(1:numFieldsIn));
insectsSS2=ynewSS(:,14*numFieldsIn+(1:numFieldsIn));
% TO VIEW TIME SERIES
%    figure;
%    plot(tnewSS,nymphsSS);
%    figure;
%    plot(tnewSS,insectsSS);
%    figure;
%    plot(tnewSS,nymphsSS2);
%    figure;
%    plot(tnewSS,insectsSS2);

% new initial conditions for epidemic spread
initValsFromSS=ynewSS(end,:);
initValsFromSS(7*numFieldsIn+invIndex)=1;    % <- disease invasion
initValsFromSS(14*numFieldsIn+invIndex)=sp2On*1; % <- insect invasion (susceptible)
initValsFromSS(3*numFieldsIn+invIndex)=initValsFromSS(3*numFieldsIn+invIndex)-sp2On*1; % <- simply reflecting invasion in loss of indiv. from WT

%%%%%% RUN second FOR EPIDEMIC WAVE
options = odeset('Events',@midScapeEvents);   % this option sets the condition for stopping the solver (i.e. when epidemic wave reaches a point on landscape)
[tnew,ynew,te1,ye,ie] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,range,propCleanSeed,cutFromWithin,is_discrim,[lh_pams eps1 epsSP2 epslscape]),[0 tMax],initValsFromSS,options);
disp('Scenario 1 epidemic run complete');
infVec=ynew(:,7*numFieldsIn+(1:numFieldsIn));
exposArrayF=ye(6*numFieldsIn+(1:numFieldsIn/2));
infArrayF=ye(7*numFieldsIn+(1:numFieldsIn/2));
expArrayF=exposArrayF/H;
incArrayF=infArrayF/H;
if sum(sum((incArrayF+expArrayF)>1.001)) % testing for errors from solver 
    disp('ERROR!!! ');
    for qq=1:numFieldsIn
        tmpfield=incArrayF(:,qq)+expArrayF(:,qq);
        inds=find(tmpfield>1.001);
        if ~isempty(inds)
            disp(['Field is ' num2str(qq)]);
            disp(['At time ' num2str(inds')]);
            disp(incArrayF(inds,qq)');
            disp(expArrayF(inds,qq)');
            return;
        end
    end
end
superAbundAF=H*(1-incArrayF-expArrayF).*(ye(3*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(4*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(:,5*numFieldsIn+(1:numFieldsIn/2)));%sp#1
superAbundBF=H*(1-incArrayF-expArrayF).*(ye(14*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(15*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(:,16*numFieldsIn+(1:numFieldsIn/2)));%sp#2
n_superAbundAF=H*(1-incArrayF-expArrayF).*(ye(0*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(1*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(:,2*numFieldsIn+(1:numFieldsIn/2)));%sp#1nymph
n_superAbundBF=H*(1-incArrayF-expArrayF).*(ye(11*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(12*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(:,13*numFieldsIn+(1:numFieldsIn/2)));%sp#2nymph
superAbundF=superAbundAF+superAbundBF; % both sp
n_superAbundF=n_superAbundAF+n_superAbundBF;


%%%% Commence output plot top panel %%%%
figure
xticmarks=0:10:100;
hold on;
subplot(4,2,3)
yyaxis left
p1=plot(1:(numFieldsIn/2),(1/H)*superAbundF./(infArrayF+1),'k');
normTo1=max((1/H)*superAbundF./(infArrayF+1));
hold on;
normTo2=max((1/H)*n_superAbundF./(infArrayF+1));
xlim([20 55])
set(gca,'ytick',0:10:200,'xtick',xticmarks,'xticklabel',[]);
ylabel('Wave-profile (Y/(I+1))');
xlabel('Field position');
yyaxis right
p3=plot(1:(numFieldsIn/2),incArrayF,'color',[0 153 212]/255); %[0.68 0.85 0.9]
ylim([0 1]);
set(gca,'ytick',0:0.2:1,'xtick',xticmarks,'xticklabel',[]);
ax1=gca;
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ax1.YAxis(1).Color = 'k';
ax1.YAxis(2).Color = [0 153 212]/255;

subplot(4,2,4)
yyaxis left
plot(1:(numFieldsIn/2),incArrayF);
ylim([0 1])
set(gca,'ytick',0:0.2:1,'xtick',xticmarks,'xticklabel',[]);
yyaxis right
plot(1:(numFieldsIn/2),(1/H)*superAbundF);
set(gca,'ytick',0:100:1000,'xtick',xticmarks,'xticklabel',[]);
ylim([0 350])
xlim([20 55])
ax1b=gca;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%% SCENARIO 2 %%%
% additional scenario specific parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVASIVE VECTOR %%%%%%%%%%%%%
eps1=1;                                                                    %%%% INVASIVE VECTOR %%%%%%%%%%%%%
epsSP2=5;                                                                  %%%% INVASIVE VECTOR %%%%%%%%%%%%%
epslscape=1;                                                               %%%% INVASIVE VECTOR %%%%%%%%%%%%%
sp2On=1;   % note SP2 switch must be on for invasive wf scenario... i.e. SP2 seeds initial density of invasive wf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVASIVE VECTOR %%%%%%%%%%%%%



%%%%%% RUN first TO insect STEADY-STATE (no pathogen!) %%%%%%%%%%%%%%%%%%%%
invIndex=1; % field location of invader
initVals=[zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) K*(1-(muu/a))*ones(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn)];
[tnewSS,ynewSS] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,range,propCleanSeed,cutFromWithin,is_discrim,[lh_pams eps1 epsSP2 epslscape]),[0 tMaxSS],initVals);
disp('Scenario 2 ss run complete');
nymphsSS=ynewSS(:,0*numFieldsIn+(1:numFieldsIn));
insectsSS=ynewSS(:,3*numFieldsIn+(1:numFieldsIn));
nymphsSS2=ynewSS(:,11*numFieldsIn+(1:numFieldsIn));
insectsSS2=ynewSS(:,14*numFieldsIn+(1:numFieldsIn));
% TO VIEW TIME SERIES (for checking steady state reached)
%    figure;
%    plot(tnewSS,nymphsSS);
%    figure;
%    plot(tnewSS,insectsSS);
%    figure;
%    plot(tnewSS,nymphsSS2);
%    figure;
%    plot(tnewSS,insectsSS2);

%%%%%% RUN second FOR EPIDEMIC WAVE
initValsFromSS=ynewSS(end,:);
initValsFromSS(7*numFieldsIn+invIndex)=1;    % <- disease invasion
initValsFromSS(14*numFieldsIn+invIndex)=sp2On*1; % <- insect invasion (susceptible)
initValsFromSS(3*numFieldsIn+invIndex)=initValsFromSS(3*numFieldsIn+invIndex)-sp2On*1; % <- simply reflecting invasion in loss of indiv. from WT 
options = odeset('Events',@midScapeEvents);    % this option sets the condition for stopping the solver (i.e. when epidemic wave reaches a point on landscape)
[tnew,ynew,te2,ye,ie] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,range,propCleanSeed,cutFromWithin,is_discrim,[lh_pams eps1 epsSP2 epslscape]),[0 tMax],initValsFromSS,options);
disp('Scenario 2 epidemic run complete');
infArrayF=ye(7*numFieldsIn+(1:numFieldsIn/2));
exposArrayF=ye(6*numFieldsIn+(1:numFieldsIn/2));
expArrayF=exposArrayF/H;
incArrayF=infArrayF/H;
if sum(sum((incArrayF+expArrayF)>1.001))    % testing for errors from solver 
    disp('ERROR!!! ');
    for qq=1:numFieldsIn/2
        tmpfield=incArrayF(qq)+expArrayF(qq);
        inds=find(tmpfield>1.001);
        if ~isempty(inds)
            disp(['Field is ' num2str(qq)]);
            disp(['At time ' num2str(inds')]);
            disp(incArrayF(qq)');
            disp(expArrayF(qq)');
            return;
        end
    end
end
superAbundAF=H*(1-incArrayF-expArrayF).*(ye(3*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(4*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(5*numFieldsIn+(1:numFieldsIn/2)));%sp#1
superAbundBF=H*(1-incArrayF-expArrayF).*(ye(14*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(15*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(16*numFieldsIn+(1:numFieldsIn/2)));%sp#2
n_superAbundAF=H*(1-incArrayF-expArrayF).*(ye(0*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(1*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(2*numFieldsIn+(1:numFieldsIn/2)));%sp#1nymph
n_superAbundBF=H*(1-incArrayF-expArrayF).*(ye(11*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(12*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(13*numFieldsIn+(1:numFieldsIn/2)));%sp#2nymph
superAbundF=superAbundAF+superAbundBF;
n_superAbundF=n_superAbundAF+n_superAbundBF;

%%%%%%%%%%%% Commence output plot second panel
xticmarks=0:10:100;
subplot(4,2,5)
yyaxis left
p1=plot(1:(numFieldsIn/2),(1/H)*superAbundF./(infArrayF+1),'k');
normTo1=max((1/H)*superAbundF./(infArrayF+1));
hold on;
normTo2=max((1/H)*n_superAbundF./(infArrayF+1));
set(gca,'ytick',0:20:200,'xtick',xticmarks,'xticklabel',[]);
ylabel('Wave-profile (Y/(I+1))')
xlabel('Field position')
yyaxis right
p3=plot(1:(numFieldsIn/2),incArrayF,'color',[0 153 212]/255);
set(gca,'ytick',0:0.2:1,'xtick',xticmarks,'xticklabel',[]);
xlim([20 55])
ylim([0 1])
ax2=gca;
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ax2.YAxis(1).Color = 'k';
ax2.YAxis(2).Color = [0 153 212]/255;

subplot(4,2,6)
yyaxis left
plot(1:(numFieldsIn/2),incArrayF);
set(gca,'ytick',0:0.2:1,'xtick',xticmarks,'xticklabel',[]);
ylim([0 1])
yyaxis right
plot(1:(numFieldsIn/2),superAbundF/H);
set(gca,'ytick',0:100:600,'xtick',xticmarks,'xticklabel',[]);
xlim([20 55])
ylim([0 350])
ax2b=gca;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%












%%% SCENARIO 3 %%%
% additional scenario specific parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
eps1=1;                                                                    %%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
epsSP2=1;                                                                  %%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
epslscape=3;                                                               %%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%
sp2On=0;  % note SP2 switch must be on for invasive wf scenario... 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIGH ENVIRONMENTAL SUITABILITY %%%%%%%%%%%%%


%%%%%% RUN first TO insect STEADY-STATE (no pathogen!) %%%%%%%%%%%%%%%%%%%%
invIndex=1; % field location of invader
initVals=[zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) K*(1-(muu/a))*ones(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn) zeros(1,numFieldsIn)];
[tnewSS,ynewSS] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,range,propCleanSeed,cutFromWithin,is_discrim,[lh_pams eps1 epsSP2 epslscape]),[0 tMaxSS],initVals);
disp('Scenario 3 ss run complete');
nymphsSS=ynewSS(:,0*numFieldsIn+(1:numFieldsIn));
insectsSS=ynewSS(:,3*numFieldsIn+(1:numFieldsIn));
nymphsSS2=ynewSS(:,11*numFieldsIn+(1:numFieldsIn));
insectsSS2=ynewSS(:,14*numFieldsIn+(1:numFieldsIn));
% TO VIEW TIME SERIES (for checking steady state reached)
%    figure;
%    plot(tnewSS,nymphsSS);
%    figure;
%    plot(tnewSS,insectsSS);
%    figure;
%    plot(tnewSS,nymphsSS2);
%    figure;
%    plot(tnewSS,insectsSS2);

%%%%%% RUN second FOR EPIDEMIC WAVE
initValsFromSS=ynewSS(end,:);
initValsFromSS(7*numFieldsIn+invIndex)=1;    % <- disease invasion
initValsFromSS(14*numFieldsIn+invIndex)=sp2On*1; % <- insect invasion (susceptible)
initValsFromSS(3*numFieldsIn+invIndex)=initValsFromSS(3*numFieldsIn+invIndex)-sp2On*1; % <- simply reflecting invasion in loss of indiv. from WT 
options = odeset('Events',@midScapeEvents);   % this option sets the condition for stopping the solver (i.e. when epidemic wave reaches a point on landscape)
[tnew,ynew,te3,ye,ie] = ode23(@(t,pops)aggreg2_wCuttings(t,pops,range,propCleanSeed,cutFromWithin,is_discrim,[lh_pams eps1 epsSP2 epslscape]),[0 tMax],initValsFromSS,options);
disp('Scenario 3 epidemic run complete');
infArrayF=ye(7*numFieldsIn+(1:numFieldsIn/2));
exposArrayF=ye(6*numFieldsIn+(1:numFieldsIn/2));
expArrayF=exposArrayF/H;
incArrayF=infArrayF/H;
if sum(sum((incArrayF+expArrayF)>1.001)) % testing for errors from solver 
    disp('ERROR!!! ');
    for qq=1:numFieldsIn/2
        tmpfield=incArrayF(qq)+expArrayF(qq);
        inds=find(tmpfield>1.001);
        if ~isempty(inds)
            disp(['Field is ' num2str(qq)]);
            disp(['At time ' num2str(inds')]);
            disp(incArrayF(qq)');
            disp(expArrayF(qq)');
            return;
        end
    end
end
superAbundAF=H*(1-incArrayF-expArrayF).*(ye(3*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(4*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(5*numFieldsIn+(1:numFieldsIn/2)));%sp#1
superAbundBF=H*(1-incArrayF-expArrayF).*(ye(14*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(15*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(16*numFieldsIn+(1:numFieldsIn/2)));%sp#2
n_superAbundAF=H*(1-incArrayF-expArrayF).*(ye(0*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(1*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(2*numFieldsIn+(1:numFieldsIn/2)));%sp#1nymph
n_superAbundBF=H*(1-incArrayF-expArrayF).*(ye(11*numFieldsIn+(1:numFieldsIn/2)))+H*expArrayF.*(ye(12*numFieldsIn+(1:numFieldsIn/2)))+H*incArrayF.*(ye(13*numFieldsIn+(1:numFieldsIn/2)));%sp#2nymph
superAbundF=superAbundAF+superAbundBF;
n_superAbundF=n_superAbundAF+n_superAbundBF;


%%%%%%%%%%%% Commence output plot third panel
xticmarks=0:10:100;
subplot(4,2,7)
yyaxis left
p1=plot(1:(numFieldsIn/2),(1/H)*superAbundF./(infArrayF+1),'k');
normTo1=max((1/H)*superAbundF./(infArrayF+1));
hold on;
normTo2=max((1/H)*n_superAbundF./(infArrayF+1));
set(gca,'ytick',0:20:200,'xtick',xticmarks);
ylabel('Wave-profile (Y/(I+1))')
xlabel('Field position')
yyaxis right
p3=plot(1:(numFieldsIn/2),incArrayF,'color',[0 153 212]/255);
set(gca,'ytick',0:0.2:1,'xtick',xticmarks);
xlim([20 55])
ylim([0 1])
ax3=gca;
set(get(get(p3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ax3.YAxis(1).Color = 'k';
ax3.YAxis(2).Color = [0 153 212]/255;

subplot(4,2,8)
yyaxis left
plot(1:(numFieldsIn/2),incArrayF);
set(gca,'ytick',0:0.2:1,'xtick',xticmarks);
ylim([0 1])
yyaxis right
plot(1:(numFieldsIn/2),superAbundF/H);
set(gca,'ytick',0:100:600,'xtick',xticmarks);
xlim([20 55])
ylim([0 350])
ax3b=gca;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%% dummy plots
subplot(4,2,1)
plot([0,1],[0,1]);
ax0=gca;
subplot(4,2,2)
plot([0,1],[0,1]);
ax0b=gca;

set(gcf, 'Units','centimeters', 'Position',[0 0 12.5 24.5])
set(gcf,'color', 'w');
set([ax0 ax0b ax1 ax2 ax3 ax1b ax2b ax3b], ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'in'     , ...
    'TickLength'  , [.01 .01] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'YGrid'       , 'off'      , ...
    'XColor'      , [.3 .3 .3], ...
    'FontSize', 12, ...
    'FontName', 'Times New Roman', ...
    'LineWidth'   , 0.2         );    %'YColor'      , [.3 .3 .3], ...
h=gcf;
h.Renderer='Painters';
orient(h,'landscape')

%saveas(gcf,'myfigure1.pdf')
