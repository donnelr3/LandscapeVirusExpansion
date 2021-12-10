function scout = aggreg2_wCuttings(t,pops,range,propCleanSeed,cutFromWithin,pams)
a=pams(1);
K=pams(2);
thet=pams(3);
H=pams(4);
muu=pams(5);
Pacq=pams(6);
Pinoc=pams(7);
plD=pams(8);
sig=pams(9);
pMigrate=pams(10);
pDisappear=pams(11);
muuN=pams(12);
plE=pams(13);
dN=pams(14);
icbn=pams(15);
eps1=pams(16);
epsSP2=pams(17);
epslscape=pams(18);

sizy=size(pops);
numFields=sizy(1)/20;
if sum(isnan(pops)>0)
    disp(['Error NANS! !!! !! !! time= ' num2str(t) ' at positions ' strtrim(cellstr(num2str(find(isnan(pops))'))')]);
    disp(['ie field #' num2str(mod(362,120))]);
    disp((pops(6*numFields+(1:numFields))/H)');
    disp((pops(7*numFields+(1:numFields))/H)');
    disp(((pops(6*numFields+(1:numFields))+pops(7*numFields+(1:numFields)))/H)');
    disp((((pops(6*numFields+(1:numFields))+pops(7*numFields+(1:numFields)))/H)>1)');
    return
end


infCuttings=0;
expCuttings=1;

% extract field specific parts (wildtype vector)
vecS_N=pops(1:numFields);                          % nymphs on S plant
vecE_N=pops(numFields+(1:numFields));              % nymphs on E plant
vecI_N=pops(2*numFields+(1:numFields));            % nymphs on I plant
vecS=pops(3*numFields+(1:numFields));              % adults on S plant
vecE=pops(4*numFields+(1:numFields));              % adults on E plant
vecI=pops(5*numFields+(1:numFields));              % adults on I plant
plExp=pops(6*numFields+(1:numFields))/H;           % prop of plants exposed
plInc=pops(7*numFields+(1:numFields))/H;           % prop of plants infected
vecSinf=pops(8*numFields+(1:numFields));           % infected vectors per S plant
vecEinf=pops(9*numFields+(1:numFields));           % infected vectors per E plant
vecIinf=pops(10*numFields+(1:numFields));          % infected vectors per I plant

% extract field specific parts (invasive vector)
vecS_N_2=pops(11*numFields+(1:numFields));         % See above comments
vecE_N_2=pops(12*numFields+(1:numFields));         
vecI_N_2=pops(13*numFields+(1:numFields));        
vecS_2=pops(14*numFields+(1:numFields));          
vecE_2=pops(15*numFields+(1:numFields));         
vecI_2=pops(16*numFields+(1:numFields));         
vecSinf_2=pops(17*numFields+(1:numFields));       
vecEinf_2=pops(18*numFields+(1:numFields));      
vecIinf_2=pops(19*numFields+(1:numFields));      

% plant carrying capacities for insect vectors not dependent on nymphs
sharedS=0*vecS_N+vecS+0*vecS_N_2+vecS_2;
sharedE=0*vecE_N+vecE+0*vecE_N_2+vecE_2;
sharedI=0+vecI+0+vecI_2;

% calculating prob. land on plant type... avoiding later singularities
disperseOnI=zeros(1,numFields);
disperseOnE=zeros(1,numFields);
disperseOnS=zeros(1,numFields);
for ii=1:numFields                               
    if plInc(ii)>1e-100, disperseOnI(ii)=(1./(H*plInc(ii))); end
    if plExp(ii)>1e-100, disperseOnE(ii)=(1./(H*plExp(ii))); end
    if (1-plExp(ii)-plInc(ii))>1e-100, disperseOnS(ii)=(1./(H*(1-plExp(ii)-plInc(ii))));end
end
probLandS=(1-plInc-plExp);
probLandI=plInc;
probLandE=plExp;

pStay=1-pMigrate;

% calculating pool of dispersing vectors within fields (to be later dispersed into plant types)
% wildtype insect vector
localInsMovemt=thet*pStay*(H*plInc.*vecI+H*plExp.*vecE+H*(1-plInc-plExp).*vecS);                       % local represents within field dispersal
migratInsMovemt=thet*(1-pDisappear)*pMigrate*(H*plInc.*vecI+H*plExp.*vecE+H*(1-plInc-plExp).*vecS);    % migrat represents between field dispersal
localINFMovemt=thet*pStay*(H*(1-plInc-plExp).*vecSinf+H*plExp.*vecEinf+H*plInc.*vecIinf);                   % localINF is within field viruliferous dispersal
migratINFMovemt=thet*(1-pDisappear)*pMigrate*(H*(1-plInc-plExp).*vecSinf+H*plExp.*vecEinf+H*plInc.*vecIinf);% migratINF is between field viruliferous dispersal
% invasive insect vector
localInsMovemt_2=thet*pStay*(H*plInc.*vecI_2+H*plExp.*vecE_2+H*(1-plInc-plExp).*vecS_2);
migratInsMovemt_2=thet*(1-pDisappear)*pMigrate*(H*plInc.*vecI_2+H*plExp.*vecE_2+H*(1-plInc-plExp).*vecS_2);
localINFMovemt_2=thet*pStay*(H*(1-plInc-plExp).*vecSinf_2+H*plExp.*vecEinf_2+H*plInc.*vecIinf_2);
migratINFMovemt_2=thet*(1-pDisappear)*pMigrate*(H*(1-plInc-plExp).*vecSinf_2+H*plExp.*vecEinf_2+H*plInc.*vecIinf_2);


% to identify incoming insects to field j take 0.5 of migratInsMovemt(j+1) and migratInsMovemt(j-1)
% WILDTYPE insect vector
migratInsMovemt_byField=zeros(1,numFields);
migratInsMovemt_byField(2:end-1)=0.5*migratInsMovemt(1:(numFields-2))+0.5*migratInsMovemt(3:numFields);
migratInsMovemt_byField(1)=0.5*migratInsMovemt(numFields)+0.5*migratInsMovemt(2);
migratInsMovemt_byField(numFields)=0.5*migratInsMovemt(numFields-1)+ 0.5*migratInsMovemt(1);
% to identify infected insects to field j take 0.5 of migratInsMovemt(j+1) and migratInsMovemt(j-1)
migratINFMovemt_byField=zeros(1,numFields);
migratINFMovemt_byField(2:end-1)=0.5*migratINFMovemt(1:(numFields-2))+0.5*migratINFMovemt(3:numFields);
migratINFMovemt_byField(1)=0.5*migratINFMovemt(numFields)+ 0.5*migratINFMovemt(2);
migratINFMovemt_byField(numFields)=0.5*migratINFMovemt(numFields-1)+ 0.5*migratINFMovemt(1);
% INVASIVE insect vector
migratInsMovemt_byField_2=zeros(1,numFields);
migratInsMovemt_byField_2(2:end-1)=0.5*migratInsMovemt_2(1:(numFields-2))+0.5*migratInsMovemt_2(3:numFields);
migratInsMovemt_byField_2(1)=0.5*migratInsMovemt_2(numFields)+0.5*migratInsMovemt_2(2);
migratInsMovemt_byField_2(numFields)=0.5*migratInsMovemt_2(numFields-1)+ 0.5*migratInsMovemt_2(1);
% to identify infected insects to field j take 0.5 of migratInsMovemt(j+1) and migratInsMovemt(j-1)
migratINFMovemt_byField_2=zeros(1,numFields);
migratINFMovemt_byField_2(2:end-1)=0.5*migratINFMovemt_2(1:(numFields-2))+0.5*migratINFMovemt_2(3:numFields);
migratINFMovemt_byField_2(1)=0.5*migratINFMovemt_2(numFields)+ 0.5*migratINFMovemt_2(2);
migratINFMovemt_byField_2(numFields)=0.5*migratINFMovemt_2(numFields-1)+ 0.5*migratINFMovemt_2(1);

%%%% MAY NEED TO ALSO WEIGHT THEM (there is some evidence that if the replacement rate of healthy plants is like monthly widespread cuttings pushes back the min)
inFlowEXP_byField=cuttingsScope2(plExp,range);   % 'discrimate' cutting is just inFlowEXP_byField
inFlowINF_byField=cuttingsScope2(plInc,range);   % 'indiscrimate' cutting is both inFlowEXP_byField and inFlowINF_byField (see eq.s)
    
scout=zeros(20*numFields,1);

% SPECIES#1 SPECIES#1 SPECIES#1 SPECIES#1 SPECIES#1 SPECIES#1 SPECIES#1 %%%
% nymphal density on healthy plants
scout(1:numFields)=a*vecS.*(1-(sharedS/(epslscape*K))).*((sharedS/(epslscape*K))<1)-(muuN+dN)*vecS_N;
% nymphal density on exposed plants
scout(1*numFields+(1:numFields))=a*vecE.*(1-(sharedE/(epslscape*K))).*((sharedE/(epslscape*K))<1)-(muuN+dN)*vecE_N;
% nymphal density on infected plants
scout(2*numFields+(1:numFields))=a*vecI.*(1-(sharedI/(epslscape*K*eps1))).*((sharedI/(epslscape*K*eps1))<1)-(muuN+dN)*vecI_N;
% adult density on healthy plants
scout(3*numFields+(1:numFields))=dN*vecS_N-(thet+muu)*vecS+(localInsMovemt+migratInsMovemt_byField').*probLandS.*disperseOnS';
% adult density on exposed plants
scout(4*numFields+(1:numFields))=dN*vecE_N-(thet+muu)*vecE+(localInsMovemt+migratInsMovemt_byField').*probLandE.*disperseOnE';
% adult density on infected plants
scout(5*numFields+(1:numFields))=dN*vecI_N-(thet+muu)*vecI+(localInsMovemt+migratInsMovemt_byField').*probLandI.*disperseOnI';

% density pathogen exposed plants                                                                     
scout(6*numFields+(1:numFields))=expCuttings*(1-propCleanSeed)*((1-cutFromWithin)*(inFlowEXP_byField)+cutFromWithin*plExp).*(plE*H*plExp+(plD+plE)*H*plInc)+Pinoc*(vecSinf+vecSinf_2)*H.*(1-plInc-plExp)-plE*H*plExp-icbn*H*plExp;  
% density pathogen infected plants      
scout(7*numFields+(1:numFields))=infCuttings*(1-propCleanSeed)*((1-cutFromWithin)*(inFlowINF_byField)+cutFromWithin*plInc).*(plE*H*plExp+(plD+plE)*H*plInc)+icbn*H*plExp-(plD+plE)*H*plInc;

% density infected vectors per S plants 
scout(8*numFields+(1:numFields))=-sig*vecSinf-thet*vecSinf-muu*vecSinf+(localINFMovemt+migratINFMovemt_byField').*probLandS.*disperseOnS';
% density infected vectors per E plants 
scout(9*numFields+(1:numFields))=-sig*vecEinf-thet*vecEinf-muu*vecEinf+(localINFMovemt+migratINFMovemt_byField').*probLandE.*disperseOnE';
% density infected vectors per I plants 
scout(10*numFields+(1:numFields))=-sig*vecIinf-thet*vecIinf-muu*vecIinf+(localINFMovemt+migratINFMovemt_byField').*probLandI.*disperseOnI'+Pacq*(vecI-vecIinf);   

% SPECIES#2 SPECIES#2 SPECIES#2 SPECIES#2 SPECIES#2 SPECIES#2 SPECIES#2 %%%
% nymphal density on healthy plants
scout(11*numFields+(1:numFields))=a*vecS_2.*(1-(sharedS/(epslscape*K*epsSP2))).*((sharedS/(epslscape*K*epsSP2))<1)-(muuN+dN)*vecS_N_2;
% nymphal density on exposed plants
scout(12*numFields+(1:numFields))=a*vecE_2.*(1-(sharedE/(epslscape*K*epsSP2))).*((sharedE/(epslscape*K*epsSP2))<1)-(muuN+dN)*vecE_N_2;
% nymphal density on infected plants
scout(13*numFields+(1:numFields))=a*vecE_2.*(1-(sharedI/(epslscape*K*epsSP2*eps1))).*((sharedI/(epslscape*K*epsSP2*eps1))<1)-(muuN+dN)*vecI_N_2;
% adult density on healthy plants
scout(14*numFields+(1:numFields))=dN*vecS_N_2-(thet+muu)*vecS_2+(localInsMovemt_2+migratInsMovemt_byField_2').*probLandS.*disperseOnS';
% adult density on exposed plants
scout(15*numFields+(1:numFields))=dN*vecE_N_2-(thet+muu)*vecE_2+(localInsMovemt_2+migratInsMovemt_byField_2').*probLandE.*disperseOnE';
% adult density on infected plants
scout(16*numFields+(1:numFields))=dN*vecI_N_2-(thet+muu)*vecI_2+(localInsMovemt_2+migratInsMovemt_byField_2').*probLandI.*disperseOnI';
% density infected vectors per S plants 
scout(17*numFields+(1:numFields))=-sig*vecSinf_2-thet*vecSinf_2-muu*vecSinf_2+(localINFMovemt_2+migratINFMovemt_byField_2').*probLandS.*disperseOnS';
% density infected vectors per E plants 
scout(18*numFields+(1:numFields))=-sig*vecEinf_2-thet*vecEinf_2-muu*vecEinf_2+(localINFMovemt_2+migratINFMovemt_byField_2').*probLandE.*disperseOnE';
% density infected vectors per I plants
scout(19*numFields+(1:numFields))=-sig*vecIinf_2-thet*vecIinf_2-muu*vecIinf_2+(localINFMovemt_2+migratINFMovemt_byField_2').*probLandI.*disperseOnI'+Pacq*(vecI_2-vecIinf_2);

end
