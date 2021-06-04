function scout = aggreg2_wCuttings(t,pops,icbn,dN,range,propCleanSeed,cutFromWithin)
global a K thet H muu Pacq Pinoc plD sig pMigrate pDisappear eps1 epsSP2 muuN plE epslscape
sizy=size(pops);
numFields=sizy(1)/20;
fn=find(pops<0);
pops(fn)=0;
%propCleanSeed=0.5;
%cutFromWithin=0.5;
infCuttings=0;
expCuttings=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extract field specific parts
vecS_N=pops(1:numFields);
vecE_N=pops(numFields+(1:numFields));
vecI_N=pops(2*numFields+(1:numFields));
vecS=pops(3*numFields+(1:numFields));
vecE=pops(4*numFields+(1:numFields));
vecI=pops(5*numFields+(1:numFields));
plExp=pops(6*numFields+(1:numFields))/H;
plInc=pops(7*numFields+(1:numFields))/H;
vecSinc=pops(8*numFields+(1:numFields));
vecEinc=pops(9*numFields+(1:numFields));
vecIinc=pops(10*numFields+(1:numFields));
% extract field specific parts (vector sp #2)
vecS_N_2=pops(11*numFields+(1:numFields));
vecE_N_2=pops(12*numFields+(1:numFields));
vecI_N_2=pops(13*numFields+(1:numFields));
vecS_2=pops(14*numFields+(1:numFields));
vecE_2=pops(15*numFields+(1:numFields));
vecI_2=pops(16*numFields+(1:numFields));
vecSinc_2=pops(17*numFields+(1:numFields));
vecEinc_2=pops(18*numFields+(1:numFields));
vecIinc_2=pops(19*numFields+(1:numFields));
sharedS=0*vecS_N+vecS+0*vecS_N_2+vecS_2; %zero weight nymphs for simplicity
sharedE=0*vecE_N+vecE+0*vecE_N_2+vecE_2;
sharedI=0*vecI_N+vecI+0*vecI_N_2+vecI_2; !
scout=zeros(20*numFields,1);
disperseOnI=zeros(1,numFields);disperseOnE=zeros(1,numFields);
for ii=1:numFields                                 % to avoid singularities
    if plInc(ii)>1e-100, disperseOnI(ii)=(1./(H*plInc(ii))); end
    if plExp(ii)>1e-100, disperseOnE(ii)=(1./(H*plExp(ii))); end
end
probLandS=(1-plInc-plExp);
probLandI=plInc;
probLandE=plExp;

pStay=1-pMigrate;

localInsMovemt=thet*pStay*(H*plInc.*vecI+H*plExp.*vecE+H*(1-plInc-plExp).*vecS);
migratInsMovemt=thet*(1-pDisappear)*pMigrate*(H*plInc.*vecI+H*plExp.*vecE+H*(1-plInc-plExp).*vecS);
localINFMovemt=thet*pStay*(H*(1-plInc-plExp).*vecSinc+H*plExp.*vecEinc+H*plInc.*vecIinc);
migratINFMovemt=thet*(1-pDisappear)*pMigrate*(H*(1-plInc-plExp).*vecSinc+H*plExp.*vecEinc+H*plInc.*vecIinc);
localInsMovemt_2=thet*pStay*(H*plInc.*vecI_2+H*plExp.*vecE_2+H*(1-plInc-plExp).*vecS_2);
migratInsMovemt_2=thet*(1-pDisappear)*pMigrate*(H*plInc.*vecI_2+H*plExp.*vecE_2+H*(1-plInc-plExp).*vecS_2);
localINFMovemt_2=thet*pStay*(H*(1-plInc-plExp).*vecSinc_2+H*plExp.*vecEinc_2+H*plInc.*vecIinc_2);
migratINFMovemt_2=thet*(1-pDisappear)*pMigrate*(H*(1-plInc-plExp).*vecSinc_2+H*plExp.*vecEinc_2+H*plInc.*vecIinc_2);


% to identify incoming insects to field j take 0.5 of migratInsMovemt(j+1) and migratInsMovemt(j-1)
migratInsMovemt_byField=zeros(1,numFields);
migratInsMovemt_byField(2:end-1)=0.5*migratInsMovemt(1:(numFields-2))+0.5*migratInsMovemt(3:numFields);
migratInsMovemt_byField(1)=0.5*migratInsMovemt(numFields)+0.5*migratInsMovemt(2);
migratInsMovemt_byField(numFields)=0.5*migratInsMovemt(numFields-1)+ 0.5*migratInsMovemt(1);
% to identify infected insects to field j take 0.5 of migratInsMovemt(j+1) and migratInsMovemt(j-1)
migratINFMovemt_byField=zeros(1,numFields);
migratINFMovemt_byField(2:end-1)=0.5*migratINFMovemt(1:(numFields-2))+0.5*migratINFMovemt(3:numFields);
migratINFMovemt_byField(1)=0.5*migratINFMovemt(numFields)+ 0.5*migratINFMovemt(2);
migratINFMovemt_byField(numFields)=0.5*migratINFMovemt(numFields-1)+ 0.5*migratINFMovemt(1);
% to identify incoming insects to field j take 0.5 of migratInsMovemt(j+1)
% and migratInsMovemt(j-1) SPECIES #2 SPECIES #2 SPECIES #2 SPECIES #2 SPECIES #2
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
inFlowINF_byField=cuttingsScope2(plInc,range);
if(infCuttings==0)
    inFlowEXP_byField=cuttingsScope2(plExp./(1-plInc),range);  % 'careful' cutting
    inFlowEXP_byField=cuttingsScope2(plExp,range);  % 'careful' cutting
    withinEXP_byField=(plExp./(1-plInc));
    withinEXP_byField=(plExp);
else
    inFlowEXP_byField=cuttingsScope2(plExp,range);   % 'indiscrimate' cutting
    withinEXP_byField=plExp;
end
    
inFlowEXP_byField(inFlowEXP_byField>1)=1;
inFlowEXP_byField(inFlowEXP_byField<0)=0;
withinEXP_byField(withinEXP_byField>1)=1;
withinEXP_byField(withinEXP_byField<0)=0;


if sum(inFlowEXP_byField>1)>0
    disp('ERROR !! inFlowEXP_byField greater than one somewhere'); 
    disp('check out frac susc'); 
    disp(plInc'); 
    disp(plExp');
    disp((1-plInc-plExp)');
end

%if t>365
%   if (flagOn~=1)
%   plot(inFlowINF_byField);
%   flagOn=1;
%   end
%end

% SPECIES#1 SPECIES#1 SPECIES#1 SPECIES#1 SPECIES#1 SPECIES#1 SPECIES#1 %%%
% nymphal abundance on healthy plants
scout(1:numFields)=a*vecS.*(1-(sharedS/(epslscape*K))).*((sharedS/(epslscape*K))<1)-(muuN+dN)*vecS_N;
% nymphal abundance on exposed plants
scout(1*numFields+(1:numFields))=a*vecE.*(1-(sharedE/(epslscape*K))).*((sharedE/(epslscape*K))<1)-(muuN+dN)*vecE_N;
% nymphal abundance on infected plants
scout(2*numFields+(1:numFields))=a*vecI.*(1-(sharedI/(epslscape*K*eps1))).*((sharedI/(epslscape*K*eps1))<1)-(muuN+dN)*vecI_N;
% abundance on healthy plants
scout(3*numFields+(1:numFields))=dN*vecS_N-(thet+muu)*vecS+(localInsMovemt+migratInsMovemt_byField').*probLandS.*(1./(H*(1-plInc-plExp)));
% abundance on exposed plants
scout(4*numFields+(1:numFields))=dN*vecE_N-(thet+muu)*vecE+(localInsMovemt+migratInsMovemt_byField').*probLandE.*disperseOnE';
% abundance on infected plants
scout(5*numFields+(1:numFields))=dN*vecI_N-(thet+muu)*vecI+(localInsMovemt+migratInsMovemt_byField').*probLandI.*disperseOnI';


% exposure in plants (full numbers)                                                                           
scout(6*numFields+(1:numFields))=expCuttings*(1-propCleanSeed)*((1-cutFromWithin)*(inFlowEXP_byField)+cutFromWithin*withinEXP_byField).*(plE*H*plExp+(plD+plE)*H*plInc)+Pinoc*(vecSinc+vecSinc_2)*H.*(1-plInc-plExp)-plE*H*plExp-icbn*H*plExp;  %%% <<<---- NOTE: First term prop to H.*(1-plInc) here because obviously not dealing with per plant unlike P_inoc
% incidence in plants (full numbers)
scout(7*numFields+(1:numFields))=infCuttings*(1-propCleanSeed)*((1-cutFromWithin)*(inFlowINF_byField)+cutFromWithin*plInc).*(plE*H*plExp+(plD+plE)*H*plInc)+icbn*H*plExp-(plD+plE)*H*plInc;

% incidence in vectors on S plants (full numbers)
scout(8*numFields+(1:numFields))=-sig*vecSinc-thet*vecSinc-muu*vecSinc+(localINFMovemt+migratINFMovemt_byField').*probLandS.*(1./(H*(1-plInc-plExp)));
% incidence in vectors on E plants (full numbers)
scout(9*numFields+(1:numFields))=-sig*vecEinc-thet*vecEinc-muu*vecEinc+(localINFMovemt+migratINFMovemt_byField').*probLandE.*disperseOnE';
% incidence in vectors on I plants (full numbers)
scout(10*numFields+(1:numFields))=-sig*vecIinc-thet*vecIinc-muu*vecIinc+(localINFMovemt+migratINFMovemt_byField').*probLandI.*disperseOnI'+Pacq*(vecI-vecIinc);   %%% <<<---- NOTE: Final term not prop to H.*plInc here because the units of infected insects on I plants is per plant
% SPECIES#2 SPECIES#2 SPECIES#2 SPECIES#2 SPECIES#2 SPECIES#2 SPECIES#2 %%%
% nymphal abundance on healthy plants
scout(11*numFields+(1:numFields))=a*vecS_2.*(1-(sharedS/(epslscape*K*epsSP2))).*((sharedS/(epslscape*K*epsSP2))<1)-(muuN+dN)*vecS_N_2;
% nymphal abundance on exposed plants
scout(12*numFields+(1:numFields))=a*vecE_2.*(1-(sharedE/(epslscape*K*epsSP2))).*((sharedE/(epslscape*K*epsSP2))<1)-(muuN+dN)*vecE_N_2;
% nymphal abundance on infected plants
scout(13*numFields+(1:numFields))=a*vecE_2.*(1-(sharedI/(epslscape*K*epsSP2*eps1))).*((sharedI/(epslscape*K*epsSP2*eps1))<1)-(muuN+dN)*vecI_N_2;
% abundance on healthy plants
scout(14*numFields+(1:numFields))=dN*vecS_N_2-(thet+muu)*vecS_2+(localInsMovemt_2+migratInsMovemt_byField_2').*probLandS.*(1./(H*(1-plInc-plExp)));
% abundance on exposed plants
scout(15*numFields+(1:numFields))=dN*vecE_N_2-(thet+muu)*vecE_2+(localInsMovemt_2+migratInsMovemt_byField_2').*probLandE.*disperseOnE';
% abundance on infected plants
scout(16*numFields+(1:numFields))=dN*vecI_N_2-(thet+muu)*vecI_2+(localInsMovemt_2+migratInsMovemt_byField_2').*probLandI.*disperseOnI';
% incidence in vectors on S plants (full numbers)
scout(17*numFields+(1:numFields))=-sig*vecSinc_2-thet*vecSinc_2-muu*vecSinc_2+(localINFMovemt_2+migratINFMovemt_byField_2').*probLandS.*(1./(H*(1-plInc-plExp)));
% incidence in vectors on E plants (full numbers)
scout(18*numFields+(1:numFields))=-sig*vecEinc_2-thet*vecEinc_2-muu*vecEinc_2+(localINFMovemt_2+migratINFMovemt_byField_2').*probLandE.*disperseOnE';
% incidence in vectors on I plants (full numbers)
scout(19*numFields+(1:numFields))=-sig*vecIinc_2-thet*vecIinc_2-muu*vecIinc_2+(localINFMovemt_2+migratINFMovemt_byField_2').*probLandI.*disperseOnI'+Pacq*(vecI_2-vecIinc_2);

end
