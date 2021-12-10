function out = cuttingsScope2(plInc,range)
numFields=length(plInc);
out=zeros(size(plInc));
if range>(numFields-1), disp('Error cuttings range broader than non-index landscape!'); return; end
for f=1:numFields
    keepTrackSum=0;
    count=0;
    for g=1:range
        ind1=(f+range+1-g);
        if(ind1<=numFields)
            keepTrackIncrem=g*plInc(ind1);
            count=count+g;
        else
            modInd=mod(ind1,numFields);
            keepTrackIncrem=g*plInc(modInd);
            count=count+g;
        end
        keepTrackSum=keepTrackSum+keepTrackIncrem;
        ind2=(f-range-1+g);
        if(ind2>0)
            keepTrackIncrem=g*plInc(ind2);
            count=count+g;
        else
            notModInd=numFields+ind2;
            keepTrackIncrem=g*plInc(notModInd);
            count=count+g;
        end
        keepTrackSum=keepTrackSum+keepTrackIncrem;
    end
    out(f)=keepTrackSum/count; 
end
end