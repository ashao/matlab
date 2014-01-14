function randn_height = make_zeroskewdist ( numtracks, numlats )


sidx=0;
eidx=0;
numheights=0;
while numheights<numtracks
    
    sidx=eidx+1;
    testdists=randn(numtracks,numlats);
    testskew=skewness(testdists,0,2);
    goodidx=find( abs(testskew)<0.001 );
    eidx=sidx+length(goodidx)-1;
    randn_height(sidx:eidx,:)=testdists(goodidx,:);
    [numheights null]=size(randn_height);
    
end

randn_height=randn_height(1:numtracks,:);



