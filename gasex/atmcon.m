function [ array ] = atmcon( array,conN,conS )

	[ntime numz numy numx]=size(array.temp);
	
	
	Nidx=find(array.y>20);
	Sidx=find(array.y<-20);
	interpidx=find(array.y<20 & array.y>-20);
	
	array.atmcon(1,Nidx,1:numx)=conN;
	array.atmcon(1,Sidx,1:numx)=conS;

	size(interp1([20 -20],[conN conS],array.y(interpidx)))
	intcon=interp1([20 -20],[conN conS],array.y(interpidx));
	size(intcon)

	plot(array.y(interpidx),intcon)	

	for i=1:length(interpidx)
	array.atmcon(1,interpidx(i),1:numx)=intcon(i);
	end
end
