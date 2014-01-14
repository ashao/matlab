function [ array ] = gasex( array )

	ocmip=load('../data/ocmip.mat');
	ocmip=ocmip.ocmip;

	[ntime numz numy numx]=size(array.temp);
	T=squeeze(array.temp(:,1,:,:));
	array.schmidt=array.schcoeff(1)-array.schcoeff(2).*T+array.schcoeff(3).*T.^2-array.schcoeff(4).*T.^3;
		
	Sk=(1-ocmip.fice).*ocmip.xkw.*array.schmidt;

	sat=ocmip.p.*array.atmcon.*squeeze(array.F(:,1,:,:));
	array.flux=Sk.*sat;

end
