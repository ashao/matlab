function [ sf6 ] = sf6_sol(  )

    datapath='~/matlab/data/';
	disp('Reading Salinity')
    salt=load([datapath 'salt.mat']);
	salt=salt.salt;
	disp('Reading Temperature')
    temp=load([datapath 'temp.mat']);
	temp=temp.temp
    solubility=load([datapath 'solubility.mat']);
    solubility=solubility.solubility 
	who

	coeffs=solubility(3).coeffs;
    [ntime numz numy numx]=size(temp.temp);

	%for time=1:ntime
		%disp([num2str(time) '/' num2str(ntime)])
		%for z=1:numz
			%disp([num2str(z) '/' num2str(numz)])
			%for y=1:numy
				
					sf6.F(:,:,:,:)=solcalc(temp.temp+273.15,salt.sal,coeffs);
				
			%end
		%end
	%end

    function F = solcalc(T,S,coeffs)
      T=290;
	S=50; 
        F= coeffs(1)+coeffs(2)*(100./T)+coeffs(3)*log(T./100)+coeffs(4)*(T./100).^2+S.*(coeffs(5)+coeffs(6)*(T./100)+coeffs(7)*(T./100).^2)
       	F=exp(F);
	disp(F) 
	pause
    end

end

