function [ Tsurf Tdeep tout coeffs ] = ...
    hw1_odefit( dt, nyears, t0, F, obs )

    options=optimset('Display','Iter');
    yearobs=1960:2000;
    data = interp1(obs.time,obs.data,yearobs);
    coeffs = fmincon( @runsol, [1.5 0.5 100 5000] ,[],[],[],[], ...
        [0.1 0.01 10 3000],[5 3 300 6000],[],options);

    function fitval = runsol( coeffs )
        
        [ Tsurf Tdeep tout ]= hw1_odeint( dt, nyears, t0, F, ...
            coeffs(1), coeffs(2),coeffs(3),coeffs(4) );
        sol=interp1(tout/365/86400,Tsurf.sol,yearobs);              
        fitval=sqrt(nansum( (sol-data).^2));
        
        
    end

end