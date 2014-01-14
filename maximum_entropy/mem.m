function [ ttd ] = mem ( prior, Cs, Cint )
%% function [ ttd ] = mem ( prior, Cs )
% Estimates the TTD at a given point using the maximum entropy method
% described in Holzer et al. 2010
%   DIMENSIONS:
%       ncon  = Number of constraints
%       ntime = Number of time measurements prior to present
%       nsource = Number of surface boundary conditions
%   INPUT:
%       prior: Guess of the ttd from the surface to the interior 
%           [ntime nsource]=size(prior);
%       Cs: Surface boundary condition of the tracer
%           [ncon ntime nsource]=size(Cs);
%       Cint: Observation at the interior (used for constraints)
%           ncon = length(Cint)
%   OUTPUT:
%       ttd: Maximum entropy TTD
%            ntime=length(ttd)

[ncon ntime nsource]=size(Cs);
options=optimset('Display','Iter','Algorithm','levenberg-marquardt');
lambda_opt = fsolve(@constraints,rand(ncon,1),options);
ttd=calcP(lambda_opt,Cs);


    function cost = constraints( lambda )
        
        Pcand=calcP(lambda, Cs);        
        cost=zeros(ncon,1);
        for con=1:ncon           
           cost(con)=sum(sum(Pcand.*squeeze(Cs(con,:,:))))-Cint(con);
        end        
        
    end

function [ Pcand ] = calcP ( lambda, Cs_scale )

    for con=1:ncon
        Cs_scale(con,:,:)=Cs_scale(con,:,:)*lambda(con);
    end    
    lambda_Cjs=squeeze(sum(Cs_scale));        
    Pcand=prior.*exp(-lambda_Cjs); 
    norm=sum(sum(Pcand));
    Pcand=Pcand/norm;        
    

end

end