function [co2,pco2,hco3,co3,dic,alk,phtotal,phfree] = f_csys(TC,S,P,flag,in1,in2,phflag,k1k2flag)
% function [co2,pco2,hco3,co3,dic,alk,phtotal,phfree] = f_csys(TC,S,P,flag,in1,in2,phflag,k1k2flag)
%
% This is a functional form of "csys.m" by Zeebe and Wolf-Gladrow. Created by
% S. Mecking 11/21/2006 following my conversion of an earlier version of
% "csys.m" with additional output arguments phfree and ptotal. See description
% of "csys.m" below. Input arguments in1 and in2 have to have units of mol/kg 
% excpet for pH (ph-scale) and pCO2 (mole fraction). Output arguments are also
% in these units even though this routine includes a conversion to umol/kg 
% (different variable names though).
%
% SM, 11/25/2006: added phflag and k1k2flag as input arguments
% SM, 12/1/2006: added Dickson and Millero refit of Mehrbach as option for 
%                k1k2flag 
%
% _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
%
% 	file: csys.m
%
%
%    	Richard E. Zeebe & Dieter A. Wolf-Gladrow
%
%               			     
%	Alfred Wegener Institute for       
%	Polar and Marine Research		  
%	P.O. Box 12 01 61		       
%	D-27515 Bremerhaven		
%	Germany	      
%	e-mail: rzeebe@awi-bremerhaven.de   wolf@awi-bremerhaven.de  	
%	www     : http://www.awi-bremerhaven.de/
%
% _/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/_/
%
% purpose: given two components of CO2 system, calculate all other.
%
%		Compare book by
%
%		Zeebe & Wolf-Gladrow, CO2 in Seawater: 
%		Equilibrium, Kinetics, Isotopes.
% updates:
%
%
%          03.01.01 k -> K
%          01.10.98 include pressure
%          08.05.95 set phflag Total or free
%
% MATLAB 5.3.1
%
% remarks: specify S, T, and P, flag-number for the two variables
%	   pH flag, choose constants (Roy or Mehrbach),
%          and values for two variables.
%
%	   the variable s is equivalent to [CO2]
%
% use:     save file in directory. open matlab. type 'csys'
%
% ---------------------------------------------------------------

%============= INPUT SECTION ===============================%
% flag = 1      h and CO2 given
% flag = 101    h and pco2 given
% flag = 2      CO2 and HCO3 given
% flag = 3      CO2 and CO3 given 
% flag = 4      CO2 and ALK given
% flag = 401    pCO2 and ALK given % SM, 11/22/06 (as in old f_csys)
% flag = 5      CO2 and DIC given
% flag = 6      h and HCO3 given
% flag = 7      h and CO3 given
% flag = 8      h and ALK given
% flag = 9      h and DIC given
% flag = 902    pH and DIC given  % SM, 11/22/06
% flag = 10     HCO3 and CO3 given
% flag = 11     HCO3 and ALK given
% flag = 12     HCO3 and DIC given
% flag = 13     CO3 and ALK given
% flag = 14     CO3 and DIC given
% flag = 15     ALK and DIC given

%============= Set flag ======================================%
% flag = 101;   % 11/22/06: passed as input argument


%--------------   pH Total or pH free  ------------%
%
%       Specify on which pH scale we are working on.
%       DIC, ALK, pCO2, CO2, CO3, and HCO3 are 
%	physical quantities which do not depend on the pH
%	scale! However, the pH does depend on the scale 
%	used. A certain set of equilibrium
%       constants belongs to a certain pH scale.
%       The scale conversion of the constants is
%       done in equic.m.
%
%       Note that from  [H+]_total > [H+]_free
%       follows           pH_total <   pH_free
%
%	The total scale is recommended.
%

if nargin<7 | isempty(phflag)
  % set default
  phflag = 0;     % 0: Total scale
                  % 1: Free scale
else
  % use value passed 
end
phflag;

%----- choose K1 and K2: Roy or Mehrbach.

if nargin<8 | isempty(k1k2flag)
  k1k2flag = 1;	% 0: Roy et al. (1993)
		% 1: Mehrbach et al (1973) as 
		%    refit by Lueker et al. (2000) on total scale.
                % 2: Mehrbach et al. (1973) as refit by Dickson and Millero 
                %    (1987) on sewater scale; assume equivalent to total scale
else
  % use value passed
end
k1k2flag;


%----- set temperature, salinity and pressure

%TC = 25.;		% deg C   % 11/22/06: passed as function argument
%S =  35.;		          % 11/22/06: passed as function argument
%P =  0.;                % bar  % 11/22/06: passed as function argument

%----  load constants
%equic;                        % 11/22/06, replaced bu funtion call
[Kh,K1,K2,Kb,Kw,Ks,Kf,Kspc,Kspa,total2free,sws2free] = f_equic(TC,S,P,phflag,k1k2flag);

bor = 1.*(416.*(S/35.))* 1.e-6;   % (mol/kg), DOE94

%---- set values for carbonate system
%---- specify two variables according to flag (see above)

% 11/22/06: passed as function arguments in1 and in2 which need to be evaluated
%           instead (see below)
%ph1   = 7.9;		%  pH
%h1    = 10^(-ph1);	%  [H+]
%s1    = 18.e-6;		%  [CO2]	(mol/kg)
%hco31 = 1750.e-6;	%  [HCO3-]	(mol/kg)
%co31  = 225.e-6;	%  [CO3--]	(mol/kg)
%alk1  = 3800.e-6;	%  ALK		(mol/kg)
%dic1  = 3900.e-6;	%  [DIC]	(mol/kg)
%pco21 = 1300.e-6;	%  pCO2		(atm)
switch flag    % 11/22/06: evaluation of input arguments
case 1      % h and s given
  h1 = in1;  s1 = in2;
case 101    % h and pco2 given
  h1 = in1;  pco21 = in2;
case 2      % s and HCO3 given
  s1 = in1;  hco31 = in2;
case 3      % s and CO3 given
  s1 = in1;  co31 = in2;
case 4      % s and ALK given
  s1 = in1;  alk1 = in2;
case 401     % pco2 and ALK given
  pco21 = in1; alk1 = in2;
  s1 = pco21*Kh;    % convert so that same calculation as in case 4 can be used
case 5      % s and DIC given
  s1 = in1;  dic1 = in2;
case 6      % h and HCO3 given
  h1 = in1;  hco31 = in2;
case 7      % h and CO3 given
  h1 = in1;  co31 = in2;
case 8      % h and ALK given
  h1 = in1;  alk1 = in2;
case 9      % h and DIC given
  h1 = in1;  dic1 = in2;
case 902      % h and DIC given
  ph1 = in1;  dic1 = in2;
  h1    = 10^(-ph1);% convert so that same calculation as in case 9 can be used 
case 10     % HCO3 and CO3 given
  hco31 = in1;  co31 = in2;
case 11     % HCO3 and ALK given
  hco31 = in1;  alk1 = in2;
case 12     % HCO3 and DIC given
  hco31 = in1;  dic1 = in2;
case 13     % CO3 and ALK given
  co31 = in1;  alk1 = in2;
case 14     % CO3 and DIC given
  co31 = in1;  dic1 = in2;
case 15     % ALK and DIC given             ===================%
  alk1 = in1;  dic1 = in2;
end
  

%----------- These values are test values, set test= 0/1 -------%
test = 0;
if test == 1;
  TC = 25.;
  S = 35.;
  P =  0.;
  equic;
  bor = (416.*(S/35.))* 1.e-6;   % (mol/kg), DOE94   
  ph1   = 8.0902;       %       pH (8.0902 Total, 8.1979 Free)
  h1    = 10^(-ph1);    %       [H+]
  s1    = 10.1304e-6;   %       [CO2]
  hco31 = 1735.9e-06;   %       [HCO3-]
  co31  = 253.9924e-6;  %       [CO3--]
  dic1  = 2000.e-6;     %       [DIC]
  alk1  = 2350.e-6;     %       Alk
  pco21 = 356.8058e-6;	%  	pCO2
end;


%============= END INPUT SECTION ===============================%
%
% rest of file contains the numerical routines.

%==================================================================%

% ----------------- case 1.) h and s given
if flag == 1
disp('flag = 1, pH and CO2 given');
h = h1;
s = s1;
dic = s*(1.+K1/h1+K1*K2/h1/h1);
hco3 = dic/(1+h1/K1+K2/h1);
co3 = dic/(1+h1/K2+h1*h1/K1/K2);
alk = s*(K1/h1+2.*K1*K2/h1/h1)+Kb*bor/(Kb+h1)+Kw/h1-h1;                       
pco2 = s/Kh;
% ----------- change units: mumol/kg
CO2     = 1.e6*s1
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
pco2    = pco2*1.e6
%alktest = (2*co31 + hco31 + bor/(1+h1/Kb) + Kw/h1 - h1)*1.e6
end;
% ----------------- case 101.) h and pco2 given
if flag == 101
disp('flag = 101, pH and pCO2 given');
h = h1;
pco2 = pco21;
s = Kh*pco2;
dic = s*(1.+K1/h1+K1*K2/h1/h1);
hco3 = dic/(1+h1/K1+K2/h1);
co3 = dic/(1+h1/K2+h1*h1/K1/K2);
alk = s*(K1/h1+2.*K1*K2/h1/h1)+Kb*bor/(Kb+h1)+Kw/h1-h1;
% ----------- change units: mumol/kg
CO2     = 1.e6*s                     
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
pco2    = pco2*1.e6
%alktest = (2*co31 + hco31 + bor/(1+h1/Kb) + Kw/h1 - h1)*1.e6
end;
% ------------ s and HCO3 given ------------------
if flag == 2
disp('flag = 2, CO2 and HCO3 given');
s = s1;
hco3 = hco31;
p3 = -hco3/K1;
p2 = s - hco3;
p1 = s*K1 - hco3*K2;
p0 = s*K1*K2;
p = [p3 p2 p1 p0];
r = roots(p);
h = max(real(r));
h*1.e12;
dic = s*(1.+K1/h+K1*K2/h/h);
%       hco3 = dic/(1+h/K1+K2/h); 
co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end
% ------------ s and CO3 given ------------------
if flag == 3
disp('flag = 3, CO2 and CO3 given');
s = s1;
co3 = co31;
p4 = -co3/K1/K2;
p3 = -co3/K2;
p2 = s-co3;
p1 = s*K1;
p0 = s*K1*K2;
p = [p4 p3 p2 p1 p0];      
r = roots(p);
h = max(real(r));
h*1.e12;
dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
%    co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end
% ------------ s and ALK given ------------------
if flag == 4 | flag == 401 % 11/22/06: added case 401
%disp('flag = 4/401, CO2 (or pCO2) and ALK given');
s = s1;
alk = alk1;
p4 = 1.; 
p3 = Kb+alk;
p2 = alk*Kb-s*K1-Kb*bor-Kw;
p1 = -s*Kb*K1-s*2.*K1*K2-Kw*Kb;
p0 = -2.*s*Kb*K1*K2;
p = [p4 p3 p2 p1 p0];
r = roots(p);
h = max(real(r));
h*1.e12;
dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
co2 = s;    % added 11/22/06 to have consistent units in function input & output
pco2 = s/Kh; % added 11/22/06 to have consistent units in function input & outp ut
%    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6;
pCO2    = s*1.e6/Kh;
HCO3    = hco3*1.e6;
CO3     = co3*1.e6;
DIC     = dic*1.e6;
ALK     = alk*1.e6;
end
% ------------ s and DIC given ------------------
if flag == 5
disp('flag = 5, CO2 and DIC given');
s = s1;
dic = dic1;
p2 = dic - s;
p1 = -s*K1;
p0 = -s*K1*K2;
p = [p2 p1 p0];
r = roots(p);
h = max(real(r));
h*1.e12;
%    dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6 
ALK     = alk*1.e6
end
% ------------ h and HCO3 given ------------------
if flag == 6
disp('flag = 6, pH and HCO3 given');
h = h1;
hco3 = hco31;
dic = hco3 * (1+h/K1+K2/h);
s = dic / (1.+K1/h+K1*K2/h/h);
h*1.e12;
%    dic = s*(1.+K1/h+K1*K2/h/h);
%    hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6        
end
% ------------ h and CO3 given ------------------
if flag == 7
disp('flag = 7, pH and CO3 given');
h = h1;
co3 = co31;
dic = co3 * (1+h/K2+h*h/K1/K2);
s = dic / (1.+K1/h+K1*K2/h/h);
h*1.e12;
%    dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
%    co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end                    
% ------------ h and ALK given ------------------
if flag == 8
disp('flag = 8, pH and ALK given');
h = h1;
alk = alk1;
s = (alk-Kw/h+h-Kb*bor/(Kb+h)) / (K1/h+2.*K1*K2/h/h);
h*1.e12;
dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
%    alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end
% ------------ h and DIC given ------------------
if flag == 9 | flag == 902 % 11/22/06: added case 902 
%disp('flag = 9/902, H+ (or pH) and DIC given');
h = h1;
dic = dic1;
s = dic / (1.+K1/h+K1*K2/h/h);
h*1.e12;
%    dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
co2 = s;    % added 11/22/06 to have consistent units in function input & output
pco2 = s/Kh; % added 11/22/06 to have consistent units in function input & outp ut
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6;
pCO2    = s*1.e6/Kh;
HCO3    = hco3*1.e6;
CO3     = co3*1.e6;
DIC     = dic*1.e6;
ALK     = alk*1.e6;
end
% ------------ HCO3 and CO3 given ------------------
if flag == 10
disp('flag = 10, HCO3 and CO3 given');
hco3 = hco31;  
co3 = co31;
p3 = -co3/K1/K2;
p2 = -co3/K2 + hco3/K1;
p1 = -co3 + hco3;
p0 = hco3*K2;
p = [p3 p2 p1 p0];
r = roots(p);
h = max(real(r));
h*1.e12;
dic = hco3 * (1+h/K1+K2/h);
s = dic / (1.+K1/h+K1*K2/h/h);
%    hco3 = dic/(1+h/K1+K2/h);
%    co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end
% ------------ HCO3 and ALK given ------------------
%
%       don't split the lines,
%       matlab does not understand it.
%
if flag == 11
disp('flag = 11, HCO3 and ALK given');
hco3 = hco31;
alk = alk1;
p5  = 1.;
p4  = alk - hco3 + K1 + Kb;
p3  = alk*(Kb+K1)-hco3*(K1+Kb+2.*K2)-Kw+K1*Kb+K1*K2-Kb*bor;
tmp = alk*(Kb*K1+K1*K2)-hco3*((Kb+2.*K2)*K1+2.*Kb*K2+K1*K2);
p2  = tmp +(-K1*Kb*bor-Kw*Kb-K1*Kw+K1*K2*Kb);
tmp = alk*Kb*K1*K2-hco3*(2.*Kb*K1*K2+K2*K1*(Kb+2.*K2));
p1  = tmp +(-K1*K2*Kb*bor-K1*Kw*Kb-K1*K2*Kw);
p0  = -hco3*2.*K2*Kb*K1*K2-K1*K2*Kw*Kb;
p   = [p5 p4 p3 p2 p1 p0];
r   = roots(p);
h   = max(real(r));
h*1.e12;         
dic = hco3 * (1+h/K1+K2/h);
s = dic / (1.+K1/h+K1*K2/h/h);
%   dic = s*(1.+K1/h+K1*K2/h/h);
%   hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end
% ------------ HCO3 and DIC given ------------------
if flag == 12
disp('flag = 12, HCO3 and DIC given');
hco3 = hco31;
dic = dic1;
p2 = hco3/K1;
p1 = hco3-dic;    
p0 = hco3*K2;
p = [p2 p1 p0];
r = roots(p);
h = min(real(r));        % min instead of max !!!!!
h*1.e12;
s = dic / (1.+K1/h+K1*K2/h/h);
dic = s*(1.+K1/h+K1*K2/h/h);
%   hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end
% ------------ CO3 and ALK given ------------------
if flag == 13
disp('flag = 13, CO3 and ALK given');   
co3 = co31;
alk = alk1;
p5  = -co3/K2+1.;
p4  = alk - co3*(K1/K2+(Kb+2.*K2)/K2) + Kb + K1;
tmp = alk*(Kb+K1)-co3*(K1+K1*(Kb+2.*K2)/K2+2.*Kb);
p3  = tmp+(-Kb*bor-Kw+K1*Kb+K1*K2);
tmp = alk*(Kb*K1+K1*K2)-co3*(K1*(Kb+2.*K2)+2.*Kb*K1);
p2  = tmp+(-Kw*Kb-K1*Kb*bor-K1*Kw+K1*K2*Kb);
tmp = alk*Kb*K1*K2-co3*2.*Kb*K1*K2-K1*Kw*Kb;
p1  = tmp+(-K1*K2*Kb*bor-K1*K2*Kw);
p0  = -K1*K2*Kw*Kb;
p   = [p5 p4 p3 p2 p1 p0];
r   = roots(p);
h   = max(real(r));
h*1.e12;
%
% NOTE : calculate dic from dic = co3*(1+h/K2+h^2/K1/K2);
%                  not from dic = hco3*(1+h/K1+K2/h);
%
dic = co3 * (1+h/K2+h^2/K1/K2);
s = dic / (1.+K1/h+K1*K2/h/h);
%   dic = s*(1.+K1/h+K1*K2/h/h);    
hco3 = dic/(1+h/K1+K2/h);
%   co3 = dic/(1+h/K2+h*h/K1/K2);
%   alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1);
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end
% ------------ CO3 and DIC given ------------------
if flag == 14
disp('flag = 14, CO3 and DIC given');
co3 = co31;
dic = dic1;
p2 = co3/K1/K2;
p1 = co3/K2;
p0 = co3-dic;
p = [p2 p1 p0];
r = roots(p);        
h = max(real(r));
h*1.e12;
s = dic / (1.+K1/h+K1*K2/h/h);
%   dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
%   co3 = dic/(1+h/K2+h*h/K1/K2);
alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1); %        6.3096e-09
CO2     = s*1.e6
pCO2    = s*1.e6/Kh
HCO3    = hco3*1.e6
CO3     = co3*1.e6
DIC     = dic*1.e6
ALK     = alk*1.e6
end
% ------------ ALK and DIC given ------------------
if flag == 15
%disp('flag = 15, ALK and DIC given');
alk = alk1;
dic = dic1;
p5  = -1.;        
p4  = -alk-Kb-K1;
p3  = dic*K1-alk*(Kb+K1)+Kb*bor+Kw-Kb*K1-K1*K2;
tmp = dic*(Kb*K1+2.*K1*K2)-alk*(Kb*K1+K1*K2)+Kb*bor*K1;
p2  = tmp+(Kw*Kb+Kw*K1-Kb*K1*K2);
tmp = 2.*dic*Kb*K1*K2-alk*Kb*K1*K2+Kb*bor*K1*K2;
p1  = tmp+(+Kw*Kb*K1+Kw*K1*K2);
p0  = Kw*Kb*K1*K2;
p   = [p5 p4 p3 p2 p1 p0];
r   = roots(p);
h   = max(real(r));
%   test = p5*h^5+p4*h^4+p3*h^3+p2*h^2+p1*h+p0
h*1.e12;
s = dic / (1.+K1/h+K1*K2/h/h);
%   dic = s*(1.+K1/h+K1*K2/h/h);
hco3 = dic/(1+h/K1+K2/h);
co3 = dic/(1+h/K2+h*h/K1/K2);
co2 = s;    % added 11/22/06 to have consistent units in function input & output
pco2 = s/Kh;
%   alk = s*(K1/h+2.*K1*K2/h/h)+Kb*bor/(Kb+h)+Kw/h-h;
% ----------- change units: mumol/kg
%h1 = 10^(-ph1); %
CO2     = s*1.e6;
pCO2    = s*1.e6/Kh;
HCO3    = hco3*1.e6;    
CO3     = co3*1.e6;
DIC     = dic*1.e6;
ALK     = alk*1.e6;

end;

% ======================================================
if phflag == 0;
%   disp('Total pH scale used.');
   phtotal = -log10(h);
   oh = Kw./h;
   phfree  = -log10(h/total2free);
end;

if phflag == 1;
%   disp('Free pH scale used.');
   phfree  = -log10(h);
   oh = Kw./h;   
   phtotal = -log10(h*total2free);
   phsws = -log10(h*sws2free);
end;

return;

