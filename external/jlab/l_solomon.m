%L_SOLOMON Script for loading Solomon Islands earthquake record

%Data downloaded from IRIS
%Incorporated Research Institutions for Seismology
%http://www.iris.edu/dms/wilber.htm
%Using the WILBER II browser

%See http://www.iris.edu/manuals/SEED_appA.htm for data description
%B=Broadband
%H=High Gain Seismometer
%Z=Vertical, N=North, E=East  

%Event: 1991/02/09 16:18:58.3
%Catalog: MHDF  Mag: 6.9  Type: MW  Contributor: HRV  
%Lat: -9.93  Lon: 159.14  Depth: 10.00
%Description: SOLOMON ISLANDS  Source: FARM

%Station: PAS - Pasadena, California, USA 
%Network: TS - TERRAscope (Southern California Seismic Network)
%Lat: 34.15 Lon: -118.17 Elev: 314.00 
%Event Name: 19910209_161858.3.farm
%Available Channels: BHE,BHN,BHZ,LHE,LHN,LHZ

%Description of data format:
%http://www.iris.edu/software/sac/
%Header information:
%http://www.iris.edu/manuals/sac/manual.html

%How to cite is at http://www.iris.edu/hq/publications/iris_citations.
%
%To cite the DMS:
%
%The facilities of the IRIS Data Management System, and specifically the
%IRIS Data Management Center, were used for access to waveform and metadata
%required in this study.
%
%The IRIS DMS is funded through the National Science Foundation and 
%specifically the GEO Directorate through the Instrumentation and 
%Facilities Program of the National Science Foundation under Cooperative 
%Agreement EAR-0552316.

%/*************************************************************************
dirname=jlab_settings('dirnames.solomon');

%I manually cut the headers and re-saved the data in the following files.
%Also, in doing so I also trimmed the last lines which were incomplete
x=load('-ascii',[dirname '/bhe.txt']);
y=load('-ascii',[dirname '/bhn.txt']);
z=load('-ascii',[dirname '/bhz.txt']);

vindex(x,y,z,1,2);

x=x-vmean(x,1);
y=y-vmean(y,1);
z=z-vmean(z,1);

%Times from headers
%E: 1991        40        16        29        50      698  
%N: 1991        40        16        28        58      397   
%Z: 1991        40        16        29         2      398  

t1e=datenum(1991,1,1)-1+40+16/24+29/60/24 + 50/60/60/24 + 698/1000/60/60/24;
t1n=datenum(1991,1,1)-1+40+16/24+28/60/24 + 58/60/60/24 + 397/1000/60/60/24;
t1z=datenum(1991,1,1)-1+40+16/24+29/60/24 +  2/60/60/24 + 398/1000/60/60/24;

datestr(t1e)  %09-Feb-1991 16:29:50
datestr(t1n)  %09-Feb-1991 16:28:58
datestr(t1z)  %09-Feb-1991 16:29:02

%t1e=t1e+652.328/24/60/60;
%t1n=t1n+600.027/24/60/60;
%t1z=t1z+604.028/24/60/60;

tx=t1e + 0.25/24/3600 *[0:length(x)-1]';
ty=t1n + 0.25/24/3600 *[0:length(y)-1]';
tz=t1z + 0.25/24/3600 *[0:length(z)-1]';

%Why does increment ssy 0.05 sec when it appears to be 0.25 sec??

num=[1:length(x)]';

iy=max(find(ty<=tx(1)));
iz=max(find(tz<=tx(1)));
vindex(ty,y,iy:length(y),1);
vindex(tz,z,iz:length(z),1);

24*3600*(tz(1)-tx(1))
24*3600*(ty(1)-tx(1))

%Very close
num=tz;

%These have different lengths!!!
vindex(num,x,y,z,1:15059,1);

x=x-vmean(x,1);
y=y-vmean(y,1);
z=z-vmean(z,1);



description='Solomon Islands Earthquake, Feb. 9 1991, recorded at Pasadena (PAS), downloaded from IRIS';
link='http://www.iris.edu/dms/wilber.htm';
creator='J. M. Lilly -- l_solomon script';
timestamp=datestr(now);
copyright='Data redistributed with permission.  The facilities of the IRIS Data Management System, and specifically the IRIS Data Management Center, were used for access to waveform and metadata required in this study.  The IRIS DMS is funded through the National Science Foundation and specifically the GEO Directorate through the Instrumentation and Facilities Program of the National Science Foundation under Cooperative Agreement EAR-0552316.'

clear solomon
matsave solomon description link creator timestamp copyright num x y z 
%\*************************************************************************




if 0


%That's as far as I got on the data I downloaded myself, try Sofia's verions

data1=load('-ascii',[dirname '/paslhz2.dat']);
data2=load('-ascii',[dirname '/paslht2.dat']);
data3=load('-ascii',[dirname '/pashlhr2.dat']);



x=[data1(:,1) data2(:,1) data3(:,1)];
t=[1:length(x)]';

matsave solomon t x
end