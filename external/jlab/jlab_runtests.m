function[]=jlab_runtests(str)
%JLAB_RUNTESTS  Run test suite for JLAB package.
%
%   JLAB_RUNTESTS runs automated tests for the JLAB package.
%
%   'jlab_runtests tests' runs all automated tests.
%   'jlab_runtests figures' makes all sample figures.
%   'jlab_runtests' with no arguments does both.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2002--2011 J.M. Lilly --- type 'help jlab_license' for details      

if nargin==0
    str='both';
end
a=now;
global BOOL_JLAB_RUNTEST
global FUNCTION_NUMBER
global FUNCTION_NUMBER_ARRAY

BOOL_JLAB_RUNTEST=[];
FUNCTION_NUMBER_ARRAY=[];

dirname=whichdir('jlab_license');
if iscell(dirname)
    dirname=dirname{1};
end
names=findfiles(dirname,'m');
for i=1:length(names)
    names{i}=names{i}(1:end-2);
end

testfailures=zeros(length(names),1);
for i=1:length(names)
    FUNCTION_NUMBER=i;
    if ~strcmp(names{i},'jlab_runtests')&&~strcmp(names{i},'Contents') %No recursion please
        fid=fopen([dirname '/' names{i} '.m'], 'r');   
        filestr=setstr(fread(fid)');
        fclose(fid);
        if strcmpi(str(1:3),'tes')||strcmpi(str(1:3),'bot')
            if ~isempty(strfind(filestr,'''--t'''))
                try
                    eval([names{i} '(''--t'');']);
                catch    
                    testfailures(i)=1;
                end
            end
        end
        if strcmpi(str(1:3),'fig')||strcmpi(str(1:3),'bot')
            if ~isempty(strfind(filestr,'''--f'''))            
                try
                    disp(['Generating figures for ' names{i} '...'])
                    eval([names{i} '(''--f'');']);
                catch    
                end
            end
        end
    end
end
close all
if strcmpi(str(1:3),'tes')||strcmpi(str(1:3),'bot') 
    disp('---------------------------------------------')
    disp(['JLAB_RUNTESTS --- ' int2str(sum(BOOL_JLAB_RUNTEST)) ' of '  int2str(length(BOOL_JLAB_RUNTEST)) ' tests passed.'])
    if sum(~BOOL_JLAB_RUNTEST)>0
          disp('    Tests in the following routines failed:')
          for i=1:length(names)
              Nfailed=length(find((FUNCTION_NUMBER_ARRAY==i)&(BOOL_JLAB_RUNTEST==0)));
              if Nfailed>0
                  disp(['          ' names{i} ', ' int2str(Nfailed) ' test(s) failed.'])
              end
          end
    end 
    if vsum(testfailures,1)>0
          disp(['    Tests in the following routines did not execute:'])
          for i=1:length(names)
              if testfailures(i)==1
                  disp(['          ' names{i} ' tests did not execute.' ])
              end
          end
    end 
end
b=now;
disp(['JLAB_RUNTESTS took ' num2str((b-a)*24*60) ' minutes.'])
