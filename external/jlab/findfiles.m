function[varargout]=findfiles(varargin)
%FINDFILES  Returns all files in a directory with a specified extension.
%
%   FILES=FINDFILES(DIRNAME,EXT) where DIRNAME is a full directory path 
%   name and EXT is a file extension returns a cell arrays FILES of all 
%   files in directory DIRNAME having extension EXT. 
%
%   The leading '.' should be omitted from EXT, and the trailing '/' 
%   should be omitted from DIRNAME.
% 
%   For example, FILES=FINDFILES('/Users/lilly/Home/data','m') returns a
%   cell array FILES of all m-files in the specified directory.
%
%   FINDFILES(DIRNAME,'') returns files having no extension.  This can
%   be used to return directory names.
%   ___________________________________________________________________
%
%   Additional options
%
%   FINDFILES(...,'include',STR1) additionally returns only those file
%   names which inlude the text STR1.
%
%   FINDFILES(...,'exclude',STR2) excludes those files names which 
%   include the text STR2.
%
%   The include and exclude options must be input after the first two
%   variables.  Both options may be input simultaneously.
%   ___________________________________________________________________
%
%   Usage: files=findfiles(pwd,'m');
%          files=findfiles(dirname,ext);    
%          files=findfiles(dirname,ext,'include','help');
%          files=findfiles(dirname,ext,'include','help','exclude','jlab');
%
%   'findfiles --t' runs a test.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2006--2009 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    findfiles_test,return
end
 
dirname=varargin{1};
ext=varargin{2};

vars=varargin(3:end);
excludestr=[];
includestr=[];
while ~isempty(vars)
    if length(vars)<2
        error('Wrong number of arguments.')
    else
        str=vars{1};
        if strcmp(str(1:3),'exc')
            excludestr=vars{2};
        elseif strcmp(str(1:3),'inc')
            includestr=vars{2};
        end
        vars=vars(3:end);
    end
end

dirstruct=dir(dirname);
files=cell(length(dirstruct),1);
for i=1:length(dirstruct)
    files{i}=dirstruct(i).name;
end

bool=false(length(files),1);
N=length(ext);
if N==0
    for i=1:length(files)
        bool(i,1)=isempty(strfind(files{i},'.'));
    end
else
    for i=1:length(files)
        if length(files{i})>N
            bool(i,1)=strcmp(files{i}(end-N:end),['.' ext]);
        end
        if ~isempty(includestr)
            bool(i,1)=bool(i,1)&~isempty(strfind(files{i}(1:end-N-1),includestr));
        end
        if ~isempty(excludestr)
            bool(i,1)=bool(i,1)&isempty(strfind(files{i}(1:end-N-1),excludestr));
        end
    end
end

index=find(bool);
if ~isempty(index)
    files=files(index);
else 
    files=[];
end

varargout{1}=files;

function[]=findfiles_test

dirname=whichdir('jlab_license');
if iscell(dirname)
    dirname=dirname{1};
end
files=findfiles(dirname,'m');
bool=0;
for i=1:length(files)
    if strcmp('jlab_license.m',files{i})
        bool=1;
    end
end
reporttest('FINDFILES found jlab_license', bool)
files=findfiles(dirname,'m','include','jlab_license');
reporttest('FINDFILES with include flag', aresame(files{1},'jlab_license.m'))

