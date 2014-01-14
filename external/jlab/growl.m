function[varargout]=growl(varargin)
%GROWL  Display text message with the Growl notification system.
%
%   GROWL(STR) displays the text message STR using Growl for Macs.
%
%   Both Growl and growlnotify must be installed on your system.  They can
%   be downloaded from http://growl.info
%   
%   Growl is available for Macs only.
%
%   GROWL with no arguments will report that a job has been completed.
%
%   By default, GROWL will also say the text message with the computer's
%   voice.  To suppress this, use GROWL(STR,'quiet');
%
%   Usage: growl
%          growl('Hello world')  
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details
 
if nargin==0
    str='Your job is done, Master';
else
    str=varargin{1};
end

str2='loud';
if nargin==2
    str2=varargin{2};
end

unix(['/usr/local/bin/growlnotify -a Matlab -m "' str '"']);
if isempty(strfind(str2,'qui'))
    unix(['say ' str]);
end
 
function[]=growl_test
 
%reporttest('GROWL',aresame())
