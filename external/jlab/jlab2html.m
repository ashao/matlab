function[varargout]=jlab2html(varargin)
%JLAB2HTML  Output JLAB help files in HTML format using M2HTML. 
%
%   JLAB2HTML  Uses the fantastic M2HTML toolbox to generate HTML files 
%   from the JLAB help files.
%
%   Then put this into m2html.css:
%
%   body {background-color:#FFFFFD;float:none;font-size:large;width:50em;
%      margin-left:auto; margin-right:auto;margin-top:0em; margin-bottom:0em; 
%      padding: 0em 0em 1em 0em;border:0em; color: black; font-family:utopia;} 
%
%   Then after this you will ftp jmlilly.net.
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details
 
cd ~
cd Desktop/Dropbox/Matlab
%m2html('mfiles','jlab', 'htmldir','doc', 'recursive','on','graph','off','template','frame', 'index','menu','source','off');
m2html('mfiles','jlab', 'htmldir','doc', 'recursive','on','graph','off', 'index','menu','source','off');


