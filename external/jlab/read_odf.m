function[varargout]=read_odf(varargin)
%READ_ODF
%
%   READ_ODF
%
%   'read_odf --t' runs a test.
%
%   Usage: []=read_odf();
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    read_odf_test,return
end
 
function[]=read_odf_test
 
%reporttest('READ_ODF',aresame())
