function[varargout]=circridge(varargin)
%CIRCRIDGE  Circular wavelet ridges
%
%   CIRCRIDGE
%
%   'circridge --t' runs a test.
%
%   Usage: []=circridge();
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2011 J.M. Lilly --- type 'help jlab_license' for details
 
if strcmp(varargin{1}, '--t')
    circridge_test,return
end
 
%   [IR,JR,XR,FR]=RIDGEWALK(W,FS) where W is a wavelet transform matrix






function[]=circridge_test
 
%reporttest('CIRCRIDGE',aresame())
