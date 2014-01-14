function[hc]=discretecolorbar(hc,cax,ci,str)
%DISCRETECOLORBAR  Plots a colorbar with discrete variation.
%
%   In filled contour plots one has discrete values of color (or of
%   shading), but the colorbar resulting from Matlab's COLORBAR command    
%   has a continuous spectrum of color.                                   
%                                                                         
%   DISCRETECOLORBAR(HC,CAX,CI) where HC is a handle to a colorbar axis   
%   (i.e. from calling HC=COLORBAR), CAX=[CMIN CMAX] are the color axis   
%   limits, and CI is a vector of contour intervals, redraws the colorbar 
%   with discrete breaks at values CI and axis limits CAX.                
%                                                                         
%   DISCRETECOLORBAR(HC,CAX,CI,STR) where STR is either 'hori' or 'vert'  
%   specifies whether the colorbar should have a horizontal or vertical   
%   orientation; STR defaults to 'vert'.                                  
%                                                                         
%   Make sure to set the color axis of the contour plot is also is set    
%   to CAX [via CAXIS(CAX)].        
%
%   Usage:  discretecolorbar(hc,cax,ci);
%           discretecolorbar(hc,cax,ci,str);
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information 
%   (C) 2001--2007 J.M. Lilly --- type 'help jlab_license' for details  
  
%discretecolorbar(hc,cax,ci,'vert')
%Matlab is shit... I should not have to do this.

if nargin==3
   str='vert';
end
%axes(hc),set(gca,'xtick',[],'ytick',[],'visible','off')
%xlabel(''),ylabel('')
%pos=get(hc,'position');
xal=get(hc,'xaxislocation');
yal=get(hc,'yaxislocation');
%hc=axes('position',pos);
set(hc,'userdata',[])

pos=get(hc,'position');
pos2=get(gca,'position');
delete(hc)
set(gca,'position',pos2);
hc=axes('position',pos);

axes(hc)
mat=(0:127)';
mat=mat/127*(cax(2)-cax(1))+cax(1);
if strcmp(str(1:4),'vert')
   contourf([0 1],mat,[mat mat],ci);
   axis(hc,[0 1 cax(1) cax(2)])
   set(hc,'ytick',ci,'xtick',[])
elseif strcmp(str(1:4),'hori')
   contourf(mat,[0 1],[mat';mat'],ci);
   axis(hc,[cax(1) cax(2) 0 1])
   set(hc,'xtick',ci,'ytick',[])
end

set(hc,'xaxislocation',xal);
set(hc,'yaxislocation',yal);
set(hc,'box','on','layer','top'),grid off
caxis(cax)

%figure,plot(1:10),hc=colorbar;discretecolorbar(hc,cax,ci,'vert');
%caxis([2.5 5.5]) 
%ci=(cax(1):.4:cax(2)]
if nargout==0,clear hc,end





