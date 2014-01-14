function[h]=dlines(varargin)
%DLINES   Add diagonal lines to a plot.
%
%   DLINES(M,STYLE) puts diagonal lines with slopes M having style
%   described by the string STYLE. 
%
%   By default, DLINES uses the lower left corner of the current axis
%   as the line origin for M>0, and the lower right corner if M<0.
%
%   DLINES(M,Y0,STYLE) instead uses Y0 as the intercept of the left-
%   hand side of the current x-axis.
%
%   (Note that the above descriptions refer to the axis sides with
%   the 'normal', rather than 'reverse', directions.)
%
%   The SYTLE strings follow the format specified in LINESTYLE.  Thus
%   STYLE='2b--' draws blue dotted lines of width 2.
%
%   STYLE is optional and defauts to STYLE='g--', a green dashed line.
%
%   H=DLINES(...) returns an array of handles H to the lines.
%
%   See also HLINES, VLINES.	
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2008 J.M. Lilly --- type 'help jlab_license' for details    

m=varargin{1};
if ischar(m)
    if strcmp(m(1:3),'--f')
        dlines_figure;return
    elseif strcmp(m(1:3),'--t')
        return
    end
end

if ~ischar(varargin{end})
    sty='g--';
else
    sty=varargin{end};
    varargin=varargin(1:end-1);
end

if length(varargin)==2
    b=varargin{2};
else
    b=[];
end


hold on

m=m(:);  
ax=axis;
[h,xo,yo,x1]=vzeros(length(m),1);

for i=1:length(m)
    if isempty(b)
        if m(i)>0
            xo(i)=ax(1);
            yo(i)=ax(3);
            x1(i)=ax(2);
        else
            xo(i)=ax(2);
            yo(i)=ax(3);
            x1(i)=ax(1);
        end
    else
       xo(i)=ax(1); 
       yo(i)=m(i).*xo(i)+b(i);
       x1(i)=ax(2);
    end
end

for i=1:length(m)
   h(i)=plot([xo(i)+sqrt(-1)*yo(i);x1(i)+sqrt(-1)*(m(i)*(x1(i)-xo(i))+yo(i))]);
end

linestyle(h,sty);
axis(ax)

if nargout==0
   clear h
end


function[]=dlines_figure
figure
subplot(2,2,1)
axis([-1 0 0 1]),dlines(-1,'b'),dlines(1,'r')
subplot(2,2,2),
axis([0 1 0 1]),dlines(-1,'b'),dlines(1,'r')
subplot(2,2,3)
axis([-1 0 -1 0]),dlines(-1,'b'),dlines(1,'r')
subplot(2,2,4)
axis([0 1 -1 0]),dlines(-1,'b'),dlines(1,'r')


