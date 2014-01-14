function[]=linestyle(varargin)
%LINESTYLE  Sets color, style, and width properties of lines.
%	    	    	    
%   LINESTYLE provides an efficient way to set common properties of
%   large groups of lines, and further allows the user to quickly
%   apply predefined sets of line colors, styles, and widths designed
%   for different purposes.
%
%   'LINESTYLE 1k-- 2.5g- m-.' sets the first line to width 1, color
%   black, and dashed style; the second to width 2.5, color green, and
%   solid style; and the third to width 1, color magenta, and dash-
%   dotted style.
%  
%   In general the format is 
%
%       LINESTYLE STR1 STR2 STR3 ...   
%
%   where each of the STRs may contain a number, specifying the width
%   of the Ith line; a letter, specifying the color; and a style
%   string.  Any two of these are optional, with a unit width solid
%   black line being the default.
%
%   By default, styles are looped if the number of STRs input is less
%   than the number of lines in the current plot. The input
%
%       LINESTYLE STR1 STR2 ... STRN +++
%
%   causes the last style input, STRN, to be repeated instead.      
%   _________________________________________________________________
%
%   Locking and unlocking 
%
%   LINESYTLE LOCK and LINESTYLE UNLOCK lock and unlock all axes in
%   the current figure.  When LOCK is on, calls to LINESTYLE or
%   LINERING are applied to all lines in the current figure.
%   _________________________________________________________________
%
%   Handle specification
%  
%   LINESTYLE -H HAN STR1 STR2 .... applies the formatting only to
%   the line handles contained in handle array HAN.
%
%   LINESTYLE(HAN,STR) also works.
%
%   HAN may also be the handle to a group of contours.
%   _________________________________________________________________
%
%   Customization
%
%   LINESTYLE permits the use of user-defined colors, such as 
%   grayscale (capital letters A--K in order of increasing darkness). 
%   See jlab_settings for more details.
%
%   LINESTYLE(NAME) or LINESTYLE NAME applies the style set NAME,
%   specified in the file JLAB_SETTINGS.  This allows one to predefine
%   useful groupings of colors, styles, and widths. Type LINESTYLE
%   with no arguments to see a list of current style sets. 
%   _________________________________________________________________
%
%   See also LINERING, JLAB_SETTINGS, FONTSIZE.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2009 J.M. Lilly --- type 'help jlab_license' for details        
  
if nargin==1
    if strcmp(varargin{1},'--t')
        return
    end
end

h=[];
linestyles=jlab_settings('linestyles');
bextendlast=0;
bloop=1;
bhandleinput=0;
na=nargin;
num=1;

nflagsin=0;
hflagin=0;

vars=varargin;

if na==0
   %if no arguments, just display available linestyle sets
   disp(' ');
   disp('Available linestyle sets:');
   disp(' ');
   disp(linestyles)
   disp(' ');
   return
end


if na==1
   %Special cases 
   if strcmp(vars{1},'lock')
        setappdata(gcf,'linestylelock',1)
	return
   elseif strcmp(vars{1},'unlock')
        setappdata(gcf,'linestylelock',0)
	return
   end
end

if ischar(vars{1})
    if length(vars{1})>1
        if strcmp(vars{1}(1:2),'-h')
           %-h flag for handle input
           hflagin=1;
           nflagsin=1;

           bhandleinput=1;
           h=vars{2};
           if ischar(h)
               eval(to_grab_from_caller(2))
               eval(['h=' varargin{2} ';']);
           end
           %h=flipud(h(:));
           na=na-2;       
           vars=vars(3:end);
        end
    end
else
    if na==2
        h=vars{1};
        bhandleinput=1;
        na=na-1;
        vars=vars(2:end);
    end
end

axes_original=gca;
if isempty(h) 
    h=flipud(linehandles(gca));
else
    %Change to axes to which the given handles belong
    hparent=get(h(1),'parent');
    if strcmp(get(hparent,'type'),'axes')
        axes(hparent)
    else
        axes(get(hparent,'parent'))
    end
end

if length(h)==1
    if strcmp(get(h,'type'),'hggroup')
        h=get(h,'children');
    end
end

if na==1
    if isfield(linestyles,vars{1})
        [widthcell,stylecell,colorcell]=linestyleapply(h,linestyles,vars{1});
        bloop=0;  %supress loop if name input
    end
end
%'hello'
%return
% 
% %/********************************************************
% %Block for applying LINESTYLE to input line handle
% if nargin==2
%   if ~ischar(vars{1});
%     h=vars{1};
%     sty=vars{2};
%     if isfield(linestyles,sty)
%        [widthcell,stylecell,colorcell]=linestyleapply(linestyles,sty);
%        bloop=0;  %supress loop if name input
%     else
%        [widthcelli,stylecelli,colorcelli]=linestyleparse(sty);    
%        widthcell{1}=widthcelli;
%        stylecell{1}=stylecelli;
%        colorcell{1}=colorcelli;
%     end
%     
%     bextendlast=1;
%     bhandleinput=1;
%     %bloop=0;  %supress loop if line handles are input
%   end
% end
% %\********************************************************
%    
% 



if isappdata(gca,'lineringpointer');
  num=getappdata(gca,'lineringpointer');   %Remember current linering
end

if bloop
   %This is specifically for when called like linestyle(h,'2g 1g 1r')
   %since normally Matlab would see this as one string, not three
   if na==1
      temp=vars{1};
      indexspaces=strfind(temp,' ');
      if ~isempty(indexspaces)
         a=[1 indexspaces+1];
         b=[indexspaces-1 length(temp)];
         vars=[];
         jj=0;
         for i=1:length(a);
             if ~isempty(temp(a(i):b(i)))
                 jj=jj+1;
                vars{i}=temp(a(i):b(i));
             end
         end
         na=length(vars);
      end
   end
   for i=1:na
     temp=vars{i};
     if strcmp(temp,'+++')  %ending in continuation symbol
         bextendlast=1;
     else
        [widthcelli,stylecelli,colorcelli]=linestyleparse(temp);
        widthcell{i}=widthcelli;
        stylecell{i}=stylecelli;
        colorcell{i}=colorcelli;
     end  
   end
end

N=length(h);
M=length(stylecell);
if bextendlast
  %Extend by looping continuing last
  for j=M+1:N
      colorcell{j}=colorcell{M}; 
      stylecell{j}=stylecell{M}; 
      widthcell{j}=widthcell{M}; 
  end
else
  %Extend by looping whole structure 
  while N>M
      colorcell((1:M)+M)=colorcell;
      stylecell((1:M)+M)=stylecell; 
      widthcell((1:M)+M)=widthcell;  
      M=length(stylecell);
  end
end

%Truncate if too long
colorcell=colorcell(1:N)';
stylecell=stylecell(1:N)';
widthcell=widthcell(1:N)';

locked=0;      %Check to see if current figures axes are locked
if isappdata(gcf,'linestylelock')
   if getappdata(gcf,'linestylelock')
        locked=1;	
   end
end

[boolline,boolpatch]=handletype(h);

if ~locked || bhandleinput    %Just apply to current axes
  linering(1*sqrt(-1));%Go to first position in linering
  set(h(boolline),{'color'},colorcell(boolline),{'linewidth'},widthcell(boolline),{'linestyle'},stylecell(boolline))
  set(h(boolpatch),{'edgecolor'},colorcell(boolpatch),{'linewidth'},widthcell(boolpatch),{'linestyle'},stylecell(boolpatch))
  linering(num*sqrt(-1));%Return to original
elseif locked && (~bhandleinput)
  h1=axeshandles(gcf);	
  for i=1:length(h1)   %Loop over all axes
       linering(1*sqrt(-1));%Go to first position in linering
       h=linehandles(h1(i));
       [boolline,boolpatch]=handletype(h);
       set(h(boolline),{'color'},colorcell(boolline),{'linewidth'},widthcell(boolline),{'linestyle'},stylecell(boolline))
       set(h(boolpatch),{'edgecolor'},colorcell(boolpatch),{'linewidth'},widthcell(boolpatch),{'linestyle'},stylecell(boolpatch))
       linering(num*sqrt(-1));%Return to original
  end
end

axes(axes_original)

function[boolline,boolpatch]=handletype(h);

htype=get(h,'type');
if ~iscell(htype)
    htypetemp=htype;
    clear htype
    htype{1}=htypetemp;
end

boolline=false(size(h));
boolpatch=false(size(h));
for j=1:length(h)
    boolline(j)= strcmp(htype{j},'line');
    boolpatch(j)= strcmp(htype{j},'patch');   
end



function[widthcell,stylecell,colorcell]=linestyleparse(str)
%Parse the linestyle input string for cell specifications

defaultwidth=1;
defaultcolor='k';
defaultstyle='-';
colors=jlab_settings('colors');
bool=true(size(str));

%/********************************************************
%linewidth
index=find(real(str)>47&real(str)<58); %find str is 0-9
if isempty(index)
   widthcell=defaultwidth;
else
   width=str2num(str(index(1):index(end)));  %encompasss '.'
   if isempty(width)
     error('Problem with linestyle string.')
   else
     widthcell=width;
     
     %remove these entries from string
     bool(index(1):index(end))=0;     
     str=str(bool);
     bool=true(size(str));
   end
end
%\********************************************************


%/********************************************************
%linecolor
index=find(real(str)>64&real(str)<123&~...
      (real(str)==real('s')|...
       real(str)==real('d')|...
       real(str)==real('v')|...
       real(str)==real('p')|...
       real(str)==real('h')|...
       real(str)==real('o')|...
       real(str)==real('x')));    ; %find str is a-Z not  s d v p h o x 
if isempty(index)
   colorcell=defaultcolor;
else
   color=str(index); 
   if length(color)>1
     error('More than one color specified.')
   else
     colorcell=getfield(colors,color);
     
     %remove these entries from string
     bool(index)=0;
     str=str(bool);
   end
end
%\********************************************************   

%/********************************************************
%linestyle
if isempty(str)
   stylecell=defaultstyle;
else
   stylecell=str;
end
%\********************************************************   

function[widthcell,stylecell,colorcell]=linestyleapply(h,linestyles,name)
%Apply contents of LINESTYLE definitions in jlab_settings to cell arrays
  
  
% %account for the fact that LINERING flips the linering
% %to put the first line on the top
% if isappdata(gca,'lineringflipped')
%    if getappdata(gca,'lineringflipped')
% 	h=flipud(h);
%    end
% end
  
colors=jlab_settings('colors');
%h=flipud(linehandles(gca));

colorx=[];
stylex=[];
widthx=[];


if ~isfield(linestyles,name)
  error('That is not a valid linestyle set name.')
else
  colorx=getfield(linestyles,name,{1});
  colorx=colorx{1};
  stylex=getfield(linestyles,name,{2});
  stylex=stylex{1};
  widthx=getfield(linestyles,name,{3});
  widthx=widthx{1};
end

%make the vectors the same length as the handles by cycling
brepeat=0;
N=length(colorx)-3;
if ~isempty(colorx)
   colorx=row2col(colorx);
   if length(colorx)>3
	if strcmp(colorx(end-2:end)','...')
	     colorx(end-2:length(h),:)=colorx(end-3,:);
	     brepeat=1;
	end
   end 		    
   while size(colorx,1)<length(h)
	 colorx=[colorx;colorx];
   end
   colorx=colorx(1:length(h),:);
end

if ~isempty(stylex)
   if size(stylex,1)==1
      stylex=row2col(stylex);
   end
   if length(stylex)==N && brepeat
      stylex(end:length(h),:)=stylex(end,:);
   else
       while size(stylex,1)<length(h)	
           stylex=[stylex;stylex];
       end
   end
   stylex=stylex(1:length(h),:);
end

if ~isempty(widthx)
   widthx=row2col(widthx);
   if length(widthx)==N && brepeat
      widthx(end:length(h),:)=widthx(end,:);
   else
       while size(widthx,1)<length(h)	
           widthx=[widthx;widthx];
       end
   end
   widthx=widthx(1:length(h),:);
end

%convert color into numeric values
temp=colorx;
clear colorx
colorx=[];
for i=1:length(temp)
    colorx(i,:)=getfield(colors,temp(i));
end


%now put into cell arrays, because that's what matlab understands
colorcell=[];
stylecell=[];
widthcell=[];

for i=1:length(h)
    if ~isempty(colorx)
	colorcell{i}=colorx(i,:);
    end    
    if ~isempty(stylex)
	stylecell{i}=stylex(i,:);
    end
    if ~isempty(widthx)
	widthcell{i}=widthx(i);
    end
end

%kludge
%h=flipud(h);
if ~isempty(widthcell),set(h,{'linewidth'},widthcell'),end
if ~isempty(stylecell),set(h,{'linestyle'},stylecell'),end
if ~isempty(colorcell),set(h,{'color'},colorcell'),end



