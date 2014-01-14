%
% JGRAPH  Fine-turning and customizing figures
%
%   See also: 
%   CONTENTS, JARRAY, JMATH, JPOLY, *JGRAPH*, JSTRINGS, JSTATS, JSIGNAL,
%   JELLIPSE, JCELL, VTOOLS, JOCEANS, JSPHERE, JSATFUN, JTRIADS, JPAPERS
%   _________________________________________________________________
%
% Specialized plotting functions
%   jimage     - Image plot with scaled cdatamapping and normal ydir.
%   jcontour   - Contouring with labelled and unlabelled contours.
%   uvplot     - Plots u and v components of velocity on the same axis.
%   stickvect  - Plots "stick vectors" for multicomponent velocity time series.
%   provec     - Generate progressive vector diagrams (simple and fancy).
%   hodograph  - Generate hodograph plots (simple and fancy).
%   fastcontour  - Lightning-fast "fake" contouring for large matrices.
%
% Graphical post-processing --- ticks, labels, and aspect ratio
%   xtick         - Sets locations of x-axis tick marks.
%   ytick         - Sets locations of y-axis tick marks.
%   ztick         - Sets locations of z-axis tick marks.
%   fixlabels     - Specify precision of axes labels.
%   ticklen       - Sets tick length of current axis. 
%   ytlpad        - Pads the ytick labels with a leading space.
%   hlines        - Add horizontal lines to a plot.
%   vlines        - Add vertical lines to a plot.
%   dlines        - Add diagonal lines to a plot.
%   latratio      - Set plot aspect ratio for latitude / longitude plot.
%
% Graphical post-processing --- extras
%   discretecolorbar - Plots a colorbar with discrete variation.
%   letterlabels  - For automatically putting letter labels on subplots.
%   timelabel     - Put month, day, or hour labels on a time axes.
%   phasecircle   - Plots little circles to indicate a phase angle.
%
% Graphical post-processing --- lines and fonts
%   linering      - Moves lines through the current line style order.  
%   yoffset       - Offsets lines in the y-direction after plotting.
%   xoffset       - Offsets lines in the x-direction after plotting.
%   linestyle     - Sets color, style, and width properties of lines.
%   fontsize      - Rapidly set title, axes, label, and text fontsizes.
%
% Graphical post-processing --- subplots
%   packcols      - Squeeze together all subplot columns of the current figure.
%   packrows      - Squeeze together all subplot rows of the current figure.
%   packboth      - Squeeze together rows and columns of the current figure.
%   letterlabels  - For automatically putting letter labels on subplots.
%
% Low-level functions
%   axeshandles	  - Returns handles to all axes children.
%   linehandes    - Finds all line and patch handles from a given set of axes.
%
% Simple graphical aliases
%   ylin          - Sets y-axis scale to linear.
%   ylog          - Sets y-axis scale to logarithmic.
%   xlin          - Sets x-axis scale to linear.
%   xlog          - Sets x-axis scale to logarithmic.
%   inticks       - Sets the 'tickdir' property of the current axis to 'in'.
%   outticks      - Sets the 'tickdir' property of the current axis to 'out'.
%   flipx         - Flips the direction of the x-axis.
%   flipy         - Flips the direction of the y-axis.
%   leftaxis      - Sets the 'yaxislocation' property of the current axis to 'left'.
%   rightaxis     - Sets the 'yaxislocation' property of the current axis to 'right'.
%   topaxis       - Sets the 'xaxislocation' property of the current axis to 'top'.
%   bottomaxis    - Sets the 'xaxislocation' property of the current axis to 'bottom'.
%   noxlabels     - Remove some or all x-axis tick mark labels.
%   noylabels     - Remove some or all y-axis tick mark labels.
%   boxon         - Sets 'box' property to 'on'.
%   boxoff        - Sets 'box' property to 'off'.
%   monthlyticks  - Set x-axis appropriate for months.
%   nocontours    - Removes contours from a CONTOURF plot.
%   axestop       - Sets the 'layer' property of the current axis to 'top'.
%   flipmap       - Flips the current colormap upside-down.
%   land          - Sets orientation to 'landscape'.
%   tall          - Sets orientation to 'tall'.
%   port          - Sets orientation to 'portrait'.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2004--2011 J.M. Lilly --- type 'help jlab_license' for details        

help Jgraph

if 0          
           %Specialized plotting functions
             jimage     %- Image plot with scaled cdatamapping and normal ydir.
             jcontour   %- Contouring with labelled and unlabelled contours.
             uvplot     %- Plots u and v components of velocity on the same axis.
             stickvect  %- Plots "stick vectors" for multicomponent velocity time series.
             provec     %- Generate progressive vector diagrams (simple and fancy).
             hodograph  %- Generate hodograph plots (simple and fancy).
             fastcontour  %- Lightning%-fast "fake" contouring for large matrices.
          
           %Graphical post%-processing %-%-%- ticks, labels, and aspect ratio
             xtick         %- Sets locations of x%-axis tick marks.
             ytick         %- Sets locations of y%-axis tick marks.
             ztick         %- Sets locations of z%-axis tick marks.
             fixlabels     %- Specify precision of axes labels.
             ticklen       %- Sets tick length of current axis. 
             ytlpad        %- Pads the ytick labels with a leading space.
             hlines        %- Add horizontal lines to a plot.
             vlines        %- Add vertical lines to a plot.
             dlines        %- Add diagonal lines to a plot.
             latratio      %- Set plot aspect ratio for latitude / longitude plot.
          
           %Graphical post%-processing %-%-%- extras
             discretecolorbar %- Plots a colorbar with discrete variation.
             letterlabels  %- For automatically putting letter labels on subplots.
             timelabel     %- Put month, day, or hour labels on a time axes.
             phasecircle   %- Plots little circles to indicate a phase angle.
          
           %Graphical post%-processing %-%-%- lines and fonts
             linering      %- Moves lines through the current line style order.  
             yoffset       %- Offsets lines in the y%-direction after plotting.
             xoffset       %- Offsets lines in the x%-direction after plotting.
             linestyle     %- Sets color, style, and width properties of lines.
             fontsize      %- Rapidly set title, axes, label, and text fontsizes.
          
           %Graphical post%-processing %-%-%- subplots
             packcols      %- Squeeze together all subplot columns of the current figure.
             packrows      %- Squeeze together all subplot rows of the current figure.
             packboth      %- Squeeze together rows and columns of the current figure.
             letterlabels  %- For automatically putting letter labels on subplots.
          
           %Low%-level functions
             axeshandles	  %- Returns handles to all axes children.
             linehandes    %- Finds all line and patch handles from a given set of axes.
          
           %Simple graphical aliases
             ylin          %- Sets y%-axis scale to linear.
             ylog          %- Sets y%-axis scale to logarithmic.
             xlin          %- Sets x%-axis scale to linear.
             xlog          %- Sets x%-axis scale to logarithmic.
             inticks       %- Sets the 'tickdir' property of the current axis to 'in'.
             outticks      %- Sets the 'tickdir' property of the current axis to 'out'.
             flipx         %- Flips the direction of the x%-axis.
             flipy         %- Flips the direction of the y%-axis.
             leftaxis      %- Sets the 'yaxislocation' property of the current axis to 'left'.
             rightaxis     %- Sets the 'yaxislocation' property of the current axis to 'right'.
             topaxis       %- Sets the 'xaxislocation' property of the current axis to 'top'.
             bottomaxis    %- Sets the 'xaxislocation' property of the current axis to 'bottom'.
             noxlabels     %- Remove some or all x%-axis tick mark labels.
             noylabels     %- Remove some or all y%-axis tick mark labels.
             boxon         %- Sets 'box' property to 'on'.
             boxoff        %- Sets 'box' property to 'off'.
             monthlyticks  %- Set x%-axis appropriate for months.
             nocontours    %- Removes contours from a CONTOURF plot.
             axestop       %- Sets the 'layer' property of the current axis to 'top'.
             flipmap       %- Flips the current colormap upside%-down.
             land          %- Sets orientation to 'landscape'.
             tall          %- Sets orientation to 'tall'.
             port          %- Sets orientation to 'portrait'.    
end

