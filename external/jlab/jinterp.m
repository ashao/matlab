function[yi]=jinterp(x,y,xi,str)
%JINTERP Matrix-matrix 1-D interpolation.
%
%   YI=JINTERP(X,Y,XI), returns the linear interpolation of Y onto XI     
%   based on the functional relationship Y(X).                            
%                                                                         
%   Unlike INTERP1, JINTERP allows X,Y, and XI to be either vectors or    
%   matrices. If more than one argument is a matrix, those matrices must  
%   be of the same size. YI is a matrix if any input argument is a        
%   matrix. All vectors should be column vectors and all matrices should  
%   have data in columns.                                                 
%                                                                         
%   Also, only data points of XI in between the maximum and minimum       
%   values of X are interpolated.                                         
%                                                                         
%   This useful, for example, in interpolating section data with          
%   nonuniform pressure levels onto standard levels.                      
%                                                                         
%   See also INTERP1.                                                      
%   __________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2008 J.M. Lilly --- type 'help jlab_license' for details    
if nargin~=4
   str='linear';
end

%convert row vectors to column vectors
if size(x,1)==1
   x=conj(x');
end
if size(y,1)==1
   y=conj(y');
end
if size(xi,1)==1
   xi=conj(xi');
end

%ensure sizes are compatible
Lx=size(x,2);
Ly=size(y,2);
Lxi=size(xi,2);
maxL=max([Lx Ly Lxi]);
bool=(Lx==1|Lx==maxL)&(Ly==1|Ly==maxL)&(Lxi==1|Lxi==maxL);
if ~bool,
   error('Arguments are not of compatible size')
end

%convert vectors to matrices
if Lx==1
   x=x*ones(1,maxL);
end
if Ly==1
   y=y*ones(1,maxL);
end
if Lxi==1
   xi=xi*ones(1,maxL);
end

yi=nan*ones(size(xi,1),maxL);

%check x for monotonicity
mdx=min(min(diff(x)));
if mdx<=0
   disp('Ensuring monotonicity of X by adding noise and sorting.')
   mx=min(min(x));
   x=x+randn(size(x))/1000/1000;
   x=sort(x,1);
end

for i=1:size(x,2)
 	colmin=min(x(isfinite(x(:,i)),i));
	colmax=max(x(isfinite(x(:,i)),i));

%	a=min(find(xi(:,i)>=colmin));
%	b=max(find(xi(:,i)<=colmax));
	index=find(xi(:,i)>=colmin&xi(:,i)<=colmax&isfinite(xi(:,i)));
	
	if ~isempty(index)>=0,
		yi(index,i)=interp1(x(:,i),y(:,i),xi(index,i),str);
	end
end




