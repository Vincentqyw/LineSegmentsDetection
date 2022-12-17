function L = sampleLine(rho,theta,m,n,pp)
%Returns locations of line samples together with the locations
%of the associated pixels in the image.
%(rho, theta) represents the line as per Fig. 1 of Tal (2012):
%rho = x*cos(theta) + y*sin(theta) <--> ax + by + c = 0 where
%a = cos(theta), b = sin(theta), c = -rho.
%pp is the principal point, in pixels.
%(m, n) is the width and height of the image in pixels.
%The returned matrix L is N x 4, where N is the number of
%samples on the line.  The 4 columns represent (x,y,xp,yp) 
%where (x,y) indicates the real-valued location of the line
%sample and (xp,yp) represents the integer-valued location
%of the associated pixel.  Note that the y index increases down.

%Translate to standard representation of line
a = cos(theta);
b = sin(theta);
c = -rho;

%Location of centre of image
x0 = pp(1);
y0 = pp(2);

%Find image bounds
left = 1 - x0;
right = m - x0;
top = 1 - y0;
bottom = n - y0;

%Identify pixels within 2 pixels of line
maxpts = round(5*sqrt(m^2+n^2));
%Pre-allocate for efficiency
xp = zeros(maxpts,1);yp = zeros(maxpts,1);
idx = 1;
%Closer to horizontal - march along x.
if theta > pi/4 && theta < 3*pi/4
    x = left:right;
    yu = floor(min(bottom,(-c+2 - x*a)/b));
    yl = ceil(max(top,(-c-2 - x*a)/b));
    for i = 1:m
        npts = max(0,yu(i)-yl(i) + 1);
        xp(idx:idx + npts-1) = x(i);
        yp(idx:idx + npts-1) = yl(i):yu(i);
        idx = idx + npts;
    end
else %Closer to vertical - march along y
    y = top:bottom;
    %if theta<pi/2, increasing rho increases x.  
    %Otherwise, decreasing rho increases x.
    s = sign(cos(theta)); 
    xu = floor(min(right,(-c+2*s - y*b)/a));
    xl = ceil(max(left,(-c-2*s - y*b)/a));
    for i = 1:n
        npts = max(0,xu(i)-xl(i) + 1);
        yp(idx:idx + npts-1) = y(i);
        xp(idx:idx + npts-1) = xl(i):xu(i);
        idx = idx + npts;
    end
end
xp = xp(1:idx-1);
yp = yp(1:idx-1);

L=[xp,yp];
% %Find orthogonal projections to line
% x = b*(b*xp-a*yp)-a*c;
% y = a*(-b*xp+a*yp)-b*c;
% 
% L = [x,y,xp,yp];
% %Sort line samples
% if theta > pi/4 & theta < 3*pi/4
%     [x,idx] = sort(x);
% else
%     [y,idx] = sort(y);
% end  
% L = L(idx,:);
% npts = size(L,1);
% npts;
% figure;
% eps = endpoints(rho,theta,left,right,top,bottom);
% plot(eps(:,1),eps(:,2),'r');
% hold on;
% plot(xp,yp,'.');
% axis ij;
% axis equal;
% axis([left,right,top,bottom]);


