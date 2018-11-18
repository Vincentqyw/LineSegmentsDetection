
% Function to plot gradient direction vectors at both pixel
%   and subpixel locations:

function[] = plot_edgelmap(img,edgemap,xzeromap,yzeromap,g1dir,...
                            newarrcol,leg1,leg2,imtitle);

[M,N]   = size(edgemap);

[I,J]   = find(edgemap==255);
K       = find(edgemap==255);

X       = J;
Y       = I;

% Set arrowlength, then define: 
%       (Uarrow,Varrow) = arrow lengths in (x,y) directions
%       (Xarrow,Yarrow) = arrow tail points 
arrowlen = 1;
Graddir = g1dir(K) + pi/2;
Uarrow  = arrowlen * cos(Graddir);
Varrow  = -arrowlen * sin(Graddir);

Xtp     = X - 0.5*Uarrow;
Ytp     = Y - 0.5*Varrow;

Xarrow  = xzeromap(K) - 0.5*Uarrow;
Yarrow  = yzeromap(K) - 0.5*Varrow;

figure; 
imagesc(img);
axis image;
colormap(gray);
hold on;
p1  = plot(0,0);
set(p1,'Color',[0 0 1]);
set(p1,'Visible','off');
p2  = plot(0,0);
set(p2,'Color',newarrcol);
set(p2,'Visible','off');
%q1 = quiver(Xtp,Ytp,Uarrow,Varrow,0); JE Jul08 - suppress pixel-localized
%edges
q2 = quiver(Xarrow,Yarrow,Uarrow,Varrow,0); %
%set(q1,'Color',[0 0 1]);
set(q2,'Color',newarrcol);
L = legend([p1,p2],leg1,leg2,-1);
set(L,'Fontsize',8);
twords = ['Edgels for Image ',imtitle];
tt  = title(twords);
set(tt,'Interpreter','none');
hold off;

return;
