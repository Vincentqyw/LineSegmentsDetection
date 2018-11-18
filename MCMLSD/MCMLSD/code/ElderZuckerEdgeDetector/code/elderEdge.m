function [edgeStruct] = elderEdge()
%Run local scale control edge detector on image (Elder & Zucker, 1998)

noise = 1; %estimated SD of pixel noise
imagefile='P1020177.jpg'; %name of image file.  Most formats are handled.  May be RGB or grayscale.
maxscale=5; %maximum scale, must be positive integer. Scales increase exponentially from 0.5 to 2^{maxscale-2}.  
edgew=10; %maximum edge blur (distance in pixels between extrema in the 2nd derivative)
conss=1; %boundary condition for scalespace computation.  0 for extension, 1 for no overlap
congrd=0; %boundary condition for derivative computation.  0 for extension, 1 for no overlap
subpixelflag=0; %generate subpixel-localized edges
%run edge detection
edgeStruct = main_edge(imagefile,maxscale,noise,edgew,conss,congrd,subpixelflag);