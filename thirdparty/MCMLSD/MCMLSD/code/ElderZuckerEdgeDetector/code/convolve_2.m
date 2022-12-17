
%##############################################################
%
% [cimg] = convolve_2(mimg,filt,bc)
% 
% Convolution of two matrices with the boundaries handled 
%   depending on the value of the boundary condition variable.
%        - 0 for extension convolution
%        - 1 for no overlap convolution 
%
%##############################################################

function[cimg] = convolve_2(mimg,filt,bc)

if (bc == 0)
    
    pad_img = pad_matrix(mimg,size(filt));
    cimg    = conv2(pad_img,filt,'same');
    cimg    = trim_matrix(cimg,size(filt));
    
else
    
    cimg = conv2(mimg,filt,'same');

    if (size(filt,1) < size(filt,2))
        k                   = floor(0.5*length(filt));
        cimg(:,1:k)         = 0;
        cimg(:,end-k+1:end) = 0;
    else
        k                   = floor(0.5*length(filt));
        cimg(1:k,:)         = 0;
        cimg(end-k+1:end,:) = 0;
    end;
    
end;

return;
        

%########################
% Function to pad matrix:
%########################
function[r] = pad_matrix(m,d)

	if (d(1) < d(2))
        k   = floor(0.5*d(2));
        ad1 = repmat(m(:,1),1,k);
        ad2 = repmat(m(:,end),1,k);
        r   = [ad1 m ad2];
	else
        k   = floor(0.5*d(1));
        ad1 = repmat(m(1,:),k,1);
        ad2 = repmat(m(end,:),k,1);
        r   = [ad1; m; ad2];
	end;

return;


%#########################
% Function to trim matrix:
%#########################
function[r] = trim_matrix(m,d)

	if (d(1) < d(2))
        k   = floor(0.5*d(2));
        r   = m(:,k+1:end-k);
	else
        k   = floor(0.5*d(1));
        r   = m(k+1:end-k,:);
	end;

return;