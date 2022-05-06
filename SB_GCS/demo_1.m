

clear 
close all;

Img=double(imread('cameraman.tif'));

%%
sigma=.8;    % scale parameter in Gaussian kernel
gfilt=fspecial('gaussian',15,sigma); % Gaussian kernel
Img_smooth=conv2(Img,gfilt,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
G=Ix.^2+Iy.^2;

% edge detector fucntion=
g=1./(1+abs(G));
u=zeros(size(Img));
u(9:end-8,9:end-8)=1;


% play with these
lambda=0.001; mu =0.0001;

u=sp_gcs(u,Img,g,lambda,mu,0.01,40); 

figure(1)
title('Final contour')
hold on
imshow(Img,[],'InitialMagnification','fit');
contour(u,[0.5 0.5],'linecolor','r');
drawnow


% Note that because the level set fucntion is regularized to L1 norm, the
% contour can be drawn at almost all levels bewteen 1 and 0 and still get
% prectically the same result.

figure(2)
title('L1 regularized level set fucntion')
surf(u,'linestyle','none')
axis tight





