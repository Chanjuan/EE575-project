function u = sp_gcs(u,f, g, lambda,mu, eps, maxiter)
%
%  Input:
%      u: level set function to be updated by level set evolution
%      f: image to be segmented
%      g: edge indicator function
%      lambda: Split Bragman  variable regularization stregth
%      mu: regularization stregth, contour level
%      eps: minimum step size
%      maxiter: maximum number of iterations
%  Output:
%      u: updated level set function after level set evolution


c1=mean(f(u>mu));
c2=mean(f(u<mu));

uprev=zeros(size(u));
b=zeros([size(u),2]);
d=zeros([size(u),2]);
ii=1;

while (norm2D(u-uprev,2)>eps) && (ii<maxiter)
    

    
%   from algorithm 3, steps: 
    % 2:
    r=(c1-f).^2-(c2-f).^2;
    % 3:
    uprev=u;
    u=Gauss_Seidel_GCS(u,r,d,b,mu,lambda);
    % 4:
    [deltau(:,:,1), deltau(:,:,2)]=gradient(u);
    d=shrink_g(deltau+b,lambda,g);
    % 5:
    b=b+deltau-d;
    % 6 & 7:
    c1=mean(f(uprev>mu));
    c2=mean(f(uprev<mu));
    
    
    ii=ii+1;
    
end
disp(['Exit in iteration number ', num2str(ii)])
end
function u=Gauss_Seidel_GCS(u,r,d,b,mu,lambda)

x=1;y=2; %for clarity

% Boundary conditions ignored

alpha=  -d(:,:,x)+b(:,:,x)-d(:,:,y)+b(:,:,y);
alpha(2:end,:)=alpha(2:end,:)+d(1:end-1,:,x)-b(1:end-1,:,x);
alpha(:,2:end)=alpha(:,2:end)+ d(:,1:end-1,y)-b(:,1:end-1,y);

beta=(-mu/lambda*r+alpha);
beta(2:end,:)=beta(2:end,:)+u(1:end-1,:);
beta(1:end-1,:)=beta(1:end-1,:)+u(2:end,:);
beta(:,2:end)=beta(:,2:end)+u(:,1:end-1);
beta(:,1:end-1)=beta(:,1:end-1)+u(:,2:end);

beta=beta*0.25;


u=max(min(beta,1),0);
end

function out=norm2D(x,n)
out=norm(x(:),n);
end

function out=shrink_g(z,lambda,g)
nz=sqrt(  sum(z.^2,3)   );

out=max(nz-lambda./g,0).*(z./nz);
out(isnan(out))=0;
end