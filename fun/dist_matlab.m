function [dist] = dist_matlab(p)
N=size(p,2);
dist=zeros(N,N);
for n=1:N
    x=p(:,n);
     m=n:N;
        y=p(:,m);
        dist(n,m)=sqrt((x(1)-y(1,:)).^2+...
                       (x(2)-y(2,:)).^2+...
                       (x(3)-y(3,:)).^2);
        dist(m,n)=dist(n,m);  
end 
end

