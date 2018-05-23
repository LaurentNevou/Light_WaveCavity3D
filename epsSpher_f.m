function[n,eps]=epsSpher_f(x,y,z,Rx,Ry,Rz,n1,n2,AbsorbingBoundaryCondition)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Optical index definition %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[X,Y,Z]=meshgrid(x,y,z);

x0=0;y0=0;z0=0;


idXYZ =  ((X-x0)/Rx).^2 + ((Y-y0)/Ry).^2 + ((Z-z0)/Rz).^2 < 1  ;  %% elipse


n = n2*idXYZ + n1*(1-idXYZ);  %% ridge optical index

%eps(Y<=0)=4+0i;                %% substrate optical index

if AbsorbingBoundaryCondition==1
    LOSS=1e-4;
    n(:,:,1)                = n(:,:,1)               + LOSS*i;
    n(:,:,end)              = n(:,:,end)             + LOSS*i;
    n(1,:,2:end-1)          = n(1,:,2:end-1)         + LOSS*i;
    n(end,:,2:end-1)        = n(end,:,2:end-1)       + LOSS*i;
    n(2:end-1,1,2:end-1)    = n(2:end-1,1,2:end-1)   + LOSS*i;
    n(2:end-1,end,2:end-1)  = n(2:end-1,end,2:end-1) + LOSS*i;
end

eps=n.^2;

break
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Name','Potential','position',[10 -50 1600 800])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,1,'fontsize',15)
hold on;grid on;

pcolor(x*1e6,z*1e6,squeeze(n(round(end/2),:,:))')

colormap(cool)
colorbar

%xlim([-1 1]*Dx/2*1e6)
%ylim([-1 1]*Dz/2*1e6)

xlabel('x (um)')
ylabel('z (um)')
title(strcat('n-xz @y=0um'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,4,'fontsize',15)
hold on;grid on;

pcolor(y*1e6,z*1e6,squeeze(n(:,round(end/2),:))') 

colormap(cool)
colorbar

%xlim([-1 1]*Dx/2*1e6)
%ylim([-1 1]*Dz/2*1e6)

xlabel('y (um)')
ylabel('z (um)')
title(strcat('n-yz @x=0um'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,2,'fontsize',15)
hold on;grid on;

idz=find(z>0);
idz=idz(1);

pcolor(x*1e6,y*1e6,squeeze(n(:,:,idz)))

colormap(cool)
colorbar

%xlim([-1 1]*Dy/2*1e6)
%ylim([-1 1]*Dy/2*1e6)

xlabel('x (um)')
ylabel('y (um)')
title(strcat('n-xy @z=',num2str(z(idz)*1e6,'%.2f'),'um'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,5,'fontsize',15)
hold on;grid on;view (-38, 20);

slice(x*1e6,y*1e6,z*1e6,n,[0],[0],[0])
xlabel('x (nm)')
ylabel('y (nm)')
zlabel('z (nm)')
title(strcat('n'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(2,3,3,'fontsize',15)
hold on;grid on;
plot(z*1e6,squeeze(n(round(end/2),round(end/2),:)) ,'b.-')
plot(z*1e6,squeeze(n(round(end/2),end,:)) ,'r.-')

xlabel('z (nm)')
ylabel('Optical index')
title(strcat('n-z @y=0um and \color{blue}x=0um ; \color{red}x=',num2str(x(end)*1e6),'um'))

