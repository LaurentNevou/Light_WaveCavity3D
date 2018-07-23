%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 23Mai2018, lne %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

c=2.99792458e8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Solving Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0_guess= c/5e-6;     %% Guess of the frequency solutions (Hz)
f0_min  = c/20e-6;    %% filter the solutions where the frequency is superior than (Hz)
f0_max  = c/2e-6;   %% filter the solutions where the frequency is inferior than (Hz)
nmodes=10;            %% number of solutions asked 

AbsorbingBoundaryCondition=0;     %% 0 or 1 (not sure it is working well...)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Optical index definition %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Three vectors (x, y and z) and one 3D-matrix n must be defined with homogeneous grid
% x, y and z [meter]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nx=13;                  %% Meshing point in x-direction
Ny=17;                  %% Meshing point in y-direction
Nz=15;                  %% Meshing point in z-direction

Dx=1E-6;                %% map X [m]
Dy=1E-6;                %% map Y [m]
Dz=1E-6;                %% map Z [m]

x = linspace(-Dx, Dx, Nx);
y = linspace(-Dy, Dy, Ny);
z = linspace(-Dz, Dz, Nz);

dx = x(2)-x(1);
dy = y(2)-y(1);
dz = z(2)-z(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose between the next 2 structures or build your own on an homogeneous grid!

n1=1; n2=3;

Lx=1.3e-6; Ly=1.2e-6; Lz=1.4e-6;
[n,eps]=epsBox_f(x,y,z,Lx,Ly,Lz,n1,n2,AbsorbingBoundaryCondition);

%Rx=0.5e-6; Ry=0.6e-6; Rz=0.8e-6;
%[n,eps]=epsSpher_f(x,y,z,Rx,Ry,Rz,n1,n2,AbsorbingBoundaryCondition);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Selection of the model %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('=======================================')

tic
if length(x)*length(y)*length(z)*3>1e4
  Nxyz=length(x)*length(y)*length(z);
  display(strcat('Warning: Take care, H=',num2str(Nxyz),'x',num2str(Nxyz),'x3 elements'))
end
[Ex,Ey,Ez,f1]=WC3D_FEM_f(x,y,z,eps,nmodes,f0_guess,f0_min,f0_max);
display(strcat('-> Finite Elements method =',num2str(toc),'sec'))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f(1:length(f1),1)=f1;
lambda(1:length(f1),1)=c./f1;


%display('=======================================')
%display('Results f(THz):')
%f*1e-12
display('=======================================')
display('Results lambda(um):')
lambda*1e6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X0fig = 50    ; Y0fig = 50;
Wfig  = 1500  ; Hfig  = 800;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(f(:,1))
    
    figure('Name',strcat('FEM method: Electrical Field-',num2str(i)),'position',[X0fig Y0fig Wfig Hfig])
    colormap(jet)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,1)
    hold on;grid on;view (-38, 20);
    
    slice(x*1e6,y*1e6,z*1e6,n,[0],[0],[0])
    contour(x*1e6,y*1e6,z*1e6,abs(n),1,'linewidth',2,'linecolor','r')
    
    PSI=abs(Ex(:,:,:,i)).^2;
    p = patch(isosurface(x*1e6,y*1e6,z*1e6,PSI,max(PSI(:))/3));
    isonormals(x*1e6,y*1e6,z*1e6,PSI, p)
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceLighting', 'gouraud')
    daspect([1,1,1])
    light ('Position', [1 1 5]);
    title(  strcat('lambda=',  num2str( lambda(i,1)*1e6,'%.2f'),'um,  Abs(Ex)'  )  )
    xlabel('x (um)')
    ylabel('y (um)')
    zlabel('z (um)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,2)
    hold on;grid on;view (-38, 20);
    
    slice(x*1e6,y*1e6,z*1e6,n,[0],[0],[0])
    
    PSI=abs(Ey(:,:,:,i)).^2;
    p = patch(isosurface(x*1e6,y*1e6,z*1e6,PSI,max(PSI(:))/3));
    isonormals(x*1e6,y*1e6,z*1e6,PSI, p)
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceLighting', 'gouraud')
    daspect([1,1,1])
    light ('Position', [1 1 5]);
    title(  strcat('lambda=',  num2str( lambda(i,1)*1e6,'%.2f'),'um,  Abs(Ey)'  )  )
    xlabel('x (um)')
    ylabel('y (um)')
    zlabel('z (um)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,3)
    hold on;grid on;view (-38, 20);
    
    slice(x*1e6,y*1e6,z*1e6,n,[0],[0],[0])
    
    PSI=abs(Ez(:,:,:,i)).^2;
    p = patch(isosurface(x*1e6,y*1e6,z*1e6,PSI,max(PSI(:))/3));
    isonormals(x*1e6,y*1e6,z*1e6,PSI, p)
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceLighting', 'gouraud')
    daspect([1,1,1])
    light ('Position', [1 1 5]);
    title(  strcat('lambda=',  num2str( lambda(i,1)*1e6,'%.2f'),'um,  Abs(Ez)'  )  )
    xlabel('x (um)')
    ylabel('y (um)')
    zlabel('z (um)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,4)
    hold on;grid on;view (-38, 20);
    
    PSI=real(Ex(:,:,:,i));
    slice(x*1e6,y*1e6,z*1e6,PSI,[0],[0],[0])
    
    title(  strcat('lambda=',  num2str( lambda(i,1)*1e6,'%.2f'),'um,  Re(Ex)'  )  )
    xlabel('x (um)')
    ylabel('y (um)')
    zlabel('z (um)')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(2,3,5)
    hold on;grid on;view (-38, 20);
    
    PSI=real(Ey(:,:,:,i));
    slice(x*1e6,y*1e6,z*1e6,PSI,[0],[0],[0])
    
    title(  strcat('lambda=',  num2str( lambda(i,1)*1e6,'%.2f'),'um,  Re(Ey)'  )  )
    xlabel('x (um)')
    ylabel('y (um)')
    zlabel('z (um)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,3,6)
    hold on;grid on;view (-38, 20);
    
    PSI=real(Ez(:,:,:,i));
    slice(x*1e6,y*1e6,z*1e6,PSI,[0],[0],[0])
    
    title(  strcat('lambda=',  num2str( lambda(i,1)*1e6,'%.2f'),'um,  Re(Ez)'  )  )
    xlabel('x (um)')
    ylabel('y (um)')
    zlabel('z (um)')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%