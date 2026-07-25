classdef piv_functions
    methods (Static)
        
        function [X,Y,Z,z,Vel_norm,Vort_X,Vort_Y,Vort_Z,Vort_norm,Qcrit,front,top,side,ImageX,ImageY,ImageZ] = HydroParam(x,X,y,z,U,V,W,sigma)
        zd = max(z)+(min(z)-max(z))/2; % recentrage de z
        z = z-zd;
        [X,Y,Z] = meshgrid(x,y,z);
        nbfiles = size(U,4);
        % Velocity gradient tensor
        dx = 2*(x(2)-x(1))*0.001; 
        dy = 2*(y(2)-y(1))*0.001;
        dz = 2*(z(2)-z(1))*0.001;
        [dudx,dvdx,dwdx,dudy,dvdy,dwdy,dudz,dvdz,dwdz] = deal(zeros(size(X,1),size(X,2),size(X,3),nbfiles));
        
        for m = 1:nbfiles
            U(:,:,:,m) = imgaussfilt3(U(:,:,:,m),sigma);
            V(:,:,:,m) = imgaussfilt3(V(:,:,:,m),sigma);
            W(:,:,:,m) = imgaussfilt3(W(:,:,:,m),sigma);
        end
        
        for m = 1:nbfiles
            for k = 1:length(z)
                for j = 2:length(y)-1
                    for i = 1:length(x)
                        dudx(i,j,k,m) = (U(i,j+1,k,m)-U(i,j-1,k,m))/dx;
                        dvdx(i,j,k,m) = (V(i,j+1,k,m)-V(i,j-1,k,m))/dx;
                        dwdx(i,j,k,m) = (W(i,j+1,k,m)-W(i,j-1,k,m))/dx;
                    end
                end
            end 
        end %dx
        
        for m = 1:nbfiles
            for k = 1:length(z)
                for j = 1:length(y)
                    for i = 2:length(x)-1
                        dudy(i,j,k,m) = (U(i+1,j,k,m)-U(i-1,j,k,m))/dy;
                        dvdy(i,j,k,m) = (V(i+1,j,k,m)-V(i-1,j,k,m))/dy;
                        dwdy(i,j,k,m) = (W(i+1,j,k,m)-W(i-1,j,k,m))/dy;
                    end
                end
            end
        end %dy
        
        for m = 1:nbfiles
            for k = 2:length(z)-1
                for j = 1:length(y)
                    for i = 1:length(x)
                        dudz(i,j,k,m) = (U(i,j,k+1,m)-U(i,j,k-1,m))/dz;
                        dvdz(i,j,k,m) = (V(i,j,k+1,m)-V(i,j,k-1,m))/dz;
                        dwdz(i,j,k,m) = (W(i,j,k+1,m)-W(i,j,k-1,m))/dz;
                    end
                end
            end
        end %dz
        
        % Norm, Vort, Qcrit
        [Vel_norm,Vort_X,Vort_Y,Vort_Z,Vort_norm,Qcrit,Omega] = deal(zeros(size(U)));
        
            for m = 1:nbfiles
                Vel_norm(:,:,:,m) = sqrt(U(:,:,:,m).^2 + V(:,:,:,m).^2 + W(:,:,:,m).^2);
                Vort_X(:,:,:,m) = dwdy(:,:,:,m)-dvdz(:,:,:,m);
                Vort_Y(:,:,:,m) = dudz(:,:,:,m)-dwdx(:,:,:,m);
                Vort_Z(:,:,:,m) = dvdx(:,:,:,m)-dudy(:,:,:,m);
                Vort_norm(:,:,:,m) = sqrt(Vort_X(:,:,:,m).^2 + Vort_Y(:,:,:,m).^2 + Vort_Z(:,:,:,m).^2);
                Om = 0.5*(Vort_X(:,:,:,m).^2+Vort_Y(:,:,:,m).^2+Vort_Z(:,:,:,m).^2);
                S = (dudx(:,:,:,m).^2+dvdy(:,:,:,m).^2+dwdz(:,:,:,m).^2)+0.5*((dudy(:,:,:,m)+dvdx(:,:,:,m)).^2+(dudz(:,:,:,m)+dwdx(:,:,:,m)).^2+(dvdz(:,:,:,m)+dwdy(:,:,:,m)).^2);
                Qcrit(:,:,:,m) = 0.5*(Om-S);
                Omega(:,:,:,m) = Om./(Om+S+20);
                
            end
        
        % Visualisation parameters
        front = [0 90]; %XY
        top = [0 0]; % XZ
        side = [90 0]; %YZ
        
        % Position de l'image du serpent au fond du volume pour les futurs plots
        ImageX = [-90 90; -90 90];
        ImageY = [90 90; -90 -90];
        ImageZ = [-52.5 -52.5; -52.5 -52.5];
        
        end
        
        
        function [RR,GG,BB] = isosurfacecolor(funct,m)
        
        % maxval = quantile(max(max(funct(:,:,:,m))),0.75);
        maxval = 30;
        ncolors = 256;
        cmap = cmocean('-balance',ncolors);
        levels = linspace(-maxval,maxval,ncolors);
        [RR,GG,BB] = deal(zeros(size(funct(:,:,:,m))));
        
        RR(:,:,:,m) = maxval<funct(:,:,:,m)*cmap(end,1);
        RR(:,:,:,m) = RR(:,:,:,m) + (funct(:,:,:,m)<-maxval)*cmap(1,1);
        GG(:,:,:,m) = maxval<funct(:,:,:,m)*cmap(end,2);
        GG(:,:,:,m) = GG(:,:,:,m) + (funct(:,:,:,m)<-maxval)*cmap(1,2);
        BB(:,:,:,m) = maxval<funct(:,:,:,m)*cmap(end,3);
        BB(:,:,:,m) = BB(:,:,:,m) + (funct(:,:,:,m)<-maxval)*cmap(1,3);
        for jj = 2:ncolors-2
            RR(:,:,:,m) = RR(:,:,:,m) + ((levels(jj)<=funct(:,:,:,m)) & (funct(:,:,:,m)<levels(jj+1)))*cmap(jj,1);
            GG(:,:,:,m) = GG(:,:,:,m) + ((levels(jj)<=funct(:,:,:,m)) & (funct(:,:,:,m)<levels(jj+1)))*cmap(jj,2);
            BB(:,:,:,m) = BB(:,:,:,m) + ((levels(jj)<=funct(:,:,:,m)) & (funct(:,:,:,m)<levels(jj+1)))*cmap(jj,3);
        end
        
        end

    end
end
