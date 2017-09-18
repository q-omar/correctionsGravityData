%GOPH 547
%Ye Sun
%Safian Omar Qureshi
%ID: 10086638
%Collaberated with Tian Yu, Younghwan Ryan Ok, Ahtisham Mohummad 

clear;
clc;
load('goph_547_lab_2_data_w2017.mat')

%% Data Preperation 
% Q1 - Q7

gz_sd=grav_survey_data(:,8); %observed gravity )in microgals) 
gz_sd(1)=0; %equating g(1) to zero

total_time=grav_survey_data(:,4); %total time of g in survey

dg_tide = grav_survey_data(:,9); %tidal variations 

p=polyfit(total_time,gz_sd,1); %finding m and b values for linear equation relationship
yval=polyval(p,total_time); %making y values 


figure %plotting gravities, tide variations and linear relation
plot(total_time,gz_sd);
title('Gravity Effect & Tidal Variations Vs. Total Time','fontweight','bold')
xlabel('Total Time (hrs)','fontweight','bold')
ylabel('Gravity Effect/Tidal Variations (microGal)','fontweight','bold')
hold on;
plot(total_time,dg_tide);
legend('Gravity Effect','Tidal Variations')
hold on;
plot(total_time,yval)
%%
%Q2

max_z=max(max(Zt)); %finding min/max of Z
min_z=min(min(Zt));

figure %plot elevation change of terrain in survey region
contourf(Xt,Yt,Zt);
h_c=colorbar; %add colorbar
cmap = flipud(colormap('parula'));
colormap(cmap);
title('Contour plot for terrain survey data','fontweight','bold');
xlabel('Easting (km)','fontweight','bold'); %adding label axis title etc
ylabel('Northing (km)','fontweight','bold');
ylabel(h_c,'Elevation (m)','fontweight','bold');

z_datum=mean(Zt(:)); %datum elevation is mean elevation for bougeur, freeair
min_z=min(min(Zt));
max_z=max(max(Zt)); %finding min/max

%Q3

x=grav_survey_data(:,5); % Extracting x,y, z data
y=grav_survey_data(:,6);
z=grav_survey_data(:,7);

x_sort=[x,y,z]; % creating 3 column matrix
[x_sort,ind_sort]=sortrows(x_sort,[1,2]); % sorting the matrix

%Q4

[x_sort,ind_uniq]=unique(x_sort,'rows','legacy'); %extrating unique entries from sorted matrix earlier

%Q5

Xg=x_sort(:,1); %storing columns x_sort as three separate column vectors Xg, Yg Zg
Yg=x_sort(:,2); 
Zg=x_sort(:,3);

Ny=sqrt(size(Xg,1)); % verifying the number of grid points in lists is square number
Nx=Ny;

%Q6
Xg = reshape(Xg,Nx,Ny); %generate X, Y, Z grids similar to grids created by meshgrid
Yg = reshape(Yg,Nx,Ny);
Zg = reshape(Zg,Nx,Ny);

Xg2 = (0:50); Yg2 = (0:50); %% verify above matixes are identical to produced below
[Xg2,Yg2] = meshgrid(Xg2,Yg2); % in the workspace, if you compare Xg to Xg_2 & Yg to Yg_2, one can see they are identical. 

%%
%Q7

g_corr=gz_sd;
g_raw = g_corr;
g_raw(1)=g_corr(1);
l=length(g_raw);

for i=2:l
    g_raw(i)= g_raw(1)+g_raw(i);
end

g_raw=g_raw(ind_sort);% extracting the sorted data
g_raw=g_raw(ind_uniq);% extracting the unique data 
g_raw=reshape(g_raw,[Nx Ny]);% turning the data into a grid 

figure
contourf(Xg,Yg,g_raw);
h_c=colorbar; %add colarbar
cmap = flipud(colormap('parula'));
colormap(cmap);
title('Raw Gravity data');
xlabel('Easting (km)', 'fontweight','bold');
ylabel('Northing (km)', 'fontweight','bold');
ylabel(h_c,'Gravity (microGal)','fontweight','bold');

%% Normal Gravity Corrections
% Q8
lat=49.1286;
ge_IGF=9.7803253359;
e_const=0.00669437999013;
k_const=0.00193185265241;
rad=lat*(pi/180);
gt=ge_IGF*((1+k_const*(sin(rad))^2)/sqrt(1-e_const*(sin(rad))^2));
gt=gt*10^8;% convert from m/s2 to microGal
g_corr(1)=g_corr(1)-gt;


% Q9
g_norm=g_corr;
l=length(g_norm);

for i=2:l
    g_norm(i)= g_norm(1)+g_norm(i);
end

g_norm=g_norm(ind_sort);% extracting the sorted data
g_norm=g_norm(ind_uniq);% extracting the unique data 
g_norm=reshape(g_norm,[Nx Ny]);% turning the data into a grid 


figure
contourf(Xg,Yg,g_norm);
h_c=colorbar; %add colarbar
cmap = flipud(colormap('parula'));
colormap(cmap);
title('Theoritical Gravity Corrected Data');
ylabel('Northing (km)', 'fontweight','bold');
ylabel(h_c,'Gravity(microGal)','fontweight','bold');
xlabel('Easting (km)', 'fontweight','bold');

%% Drift Corrections
% Q10 - Q11: doing tidal drift correction and applying it as g_corr 

g_corr = g_corr - dg_tide; g_tide = g_corr;

l=length(g_tide);

for i=2:l
    g_tide(i)= g_tide(1)+g_tide(i);
end

g_tide=g_tide(ind_sort);% extracting the sorted data
g_tide=g_tide(ind_uniq);% extracting the unique data 
g_tide=reshape(g_tide,[Nx Ny]);% turning the data into a grid 

figure
contourf(Xg,Yg,g_tide);
h_c=colorbar; %add colarbar
cmap = flipud(colormap('parula'));
colormap(cmap);
title('Tidal Drift Corrected Data');
xlabel('Easting (km)', 'fontweight','bold');
ylabel('Northing (km)', 'fontweight','bold');
ylabel(h_c,'Gravity(microGal)','fontweight','bold');

%%
%Q12

i1=1; % beginning of 1st loop

for i2=1:size(gz_sd) % end of current loop
    
    if grav_survey_data(i2,1)==2
        drift=g_corr(i2);% part a:
        dt=total_time(i2)-total_time(i1);% time difference
        drift_rate=drift/dt;%part b: drift rate
        g_corr(i1+1:i2-1)=g_corr(i1+1:i2-1)-drift_rate*(total_time(i1+1:i2-1)-total_time(i1));%part c: correcting for instrument drift
        g_corr(i2:end)=g_corr(i2:end)-drift;
    end
    
    if i2>1 && grav_survey_data(i2,1)==1
        i1=i2; %save index for beginning of next loop
        drift=g_corr(i1);
        g_corr(i1:end)=g_corr(i1:end)-drift;
    end
    
end

% Q13 % If one opens the workspace for g_corr they can see that g_corr==0
% for all points where base_point==1 or 2.
l=length(g_corr);
for i=2:l
    g_corr(i)= g_corr(1)+g_corr(i);
end
   
figure
plot(total_time,g_corr);
xlabel('Total Time (hrs)','fontweight','bold')
ylabel('Corrected Gravity Effect (microGal)','fontweight','bold')
title('Corrected Gravity Effect Vs. Total Time','fontweight','bold')
    
%%        
% Q14:
g_corr=g_corr(ind_sort);% extracting the sorted data
g_corr=g_corr(ind_uniq);% extracting the unique data
g_corr=reshape(g_corr, [Nx Ny]);

g_drift=g_corr;

figure
contourf(Xg,Yg,g_drift);
h_c = colorbar;
cmap = flipud(colormap('parula'));
colormap(cmap);
title('All time dependent Corrections Applied');
ylabel('Northing (km)', 'fontweight','bold');
ylabel(h_c,'Gravity(microGal)','fontweight','bold');
xlabel('Easting (km)', 'fontweight','bold');

%% Elevation & Terrain Corrections
% Q15: Free air correction

for i=1:numel(g_corr)
    dz(i)=Zg(i)-z_datum;
    dg_FA(i)=-0.3086*1000*dz(i); %x1000 since our data is in km and formula for meter
    g_corr(i)=g_corr(i)+dg_FA(i);
end

% Q16: plotting free air correct data
g_FA=g_corr;

figure
contourf(Xg,Yg,g_FA);
h_c = colorbar;
cmap = flipud(colormap('parula'));
colormap(cmap);
title('Free Air Corrected Data');
ylabel('Northing (km)', 'fontweight','bold');
ylabel(h_c,'Gravity(microGal)','fontweight','bold');
xlabel('Easting (km)', 'fontweight','bold');

% Q17: Bouguer plate correction
for i=1:numel(g_corr)
    dz(i)=Zg(i)-z_datum;
    rho=2.65;
    dg_BP(i)=0.04193*rho*1000*dz(i); %x1000 since our data is in km and formula for meter
    g_corr(i)=g_corr(i)+dg_BP(i);
end

% Q18: Plotting Bouguer plate correction
g_elev=g_corr;

figure
contourf(Xg,Yg,g_elev);
h_c=colorbar; %add colarbar
cmap = flipud(colormap('parula'));
colormap(cmap);
title('Elevation Corrections Applied');
ylabel('Northing (km)', 'fontweight','bold');
ylabel(h_c,'Gravity(microGal)','fontweight','bold');
xlabel('Easting (km)', 'fontweight','bold');

%% 
%Q19

dg_terr=zeros(size(g_corr),'like',g_corr);% initilizing terrain correction array to zeros and same dimensions as g_corr

dx=1000;
dy=dx;
dA=dx^2;
n=size(dg_terr,1);
n2=size(Zt,1);

for i=1:n
    for j=1:n
        xi=[Xg(i,j)*1000,Yg(i,j)*1000,z_datum]; %part a
        for i2=1:n2
            for j2=1:n2
                xm=[Xt(i2,j2)*1000,Yt(i2,j2)*1000,(Zt(i2,j2)-z_datum)*0.5];% part b:i
                dm=(Zt(i2,j2)-z_datum)*rho*1000*dA;% part b:ii
                grav_eff=grav_eff_point(xi,xm,dm);%part b:iii
                dg_terr(i,j)=dg_terr(i,j)+abs(grav_eff); % part b:iv partial
            end
        end
    end
end

% part b:iv partial increment current value by absolute gravity effect
inc_dg_terr=abs(dg_terr)*10^8;
dg_terr=dg_terr+inc_dg_terr;

% Q20:
g_corr=g_corr+dg_terr;
g_terr=g_corr;

figure
contourf(Xg,Yg,g_terr);
h_c=colorbar; %add colarbar
cmap = flipud(colormap('parula'));
colormap(cmap);
title('Terrain Corrections Applied');
ylabel('Northing (km)', 'fontweight','bold');
ylabel(h_c,'Gravity(microGal)','fontweight','bold');
xlabel('Easting (km)', 'fontweight','bold');

%%
%Q21: eliminate regional variations
dx=1000; %recreated these variables here because it takes a long time to run terrain corrections, needed for derivatives 
dy=dx;
N=size(g_corr,1);
M=N;
tot_g_corr=sum(sum(g_corr));
dg_rgn1=tot_g_corr/(N*M);

g_corr=g_corr-dg_rgn1;
g_anom=g_corr;

figure
contourf(Xg,Yg,g_anom);
h_c=colorbar; %add colarbar
cmap = flipud(colormap('parula'));
colormap(cmap);
title('Regional Corrections Applied');
ylabel('Northing (km)', 'fontweight','bold');
ylabel(h_c,'Gravity(microGal)','fontweight','bold');
xlabel('Easting (km)', 'fontweight','bold');

%% 1st and 2nd derivatives
%for this part I am not too confident of my answers, the figures produced
%dont seem to be correct. Worked with Ryan Younghwan Ok and Tian Yu 

% finding the lateral extent of the anomalies

dg_dx = zeros(Nx); dg_dy = dg_dx; dg_dx_2 = dg_dy; dg_dy_2 = dg_dx_2;

g_z = g_corr;

for i = 2:Nx-1
    for j = 2:Nx-1
        dg_dx(i,j) = (g_z(i,j)-g_z(i,j-1))/dx;
        dg_dy(i,j) = (g_z(i,j)-g_z(i-1,j))/dy;
    end
end

% find second derivative in x and y directions to use Laplace's equation to
% find second vertical derivative

for i = 2:Nx-1
    for j = 2:Nx-1
        dg_dy_2(i,j) = (g_z(i+1,j)-2*g_z(i,j)+g_z(i-1,j))/dy^2;
        dg_dx_2(i,j) = (g_z(i,j+1)-2*g_z(i,j)+g_z(i,j-1))/dx^2;
    end
end

% Enhance near surface effects using the second vertical derivative
dg_dz_2 = -(dg_dx_2 + dg_dy_2);

mindg_dz = min(min(dg_dz_2));
maxdg_dz = max(max(dg_dz_2));


mindg_dx = min(min(dg_dx));
maxdg_dx = max(max(dg_dx));

mindg_dy = min(min(dg_dy));
maxdg_dy = max(max(dg_dy));

figure(12);
contourf(Xg,Yg,dg_dx);
h = colorbar;
cmap = flipud(colormap('parula'));
colormap(cmap);
title('1^{st} derivative of gravity in x direction');
xlabel('Easting (km)'); ylabel('Northing (km)'); ylabel(h,'\muGal/m');
caxis([mindg_dx,maxdg_dx]);
axis equal; 
figure(13);
contourf(Xg,Yg,dg_dy);
h = colorbar;
cmap = flipud(colormap('parula'));
colormap(cmap);
title('1^{st} derivative of gravity in y direction');
xlabel('Easting (km)'); ylabel('Northing (km)'); ylabel(h,'\muGal/m');
caxis([mindg_dy,maxdg_dy]);
axis equal;

figure(14);
contourf(Xg,Yg,dg_dz_2);
h = colorbar;
cmap = flipud(colormap('parula'));
colormap(cmap);
title('2^{nd} derivative of gravity in vertical direction');
xlabel('Easting (km)'); ylabel('Northing (km)'); ylabel(h,'\muGal/m^2');
caxis([mindg_dz,maxdg_dz]);
axis equal;



%% Excess mass calculation
g_z = g_corr * 10^-8;
G=6.674*10^(-11)
avg_g = 0;

for j = 6:26
    for i = 26:41
        avg_g = avg_g + (g_z(i,j)+g_z(i+1,j)+g_z(i+1,j+1)+g_z(i,j+1))/4;
    end
end

Mass = dA * avg_g / (2*pi*G); disp(['Mass of  anomaly ' num2str(Mass) ,' kg']);