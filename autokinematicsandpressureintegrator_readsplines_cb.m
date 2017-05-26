%May 21, 2014

%This script uses a spline file from the auto foil kinematics vi to
%generate the interface along which forces and torques due to pressure are
%calculated.  Pressure field data comes from the Queen2 pressure field
%calculator.


%Note that the cbrewer function must be loaded into MatLab prior to using
%this script



clear all                                                     

close all                                                       %close open figures

j=0;                                                         %initialize index for movie frame storage

first = 1;                                                  % first image and pressure file number

last = 1981;                                                 % last image and pressure file number

increment = 10;                                              % increment between file numbers



numformat = '%05d';                                         % format of numbers in file name. '%05d' is five digit format with leading zeros

delimiter = ',';                                             % delimiter between columns in pressure file

headerlines = 0;                                            % number of header lines in pressure file


spldata = dlmread('H:\Nonunif_performance_PIV\Pressure code analysis\dynamic\3_3_heave_1.5Hz_1.5cm_0.3ms\results_3cyc\splines_frm1-2000_incr10.txt');          %read in spline data file


pressfilestub = 'H:\Nonunif_performance_PIV\Pressure code analysis\dynamic\3_3_heave_1.5Hz_1.5cm_0.3ms\results_3cyc\queen2_rodsim_dT10ms_';                                  % string prefix for pressure files

mov=VideoReader('H:\Nonunif_performance_PIV\Data_videos\dynamic\3_3_heave_1.5Hz_1.5cm_0.3ms.avi');          % read in movie
c1=1;                                                       %initialize index for spline file x-data column count
c2=2;                                                       %initialize index for spline file y-data column count
fr=first;                                                   %set index for movie frame read
frstep=10;


mycolors=cbrewer('div','RdBu',128);                       %Get ColorBrewer palette
fci = 0;
for row = 0:1:127
    fci = fci+1;
    flippedcolors(fci,:) = mycolors(128-row,:);
end



for filenum = first:increment:last-increment                          % loop through files                                                 
clear A                                                     %clear last iteration's image to start fresh
clear X                                                     %clear last iteration's image
    clf                                                     % clear current figure

    filenum                                                 % display current file number
press = dlmread([pressfilestub num2str(filenum,numformat) '.dat'],delimiter,headerlines,0); % read in pressure file

    xdata=spldata(2:end,c1);       %Get x data for the spline
    xdata=xdata/1000.;          %Convert from mm to m to match the pressure field data
    xdata=min(press(:,1))+xdata;    %Scale to match the pressure field
    ydata=spldata(2:end,c2);       %Get y data for the spline
    ydata=ydata/1000.;          %Convert from mm to m to match the pressure field data
    ydata=max(press(:,2))-ydata-2/1000;         %Scale to match the pressure field and center spline on foil
    
im=flipud(read(mov,spldata(1,c1)));    %Get current image, flip upside-down to move origin from upper left to lower left

im = imfill(im);                                                          % even out image texture



% Remove image boundaries

%im(1:120,:) = 0;

%im(700:end,:) = 0;

%%%%%


figure(1)

plot(xdata,ydata,'LineWidth',2,'Color','white')                    %plot spline in white
set(gca,'Color',[0 0 0]);                       %make background black
axis([min(press(:,1)) max(press(:,1)) min(press(:,2)) max(press(:,2))]);
set(gca,'YDir','normal')                        %make sure y-axis runs in correct direction (positive up) 
X=getframe;                                     %get the figure
A=X.cdata;                                      %make it an image




                                                                %im2BW converts from greyscale to BW; the .## is the threshold to make black versus white
BWs = bwareaopen(im2bw(A,.7),50);                                      % remove nything that isn't the foil

                                                                          
                                                                          %want boundary outside the "not a number" region
se90 = strel('line', 15, 90);                                           % kernel for padding vertical boundaries of tracked object

se0 = strel('line', 15, 0);                                             % kernel for padding horizontal boundaries of tracked object



BWsfinal = imdilate(BWs, [se90 se0]);                                   % pad objects in binary image
BWsfinal=flipud(BWsfinal);                                              %flip binary image so that the origin is in lower left corner




                                                                %im2BW converts from greyscale to BW; the .## is the threshold to make black versus white
BWs = bwareaopen(im2bw(A,.7),200);                                      % remove nything that isn't the foil

                                                                          
                                                                          %want boundary outside the "not a number" region
se90 = strel('line', 15, 90);                                           % kernel for padding vertical boundaries of tracked object

se0 = strel('line', 15, 0);                                             % kernel for padding horizontal boundaries of tracked object



BWsfinal = imdilate(BWs, [se90 se0]);                                   % pad objects in binary image
BWsfinal=flipud(BWsfinal);                                              %flip binary image so that the origin is in lower left corner











%Turn these lines on only to check the conversion to a binary image - keep OFF when running the code to completion

% figure(2)
% imshow(BWsfinal)
% return








B = bwboundaries(BWsfinal);                                             % find boundaries of remaining objects in binary image



[max_Bsize, max_Bindex] = max(cellfun('size', B, 1));                   % find largest object in binary image (assumed to be object of interest)



boundraw = B{max_Bindex};                                               % extract coordinates of object

boundxraw = boundraw(1:length(boundraw)/200:end,2);                     % create 200-point boundary coordinates

boundyraw = boundraw(1:length(boundraw)/200:end,1);                     % create 200-point boundary coordinates



boundx = (boundxraw./size(A,2))*range(press(:,1))+min(press(:,1));      % convert boundary coordinates from pixels to physical units

boundy = (boundyraw./size(A,1))*range(press(:,2))+min(press(:,2));      % convert boundary coordinates from pixels to physical units



boundx = smooth(boundx);

boundy = smooth(boundy);



surfdx = boundx - circshift(boundx,1);                                  % compute x-spacing between boundary coordinates

surfdy = boundy - circshift(boundy,1);                                  % compute y-spacing between boundary coordinates

   

surfnormx = surfdy;                                                     % compute x-component of surface normal

surfnormy = -surfdx;                                                    % compute y-component of surface normal

   

surfunitnormx = surfnormx./(sqrt(surfnormx.^2 + surfnormy.^2));         % compute x-component of unit surface normal

surfunitnormy = surfnormy./(sqrt(surfnormx.^2 + surfnormy.^2));         % compute y-component of unit surface normal

   

surfposx = (boundx + circshift(boundx,1))./2;                           % compute x-midpoint of each object surface facet

surfposy = (boundy + circshift(boundy,1))./2;                           % compute y-midpoint of each object surface facet




%The following section is for pressure force/torque calculation (uses the queen2 pressure field to calculate forces/torques along the interface made above)


surfpress = griddata(press(:,1),press(:,2),press(:,7),surfposx,surfposy,'nearest');       % compute pressure at midpoint of each surface facet



time(((fr-first)/frstep)+1) = (fr)*0.001;



forcex(((fr-first)/frstep)+1) = nansum(sqrt(surfdx.^2+surfdy.^2).*-surfunitnormx.*surfpress);    % compute x-force due to pressure and store in temporal record

forcey(((fr-first)/frstep)+1) = nansum(sqrt(surfdx.^2+surfdy.^2).*-surfunitnormy.*surfpress);    % compute y-force due to pressure and store in temporal record



%continue

[geom, iner, cpmo] = polygeom(surfposx, surfposy);                      % compute geometric parameters of object

   

surfcentx(((fr-first)/frstep)+1) = geom(2);                     % compute x-centroid and store in temporal record

surfcenty(((fr-first)/frstep)+1) = geom(3);                     % compute y-centroid and store in temporal record

 

axisang1x = cos(cpmo(2));                                            % compute first principal axis

axisang1y = sin(cpmo(2));

  

axisang2x = cos(cpmo(4));                                            % compute second principal axis

axisang2y = sin(cpmo(4));

 

bodyangle(((fr-first)/frstep)+1) = cpmo(2)+pi;                  % compute body angle and store in temporal record



%%%%% compute radial distance from each surface point to axis [USED FOR AXISYMMETRIC BODIES]

        a = ([(surfcentx(((fr-first)/frstep)+1)+axisang2x)*ones(size(surfposx)),(surfcenty(((fr-first)/frstep)+1)+axisang2y)*ones(size(surfposx)),0*ones(size(surfposx))] - [surfcentx(((fr-first)/frstep)+1)*ones(size(surfposx)),surfcenty(((fr-first)/frstep)+1)*ones(size(surfposx)),0*ones(size(surfposx))]);

        b = [surfposx,surfposy,0*ones(size(surfposx))] - [surfcentx(((fr-first)/frstep)+1)*ones(size(surfposx)),surfcenty(((fr-first)/frstep)+1)*ones(size(surfposx)),0*ones(size(surfposx))];

        crossab = cross(a,b,2);

        R = abs(crossab(:,3)) ./ hypot(a(:,1),a(:,2));

%%%%% 

   





[surfmin minindex] = min(boundx);                       %Get the minimum and maximum x values in the boundary
[surfmax maxindex] = max(boundx);

%Add the leading edge edge to list of leading edge coordinates
surfheadx(((fr-first)/frstep)+1) = min(boundx);

surfheady(((fr-first)/frstep)+1) = boundy(minindex);




%Add the trailing edge to list of trailing edge coordinates
surftailx(((fr-first)/frstep)+1) = max(boundx);

surftaily(((fr-first)/frstep)+1) = boundy(maxindex);



%%%%% compute moment arm for torque calculation

armx = surfposx-surfheadx(((fr-first)/frstep)+1);   

army = surfposy-surfheady(((fr-first)/frstep)+1);

arm = cat(2, armx, army, zeros(size(armx)));                            % moment arm



surfunitnorm = cat(2, surfunitnormx, surfunitnormy, zeros(size(armx))); % unit normal

darea = sqrt(surfdx.^2 + surfdy.^2);                                    % area (per unit depth) of each surface facet

  

dtorque = -surfpress.*(cross(arm, surfunitnorm, 2)*[0;0;1]).*darea;     % calculate torque on each surface facet

  

torque(((fr-first)/frstep)+1) = nansum(dtorque);                   % calculate net torque on body and store in temporal record    

%continue


                                          




%Get the blanked region and make it black
figure(4)
contourf(reshape(press(:,1),128,128),reshape(press(:,2),128,128),reshape(press(:,7),128,128),24,'LineColor','none')     %Display the pressure field
axis([min(press(:,1)) max(press(:,1)) min(press(:,2)) max(press(:,2))])
caxis([-80 1000])
colormap('gray')                                %Make field greyscale
X=getframe;                                     %get the figure
D=X.cdata;                                      %make it an image

BWs = im2bw(D,0.99);             %Get the white region (blanked foil)
BWs = flipud(BWs);                              %Move the origin to the right place

se90 = strel('line', 3, 90);                                           % kernel for padding vertical boundaries of tracked object
se0 = strel('line', 3, 0);                                             % kernel for padding horizontal boundaries of tracked object
BWs = imdilate(BWs, [se90 se0]);                                        %Dilate the foil based on padding kernels above

B = bwboundaries(BWs);                                             % find boundaries of remaining objects in binary image
[max_Bsize, max_Bindex] = max(cellfun('size', B, 1));                   % find largest object in binary image (the foil)


braw = B{max_Bindex};                                               % extract coordinates of object
bxraw = braw(1:length(braw)/200:end,2);                     % create 200-point boundary coordinates
byraw = braw(1:length(braw)/200:end,1);                     % create 200-point boundary coordinates

bx = (bxraw./size(BWs,2))*range(press(:,1))+min(press(:,1))+0.001;      % convert boundary coordinates from pixels to physical units
by = (byraw./size(BWs,1))*range(press(:,1))+min(press(:,2))-0.001;      % convert boundary coordinates from pixels to physical units

bx = smooth(bx);
by = smooth(by);







%Nondimensionalize pressure.
cp = press(:,7)/(1000*0.3*0.3);


%Plot the pressure field in ColorBrewer palette and the blanked foil in black
figure(6)
contourf(reshape(press(:,1),128,128),reshape(press(:,2),128,128),reshape(press(:,7),128,128),24,'LineColor','none')     %Plot pressure field
%contourf(reshape(press(:,1),128,128),reshape(press(:,2),128,128),reshape(cp,128,128),24,'LineColor','none')     %Plot Cp field
colormap(flippedcolors)      %Set color palette to ColorBrewer

hold on

fill(bx,by,[0 0 0]);            %Plot the foil
c = colorbar;
caxis([-0.8 0.8])
hold off

xlabel('position [m]')
ylabel('position [m]')
c.Label.String = 'Pressure [Pa]';
%c.Label.String = 'C_{p}';

%save out tiff stack of images
%print('H:\Nonunif_performance_PIV\Paper_Pcode_val\Figures\sample_Pfields\3_3_heave_1.5Hz_1.5cm_0.3ms_frm841.tif', '-dtiffn','-r600');
return                                                                       %remove return when ready to run


%Optionally write interface for calculation to dat files - typically off
%dlmwrite(['foilinterface_thick_' num2str(filenum,numformat) '.dat'], cat(1,zeros(3,2),cat(2,boundx,boundy)), ' ');


j=j+1;                                                                    % increment index for interface movie
M(j) = getframe(gcf);                                             % store frame for movie

c1=c1+2;%*frstep;                                                                %increment indices for x and y data columns in the spline file
c2=c2+2;%*frstep;

fr=fr+frstep;                                                          %Increment frame number (for incr2 PIV)
%pause(0.5)                                                             %pause to see output - comment out to run code faster

end

xlswrite('C:\Users\Lauder Lab\Desktop\xforce_tail.xlsx',forcex)
xlswrite('C:\Users\Lauder Lab\Desktop\yforce_tail.xlsx',forcey)
xlswrite('C:\Users\Lauder Lab\Desktop\torque_tail.xlsx',torque)
%movie2avi(M,'F:\Nonunif_performance_PIV\Pressure code analysis\tail_mid\tail_0angle_1.5Hz_1.5cm_0.3ms_mid\results_3cyc\queen2_rodsim_dT10ms_cb.avi')