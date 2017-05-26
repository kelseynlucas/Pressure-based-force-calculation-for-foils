%May 21, 2016

%This script uses a spline file from the auto foil kinematics vi to
%generate the interface along which forces and torques due to pressure are
%calculated.  Pressure field data comes from the Queen2 pressure field
%calculator.  This version is for use when there are two foil
%segments visible.


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


spldata = dlmread('F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\splines_frm1-2000_incr10.txt');          %read in spline data file


pressfilestub = 'F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\queen2_rodsim_dT10ms_';                                  % string prefix for pressure files

mov=VideoReader('F:\Nonunif_performance_PIV\Data_videos\tail\tail_0angle_1.5Hz_1.5cm_0.3ms_gap.avi');          % read in movie
c1=1;                                                       %initialize index for spline file x-data column count
c2=2;                                                       %initialize index for spline file y-data column count
fr=first;                                                   %set index for movie frame read
frstep=10;




%Make equally spaced x coordinates for the height calculation
xs = linspace(0,0.185,200);
xs = xs';

xcurr = xs(1);
htspl =[];
d=0;

m1=35/129.8;  %Some constants (slope) from formulas describing tail-shaped foil outline
m2=35/55.2;

%Build height dictionary by calculating distance between lines describing
%the tail-shaped foil outline
%For every x in xs,
for i = 1:1:200,

        
    if xs(i,1) < 0.1298,
        h = -2*m1*xs(i,1)+0.090;
        h
    end
        
    %If that distance is more than the distance to the minimum peduncle
    %depth,
    if xs(i,1) >= 0.1298
        h = 2*m2*xs(i,1)-2*0.185*m2+0.090;
    end


    %Collect the height data
    htspl = [htspl; h];
        
        
end








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

%Break the xdata into front and back components, where the front data
%contains the front segment of the foil, etc.
flag = 0;
xfdata = [];
yfdata = [];
xbdata = [];
ybdata = [];
data_index = 1;

%Until hit the gap (indicated by a large jump in x), put the data in the
%front vector.  Then stop.
while flag < 1
    xfdata = [xfdata, xdata(data_index)];
    yfdata = [yfdata, ydata(data_index)];    
    old_index = data_index;
    data_index = data_index + 1;
    
    if abs(xdata(data_index) - xdata(old_index)) > 0.0018
        flag = 1;
    end
    
end

%After the gap, put the data in the back vector.
while data_index <= size(xdata,1)
    xbdata = [xbdata, xdata(data_index)];
    ybdata = [ybdata, ydata(data_index)];
    data_index = data_index + 1;
end


%Plot the foil as two pieces
figure(1)

plot(xfdata,yfdata,'LineWidth',2,'Color','white')                    %plot front spline in white
hold on
plot(xbdata,ybdata,'LineWidth',2,'Color','white')                    %plot front spline in white
set(gca,'Color',[0 0 0]);                       %make background black
axis([min(press(:,1)) max(press(:,1)) min(press(:,2)) max(press(:,2))]);        %scale axes to match those in the velocity file
set(gca,'YDir','normal')                        %make sure y-axis runs in correct direction (positive up)
X=getframe;                                     %get the figure
A=X.cdata;                                      %make it an image
hold off



                                                                %im2BW converts from greyscale to BW; the .## is the threshold to make black versus white
BWs = bwareaopen(im2bw(A,.7),50);                                      % remove nything that isn't the foil

                                                                          
                                                                          %want boundary outside the "not a number" region
se90 = strel('line', 15, 90);                                           % kernel for padding vertical boundaries of tracked object

se0 = strel('line', 15, 0);                                             % kernel for padding horizontal boundaries of tracked object



BWsfinal = imdilate(BWs, [se90 se0]);                                   % pad objects in binary image
BWsfinal=flipud(BWsfinal);                                              %flip binary image so that the origin is in lower left corner


%Turn these lines on only to check the conversion to a binary image - keep OFF when running the code to completion

%figure(2)
%imshow(BWsfinal)
%return



B = bwboundaries(BWsfinal);                                             % find boundaries of remaining objects in binary image

%Get the front part of the foil (larger piece)

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




%Repeat above for second foil bit

[max_Bsize2, max_Bindex2] = min(cellfun('size', B, 1));                   % find 2nd largest object in binary image (assumed to be object of interest)



boundraw2 = B{max_Bindex2};                                               % extract coordinates of object

boundxraw2 = boundraw2(1:length(boundraw2)/75:end,2);                     % create 75-point boundary coordinates

boundyraw2 = boundraw2(1:length(boundraw2)/75:end,1);                     % create 75-point boundary coordinates



boundx2 = (boundxraw2./size(A,2))*range(press(:,1))+min(press(:,1));      % convert boundary coordinates from pixels to physical units

boundy2 = (boundyraw2./size(A,1))*range(press(:,2))+min(press(:,2));      % convert boundary coordinates from pixels to physical units



boundx2 = smooth(boundx2);

boundy2 = smooth(boundy2);



surfdx2 = boundx2 - circshift(boundx2,1);                                  % compute x-spacing between boundary coordinates

surfdy2 = boundy2 - circshift(boundy2,1);                                  % compute y-spacing between boundary coordinates

    

surfnormx2 = surfdy2;                                                     % compute x-component of surface normal

surfnormy2 = -surfdx2;                                                    % compute y-component of surface normal

   

surfunitnormx2 = surfnormx2./(sqrt(surfnormx2.^2 + surfnormy2.^2));         % compute x-component of unit surface normal

surfunitnormy2 = surfnormy2./(sqrt(surfnormx2.^2 + surfnormy2.^2));         % compute y-component of unit surface normal

   

surfposx2 = (boundx2 + circshift(boundx2,1))./2;                           % compute x-midpoint of each object surface facet

surfposy2 = (boundy2 + circshift(boundy2,1))./2;                           % compute y-midpoint of each object surface facet



%Put the boundaries info for the two objects together
boundx = [boundx; boundx2];
boundy = [boundy; boundy2];
surfdx = [surfdx; surfdx2];
surfdy = [surfdy; surfdy2];
surfnormx = [surfnormx; surfnormx2];
surfnormy = [surfnormy; surfnormy2];
surfunitnormx = [surfunitnormx; surfunitnormx2];
surfunitnormy = [surfunitnormy; surfunitnormy2];
surfposx = [surfposx; surfposx2];
surfposy = [surfposy; surfposy2];






%Get the x-value of the first point actually on the foil (accounts for
%boundary dilation)
xMinInt = min(xdata);
%Shift height spline x's to match current axes
posshifted = xs + xMinInt;
%Initialize storage for the height data
hdata =[];

%Look up heights along the boundary of foil
for i = 1:1:length(boundx),
   if boundx(i) <= xMinInt,
       hdata(i) = htspl(1);
   elseif boundx(i) > xMinInt & boundx(i) < max(posshifted),
       [val idx] = min(abs(posshifted-boundx(i)));
       hdata(i) = htspl(idx);
   else,
       hdata(i) = htspl(200);
   end
    
end
hdata = hdata';


%Use the height dictionary to estimate the height at each surface facet's
%midpoint
heights = griddata(boundx,boundy,hdata,surfposx,surfposy,'nearest');


%return


%The following section is for pressure force/torque calculation - (uses the queen2 pressure field to calculate forces/torques along the interface made above)


surfpress = griddata(press(:,1),press(:,2),press(:,7),surfposx,surfposy,'nearest');       % compute pressure at midpoint of each surface facet



time(((fr-first)/frstep)+1) = (fr)*0.001;



forcex(((fr-first)/frstep)+1) = nansum(sqrt(surfdx.^2+surfdy.^2).*-surfunitnormx.*surfpress.*heights);    % compute x-force due to pressure and store in temporal record

forcey(((fr-first)/frstep)+1) = nansum(sqrt(surfdx.^2+surfdy.^2).*-surfunitnormy.*surfpress.*heights);    % compute y-force due to pressure and store in temporal record



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

  

dtorque = -surfpress.*(cross(arm, surfunitnorm, 2)*[0;0;1]).*darea.*heights;     % calculate torque on each surface facet

  

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
bxraw = braw(1:length(braw)/125:end,2);                     % create 125-point boundary coordinates
byraw = braw(1:length(braw)/125:end,1);                     % create 125-point boundary coordinates

bx = (bxraw./size(BWs,2))*range(press(:,1))+min(press(:,1))+0.001;      % convert boundary coordinates from pixels to physical units
by = (byraw./size(BWs,1))*range(press(:,1))+min(press(:,2))-0.001;      % convert boundary coordinates from pixels to physical units

bx = smooth(bx);
by = smooth(by);

%repeat for 2nd foil segment
[max_Bsize2, max_Bindex2] = min(cellfun('size', B, 1));                   % find largest object in binary image (the foil)


braw2 = B{max_Bindex2};                                               % extract coordinates of object
bxraw2 = braw2(1:length(braw2)/75:end,2);                     % create 75-point boundary coordinates
byraw2 = braw2(1:length(braw2)/75:end,1);                     % create 75-point boundary coordinates

bx2 = (bxraw2./size(BWs,2))*range(press(:,1))+min(press(:,1))+0.001;      % convert boundary coordinates from pixels to physical units
by2 = (byraw2./size(BWs,1))*range(press(:,1))+min(press(:,2))-0.001;      % convert boundary coordinates from pixels to physical units

bx2 = smooth(bx2);
by2 = smooth(by2);



%Plot the pressure field in ColorBrewer palette and the blanked foil in black
figure(6)
contourf(reshape(press(:,1),128,128),reshape(press(:,2),128,128),reshape(press(:,7),128,128),24,'LineColor','none')     %Plot pressure field
colormap(flippedcolors)      %Set color palette to ColorBrewer

hold on

fill(bx,by,[0 0 0]);            %Plot the foil
fill(bx2,by2,[0 0 0]);            %Plot the rest of the foil
c = colorbar;
caxis([-80 80])
ylabel(c, 'Pressure [Pa]')
hold off

xlabel('position [m]')
ylabel('position [m]')




%return                                                                       %remove return when ready to run


%Optionally write interface for calculation to dat files - typically off
%dlmwrite(['foilinterface_thick_' num2str(filenum,numformat) '.dat'], cat(1,zeros(3,2),cat(2,boundx,boundy)), ' ');


j=j+1;                                                                    % increment index for interface movie
M(j) = getframe(gcf);                                             % store frame for movie

c1=c1+2;%*frstep;                                                                %increment indices for x and y data columns in the spline file
c2=c2+2;%*frstep;

fr=fr+frstep;                                                          %Increment frame number (for incr2 PIV)
%pause(0.5)                                                             %pause to see output - comment out to run code faster

end

xlswrite('F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\queen2_rodsim_dT10ms_xforce_fixed.xlsx',forcex)
xlswrite('F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\queen2_rodsim_dT10ms_yforce_fixed.xlsx',forcey)
xlswrite('F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\queen2_rodsim_dT10ms_torque_fixed.xlsx',torque)
movie2avi(M,'F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\queen2_rodsim_dT10ms_cb_fixed.avi')