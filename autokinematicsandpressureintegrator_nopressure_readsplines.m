%May 21, 2014

%This script uses a spline file from the auto foil kinematics vi to
%generate the boundary interface used for blanking in the Queen2 pressure
%field calculator.




clear all                                                     
close all                                                    %close open figures

j=0;                                                         %initialize index for movie frame storage

first = 1;                                                  % first image and pressure file number (note that, in this case, pressure is actually a velocity file)

last = 2000;                                                 % last image and pressure file number

increment = 10;                                              % increment between file numbers



numformat = '%05d';                                         % format of numbers in file name. '%05d' is five digit format with leading zeros

delimiter = '\t';                                             % delimiter between columns in pressure file

headerlines = 1;                                            % number of header lines in pressure file


spldata = dlmread('F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\splines_frm1-2000_incr10.txt');          %read in spline data file


pressfilestub = 'F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\full_stack\B';             % string prefix for *velocity* files

mov=VideoReader('F:\Nonunif_performance_PIV\Data_videos\tail\tail_0angle_1.5Hz_1.5cm_0.3ms_gap.avi');          % read in movie
c1=1;                                                       %initialize index for spline file x-data column count
c2=2;                                                       %initialize index for spline file y-data column count
fr=1;                                                   %set index for movie frame read
frstep=10;

for filenum = first:increment:last                          % loop through files                                                 
clear A                                                     %clear last iteration's image to start fresh
clear X                                                     %clear last iteration's image
    clf                                                     % clear current figure

    filenum                                                 % display current file number
    
press = dlmread([pressfilestub num2str(filenum,numformat) '.dat'],delimiter,headerlines,0); % read in velocity file

    xdata=spldata(2:end,c1);       %Get x data for the spline
    xdata=min(press(:,1))+xdata;    %Scale data to match velocity file
    ydata=spldata(2:end,c2);       %Get y data for the spline
    ydata=max(press(:,2))-ydata-2;      %Scale data to match velocity file and center spline on foil
    
im=flipud(read(mov,spldata(1,c1)));    %Get current image, flip upside-down to move origin from upper left to lower left

im = imfill(im);                                                          % even out image texture



% Remove image boundaries

%im(1:200,:) = 0;                                                         %45+47 specific to dataset, which gets rid of bright walls

%im(700:end,:) = 0;





figure(1)

plot(xdata,ydata,'LineWidth',2,'Color','white')                    %plot spline in white
set(gca,'Color',[0 0 0]);                       %make background black
axis([min(press(:,1)) max(press(:,1)) min(press(:,2)) max(press(:,2))]);        %scale axes to match those in the velocity file
set(gca,'YDir','normal')                        %make sure y-axis runs in correct direction (positive up)
X=getframe;                                     %get the figure
A=X.cdata;                                      %make it an image

%return


                                                             %im2BW converts from greyscale to BW; the .## is the threshold to make black versus white
BWs = bwareaopen(im2bw(A,.6),200);                                      % remove anything that isn't the foil

BWs(1:5,:) = 0;                                                         %45+47 specific to dataset, which gets rid of bright walls

BWs(425:end,:) = 0;

BWs(:,1:5) = 0;
BWs(:,325:end) = 0;

se90 = strel('line', 10, 90);                                           % kernel for padding vertical boundaries of tracked object

se0 = strel('line', 10, 0);                                             % kernel for padding horizontal boundaries of tracked object



BWsfinal = imdilate(BWs, [se90 se0]);                                   % pad objects in binary image; smooths out tracked object
BWsfinal=flipud(BWsfinal);                                              %flip binary image so that the origin is in lower left corner


%Turn these lines on only to check the conversion to a binary image - keep OFF when running the code to completion

%figure(2)
%imshow(BWsfinal)
%axis on
%return



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







%Make boundaries:

figure(2)
imshow(im,'XData',[min(press(:,1)) max(press(:,1))],'YData',[min(press(:,2)) max(press(:,2))]);     %Display original foil image, scaled to match axes in the velocity file

axis on
set(gca,'YDir','normal')            %make sure y-axis is oriented correctly (positive up)

hold on
plot(boundx,boundy,'.')                                                 % plot object boundary

quiver(surfposx, surfposy, surfunitnormx, surfunitnormy,1,'y')
hold off








return                                                                       %remove return when ready to run



%Write boundaries to dat file
dlmwrite(['F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\interface_rodsim_' num2str(filenum,numformat) '.dat'], cat(1,zeros(1,2),cat(2,boundx,boundy)), '\t');


j=j+1;                                                                    % increment index for interface movie
M(j) = getframe(gcf);                                             % store frame for movie

c1=c1+2%*frstep;                                                                %increment indices for x and y data columns in the spline file
c2=c2+2%*frstep;
fr=fr+frstep;
%pause(0.5)                                                             %pause to see output - comment out to run code faster



end

movie2avi(M, 'F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_0angle_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\interface_rodsim_incr10.avi')