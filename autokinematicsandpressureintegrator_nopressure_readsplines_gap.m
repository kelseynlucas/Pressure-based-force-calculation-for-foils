%May 21, 2016

%This script uses a spline file from the auto foil kinematics vi to
%generate the boundary interface used for blanking in the Queen2 pressure
%field calculator.  This version is for use when there are two foil
%segments visible.




clear all                                                     
close all                                                    %close open figures

j=0;                                                         %initialize index for movie frame storage

first = 1;                                                  % first image and pressure file number (note that, in this case, pressure is actually a velocity file)

last = 2000;                                                 % last image and pressure file number

increment = 10;                                              % increment between file numbers



numformat = '%05d';                                         % format of numbers in file name. '%05d' is five digit format with leading zeros

delimiter = '\t';                                             % delimiter between columns in pressure file

headerlines = 1;                                            % number of header lines in pressure file


spldata = dlmread('F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_heave_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\splines_frm1-2000_incr10.txt');          %read in spline data file


pressfilestub = 'F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_heave_1.5Hz_1.5cm_0.3ms_gap\full_stack\B';             % string prefix for *velocity* files

mov=VideoReader('F:\Nonunif_performance_PIV\Data_videos\tail\tail_heave_1.5Hz_1.5cm_0.3ms_gap.avi');          % read in movie
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

%Split data into front "f" and back "b" data lists
%Marker for when we've moved from front to back
flag = 0;
%Initialize storage
xfdata = [];
yfdata = [];
xbdata = [];
ybdata = [];
%set up an index
data_index = 1;

%While we are still on the front segment,
while flag < 1
    %Add outline points to the front list
    xfdata = [xfdata, xdata(data_index)];
    yfdata = [yfdata, ydata(data_index)];    
    old_index = data_index;
    data_index = data_index + 1;
    
    %If the gap between current x data point and next data point is really big, we're in the gap
    %between foil segments.
    if abs(xdata(data_index) - xdata(old_index)) > 1.8
        flag = 1;
    end
    
end

%While we're looking at the back segment, add points to the back list
while data_index <= size(xdata,1)
    xbdata = [xbdata, xdata(data_index)];
    ybdata = [ybdata, ydata(data_index)];
    data_index = data_index + 1;
end


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

%return


                                                             %im2BW converts from greyscale to BW; the .## is the threshold to make black versus white
BWs = bwareaopen(im2bw(A,.6),50);                                      % remove anything that isn't the foil

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

boundxraw = boundraw(1:length(boundraw)/125:end,2);                     % create 125-point boundary coordinates

boundyraw = boundraw(1:length(boundraw)/125:end,1);                     % create 125-point boundary coordinates



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



%Here, add in a thin area connecting the front and back foil pieces so that
%the blanking can perform correctly (assumes one object, but if the thin
%area is small enough to not enclose any points on the velocity vector
%field grid, it would be invisible to queen2).


%Find the trailing edge x and y coordinates of the front portion, and the
%row these coordinates occur in.
[TEfx, front_index] = max(boundx);
TEfy = boundy(front_index);

%Assume the first point in the points list for the back foil portion is the
%leading edge.  Extract the coordinates.
LEbx = boundx2(1);
LEby = boundy2(1);

%Find the trailing edge of the back portion and note the row the
%coordinates are in.
[TEbx, back_index] = max(boundx2);
TEby = boundy2(back_index);


%Initialize storage for the blanking boundaries master points list.
blankx = [];
blanky = [];

%Initialize an index
currin = 1;

%While the index is less the the index of the front portion's trailing
%edge,
while currin <= front_index
    %Add the boundary points to the master list
    blankx = [blankx; boundx(currin)];
    blanky = [blanky; boundy(currin)];
    
    %Go on to the next row
    currin = currin + 1;
end


%Determine if the leading edge of the back portion is higher or lower than
%the trailing edge of the front portion.
if LEby-TEfy <=0
    pathud = 0;
else
    pathud = 1;
end


%If we're on the first velocity field, extract the grid spacing
%information.
if filenum == first
    velgridx = unique(press(:,1));
    velgridy = unique(press(:,2));
    xgridspacing = abs(velgridx(2) - velgridx(1));
    ygridspacing = abs(velgridy(2) - velgridy(1));
    
end

%Initialize an index
i = 1;

%Check to make sure that the leading/trailing edge coordinates aren't very
%close to a velocity grid point.  If they are, displace them slightly.
%Just ensures that the thin area won't ever contain velocity vectors due to
%a rounding mistake.

while i <= size(velgridx,1)
    testTEx = abs(velgridx(i)-TEfx);
    
    if testTEx < 0.15*xgridspacing
        TEfx = TEfx + 0.2*xgridspacing;
    end
        
    i=i+1;
    
end

i=1;
while i <= size(velgridx,1)
    testTEy = abs(velgridy(i)-TEfy);
    
    if testTEy < 0.15*ygridspacing
        TEfy = TEfy + 0.2*ygridspacing;
    end
        
    i=i+1;
    
    
end
i=1;
while i <= size(velgridx,1)
    testLEx = abs(velgridx(i)-LEbx);
    
    if testLEx < 0.15*xgridspacing
        LEbx = LEbx + 0.2*xgridspacing;
    end
    
    i=i+1;
    
    
end
i=1;
while i <= size(velgridx,1)
    testLEy = abs(velgridy(i)-LEby);
    
    if testLEy < 0.15*ygridspacing
        LEby = LEby + 0.2*ygridspacing;
    end
        
    i=i+1;
    
    
end



%Initialize storage for the thin area's boundary coordinates.  Because the
%boundaries are read in a counterclockwise fashion, we need a list of
%points for moving left to right along one side of the area, and a list
%of points for moving right to left back on the other side of the area.
xpts_l = [];
ypts_l = [];
xpts_r = [];
ypts_r = [];

%Set the starting values for x and y coordinates on the area's boundaries.

%Will travel horizontally until we reach the leading edge x coordinate of
%the back portion, and then vertically up or down until we reach the
%leading edge of the back portion.  This way, we make sure that we always
%are moving between grid points.

xcurr = TEfx + xgridspacing;
ycurr = TEfy;

%Give the area thickness
othery = ycurr+0.1*ygridspacing;

%Until we reach the leading edge of the back portion, step along the
%x-axis, keep y constant.  Add points to the end of one list set to get
%left-to-right travel, and add points to the beginning of the other list so
%this list, when complete, is a list of points traveling
%right-to-left.

while xcurr < LEbx
    xpts_l = [xpts_l; xcurr];
    ypts_l = [ypts_l; ycurr];
    
    xpts_r = [xcurr; xpts_r];
    
    ypts_r = [othery; ypts_r];    
    xcurr = xcurr + xgridspacing;
    
end

%In the loop above, we end with xcurr being one step past the leading edge 
%of the back portion.  So, take a step back.
xcurr = xcurr - xgridspacing;

%Since we'll travel vertically, give the area thickness using the
%x-dimension.
otherx = xcurr+0.1*xgridspacing;

%Add in a point to turn the corner on the way back (accounts for the area's
%thickness)
xpts_r = [otherx; xpts_r];
ypts_r = [othery; ypts_r];

%Do this version if traveling down.
if pathud == 0
    
    %As long as ycurr is bigger than the leading edge y coordinate (we're
    %above the leading edge), take steps down and add points.
    while ycurr > LEby
       
        xpts_l = [xpts_l; xcurr];
        ypts_l = [ypts_l; ycurr];
    
        xpts_r = [otherx; xpts_r];
        ypts_r = [ycurr; ypts_r];
    
        ycurr = ycurr - ygridspacing;
    
    end
    
end

%Do this version if traveling up.
if pathud == 1
    
    %As long as ycurr is smaller than the leading edge y coordinate (we're
    %below the leading edge), take steps up and add points.
    while ycurr < LEby
       
        xpts_l = [xpts_l; xcurr];
        ypts_l = [ypts_l; ycurr];
    
        xpts_r = [otherx; xpts_r];
        ypts_r = [ycurr; ypts_r];
    
        ycurr = ycurr + ygridspacing;
    
    end
    
end


%Add the left-to-right list to the master points list
blankx = [blankx; xpts_l];
blanky = [blanky; ypts_l];


%Add every point in the back part of the foil to the master list
blankx = [blankx; boundx2];
blanky = [blanky; boundy2];

%Add the right-to-left thin area list to the master list
blankx = [blankx; xpts_r];
blanky = [blanky; ypts_r];

%Add all the points after the row where the trailing edge is to the master
%list (this part of boundx&boundy indicate right-to-left travel along the
%front portion of the foil)
currin = front_index;

while currin <= size(boundx,1)
    
    blankx = [blankx; boundx(currin)];
    blanky = [blanky; boundy(currin)];
    currin = currin + 1;
end




%Put the boundaries info for the two objects together for plotting purposes
boundplotx = [boundx; boundx2];
boundploty = [boundy; boundy2];







%Make boundaries:

figure(2)
imshow(im,'XData',[min(press(:,1)) max(press(:,1))],'YData',[min(press(:,2)) max(press(:,2))]);     %Display original foil image, scaled to match axes in the velocity file

axis on
set(gca,'YDir','normal')            %make sure y-axis is oriented correctly (positive up)

hold on
plot(boundplotx,boundploty,'.')                                                 % plot object boundary

quiver(surfposx, surfposy, surfunitnormx, surfunitnormy,1,'y')
quiver(surfposx2, surfposy2, surfunitnormx2, surfunitnormy2,1,'y')
hold off







%return                                                                       %remove return when ready to run



%Write boundaries to dat file
dlmwrite(['F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_heave_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\interface_rodsim_' num2str(filenum,numformat) '.dat'], cat(1,zeros(1,2),cat(2,blankx,blanky)), '\t');


j=j+1;                                                                    % increment index for interface movie
M(j) = getframe(gcf);                                             % store frame for movie

c1=c1+2;%*frstep;                                                                %increment indices for x and y data columns in the spline file
c2=c2+2;%*frstep;
fr=fr+frstep;
%pause(0.5)                                                             %pause to see output - comment out to run code faster



end

movie2avi(M, 'F:\Nonunif_performance_PIV\Pressure code analysis\tail_gap\tail_heave_1.5Hz_1.5cm_0.3ms_gap\results_3cyc\interface_rodsim_incr10.avi')