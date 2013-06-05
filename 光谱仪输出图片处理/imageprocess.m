function new_curve = imageprocess(filename,lamda_start,lamda_end,power_start,power_end,x,y)
%this function extracts data from image acquired by the OSA
%
% USAGE:
%
%new_curve = imageprocess(filename,lamda_start,lamda_end,power_start,power_end,x,y)
%
% INPUT:
%
% filename        filename of the image
% lamda_start     start wavelength in the image
% lamda_end       end wavelength in the image
%
%
% x               the x position of the leftdown dot of the sweeping grid
% y               the y position of the leftdown dot of the sweeping grid
original_image = imread(filename);
new_curve = zeros(1000,2);
d_lamda = (lamda_end-lamda_start)/500;
d_power = (power_end-power_start)/344;
dots = 1;
for column=x:size(original_image,2)
   for row=y-344:size(original_image,1)
      if original_image(row,column,1)>=100 && original_image(row,column,2)==0 ...
              && abs(original_image(row,column,3)-original_image(row,column,1))<=10
          new_curve(dots,1) = lamda_start + d_lamda*(column-x-1);
          new_curve(dots,2) = power_start + (y-row)*d_power;
          dots=dots+1;
      end
   end
end
new_curve = new_curve(any(new_curve,2),:);

