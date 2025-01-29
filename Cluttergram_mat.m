function [f_x,ti,RR] = Cluttergram_mat(resolution,mola_res,m,M)


% resolution: Index-step of the orbiter tranjectory
% mola_res: Index-step of the MOLA (topography) 
% m: 
% M: x,y,z points of topography
% m: 

% f_x: Distance of the orbiter with respect to the first location of the
% orbiter in the survey
% ti: Time axis
% RR: 2D matrix with cluttergram






% Mars radius
R=3389.5*1000;
% Time step of the data
time_step=(4.6009e-07)/2;

% Transform the latitude,longitude to x,y,z
x = R*cosd(M(:,2)).*cosd(M(:,1));
y = R*cosd(M(:,2)).*sind(M(:,1));
z = R*sind(M(:,2));

% Find the minimum and maximum longitude and latitude
Mi_La=min(M(:,2));
Ma_La=max(M(:,2));
Mi_Lo=min(M(:,1));
Ma_Lo=max(M(:,1));


% Calculate the distance from the center of Marx
d2=R+M(:,3);    

% Scale for plotting purpolses
lamda=d2./R;
x=lamda.*x;
y=lamda.*y;
z=lamda.*z;    
disp('Displaying MOLA and SHARAD Orbit')
tri = delaunay(x,y);
[r,c] = size(tri);
lp=M(:,3);
h = trisurf(tri, x, y, z,lp);shading interp;
colormap('hot')
xlabel('x-axis (m)');
ylabel('y-axis (m)');
zlabel('z-axis (m)');
axis equal
axis vis3d
view(135,-6)
axis off


    

% Open the coordinates of the oribter (longtidue, magnitude, altitude)
ii=0;
for i=1:1:length(m.data(:,1));
    if m.data(i,2)>Mi_Lo & m.data(i,2)<Ma_Lo & m.data(i,1)>Mi_La & m.data(i,1)<Ma_La;
        ii=ii+1;
        La(ii)=m.data(i,1);
        Lo(ii)=m.data(i,2);
        To(ii)=mean(m.data(:,4))*1000;
    end
end

% Transform to x,y,z
x_a1 = To.*cosd(La).*cosd(Lo);
y_a1 = To.*cosd(La).*sind(Lo);
z_a1 = To.*sind(La);

% Reduce the resolution of the measurements
xp=1:length(x_a1);
nxp=1:resolution:length(x_a1);
x_a=interp1(xp,x_a1,nxp,'spline');
y_a=interp1(xp,y_a1,nxp,'spline');
z_a=interp1(xp,z_a1,nxp,'spline');

% This is for plotting the projection of the satelite orbiter to the
% surface
x_aa = R.*cosd(La).*cosd(Lo);
y_aa = R.*cosd(La).*sind(Lo);
z_aa = R.*sind(La);


% We need this part for plotting the cluttergram later on; 
% f_x is the distance of the orbiter from the begining of the survey
for i=2:length(x_a);
    b1=[x_a(i), y_a(i), z_a(i)];
    b2=[x_a(1), y_a(1), z_a(1)];
    f_x(i)=norm(b2-b1);
end
f_x=[0,f_x];

% Plotting the position of the orbiter and its projection to the surface
hold on;
plot3(x_a,y_a,z_a,'k-','linewidth',2);hold on
plot3(x_aa,y_aa,z_aa,'g-','linewidth',2);axis equal
legend('MOLA Map','Actual Orbit', 'Projected Orbit')



figure;




% Initialise the cluttergram
B=zeros(10000,length(x_a));
r=[x,y,z];
ptCloud = pointCloud(r);
% Estimating the normal vectors from a point cloud.
normals = pcnormals(ptCloud);

reverseStr='';
for i=1:length(x_a)-1;
   
   msgg = ['Calculating : ', num2str(round(100*i/length(x_a))), '%%'];
   msg= sprintf('%s ',msgg);
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msgg));

   c1=[x_a(i),y_a(i),z_a(i)];
   cc2=[x_a(i+1),y_a(i+1),z_a(i+1)];
   
   % Vector n1 is tangential to the orbiter tranjectory
   n1=(cc2-c1);
   
   % cross(n1,c1) is the unit vector perpedicular to the E-plane
   % We need that to incorporate the radiation pattern of the antenna in
   % the cluttergram
   dipole=cross(n1,c1)./norm(cross(n1,c1));
    
   ii=0;
   % Go through all the points in the MOLA topography
   for j= 1:mola_res:length(x);
      c2=[x(j),y(j),z(j)];
      % Calcualte the distance between the orbiter position and a given
      % x,y,z position in MOLA
      di=norm(c2-c1);
      
      % Calculate the two-ways travel time 
      s=(2*di/(2.99*10^8));
      % If the MOLA position is too far away, then omit this measurement
      if s<0.00195;
      
   
          vec=c1-c2;
          vec2=normals(j,:);
   
          % Calculate the angle of the ray path to the plane x,y,z
          degr=abs(acosd(dot(vec, vec2) / (norm(vec) * norm(vec2))));
          degr = degr+90;
          degr=abs(180-degr); 
         
          degr = deg2rad(degr);
   
          % Calculate the directivity of the dipole antenna
          proj_n1=((vec*n1')/(norm(n1)^2))*n1;
          pro=vec-proj_n1;
          theta=acosd((pro*dipole')/(norm(pro)*norm(dipole)));
          Par=sind(theta).^3;
          
          B(round(s/time_step) + 1,i)= B(round(s/time_step) + 1,i)+(degr^4)*Par/((di)^4);
       end
   end
        
end

% Migrate the data to collapse hyperbolae
[migRF zz] = ezfkmig(B,time_step,mean(diff(f_x(3:end))),1);

% Process data for plotting purposes
RR=imadjust(migRF);

subplot(131),imagesc(f_x,2*zz/(3*10^8),B);
title("Un-focused raw data")
axis([f_x(1),f_x(end),1.88*10^-3, 1.96*10^-3]);
xlabel('Distance (m)');
ylabel('Time (s)');
colormap('bone')

subplot(132),imagesc(f_x,2*zz/(3*10^8),RR);
title("Focused")

axis([f_x(1),f_x(end),1.88*10^-3, 1.96*10^-3]);
xlabel('Distance (m)');
ylabel('Time (s)');
colormap('bone')
subplot(133),imagesc(f_x,2*zz/(3*10^8),RR.^2);
title("Focused and squared")
axis([f_x(1),f_x(end),1.88*10^-3, 1.96*10^-3]);
xlabel('Distance (m)');
ylabel('Time (s)');
colormap('bone')
ti = 2*zz/(3*10^8);

end

