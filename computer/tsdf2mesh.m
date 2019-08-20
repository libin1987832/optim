% ---------------------------------------------------------
% Loads a TSDF voxel grid from a binary file (tsdf.bin) and
% creates a mesh (saved to mesh.ply), which can be viewed
% with a 3D viewer like Meshlab.
%
% Author: Andy Zeng, Princeton University, 2016
% ---------------------------------------------------------

% Load TSDF voxel grid from binary file
%  x=[3,10,20,35];
% y1=[0.823,0.833,0.8698,0.8514];
% y2=[0.9833,0.9737,0.9912,0.9847];
% y3=[0.7,0.5789,0.57017,0.504];
% plot([x;x;x]',[y1;y2;y3]')
% legend('raw points','RANSAC','the area filter')
% y1=[0.640625,0.032,0.0247,0.002437];
% y2=[0.46875,0.07037,0.0222,0.007784];
% y3=[0.46875,0.07037,0.0222,0.007427];
% plot([x;x;x]',[y1;y2;y3]')
% legend('raw points','RANSAC','the area filter')
% return;




clear;
fid = fopen('F://tsdf//test_app.bin','rb');
tsdfHeader = fread(fid,8,'single');
voxelGridDim = tsdfHeader(1:3);
voxelGridOrigin = tsdfHeader(4:6);
voxelSize = tsdfHeader(7);
truncMargin = tsdfHeader(8);

tsdf_r = fread(fid,voxelGridDim(1)*voxelGridDim(2)*voxelGridDim(3),'single');
fclose(fid);




fid = fopen('F://tsdf//test_app_color.bin','rb');
tsdfHeader_c = fread(fid,8,'single');
voxelGridDim_c = tsdfHeader_c(1:3);
voxelGridOrigin_c = tsdfHeader_c(4:6);
voxelSize_c = tsdfHeader_c(7);
truncMargin_c = tsdfHeader_c(8);
tsdf_c = fread(fid,voxelGridDim_c(1)*voxelGridDim_c(2)*voxelGridDim_c(3),'int32');
tsdf_c=int32(tsdf_c);
fclose(fid);
%tsdf=tsdf(1:2:end);
% Convert from TSDF to mesh  
% tsdf=ones([voxelGridDim(1),voxelGridDim(2),voxelGridDim(3)]);
% [m,n]=size(tsdf_r);
% tsdf_index = reshape(tsdf_r,[4,m/4]);
% for i=1:m/4
%     tsdf(tsdf_index(1,i)+1,tsdf_index(2,i)+1,tsdf_index(3,i)+1)=tsdf_index(4,i);
% end








 tsdf = reshape(tsdf_r,[voxelGridDim_c(1),voxelGridDim_c(2),voxelGridDim_c(3)]);
 tsdf_cc = reshape(tsdf_c,[voxelGridDim_c(1),voxelGridDim_c(2),voxelGridDim_c(3)]);
%  tsdf_cc(:,:,:)=8454016;
 color=tsdf_cc;
 [nY,nX,nZ] = size(tsdf); 
 [X,Y,Z] = meshgrid(1:nX,1:nY,1:nZ); 
 [faces_1,points_1,col_1] = MarchingCubes(X,Y,Z,tsdf,0,tsdf_cc);
 faces=faces_1';
  points=points_1';
  col=col_1';
 uint8data=uint8(repmat([0;0;0],1,size(points,2)));
 for i=1:size(points,2)
%     uint8data(1,i) = bitand( bitshift(col(i),0), 255);
%     uint8data(2,i) = bitand( bitshift(col(i), -8), 255);
%     uint8data(3,i) = bitand( bitshift(col(i), -16), 255);
    int_x = floor(points(1,i));minEven_x = int_x+mod(int_x,2); 
    int_y = floor(points(2,i));minEven_y = int_y+mod(int_y,2); 
    int_z = floor(points(3,i));minEven_z = int_z+mod(int_z,2); 
    col(i)=0;
    for x_i=-2:2
        for y_i=-2:2
            for z_i=-2:2
                if minEven_x+x_i>0 & minEven_y+y_i>0 & minEven_z+z_i>0 & tsdf_cc(minEven_x+x_i,minEven_y+y_i,minEven_z+z_i)>0
                    col(i)=tsdf_cc(minEven_x+x_i,minEven_y+y_i,minEven_z+z_i);
                end
            end
        end
    end
    uint8data(1,i) = bitand( bitshift(col(i),0), 255);
    uint8data(2,i) = bitand( bitshift(col(i), -8), 255);
    uint8data(3,i) = bitand( bitshift(col(i), -16), 255);
 end
%  return
% [x,y,z] = meshgrid(-2:4/voxelGridDim(1):2-4/voxelGridDim(1),-2:4/voxelGridDim(2):2-4/voxelGridDim(2),-2:4/voxelGridDim(3):2-4/voxelGridDim(3)); 
% [F,V] = MarchingCubes(x,y,z,tsdf,0);
% xlabel('x');ylabel('y');zlabel('z'); 
% grid on; view([1,1,1]); axis equal on; camlight; lighting gouraud
% returnr
% fv = isosurface(tsdf,0);
% points = fv.vertices';
% faces = fv.faces';


% Set mesh color (light blue)
%color = uint8(repmat([175;198;233],1,size(points,2)));
% uint8data=uint8(repmat([175;198;233],1,size(points,2)));
% for i=1:size(points,2)
% x=points(1,i);
% y=points(2,i);
% z=points(3,i);
% dx=fix(x);Cx=fix(x+1);dz=fix(z);Cz=fix(z+1);dy=fix(y);Cy=fix(y+1);
% if Cx>255
%     dx=255;Cx=255;
% end
% if Cy>255
%     dy=255;Cy=255;
% end
% if Cz>255
%     dz=255;Cz=255;
% end
% if color(dx,dy,Cz)>0 
%     dz=Cz;
% elseif color(dx,Cy,dz)>0
%     dy=Cy;
% elseif color(dx,Cy,Cz)>0;
%     dy=Cy;dz=Cz;
%  elseif color(Cx,dy,dz)>0
%     dx=Cx;
% elseif color(Cx,Cy,dz)>0;
%      dx=Cx;dy=Cy;
% elseif color(Cx,dy,Cz)>0
%      dx=Cx;dz=Cz;
% elseif color(Cx,Cy,Cz)>0;
%      dx=Cx;dy=Cy;dz=Cz;
% end
% if color(dx,dy,dz)>0    
% uint8data(3,i) = bitand( bitshift(color(dx,dy,dz),0), 255);
% uint8data(2,i) = bitand( bitshift(color(dx,dy,dz), -8), 255);
% uint8data(1,i) = bitand( bitshift(color(dx,dy,dz), -16), 255);
% end
% end


% Transform mesh from voxel coordinates to camera coordinates
meshPoints(1,:) = voxelGridOrigin(1) + points(2,:)*voxelSize_c; % x y axes are swapped from isosurface
meshPoints(2,:) = voxelGridOrigin(2) + points(1,:)*voxelSize_c;
meshPoints(3,:) = voxelGridOrigin(3) + points(3,:)*voxelSize_c;

% Write header for mesh file
data = reshape(typecast(reshape(single(meshPoints),1,[]),'uint8'),3*4,[]);
% data = [data; color];
data = [data; uint8data];
fid = fopen('F://tsdf//mesh2.ply','w');
fprintf (fid, 'ply\n');
fprintf (fid, 'format binary_little_endian 1.0\n');
fprintf (fid, 'element vertex %d\n', size(data,2));
fprintf (fid, 'property float x\n');
fprintf (fid, 'property float y\n');
fprintf (fid, 'property float z\n');
fprintf (fid, 'property uchar red\n');
fprintf (fid, 'property uchar green\n');
fprintf (fid, 'property uchar blue\n');
fprintf (fid, 'element face %d\n', size(faces,2));
fprintf (fid, 'property list uchar int vertex_index\n');
fprintf (fid, 'end_header\n');

% Write vertices
fwrite(fid, data,'uint8');

% Write faces
faces = faces([3 2 1],:); % reverse the order to get a better normal    
faces_data = int32(faces-1);
faces_data = reshape(typecast(reshape(faces_data,1,[]),'uint8'),3*4,[]);
faces_data = [uint32(ones(1,size(faces,2))*3); faces_data];        
fwrite(fid, faces_data,'uint8');

fclose(fid);