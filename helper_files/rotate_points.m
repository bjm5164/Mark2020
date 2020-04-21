function [rp] = rotate_points(points,theta,axis)
load Neuropil_Mesh_Object.mat
if nargin == 2
    axis = input('Which Axis? X:1 Y:2 Z:3')

    else
    end

theta_r = deg2rad(theta);

if axis == 1
    R = makehgtform('xrotate',theta_r);
elseif axis == 2
    R = makehgtform('yrotate',theta_r);
elseif axis == 3
    R = makehgtform('zrotate',theta_r);
else
    error('Wrong Axis')
end

R = R(1:3,1:3);
s = points';

axis_lims_MLDV = [min(NPM.vertices(:,1)), max(NPM.vertices(:,1)),min(NPM.vertices(:,2)), max(NPM.vertices(:,2))];
axis_lims_APDV = [min(NPM.vertices(:,3)), max(NPM.vertices(:,3)),min(NPM.vertices(:,1)), max(NPM.vertices(:,1))];
center.x = (axis_lims_MLDV(2) + axis_lims_MLDV(1))/2 ;
center.y = (axis_lims_MLDV(4)+axis_lims_MLDV(3))/2;
center.z = (axis_lims_APDV(4)+axis_lims_MLDV(3))/2;

center_coords = repmat([center.x,center.y,center.z]',1,length(s(1,:)));
so = R*(s-center_coords) + center_coords;

rp = so';

end