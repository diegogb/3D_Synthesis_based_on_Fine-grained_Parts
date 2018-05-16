%%%%%
% NOTE:
%%%
% THIS CODE IS USED UNDER PERMISSION OF Noa Fish (and Oliver van Kaick)
%%%
% This function is used for the project: "3D Shape Synthesis based on Fine-Grained segments",
% Diego Gonzalez's Master's Thesis - Supervisor: Oliver van Kaick, (Carleton University, 2016 - 2017)
%
% References:
%   (a) Fish, N., Averkiou, M., van Kaick, O., et al. "Meta-representation of shape families". (2014)
%%%%%
%--------------------------------------------------------------------------
function part_rep = calc_axis_dir(part_rep)

axes = part_rep.axes;
exts = axes(4,:)';
cen = part_rep.center;
new_axes = axes;
%new_origin = part_rep.origin;
ref_vec = [1;1;1];
ref_vec = ref_vec / norm(ref_vec);

for i=1:3
    ax = axes(1:3,i);
    ang1 = vec_angle(ref_vec, ax);
    ang2 = vec_angle(ref_vec, -ax);
    
    if(ang2 < ang1)
        ax = -ax;
        new_axes(1:3,i) = ax;
    end
end

new_origin = cen;
new_origin = new_origin - 0.5*exts(1)*new_axes(1:3,1);
new_origin = new_origin - 0.5*exts(2)*new_axes(1:3,2);
new_origin = new_origin - 0.5*exts(3)*new_axes(1:3,3);

part_rep.axes = new_axes;
part_rep.origin = new_origin;
