function [Out] = reshape_3dbin(name,nx,ny)

file = name;
fid = fopen(file,'r');
bin = fread(fid,'double');
t_max = length(bin)/nx/ny;
Out = zeros(nx,ny,t_max);

for t = 1:t_max
    for j = 1:ny
        for i = 1:nx
            Out(i,j,t) = bin(i+(j-1)*nx+(t-1)*nx*ny);
        end
    end
end