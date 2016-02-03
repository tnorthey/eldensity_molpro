function [] = mat2dx_simple(densityMatrix,outfile,minX,minY,minZ,dX)
%MAT2DX_SIMPLE 

input.densityMatrix=densityMatrix;
input.outfile=outfile;
input.minX=minX;
input.minY=minY;
input.minZ=minZ;
input.voxelLength=dX;

mat2dx(input)

end

