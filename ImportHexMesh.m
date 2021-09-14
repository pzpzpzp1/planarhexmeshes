function mesh = ImportHexMesh(meshFile)
if nargin==0
    meshFile = 'meshes/bone.mesh';
end

[~, ~, ext] = fileparts(meshFile);

if strcmp(ext, '.mesh')
    fid = fopen(meshFile,'r');
    textscan(fid, 'MeshVersionFormatted %d');
    textscan(fid, 'Dimension %d', 'Whitespace', ' \t\n', 'MultipleDelimsAsOne', true);
    
    textscan(fid, 'Vertices %d', 'Whitespace', ' \t\n', 'MultipleDelimsAsOne', true);
    verts = textscan(fid, '%f %f %f %*d', 'CollectOutput', true, 'MultipleDelimsAsOne', true);
    
    textscan(fid, 'Quadrilaterals %d', 'Whitespace', ' \t\n');
    textscan(fid, '%d %d %d %d %*d', 'CollectOutput', true, 'MultipleDelimsAsOne', true);
    
%     textscan(fid, 'Edges %d', 'Whitespace', ' \t\n');
%     textscan(fid, '%d %d %*d', 'CollectOutput', true, 'MultipleDelimsAsOne', true);
    
    textscan(fid, 'Hexahedra %d', 'Whitespace', ' \t\n');
    hexes = textscan(fid, '%d %d %d %d %d %d %d %d %*d', 'CollectOutput', true, 'MultipleDelimsAsOne', true);
    
    fclose(fid);

    verts = verts{:};
    hexes = double(hexes{:});
else
    error('unknown type');
end
mesh.points = verts;
mesh.cells = hexes;

end