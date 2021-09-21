function flipmesh_on_harddrive(filename)
    % filename = 'meshes/sing1.vtk'; %  DO NOT RUN WITHOUT ARGS
        
    mesh = load_vtk(filename);
    mesh.cells = mesh.cells(:,[5 6 7 8 1 2 3 4]);
    save_vtk(mesh,filename);

end