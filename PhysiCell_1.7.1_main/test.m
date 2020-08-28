cd output
MCDS=read_MultiCellDS_xml('output00000024.xml');
k = find( MCDS.mesh.Z_coordinates == 16 );
figure(1)
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(1).data(:,:,k),20) ;
caxis([ 0 1 ])
    axis image
    colorbar
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(1).name , ...
        MCDS.continuum_variables(1).units , ...
        MCDS.metadata.current_time , ...
        MCDS.metadata.time_units, ...
        MCDS.mesh.Z_coordinates(k), ...
        MCDS.metadata.spatial_units ) );

    
k = find( MCDS.mesh.Z_coordinates == 16 );
figure(2)
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(2).data(:,:,k),20) ;
caxis([ 0 1 ])
    axis image
    colorbar
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(2).name , ...
        "a.u", ...
        MCDS.metadata.current_time , ...
        MCDS.metadata.time_units, ...
        MCDS.mesh.Z_coordinates(k), ...
        MCDS.metadata.spatial_units ) );
    
    
k = find( MCDS.mesh.Z_coordinates == 16 );
figure(3)
contourf( MCDS.mesh.X(:,:,k), MCDS.mesh.Y(:,:,k), MCDS.continuum_variables(3).data(:,:,k),20) ;
caxis([ 0 1 ])
    axis image
    colorbar
    xlabel( sprintf( 'x (%s)' , MCDS.metadata.spatial_units) );
    ylabel( sprintf( 'y (%s)' , MCDS.metadata.spatial_units) );

    title( sprintf('%s (%s) at t = %3.2f %s, z = %3.2f %s', MCDS.continuum_variables(3).name , ...
        "a.u", ...
        MCDS.metadata.current_time , ...
        MCDS.metadata.time_units, ...
        MCDS.mesh.Z_coordinates(k), ...
        MCDS.metadata.spatial_units ) );
cd ..