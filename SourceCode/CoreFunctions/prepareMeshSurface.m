function structSurface = prepareMeshSurface( structMeshPreparation, doPlot) 
    warning('off', 'MATLAB:MKDIR:DirectoryExists');

    % Load the mesh
    meshIn = stlreadNonCleaning( structMeshPreparation.fileNameMesh );
  
    % Clean some issues with the mesh, so it behaves nicer.
    meshClean = surfaceMesh(meshIn.vertices,meshIn.faces);
    removeDefects(meshClean,'duplicate-vertices');
    removeDefects(meshClean,'duplicate-faces');
    removeDefects(meshClean,'nonmanifold-edges');
    surfaceVerticesClean = meshClean.Vertices;
    
    indexesCleanToRawMesh = meshClean.Faces';
    indexesCleanToRawMesh = indexesCleanToRawMesh(:);
    surfaceVertices = surfaceVerticesClean( indexesCleanToRawMesh, : );
    surfaceFaces = reshape( 1:length(indexesCleanToRawMesh), 3, length(indexesCleanToRawMesh)/3 )';

    % Remove the mean
    surfaceVertices = surfaceVertices - mean( surfaceVertices );
    
    % Rotate and position the mesh
    rotMatMesh = rotz(structMeshPreparation.orientation(3) ) *  roty(structMeshPreparation.orientation(2) ) * rotx(structMeshPreparation.orientation(1) );
    surfaceVertices = surfaceVertices * rotMatMesh;
    surfaceVertices = surfaceVertices * structMeshPreparation.vertexScaling;
    surfaceVertices = surfaceVertices + structMeshPreparation.position;

    % Save the mesh in a local struct for normal calculation
    localStructSurface = struct();
    localStructSurface.faces = surfaceFaces;
    localStructSurface.vertices = surfaceVertices;
    surfaceNormals = patchnormals( localStructSurface );
    if( structMeshPreparation.FLIPNORMALS==1 )
        for cntNormal = 1 : size( surfaceNormals, 1 )
            surfaceNormals( cntNormal, : ) = -surfaceNormals( cntNormal, : ) / norm(  surfaceNormals( cntNormal, : ) );
        end
    else
        for cntNormal = 1 : size( surfaceNormals, 1 )
            surfaceNormals( cntNormal, : ) = surfaceNormals( cntNormal, : ) / norm(  surfaceNormals( cntNormal, : ) );
        end
    end

    
    % Setup the frequency response of the reflectivity
    nFreqsSim = length( structMeshPreparation.vecFreqSim );
    

    FV = struct();
    FV.faces = meshClean.Faces;
    FV.vertices = meshClean.Vertices;
    hashOfMeshIn = DataHash( meshIn );
    mkdir( 'DataCalculated/CurvatureStorage/' );
    filenameCurvatureStruct = [ 'DataCalculated/CurvatureStorage/' hashOfMeshIn '.mat' ];
    
    if( structMeshPreparation.precomputeCurvature == 1 )
        if( exist( filenameCurvatureStruct, 'file' ) == 2 )
            load( filenameCurvatureStruct )
        else
            [PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]= GetCurvatures( FV ,0);
            structCurvature = struct();
            structCurvature.PrincipalCurvatures = PrincipalCurvatures;
            structCurvature.PrincipalDir1 = PrincipalDir1;
            structCurvature.PrincipalDir2 = PrincipalDir2;
            structCurvature.FaceCMatrix = FaceCMatrix;
            structCurvature.VertexCMatrix = VertexCMatrix;
            structCurvature.Cmagnitude = Cmagnitude;

            save( filenameCurvatureStruct, 'structCurvature' );
        end
    else
        [PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]= GetCurvatures( FV ,0);
        
        structCurvature = struct();
        structCurvature.PrincipalCurvatures = PrincipalCurvatures;
        structCurvature.PrincipalDir1 = PrincipalDir1;
        structCurvature.PrincipalDir2 = PrincipalDir2;
        structCurvature.FaceCMatrix = FaceCMatrix;
        structCurvature.VertexCMatrix = VertexCMatrix;
        structCurvature.Cmagnitude = Cmagnitude;
    end


    

    curvatureMagnitudeClean = sqrt( sum( structCurvature.PrincipalCurvatures.^2) );
    curvatureMagnitude = curvatureMagnitudeClean( indexesCleanToRawMesh );
    
    if( doPlot == 1 )
        figure; hist(curvatureMagnitudeClean(:),100);
    end

    surfaceBRDF = zeros( size( surfaceVertices, 1), nFreqsSim );

    % A slope of 1 feels "natural" but is way to slow in material transition. It should be 8 or 10 or so. Therefore, we add in this
    % factor.
    sigmoidSlopeMultiplier = 8;

    for cntFreq = 1 : nFreqsSim
        BRDFTransitionPosition = structMeshPreparation.BRDFTransitionPosition;
        BRDFTransitionSlope = structMeshPreparation.BRDFTransitionSlope;
        BRDFExponentSpecular = structMeshPreparation.BRDFExponentSpecular( cntFreq );
        BRDFExponentDiffractive = structMeshPreparation.BRDFExponentDiffractive( cntFreq );
        surfaceBRDF( :, cntFreq ) = mixSigmoidProperties( curvatureMagnitude, sigmoidSlopeMultiplier*BRDFTransitionSlope, BRDFTransitionPosition, BRDFExponentDiffractive, BRDFExponentSpecular );
    end

    surfaceMaterial =  zeros( size( surfaceVertices, 1), nFreqsSim );
    for cntFreq = 1 : nFreqsSim
        materialTransitionPosition =structMeshPreparation.materialSTransitionPosition;
        materialTransitionSlope = structMeshPreparation.materialSTransitionSlope;
        materialStrengthSpecular = structMeshPreparation.materialStrengthSpecular( cntFreq );
        materialStrengthDiffrative = structMeshPreparation.materialStrengthDiffractive( cntFreq );
        surfaceMaterial( :, cntFreq ) = mixSigmoidProperties( curvatureMagnitude, sigmoidSlopeMultiplier*materialTransitionSlope, materialTransitionPosition, materialStrengthDiffrative, materialStrengthSpecular );
    end

    FVPlot = struct();
    FVPlot.faces = surfaceFaces;
    FVPlot.vertices = surfaceVertices;

    if( doPlot == 1 )
        figure(132);
        clf
            subplot(2,2,1)
                colormap cool
                mesh_h=patch(FVPlot,'FaceVertexCdata',surfaceBRDF(:,1),'edgecolor', 'interp', 'facecolor','interp','EdgeAlpha',0.2);
                %set some visualization properties
                set(mesh_h,'ambientstrength',0.35);
                % axis off
                view([-45,35.2]);
                camlight();
                lighting phong
                colorbar();
                axis equal
                view( [47.0033   -4.3524]);
            subplot(2,2,2)
                colormap cool
                mesh_h=patch(FVPlot,'FaceVertexCdata',surfaceBRDF(:,14),'edgecolor', 'interp', 'facecolor','interp','EdgeAlpha',0.2);
                %set some visualization properties
                set(mesh_h,'ambientstrength',0.35);
                % axis off
                view([-45,35.2]);
                camlight();
                lighting phong
                colorbar();
                axis equal
                view( [47.0033   -4.3524]);            
           subplot(2,2,3)
                colormap cool
                mesh_h=patch(FVPlot,'FaceVertexCdata',surfaceMaterial(:,1),'edgecolor', 'interp', 'facecolor','interp','EdgeAlpha',0.2);
                %set some visualization properties
                set(mesh_h,'ambientstrength',0.35);
                % axis off
                view([-45,35.2]);
                camlight();
                lighting phong
                colorbar();
                axis equal
                view( [47.0033   -4.3524]);
            subplot(2,2,4)
                colormap cool
                mesh_h=patch(FVPlot,'FaceVertexCdata',surfaceMaterial(:,14),'edgecolor', 'interp', 'facecolor','interp','EdgeAlpha',0.2);
                %set some visualization properties
                set(mesh_h,'ambientstrength',0.35);
                % axis off
                view([-45,35.2]);
                camlight();
                lighting phong
                colorbar();
                axis equal
                view( [47.0033   -4.3524]);                   

    end
    % And now we store the surface mesh into the struct.
    structSurface = struct();
    structSurface.surfaceVertices = surfaceVertices;
    structSurface.surfaceFaces = surfaceFaces;
    structSurface.surfaceNormals = surfaceNormals;
    structSurface.surfaceMaterial = surfaceMaterial;
    structSurface.BRDF = surfaceBRDF;
    structSurface.vecFreqSim = structMeshPreparation.vecFreqSim;
    structSurface.structCurvature = structCurvature;
    structSurface.structMeshPreparation = structMeshPreparation;

    if( doPlot == 1 )
            nNormals = size( structSurface.surfaceNormals, 1 );
    normalSpacing = round( nNormals / 1000);
    idxNormalsPlot = 1 : normalSpacing : nNormals;
    figure(1212)       
        clf
        hp = patch('faces', structSurface.surfaceFaces, 'vertices', structSurface.surfaceVertices, 'FaceColor', [ 1 0 0], 'EdgeAlpha', 0.4); 
        % axis equal;
        xlabel( 'X-Axis' )
        ylabel( 'Y-Axis' )
        zlabel( 'Z-Axis' )
        grid on

        hold on
        for cntNormal = 1 : length( idxNormalsPlot )
            curIdxPlot = idxNormalsPlot( cntNormal );
            curPointOrigin = structSurface.surfaceVertices( curIdxPlot, : );
            curPointNormal = curPointOrigin + structSurface.surfaceNormals( curIdxPlot, : ) * 0.02;
            curPointsPlot = [ curPointOrigin ; curPointNormal ];
            % plot3( curPointOrigin(:,1), curPointOrigin(:,2), curPointOrigin(:,3), 'g.' );
            plot3( curPointsPlot(:,1), curPointsPlot(:,2), curPointsPlot(:,3), 'g' );
        end
        hold off;
        axis equal
        set( gca, 'view', [21.6026   19.6988]);
    end


end