
function pointsArray = generateCircularArray( distanceMicrophones, radiusArray, doPlot )

    xall=[]; yall=[];
    dy = sqrt(3)/2 * distanceMicrophones;
    ny = floor(radiusArray/dy);
    for i=-ny:ny
        y = dy*i;
        if rem(i, 2)==0
            nx = floor((sqrt(radiusArray^2 - y^2))/distanceMicrophones);
            x = (-nx:nx)'*distanceMicrophones;
        else
            nx = floor((sqrt(radiusArray^2 - y^2)-distanceMicrophones/2)/distanceMicrophones);
            x = (-nx-0.5:nx+0.5)'*distanceMicrophones;
        end
        xall = [xall; x];
        yall = [yall; y*ones(size(x))];
    end

    if( doPlot == 1 )    
        figure;
            plot(xall(:), yall(:), '.');
            hold on
                theta = 0:360;
                plot(radiusArray*cosd(theta), radiusArray*sind(theta), 'r')
            hold off
            axis equal
    end    

    pointsArray = [ xall(:), yall(:) ];