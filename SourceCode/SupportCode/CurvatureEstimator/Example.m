%% This script is an example of how to implement the functions for estimating curvature and their derivatives 
%Author: Itzik Ben Shabat
%Last Update: July 2014

%%
%Clear all variables, close all windows and clear command window 
clear all
close all
clc
%% Generate the example triangle mesh example
    [ surfaceFaces, surfaceVertices ] = stlread("../..//Models/leaf2_dragonfly2.stl");
    FV = struct();
    FV.faces = surfaceFaces;
    FV.vertices = surfaceVertices*1000;

%% calcualte curvatures


mesh = surfaceMesh(surfaceVertices,surfaceFaces);
removeDefects(mesh,'nonmanifold-edges');
removeDefects(mesh,'duplicate-faces');
removeDefects(mesh,'duplicate-vertices');

    FV = struct();
    FV.faces = mesh.Faces;
    FV.vertices = mesh.Vertices;

getderivatives=0;
[PrincipalCurvatures,PrincipalDir1,PrincipalDir2,FaceCMatrix,VertexCMatrix,Cmagnitude]= GetCurvatures( FV ,1);
% 
PrincipalCurvatures(PrincipalCurvatures>10 ) = 10;
PrincipalCurvatures(PrincipalCurvatures<-10 ) = -10;


GausianCurvature=PrincipalCurvatures(1,:).*PrincipalCurvatures(2,:);
%% Draw the mesh to the screen 
figure('name','Triangle Mesh Curvature Example','numbertitle','off','color','w');
colormap cool

threshCurvatureDiffraction = 20;
% diffractionVertices =  scaled_shifted_sigmoid(Cmagnitude, 0.1, 10 );

diffractionVertices =   mixSigmoidProperties( Cmagnitude, 0.1, 10, 10, 50 );
% diffractionVertices( Cmagnitude > threshCurvatureDiffraction ) = 100;
% diffractionVertices( Cmagnitude <= threshCurvatureDiffraction ) = 0;
mesh_h=patch(FV,'FaceVertexCdata',diffractionVertices,'facecolor','interp','edgecolor','interp','EdgeAlpha',0.2);


%set some visualization properties
set(mesh_h,'ambientstrength',0.35);
axis off
view([-45,35.2]);
camlight();
lighting phong
colorbar();
axis equal

%%

%%
function y = scaled_shifted_sigmoid(x, a, b)
    y = 1 ./ (1 + exp(-a * (x - b)));
end

