clear;clc;
matrix = randi( 100, [1000 1000] );
matrix2 = matrix;
matrix3 = matrix;

rows2remove=[1 2 999 1000];
cols2remove=[1 2 999 1000];


tic
matrix2 = matrix2( 3:998, 3:998 );
toc

tic
matrix(:,cols2remove)=[];
matrix(rows2remove,:)=[];
toc