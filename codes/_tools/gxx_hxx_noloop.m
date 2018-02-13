%GXX_HXX_NOLOOP.M
%[gxx,hxx] = gxx_hxx_noloop(fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,hx,gx) 
%finds the 3-dimensional arrays gxx and hxx necessary to compute the 2nd order approximation 
%to the decision rules of a DSGE model of the form E_tf(yp,y,xp,x)=0, with solution 
%xp=h(x,sigma) + sigma * eta * ep and y=g(x,sigma). For more details, see  
%``Solving Dynamic General Equilibrium Models Using a Second-Order Approximation to the Policy 
%Function,'' by Stephanie Schmitt-Grohe and Martin Uribe, JEDC, January 2004, p. 755-775. 
%
%INPUTS: First and second derivatives of f and first-order approximation to the functions g and 
%h: fx, fxp, fy, fyp, fypyp, fypy, fypxp, fypx, fyyp, fyy, fyxp, fyx, fxpyp, fxpy, fxpxp, fxpx, 
%fxyp, fxy, fxxp, fxx, hx, gx
%
%OUTPUTS: Second-order derivatives of the functions g and h with respect to x, evaluated 
%at (x,sigma)=(xbar,0), where xbar=h(xbar,0). That is, hxx gxx
%
%gxx_hxx_noloop.m is faster than gxx_hxx.m, but runs more often into memory problems 
%
%We follow Paul Klein (2005) in using the Hessian definition of Magnus and Neuberger to 
%find gxx and hxx
%
%(c) Stephanie Schmitt-Grohe and Martin Uribe
%
%Date July 11, 2005
function [gxx,hxx] = gxx_hxx_noloop(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx)
nfypyp=hessian(nfypyp);
nfypy=hessian(nfypy);
nfypxp=hessian(nfypxp);
nfypx=hessian(nfypx);

nfyyp=hessian(nfyyp);
nfyy =hessian(nfyy);
nfyxp=hessian(nfyxp);
nfyx =hessian(nfyx);

nfxpyp=hessian(nfxpyp);
nfxpy=hessian(nfxpy);
nfxpxp=hessian(nfxpxp);
nfxpx=hessian(nfxpx);

nfxyp=hessian(nfxyp);
nfxy=hessian(nfxy);
nfxxp=hessian(nfxxp);
nfxx=hessian(nfxx);


nx=size(nfx,2);
ny = size(nfy,2);
n = nx + ny;
 
A = kron(eye(n),hx'*gx')*(nfypyp*gx*hx + nfypxp*hx + nfypy*gx + nfypx) + ...
    kron(eye(n),gx')    *(nfyyp *gx*hx + nfyxp *hx + nfyy *gx + nfyx ) + ...    
    kron(eye(n),hx')    *(nfxpyp*gx*hx + nfxpxp*hx + nfxpy*gx + nfxpx) + ...
                         (nfxyp *gx*hx + nfxxp *hx + nfxy *gx + nfxx);   
B = kron(nfyp, hx') ; 
C = kron(nfy, eye(nx));
D = kron(nfyp*gx, eye(nx))+ kron(nfxp, eye(nx)); 

Qq = -[ kron(hx',B)+kron(eye(nx),C), kron(eye(nx), D)] \ (A(:));

gxx = permute(reshape(Qq(1:nx^2*ny),nx,ny,nx),[2 1 3]); %this reshaping is done so that we get the output in the same order as we used to
hxx = permute(reshape(Qq(nx^2*ny+1:end),nx,nx,nx),[2 1 3]);%reshape to get the same as in the previous version version of gxx_hxx.m

function x = hessian(x)
[n1,n2,n3]=size(x);
x=reshape(permute(x,[ 2 1 3]), n1*n2, n3);