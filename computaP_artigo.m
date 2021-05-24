function q = computaP_artigo(a,p0,kE,kI,phoI,phoE,cI,cE,lmax)
%%%% plota campo de pressão conforme figura 4 do artigo Matheus/prof Antonio  "Acoustic scattering and forces on a arbitrrily sized fluid sphere"

Pinc = @(z) p0*exp(i*kE*z); %onda incidente na direção z
Pint_ant = @(r,theta,phi) 0; %valor de inicialização
Psca_ant = @(r,theta,phi) 0; %valor de inicialização


for l = 1:lmax
  for m = -l:+1

    Y = @(theta,phi) (theta>=0).*(theta<=pi).*(phi>=0).*(phi<=2*pi).*( sph_harmonic(l, m, theta, phi) );

    j = @(l,x) sqrt(pi/2)*besselj(l+0.5,x)./sqrt(x);
    y = @(l,x) sqrt(pi/2)*bessely(l+0.5,x)./sqrt(x);
    h = @(l,x) j(l,x) + i*y(l,x);

    g = @(r,theta,phi) Pinc(r.*cos(theta)).*Y(theta,phi)'.*sin(theta);
    %G = @(l,m,r) integral2(g,0,r,0,pi,0,2*pi)./(p0*j(l,kE*r));
    G = @(l,m) ( sqrt(4*pi*(2*l+1))*(i^l)*delta(m,0) ); 
    
    alfa = phoI*cI/(phoE*cE);

    lambE = 2*pi/kE; lambI = 2*pi/kI;
    dr = min([lambE lambI a])/20;
    dx = kE*dr;

    hlin = @(l,x) ( h(l,x+dx) - h(l,x) )./dx ;
    jlin = @(l,x) ( j(l,x+dx) - j(l,x) )./dx ;

    xE = kE*a;
    xI = kI*a;
    A = @(l,m) ( ( j(l,xE).*hlin(l,xE) - jlin(l,xE).*h(l,xE) )./( j(l,xI).*hlin(l,xI) - (1/alfa)*jlin(l,xI).*h(l,xI) ) ).*G(l,m);   %,a) ;
    C = @(l,m) ( ( alfa*j(l,xE).*hlin(l,xE) - jlin(l,xE).*h(l,xE) )./( j(l,xI).*hlin(l,xI) - alfa*jlin(l,xI).*h(l,xI) ) ).*G(l,m);   %,a) ;

    Pint = @(r,theta,phi) Pint_ant + p0*A(l,m).*j(l,kI*r).*Y(theta,phi);
    Pint_ant =  @(r,theta,phi) Pint(r,theta,phi);
    
    Psca = @(r,theta,phi) Psca_ant(r,theta,phi) + p0*C(l,m).*h(l,kI*r).*Y(theta,phi);
    Psca_ant = @(r,theta,phi) Psca(r,theta,phi);  
    
  end
end

k = max([kE kI]);
xmin = -4*k*a; xmax = 8*k*a; 
ymin = - 3*k*a; ymax = 3*k*a;
dx = k*dr; dy = dx;

x = xmin:dx:xmax; y = ymin:dy:ymax;

[X,Y] = meshgrid(x,y);
R = sqrt(X.^2 + Y.^2); THETA = zeros(size(R));; PHI = atan(Y./X);

P = Pinc(R) + Pint(R,THETA,PHI) + Psca(R,THETA,PHI);

surf(X,Y,P); colorbar

q = 1;

end

 


