
function DCDT = Diffusion_1D_4Cmpts(t2,C,x,D_I,D_F,D_E,D_S,k_F,k_B, phi_IF,phi_FE,phi_ES, nI,nF,nE,nS, A_F, A_I)

%NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
%AND FORWARD.
dx = x(2)-x(1);
dCdt = zeros(1,length(x));
dCdx = zeros(1,length(x));
d2Cdx2 = zeros(1,length(x));

  C_intf_IFa = (A_I.*D_I.*C(nI)+(A_F*D_F.*C(nI+1)))./((A_F*D_F.*phi_IF)+A_I*D_I); %Right before IVR/epithelium interface
  C_intf_IFb = (A_I*D_I.*C(nI)+(A_F*D_F.*C(nI+1)))./((A_I*D_I./phi_IF)+A_F*D_F); %Right after IVR/epithelium interface
%C_intf_IFa = (D_I.*C(nI)+(D_F.*C(nI+1)))./((D_F.*phi_IF)+D_I); %Right before IVR/epithelium interface
%C_intf_IFb = (D_I.*C(nI)+(D_F.*C(nI+1)))./((D_I./phi_IF)+D_F); %Right after IVR/epithelium interface

 C_intf_FEa = (D_F.*C(nI+nF)+(D_E.*C(nI+nF+1)))./((D_E.*phi_FE)+D_F); %Right before epithelium/stroma interface
 C_intf_FEb = (D_F.*C(nI+nF)+(D_E.*C(nI+nF+1)))./((D_F./phi_FE)+D_E); %Right after epithelium/stroma interface

C_intf_ESa = (D_E.*C(nI+nF+nE)+(D_S.*C(nI+nF+nE+1)))./((D_S.*phi_ES)+D_E); %Right before epithelium/stroma interface
C_intf_ESb = (D_E.*C(nI+nF+nE)+(D_S.*C(nI+nF+nE+1)))./((D_E./phi_ES)+D_S); %Right after epithelium/stroma interface


for i = 1 %B.C: zero flux (dcdx=0)
    dCdt(i) = 4*D_I.*(C(i+1)-C(i))./(dx.^2);
end
for i = 2:(nI-1) %in IVR
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dx.^2);
    dCdt(i) = ((D_I/x(i)).*dCdx(i)) + D_I.*d2Cdx2(i);
   %dCdt(i) = D_I.*d2Cdx2(i);
end
for i = nI %right before interface IVR-fluid
    dCdx(i) = (C_intf_IFa - C(i))./dx;
    d2Cdx2(i) = (C(i-1)-2.*C(i)+C_intf_IFa)./(dx.^2);
    dCdt(i) = (((D_I)/x(i)).*dCdx(i)) + D_I.*d2Cdx2(i);
    %dCdt(i) = D_I.*d2Cdx2(i);

end
for i = nI+1 %right after interface IVR-fluid
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = ((C_intf_IFb)-2.*C(i)+C(i+1))./(dx.^2);
    dCdt(i) = D_F.*d2Cdx2(i) - k_F.*C(i);
end
for i = (nI+2):(nI+nF)-1 %in fluid
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dx.^2);
    dCdt(i) = D_F.*d2Cdx2(i)- k_F.*C(i);
end
for i = nI+nF % right before fluid-epithelium interface
    dCdx(i) = (C_intf_FEa - C(i))./dx;
    d2Cdx2(i) = (C(i-1)-2.*C(i)+C_intf_FEa)./(dx.^2);
    dCdt(i) = D_F.*d2Cdx2(i)- k_F.*C(i);
end
for i = nI+nF+1 % right after fluid-epithelium interface
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = ((C_intf_FEb)-2.*C(i)+C(i+1))./(dx.^2);
    dCdt(i) = D_E.*d2Cdx2(i);
end

for i = nI+nF+2:nI+nF+nE-1 % in epithelium
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dx.^2);  
    dCdt(i) = D_E.*d2Cdx2(i);
end


for i = nI+nF+nE % right before epithelium-stroma interface
    dCdx(i) = (C_intf_ESa - C(i))./dx;
    d2Cdx2(i) = (C(i-1)-2.*C(i)+C_intf_ESa)./(dx.^2);
    dCdt(i) = D_E.*d2Cdx2(i);
end
for i = nI+nF+nE+1 % right after epithelium-stroma interface
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = ((C_intf_ESb)-2.*C(i)+C(i+1))./(dx.^2);
    dCdt(i) = D_S.*d2Cdx2(i)- k_B.*C(i);
end
for i = nI+nF+nE+2:nI+nF+nE+nS-1 % in stroma
    dCdx(i) = (C(i+1)-C(i))./dx;
    d2Cdx2(i) = (C(i+1)-2.*C(i)+C(i-1))./(dx.^2);  
    dCdt(i) = D_S.*d2Cdx2(i)- k_B.*C(i);
end

for i = (nI+nF+nE+nS) % end of stroma
    %dCdt(i) = D_S.*(-2*C(i)+C(i-1))./(dx.^2) - k_B.*C(i); %BC zero conc
    dCdt(i) = 2*D_S*(C(i-1)-C(i))./(dx.^2) - k_B.*C(i); % BC zero flux
end

DCDT = dCdt';

end