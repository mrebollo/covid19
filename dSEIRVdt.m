function dYdt = dSEIRVdt (SEIRV)
% function dYdt = dSEIRVdt (t,SEIRV) % modelo SEIRV que usa constantes globales
%%
global Lambda betaE betaI betaV mu alfa w gamma xi1 xi2 sigma;
S=SEIRV(1); E=SEIRV(2); I=SEIRV(3); R=SEIRV(4); V=SEIRV(5);
%%
dSdt = Lambda-betaE(E)*S*E-betaI(I)*S*I-betaV(V)*S*V-mu*S;
dEdt = betaE(E)*S*E+betaI(I)*S*I+betaV(V)*S*V-(alfa+mu)*E;
dIdt = alfa*E-(w+gamma+mu)*I;
dRdt = gamma*I-mu*R;
dVdt = xi1*E+xi2*I-sigma*V;
%%
dYdt = [dSdt ; dEdt; dIdt; dRdt; dVdt ]
end