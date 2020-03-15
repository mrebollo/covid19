% componentes del vector y = SEIRV
%S=1;E=2;I=3;R=4;V=5;
%dSEIRVdt = @(t,y) [
% dS/dt =
%Lambda-betaE(y(E))*y(S)*y(E)-betaI(y(I))*y(S)*y(I)-betaV(y(V))*y(S)*y(V)-mu*y(S)
% dE/dt =
%betaE(y(E))*y(S)*y(E)+betaI(y(I))*y(S)*y(I)+betaV(y(V))*y(S)*y(V)-(alfa+mu)*y(E)
% dI/dt =
%alfa*y(E)-(w+gamma+mu)*y(I)
% dR/dt =
%gamma*y(I)-mu*y(R)
% dV/dt =
%xi1*y(E)+xi2*y(I)-sigma*y(V)
%];

% parametros del modelo SEIRV http://doi.org/10.3934/mbe.2020148
Lambda = 271.23 ;% flujo de poblacion hacia Wuhan (por dia)
betaE0 = 3.11e-8;% tasa de transmision entre S y E (/personas/dia)
betaI0 = 0.62e-8;% tasa de transmision entre S e I (/personas/dia)
betaV0 = 1.03e-8;% tasa de transmision entre S y V (/personas/dia)
c = 1.01e-4;% coeficiente de ajuste de la transmision
mu = 3.01e-5;% tasa natural de fallecimientos (por día)
alfa = 1/7 ;% 1/alfa = periodo de incubacion (/dias)
w = 0.01 ;% tasa de fallceimientos por la infeccion (por dia)
gamma = 1/15 ;% tasa de recuperacion de la infeccion (por dia)
sigma = 1 ;% tasa de eliminacion del virus del ambiente (por dia)
xi1 = 2.30 ;% dispersion del virus por infectados asintomaticos
xi2 = 0 ;% y por infectados aislados (por persona por dia por ml)

% funciones de contagio (positivas y decrecientes con su argumento)
betaE = @(E) betaE0/(1+c*E); % para expuestos
betaI = @(I) betaI0/(1+c*I); % para infectados
betaV = @(V) betaV0/(1+c*V); % para el ambiente

% condicion inicial
S0 = 8998505 ;% susceptibles de infectar
E0 = 1000 ;% expuestos (infectados asintomaticos)
I0 = 475 ;% infectados (con sintomas y aislados)
R0 = 10 ;% recuperados
V0 = 10000 ;% concentracion de virus en el ambiente


% resolucion numerica
SEIRV0 = [S0; E0; I0; R0; V0];
[t,SEIRV] = ode23t(dSEIRVdt,[0,3000],SEIRV0);



% dibuja el plano de fase
I0v = [ 5 10 15 15 15 15 15 10 6 3 0 0 0 0 ]*1e4;
E0v = [ 0 0 0 1 2 3 4 4 4 4 3.5 2.5 1.5 0.5 ]*1e4;
clf; box on;
for iv=1:length(I0v)
SEIRV0 = [S0; E0v(iv); I0v(iv); R0; V0];
[t,SEIRV] = ode23t(dSEIRVdt,[0,50000],SEIRV0);
subplot(2,2,1); plot(SEIRV(:,E),SEIRV(:,I),'-'); hold on;
subplot(2,2,2); semilogx(t,SEIRV(:,E),t,SEIRV(:,I)); hold on;
end;
subplot(2,2,1); xlabel('E'); ylabel('I'); title('phase portrait');
subplot(2,2,2); xlabel('t'); title('E, I'); axis([0.1 5e4 0 1.5e5])

% calcula punto de equilibrio
SEIRVstar = fsolve(@(y) dSEIRVdt(0,y), SEIRV0)
SEIRVstar([E, I])







% probar toda la odesuite
%odesuite = {"ode23","ode45","ode15s","ode23s","ode23t","ode23tb"};
%for ii=1:length(odesuite)
%time = cputime;
%eval(strcat('[t,SEIRV]=',odesuite{ii},'(dSEIRVdt,[0,30000],SEIRV0);'));
%cost(ii) = cputime-time,
%end


