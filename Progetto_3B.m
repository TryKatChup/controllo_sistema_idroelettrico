%Progetto Tipologia 3B
%%%%%%%%
%   Sistema controllo centrale idroelettrica
%       Bernardi Daniel
%       Chichifoi Karina
%       Ivan Andrei Daniel
%       Pizzini Cavagna Hiari
%%%%%%%%

%% Definizione sistema
%x_dot=f(x,u)

%x_dot_1=0;
%x_dot_2=-C_d*u*x_2*abs(x_2) -R_0*x_2*abs(x_2) +x_1;

%y=-eta*x_1*x_2

%%Specifiche del progetto:
    % errore a regime nullo con riferimento a gradino w(t)=1(t)*W
    % Per la robustezza, deve garantire 
    %   - Mf>=45°
    % Sovraelongazione accettabile 
    %   - S%<=5%
    % Tempo di assestamento al 5%
    %   - T_a5<=0.15 [s]
    % Il rumore deve essere abbattuto di almeno B_n=30 volte, 
    %   caratteristiche rumore: da omega_n>1500 [rad/s], ampiezza A_n=0.05;

%y_ref uscita di riferimento con notazione W
W=40;
A_n=0.05;
y_ref=40;

% Def. parametri del sistema
C_d=(2.2*pi^2);
R_0=35;
eta=0.65;
omega_n=1500;
%x_1, x_2 punti di equilibrio
x_1=10;
x_2=5;
y_notlin=-eta*x_1*x_2;
u_=((x_1/(x_2*abs(x_2))- R_0) / C_d); u_
%Matrici punto di equilibrio (x_1,x_2)=(10,5);
A= [0   ,   0;
    1   ,   (-C_d*u_-R_0)*abs(x_2)*2];

B=[0;
   -C_d*x_2*abs(x_2)];

C=[-eta*x_2 , -eta*x_1];

D=0;

%% Funzione di trasferimento
s=tf('s');
[N,D]=ss2tf(A,B,C,D);
G=tf(N,D);
zpk(G)
%Definizione dell'intervallo di frequenze del diagramma di Bode
omega_plot_min=10^(-2);
omega_plot_max=10^5;

%Il diagramma di Bode presenta anche le limitazioni per il disturbo
%di misura (zona gialla)
figure();
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-29.9,-29.9,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
hold on;
[Mag,phase,w]=bode(G,{omega_plot_min,omega_plot_max});
margin(Mag,phase,w);
grid on;
hold off;
title("Funzione di trasferimento iniziale");
%Come si nota G(s) attraversa la zona proibita 

%Check del sistema in closed loop senza nessun controllo
% F è la funzione di sensitività complementare
% F = R(s)*G(S)/(1 + R(s)*G(s)), oppure
% F = L(s)/(1 + L(s)), ovvero L(s) in retroazione
%Qui riportata senza R(s), ancora da elaborare
figure();
hold on;
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-29.9,-29.9,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
margin(Mag,phase,w);
grid on;

%Definizione di F
F=G/(1+G);
margin(F);
legend('G','F')
title("Funzione di trasferimento ad anello aperto e chiuso");
hold off;

figure();
step(F);

%% Shaping L(s)
%L(s) funzione di trasferimento ad anello aperto
%L(s)=G(s)*R(s)
%R_s(s) REGOLATORE STATICO

%% Prestazioni statiche
%Regolazione dell'errore a regime a 0
%Utilizzo semplice regolatore statico con un polo e guadagno statico libero
%Il polo è necessario poichè la G(s) di partenza NON ha un polo
%nell'origine, così facendo non dovremo preoccuparci del valore di W
%Da pensare anche al guadagno statico mu_s, qui = 1, ma utile per ulteriori
%vincoli
mu_s=1;
R_s=mu_s/s;

%Otteniamo nuova funzione di trasferimento, in serie con il regolatore 
G_=R_s*G;

figure();
hold on;
[Mag,phase,w]=bode(G_,{omega_plot_min,omega_plot_max});
%Zona proibita
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[- 29.9,-29.9,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);

margin(Mag,phase,w);
margin(G);
grid on;
legend('G(s)*R_s(s)','G');
title("Funzione trasferimento con regolatore statico");
hold off;
%Frequenza di taglio a 59.3
%Il polo aggiunto risolve il problema dell'errore a regime e pone la G
%al di fuori della zona gialla

%Check del sistema in closed loop con l'aggiunta del controllo
F=G_/(1+G_);
zpk(F)

figure();
step(F);
%Si hanno delle osciallazioni a causa dei polo c.c.
%I tempi di assestamento sono ancora alti, intorno ai 3 secondi 
% -->0.15 obiettivo, prossima sezione

figure();
margin(G_);
hold on;
margin(F);
legend('G(s)*R_s(s)','F(s)');
title("Funzione di trasferimento con regolatore ad anello aperto e chiuso");
%Il fatto che il diagramma di Bode di F sia a 0 nelle frequenze basse va benissimo
%in quanto è l'intervallo che va a incidere sul riferimento W, non
%modificandolo
%Il picco di F(s) è detto PICCO DI RISONANZA
%Intorno a quelle frequenze, l'ingresso viene enormemente amplificato 
%Dopo omega_n si abbassa per attenuare il disturbo n ad alte frequenze

%NB! Abbiamo già un unica frequenza di attraversamento, il che significa
%che rientriamo nei criteri di Bode

%% Prestazioni dinamiche
%Regolazione sovraelongazione, tempo assestamento 
%I due valori sono legati alla frequenza di taglio t.c.
% omega_c >= 300 /(T* * M_f)
%Il valore 300 corrisponde ad una costante del T di assestamento al 5%
%T* è il nostro obbiettivo 0.15, M_f a 45°
%La frequenza di taglio attuale di G_ è 59.3 (ottenibile da grafico precendete)
%Calcolo xi = M_f * 100, girando la formula S% = 100* exp(-pi*xi(sqrt(1-xi^2),
%considerando S% <= 5%;
xi=sqrt(log(0.05)^2/(pi^2+log(0.05)^2));
%xi = 0.69
M_f=xi*100;
%M_f = 69°, molto più alto della consegna 45°, ma terremo questo valore per
%garantire le migliori performance ottenibili
%La frequenza di taglio minima omega_c_min sarà quindi:
omega_c_min=300/(0.15 * M_f);
%omega_c_min = 28.98 
omega_c_max=100;

figure();
hold on;
[Mag,phase,w]=bode(G_,{omega_plot_min,omega_plot_max});
patch([omega_n,omega_plot_max,omega_plot_max,omega_n],[-29.9,-29.9,100,100],'y','FaceAlpha',0.3,'EdgeAlpha',0);
legend('Zona disturbo n');
margin(Mag,phase,w);
patch([omega_c_min,omega_n,omega_n,omega_c_min],[-180+M_f,-180+M_f,-270,-270],'g','FaceAlpha',0.3,'EdgeAlpha',0); 
legend('G(s)*R_s(s)','X Margine fase');

%Determino il valore del guadagno del REGOLATORE DINAMICO R_d(s)
%mu_d è il guadagno del regolatore dinamico, ottenuto grazie alla funz
%evalfr -> evaluates frequency che determina il modulo di G_ nella
%omega_c_min. 
%NB! 1i * omega_c_min + jomega ((jw))
%poichè evalfr ha argomenti (sistema,frequency) 

omega_c_star=omega_c_min;
M_star=abs(evalfr(G_,1i*omega_c_star));
phi_star=(M_f)-180-angle(evalfr(G_,1i*omega_c_star))*180/pi;
tau=(M_star-cos(phi_star*pi/180))/(omega_c_star*sin(phi_star*pi/180));
alpha=(cos(phi_star*pi/180)-1/M_star)/(omega_c_star*sin(phi_star*pi/180));

R_d=(1+tau*s)/(1+alpha*tau*s);
%Il regolatore dinamico così progettato non rispetta i requisiti di 
%riduzione del disturbo n.
%Il guadagno statico non è bloccato da R_s, quindi posso sfruttarlo
%per ridurre l'ampiezza per uscire dalla zona gialla
%mu_d=0.30;
mu_d_val=mag2db(abs(evalfr(G_,1i*omega_c_min)));
mu_d=db2mag(-mu_d_val); 
R_d=mu_d*R_d;

%Posso ora ricavare il regolatore completo, moltiplicando lo statico per il
%dinamico
R=R_d*R_s;
[NR,DR] = tfdata(R);

%Ottengo ora L(s), data dalla dall'applicazione del regolatore finale a G
%iniziale
L=R*G;
F=L/(1+L);

margin(F);
margin(L);
%patch([omega_c_min,omega_c_max,omega_c_max,omega_c_min],[-180+M_f,-180+M_f,-270,-270],'g','FaceAlpha',0.2,'EdgeAlpha',0); 
legend('G(s)*R_s(s)','X Margine fase', 'F', 'L');
hold off;
%Il fatto che |F(s)| sia = 0       per omega << omega_c 
%                        = |L(s)|  per omega >> omega_c 
%conferma che sia giusto

figure();
step(F);
%I requisiti di sovraelongazione e tempo di assestamento al 5% sono 
%soddisfatti sia per T_a5(0.15) che per T_a0(0.04)
%Il margine di fase è rispettato poichè nella nuova frequenza di taglio
%omega_c non si è nella zona proibita

%Ampiezza L in omega_n < -30 -> requisito disturbo n rispettato
mag2db(abs(evalfr(L,omega_n)))

return;
