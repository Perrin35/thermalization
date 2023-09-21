OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(-0.49931061) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2274269) q[0];
sx q[0];
rz(-1.6343071) q[0];
sx q[0];
rz(-2.8660197) q[0];
x q[1];
rz(-2.2917065) q[2];
sx q[2];
rz(-1.4404785) q[2];
sx q[2];
rz(-0.40638129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5159113) q[1];
sx q[1];
rz(-1.4206919) q[1];
sx q[1];
rz(2.4019269) q[1];
rz(-pi) q[2];
rz(0.026704549) q[3];
sx q[3];
rz(-1.4674125) q[3];
sx q[3];
rz(-1.0552693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.22380655) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4734128) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(-1.0043253) q[0];
rz(-1.5197808) q[1];
sx q[1];
rz(-2.2147949) q[1];
sx q[1];
rz(1.0027592) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28264499) q[0];
sx q[0];
rz(-2.1574321) q[0];
sx q[0];
rz(-1.0413917) q[0];
rz(-pi) q[1];
rz(2.9687256) q[2];
sx q[2];
rz(-1.3860518) q[2];
sx q[2];
rz(-2.6181521) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1042852) q[1];
sx q[1];
rz(-2.1530495) q[1];
sx q[1];
rz(-2.9557455) q[1];
rz(-pi) q[2];
rz(0.59870436) q[3];
sx q[3];
rz(-1.2748534) q[3];
sx q[3];
rz(-0.89950365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(-2.7271467) q[3];
sx q[3];
rz(-0.60025418) q[3];
sx q[3];
rz(-3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8931483) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-2.3625968) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(0.88358203) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1012264) q[0];
sx q[0];
rz(-1.7352312) q[0];
sx q[0];
rz(2.7357487) q[0];
rz(0.29404624) q[2];
sx q[2];
rz(-1.9532734) q[2];
sx q[2];
rz(1.2750212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.59144943) q[1];
sx q[1];
rz(-1.5482229) q[1];
sx q[1];
rz(-1.8557465) q[1];
rz(0.65914764) q[3];
sx q[3];
rz(-1.4775606) q[3];
sx q[3];
rz(-2.6574709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8399923) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(2.7123614) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-2.0638549) q[3];
sx q[3];
rz(-2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(-2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(0.5501737) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4470091) q[0];
sx q[0];
rz(-0.83568474) q[0];
sx q[0];
rz(-1.4674594) q[0];
rz(-pi) q[1];
rz(0.29454622) q[2];
sx q[2];
rz(-1.4809594) q[2];
sx q[2];
rz(1.3603269) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95851129) q[1];
sx q[1];
rz(-1.4554879) q[1];
sx q[1];
rz(2.1469841) q[1];
rz(-pi) q[2];
rz(0.22927852) q[3];
sx q[3];
rz(-1.6212654) q[3];
sx q[3];
rz(-0.63427395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.52175534) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(-0.27553976) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(-0.47079852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(2.450401) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(2.2263288) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4782151) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(-0.94232725) q[0];
x q[1];
rz(1.615633) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(0.63477883) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3770959) q[1];
sx q[1];
rz(-0.15863523) q[1];
sx q[1];
rz(-1.7736969) q[1];
rz(-0.74113412) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(-2.4853064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.22121945) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(0.4635703) q[2];
rz(2.5772337) q[3];
sx q[3];
rz(-2.1488991) q[3];
sx q[3];
rz(0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1148465) q[0];
sx q[0];
rz(-2.5214054) q[0];
sx q[0];
rz(2.0625431) q[0];
rz(-0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(-2.7450096) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81178938) q[0];
sx q[0];
rz(-1.3178696) q[0];
sx q[0];
rz(-2.0407709) q[0];
x q[1];
rz(-0.47607143) q[2];
sx q[2];
rz(-1.8487612) q[2];
sx q[2];
rz(-1.4897886) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.086416883) q[1];
sx q[1];
rz(-0.91311087) q[1];
sx q[1];
rz(-0.1187101) q[1];
rz(-pi) q[2];
rz(-0.58432213) q[3];
sx q[3];
rz(-1.366368) q[3];
sx q[3];
rz(1.0014597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98809272) q[2];
sx q[2];
rz(-0.92553878) q[2];
sx q[2];
rz(-1.5863824) q[2];
rz(-1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(2.356785) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.3762459) q[1];
sx q[1];
rz(-2.0152337) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61486812) q[0];
sx q[0];
rz(-1.65979) q[0];
sx q[0];
rz(2.8429549) q[0];
rz(-pi) q[1];
rz(0.14234219) q[2];
sx q[2];
rz(-1.0668313) q[2];
sx q[2];
rz(0.0056643639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3583402) q[1];
sx q[1];
rz(-2.0166774) q[1];
sx q[1];
rz(-1.7794442) q[1];
x q[2];
rz(-3.0464826) q[3];
sx q[3];
rz(-1.9169807) q[3];
sx q[3];
rz(2.9444875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9324947) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-2.9879925) q[2];
rz(-2.8364654) q[3];
sx q[3];
rz(-1.1161476) q[3];
sx q[3];
rz(-1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-3.0287108) q[0];
rz(-2.1408634) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(0.032756068) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8772802) q[0];
sx q[0];
rz(-0.41752975) q[0];
sx q[0];
rz(0.89950048) q[0];
x q[1];
rz(2.6997386) q[2];
sx q[2];
rz(-1.0519069) q[2];
sx q[2];
rz(-1.9372802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0667463) q[1];
sx q[1];
rz(-0.42145887) q[1];
sx q[1];
rz(1.4261817) q[1];
x q[2];
rz(1.350358) q[3];
sx q[3];
rz(-2.7535097) q[3];
sx q[3];
rz(-0.9006075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1774566) q[2];
sx q[2];
rz(-2.5033341) q[2];
sx q[2];
rz(-0.62210554) q[2];
rz(-1.1671676) q[3];
sx q[3];
rz(-2.0303576) q[3];
sx q[3];
rz(-2.7511403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0416097) q[0];
sx q[0];
rz(-2.6384625) q[0];
sx q[0];
rz(1.6149678) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(2.656235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.582726) q[0];
sx q[0];
rz(-0.14478806) q[0];
sx q[0];
rz(-2.7607714) q[0];
x q[1];
rz(-0.59818563) q[2];
sx q[2];
rz(-1.7028114) q[2];
sx q[2];
rz(-0.26192947) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.34708729) q[1];
sx q[1];
rz(-0.80650389) q[1];
sx q[1];
rz(1.1361213) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5076809) q[3];
sx q[3];
rz(-2.4571107) q[3];
sx q[3];
rz(2.4661635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1099403) q[2];
sx q[2];
rz(-1.3039219) q[2];
sx q[2];
rz(-2.9821441) q[2];
rz(-1.4032646) q[3];
sx q[3];
rz(-2.7519029) q[3];
sx q[3];
rz(0.015550912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(-1.4341226) q[0];
rz(1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(2.5433345) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82829581) q[0];
sx q[0];
rz(-1.3643571) q[0];
sx q[0];
rz(2.9199827) q[0];
rz(-pi) q[1];
rz(-1.7190476) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(1.460618) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.977539) q[1];
sx q[1];
rz(-2.7643449) q[1];
sx q[1];
rz(-1.3110728) q[1];
x q[2];
rz(-0.80501276) q[3];
sx q[3];
rz(-2.9152438) q[3];
sx q[3];
rz(-1.5387907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(2.5478798) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(-1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3020637) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.2364173) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-1.6937704) q[2];
sx q[2];
rz(-2.1868144) q[2];
sx q[2];
rz(-1.1745324) q[2];
rz(-0.39246172) q[3];
sx q[3];
rz(-0.93467181) q[3];
sx q[3];
rz(-1.296464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
