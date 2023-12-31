OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3061476) q[0];
sx q[0];
rz(-2.4581576) q[0];
sx q[0];
rz(-0.47877065) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(5.1217084) q[1];
sx q[1];
rz(6.9245467) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35533479) q[0];
sx q[0];
rz(-1.1874275) q[0];
sx q[0];
rz(-1.8281561) q[0];
x q[1];
rz(0.85513656) q[2];
sx q[2];
rz(-2.5980066) q[2];
sx q[2];
rz(0.91993514) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8909059) q[1];
sx q[1];
rz(-2.2911934) q[1];
sx q[1];
rz(2.5198063) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5121243) q[3];
sx q[3];
rz(-1.2347722) q[3];
sx q[3];
rz(2.0365086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.59149867) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(-2.6634898) q[2];
rz(-1.452662) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(12/(7*pi)) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96145445) q[0];
sx q[0];
rz(-0.77471662) q[0];
sx q[0];
rz(2.1226728) q[0];
rz(1.6628751) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(0.63308024) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79257747) q[0];
sx q[0];
rz(-1.4903869) q[0];
sx q[0];
rz(0.76064674) q[0];
rz(-pi) q[1];
x q[1];
rz(0.622153) q[2];
sx q[2];
rz(-0.85217798) q[2];
sx q[2];
rz(-2.6785786) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.1795579) q[1];
sx q[1];
rz(-1.1711367) q[1];
sx q[1];
rz(2.7707997) q[1];
rz(-0.88926104) q[3];
sx q[3];
rz(-1.325843) q[3];
sx q[3];
rz(-2.3219061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7818266) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(0.11745545) q[2];
rz(0.30101267) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(-1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2837219) q[0];
sx q[0];
rz(-0.71325934) q[0];
sx q[0];
rz(0.088407956) q[0];
rz(1.1075426) q[1];
sx q[1];
rz(-0.91845599) q[1];
sx q[1];
rz(-2.450768) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72969499) q[0];
sx q[0];
rz(-1.6459961) q[0];
sx q[0];
rz(-2.9883283) q[0];
x q[1];
rz(-0.76491852) q[2];
sx q[2];
rz(-1.3996482) q[2];
sx q[2];
rz(-0.82676065) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1580116) q[1];
sx q[1];
rz(-2.9576655) q[1];
sx q[1];
rz(1.7137101) q[1];
rz(-pi) q[2];
rz(1.2191804) q[3];
sx q[3];
rz(-1.9295441) q[3];
sx q[3];
rz(0.84695942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2174125) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(-1.4364093) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(-1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0687662) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(2.6191214) q[0];
rz(-0.32896313) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(2.5879588) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41933665) q[0];
sx q[0];
rz(-2.0291078) q[0];
sx q[0];
rz(2.3769747) q[0];
rz(-pi) q[1];
rz(-2.3217818) q[2];
sx q[2];
rz(-0.95633436) q[2];
sx q[2];
rz(-2.7053506) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.5987293) q[1];
sx q[1];
rz(-2.3485564) q[1];
sx q[1];
rz(1.7584156) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5418566) q[3];
sx q[3];
rz(-0.77416285) q[3];
sx q[3];
rz(-0.72284568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5840977) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(-2.6468357) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5938063) q[3];
sx q[3];
rz(-1.2899227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49098) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(-0.96013534) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(2.9575612) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83443123) q[0];
sx q[0];
rz(-2.8330028) q[0];
sx q[0];
rz(0.35085268) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8475624) q[2];
sx q[2];
rz(-1.9230611) q[2];
sx q[2];
rz(2.2500452) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.80464333) q[1];
sx q[1];
rz(-0.78362432) q[1];
sx q[1];
rz(-0.6964535) q[1];
rz(-pi) q[2];
rz(-1.3695413) q[3];
sx q[3];
rz(-2.2866837) q[3];
sx q[3];
rz(-1.9191238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.90157834) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(-2.0007755) q[2];
rz(2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(2.0146577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7891156) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(-2.254803) q[0];
rz(2.9011762) q[1];
sx q[1];
rz(-1.0188894) q[1];
sx q[1];
rz(-0.14850798) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37489031) q[0];
sx q[0];
rz(-2.3935211) q[0];
sx q[0];
rz(-0.23776777) q[0];
rz(-pi) q[1];
x q[1];
rz(2.912942) q[2];
sx q[2];
rz(-1.6778523) q[2];
sx q[2];
rz(0.46496898) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.099523274) q[1];
sx q[1];
rz(-1.6798786) q[1];
sx q[1];
rz(-1.3101577) q[1];
rz(-pi) q[2];
rz(2.4466483) q[3];
sx q[3];
rz(-0.57999014) q[3];
sx q[3];
rz(-2.490173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.85577661) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(-0.4883858) q[2];
rz(-2.6521902) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(1.874812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4483036) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(2.3235902) q[0];
rz(-1.2524293) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(-2.0163527) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97911191) q[0];
sx q[0];
rz(-1.1971803) q[0];
sx q[0];
rz(-3.0309249) q[0];
rz(-pi) q[1];
rz(-0.15525012) q[2];
sx q[2];
rz(-1.2755738) q[2];
sx q[2];
rz(-2.070602) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.352467) q[1];
sx q[1];
rz(-1.8450292) q[1];
sx q[1];
rz(-1.2624361) q[1];
rz(-pi) q[2];
rz(-0.66594395) q[3];
sx q[3];
rz(-2.0332554) q[3];
sx q[3];
rz(2.8921814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0730878) q[2];
sx q[2];
rz(-0.98098522) q[2];
sx q[2];
rz(-0.64458624) q[2];
rz(1.6623496) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(-1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709764) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(1.6059426) q[0];
rz(1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(-2.1309526) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3686417) q[0];
sx q[0];
rz(-0.7612345) q[0];
sx q[0];
rz(1.3377454) q[0];
rz(-0.29823093) q[2];
sx q[2];
rz(-1.094162) q[2];
sx q[2];
rz(2.6963866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1809363) q[1];
sx q[1];
rz(-1.8005383) q[1];
sx q[1];
rz(0.22682637) q[1];
rz(-pi) q[2];
rz(-1.478785) q[3];
sx q[3];
rz(-1.9698471) q[3];
sx q[3];
rz(-1.9012746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(2.4475205) q[2];
rz(2.5726035) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(-0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086031832) q[0];
sx q[0];
rz(-1.1431575) q[0];
sx q[0];
rz(0.49474299) q[0];
rz(-0.61839473) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(0.075597413) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6583017) q[0];
sx q[0];
rz(-2.1245983) q[0];
sx q[0];
rz(-2.7730586) q[0];
rz(-pi) q[1];
rz(-1.7463023) q[2];
sx q[2];
rz(-1.5117466) q[2];
sx q[2];
rz(0.67948558) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1739993) q[1];
sx q[1];
rz(-0.6951957) q[1];
sx q[1];
rz(-2.980152) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5209386) q[3];
sx q[3];
rz(-1.7194347) q[3];
sx q[3];
rz(1.8295446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8081234) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(-1.3405651) q[2];
rz(0.30424413) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(-2.5568967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8463523) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(-2.9647968) q[0];
rz(-1.2416174) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(0.44100824) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.037017578) q[0];
sx q[0];
rz(-1.3694166) q[0];
sx q[0];
rz(-1.3637533) q[0];
x q[1];
rz(0.44956019) q[2];
sx q[2];
rz(-2.9712147) q[2];
sx q[2];
rz(2.5296488) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8969352) q[1];
sx q[1];
rz(-1.1420297) q[1];
sx q[1];
rz(-2.7662617) q[1];
rz(-pi) q[2];
rz(-0.62426626) q[3];
sx q[3];
rz(-1.3253951) q[3];
sx q[3];
rz(-2.0516968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.56197721) q[2];
sx q[2];
rz(-0.56869555) q[2];
sx q[2];
rz(2.8397172) q[2];
rz(-2.2484696) q[3];
sx q[3];
rz(-1.2877269) q[3];
sx q[3];
rz(2.9368403) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72398913) q[0];
sx q[0];
rz(-1.8483193) q[0];
sx q[0];
rz(1.666477) q[0];
rz(-3.1148615) q[1];
sx q[1];
rz(-1.4550799) q[1];
sx q[1];
rz(1.4310238) q[1];
rz(2.0886291) q[2];
sx q[2];
rz(-0.79635194) q[2];
sx q[2];
rz(1.2780381) q[2];
rz(0.55587739) q[3];
sx q[3];
rz(-2.0049958) q[3];
sx q[3];
rz(-0.47331664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
