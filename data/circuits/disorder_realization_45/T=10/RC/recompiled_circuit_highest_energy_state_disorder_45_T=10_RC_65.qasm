OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2019914) q[0];
sx q[0];
rz(-2.658598) q[0];
sx q[0];
rz(1.510409) q[0];
rz(3.0022439) q[1];
sx q[1];
rz(-2.5820093) q[1];
sx q[1];
rz(-2.411627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48929991) q[0];
sx q[0];
rz(-2.518939) q[0];
sx q[0];
rz(-1.8239418) q[0];
rz(-2.5467186) q[2];
sx q[2];
rz(-1.7584137) q[2];
sx q[2];
rz(-1.7930195) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12355676) q[1];
sx q[1];
rz(-0.58096051) q[1];
sx q[1];
rz(0.058079795) q[1];
x q[2];
rz(-1.8584941) q[3];
sx q[3];
rz(-1.3511063) q[3];
sx q[3];
rz(0.9385329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1175179) q[2];
sx q[2];
rz(-0.2710318) q[2];
sx q[2];
rz(-2.3226341) q[2];
rz(3.135318) q[3];
sx q[3];
rz(-1.9082021) q[3];
sx q[3];
rz(2.1591469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(2.5421802) q[0];
rz(0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(-0.74554602) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4670938) q[0];
sx q[0];
rz(-1.8601396) q[0];
sx q[0];
rz(1.4379005) q[0];
rz(-pi) q[1];
rz(1.2471871) q[2];
sx q[2];
rz(-0.83728803) q[2];
sx q[2];
rz(0.74786438) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88953253) q[1];
sx q[1];
rz(-0.65450689) q[1];
sx q[1];
rz(-0.26328523) q[1];
rz(-pi) q[2];
rz(2.638444) q[3];
sx q[3];
rz(-0.9434349) q[3];
sx q[3];
rz(-0.37093758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(1.4073184) q[2];
rz(-0.84434858) q[3];
sx q[3];
rz(-2.307939) q[3];
sx q[3];
rz(-2.1048529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5832962) q[0];
sx q[0];
rz(-1.946227) q[0];
sx q[0];
rz(2.1287647) q[0];
rz(-0.60802513) q[1];
sx q[1];
rz(-1.5826179) q[1];
sx q[1];
rz(-1.8720522) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69508775) q[0];
sx q[0];
rz(-1.9534612) q[0];
sx q[0];
rz(-1.0951359) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6744587) q[2];
sx q[2];
rz(-0.67358635) q[2];
sx q[2];
rz(-1.4530593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.041965466) q[1];
sx q[1];
rz(-1.3746975) q[1];
sx q[1];
rz(2.4134497) q[1];
rz(1.2392912) q[3];
sx q[3];
rz(-1.1024144) q[3];
sx q[3];
rz(1.4862332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.515392) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(-1.8499648) q[2];
rz(-3.1332704) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(0.21361175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13852791) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(-1.7648765) q[0];
rz(1.5273013) q[1];
sx q[1];
rz(-1.2034028) q[1];
sx q[1];
rz(0.52070224) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5138869) q[0];
sx q[0];
rz(-1.2180274) q[0];
sx q[0];
rz(-1.1350495) q[0];
rz(-pi) q[1];
rz(2.7815656) q[2];
sx q[2];
rz(-1.0485149) q[2];
sx q[2];
rz(1.6619267) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5508073) q[1];
sx q[1];
rz(-1.1414764) q[1];
sx q[1];
rz(-0.6145668) q[1];
rz(-pi) q[2];
rz(-2.6163231) q[3];
sx q[3];
rz(-2.161918) q[3];
sx q[3];
rz(-2.1256465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.5248096) q[2];
sx q[2];
rz(-0.7889792) q[2];
sx q[2];
rz(1.3516124) q[2];
rz(1.0194408) q[3];
sx q[3];
rz(-2.5717058) q[3];
sx q[3];
rz(1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.549642) q[0];
sx q[0];
rz(-1.4627946) q[0];
sx q[0];
rz(-0.098966448) q[0];
rz(2.4348266) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(-0.4695355) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5564559) q[0];
sx q[0];
rz(-3.0017108) q[0];
sx q[0];
rz(1.3993457) q[0];
rz(-0.97909285) q[2];
sx q[2];
rz(-0.62034494) q[2];
sx q[2];
rz(1.1466591) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5597804) q[1];
sx q[1];
rz(-0.99039927) q[1];
sx q[1];
rz(2.0349166) q[1];
rz(-1.4245701) q[3];
sx q[3];
rz(-1.0195882) q[3];
sx q[3];
rz(2.220038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2609451) q[2];
sx q[2];
rz(-1.8069043) q[2];
sx q[2];
rz(-1.8355231) q[2];
rz(-2.7169054) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(2.2556321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11923085) q[0];
sx q[0];
rz(-0.62218085) q[0];
sx q[0];
rz(1.8192044) q[0];
rz(-2.0139096) q[1];
sx q[1];
rz(-1.6219982) q[1];
sx q[1];
rz(1.3816396) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36395006) q[0];
sx q[0];
rz(-2.3780883) q[0];
sx q[0];
rz(-1.6958773) q[0];
rz(-pi) q[1];
rz(-0.83397978) q[2];
sx q[2];
rz(-0.72089855) q[2];
sx q[2];
rz(-2.3490459) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0387011) q[1];
sx q[1];
rz(-2.3159317) q[1];
sx q[1];
rz(-1.3172512) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30140169) q[3];
sx q[3];
rz(-1.5001138) q[3];
sx q[3];
rz(1.7168728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7606925) q[2];
sx q[2];
rz(-1.5194632) q[2];
sx q[2];
rz(-0.48119989) q[2];
rz(-0.5091269) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5803489) q[0];
sx q[0];
rz(-2.6600397) q[0];
sx q[0];
rz(0.089381889) q[0];
rz(-2.0629758) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(2.1481029) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2334796) q[0];
sx q[0];
rz(-2.3589239) q[0];
sx q[0];
rz(3.0872869) q[0];
rz(-2.2682796) q[2];
sx q[2];
rz(-0.9160348) q[2];
sx q[2];
rz(2.8755434) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.76748874) q[1];
sx q[1];
rz(-1.933177) q[1];
sx q[1];
rz(0.15717536) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0519846) q[3];
sx q[3];
rz(-2.1362274) q[3];
sx q[3];
rz(2.4861002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.79954687) q[2];
sx q[2];
rz(-2.5265103) q[2];
sx q[2];
rz(-2.7395524) q[2];
rz(-2.0461931) q[3];
sx q[3];
rz(-1.9270555) q[3];
sx q[3];
rz(-0.57797617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47356975) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(1.8079669) q[0];
rz(2.2857621) q[1];
sx q[1];
rz(-1.7355093) q[1];
sx q[1];
rz(2.2241101) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8569023) q[0];
sx q[0];
rz(-1.7698341) q[0];
sx q[0];
rz(-3.0375098) q[0];
rz(-pi) q[1];
rz(3.0227674) q[2];
sx q[2];
rz(-2.1911494) q[2];
sx q[2];
rz(-2.9147749) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5405824) q[1];
sx q[1];
rz(-1.9153908) q[1];
sx q[1];
rz(0.2356727) q[1];
x q[2];
rz(-2.4948984) q[3];
sx q[3];
rz(-1.1521253) q[3];
sx q[3];
rz(-0.44548098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.63460073) q[2];
sx q[2];
rz(-1.4054106) q[2];
sx q[2];
rz(-1.290192) q[2];
rz(2.2318132) q[3];
sx q[3];
rz(-1.6981643) q[3];
sx q[3];
rz(3.1006052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9057587) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(-0.27715096) q[0];
rz(1.9174891) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(-1.7701497) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0540028) q[0];
sx q[0];
rz(-1.3773019) q[0];
sx q[0];
rz(-0.67275472) q[0];
x q[1];
rz(-0.65176378) q[2];
sx q[2];
rz(-2.4874176) q[2];
sx q[2];
rz(-0.68133611) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7595948) q[1];
sx q[1];
rz(-1.7773668) q[1];
sx q[1];
rz(0.91464154) q[1];
rz(-pi) q[2];
rz(0.84884642) q[3];
sx q[3];
rz(-1.4392142) q[3];
sx q[3];
rz(2.2717486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93854967) q[2];
sx q[2];
rz(-2.2800192) q[2];
sx q[2];
rz(2.4533563) q[2];
rz(-2.8271683) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(-1.2615874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37453434) q[0];
sx q[0];
rz(-2.2735167) q[0];
sx q[0];
rz(2.5700997) q[0];
rz(-0.68069619) q[1];
sx q[1];
rz(-1.1704159) q[1];
sx q[1];
rz(-2.4933955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8342736) q[0];
sx q[0];
rz(-1.5569485) q[0];
sx q[0];
rz(-1.5939006) q[0];
x q[1];
rz(0.83309116) q[2];
sx q[2];
rz(-2.0948185) q[2];
sx q[2];
rz(-1.810871) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.98013377) q[1];
sx q[1];
rz(-0.69936434) q[1];
sx q[1];
rz(-3.1070263) q[1];
rz(-pi) q[2];
rz(0.33721029) q[3];
sx q[3];
rz(-2.5993532) q[3];
sx q[3];
rz(2.6553939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8615243) q[2];
sx q[2];
rz(-1.0736977) q[2];
sx q[2];
rz(1.466922) q[2];
rz(0.17659771) q[3];
sx q[3];
rz(-2.4956775) q[3];
sx q[3];
rz(1.2944029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8661154) q[0];
sx q[0];
rz(-1.7452411) q[0];
sx q[0];
rz(-1.2757975) q[0];
rz(2.7571309) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(-1.8865042) q[2];
sx q[2];
rz(-1.2366625) q[2];
sx q[2];
rz(-1.2319215) q[2];
rz(-2.3220358) q[3];
sx q[3];
rz(-2.4278276) q[3];
sx q[3];
rz(0.62558382) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
