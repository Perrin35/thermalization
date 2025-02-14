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
rz(3.6245873) q[0];
sx q[0];
rz(10.935187) q[0];
rz(-0.13934879) q[1];
sx q[1];
rz(-0.55958334) q[1];
sx q[1];
rz(-0.72996563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18099526) q[0];
sx q[0];
rz(-0.9708403) q[0];
sx q[0];
rz(-2.9636895) q[0];
rz(-pi) q[1];
rz(0.59487409) q[2];
sx q[2];
rz(-1.7584137) q[2];
sx q[2];
rz(-1.7930195) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19300191) q[1];
sx q[1];
rz(-2.1506502) q[1];
sx q[1];
rz(1.5327044) q[1];
x q[2];
rz(2.9128051) q[3];
sx q[3];
rz(-1.8513894) q[3];
sx q[3];
rz(-0.56786637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.0240747) q[2];
sx q[2];
rz(-0.2710318) q[2];
sx q[2];
rz(-2.3226341) q[2];
rz(-0.0062746127) q[3];
sx q[3];
rz(-1.9082021) q[3];
sx q[3];
rz(2.1591469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5821424) q[0];
sx q[0];
rz(-1.5511976) q[0];
sx q[0];
rz(-0.5994125) q[0];
rz(0.86743152) q[1];
sx q[1];
rz(-2.1116833) q[1];
sx q[1];
rz(-0.74554602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6744989) q[0];
sx q[0];
rz(-1.8601396) q[0];
sx q[0];
rz(1.4379005) q[0];
rz(-0.33919427) q[2];
sx q[2];
rz(-0.789398) q[2];
sx q[2];
rz(-1.2123002) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.47059575) q[1];
sx q[1];
rz(-1.4116916) q[1];
sx q[1];
rz(-2.5039423) q[1];
rz(-pi) q[2];
x q[2];
rz(0.50314869) q[3];
sx q[3];
rz(-2.1981578) q[3];
sx q[3];
rz(-0.37093758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.6856689) q[2];
sx q[2];
rz(-0.29752877) q[2];
sx q[2];
rz(1.4073184) q[2];
rz(2.2972441) q[3];
sx q[3];
rz(-2.307939) q[3];
sx q[3];
rz(1.0367397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5832962) q[0];
sx q[0];
rz(-1.1953657) q[0];
sx q[0];
rz(-1.012828) q[0];
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
rz(2.0464567) q[0];
rz(-3.059194) q[2];
sx q[2];
rz(-2.2401056) q[2];
sx q[2];
rz(-1.5853887) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.041965466) q[1];
sx q[1];
rz(-1.3746975) q[1];
sx q[1];
rz(0.72814299) q[1];
rz(-0.49130398) q[3];
sx q[3];
rz(-1.8654556) q[3];
sx q[3];
rz(0.06959411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6262007) q[2];
sx q[2];
rz(-1.0582558) q[2];
sx q[2];
rz(-1.8499648) q[2];
rz(0.0083222566) q[3];
sx q[3];
rz(-0.95299995) q[3];
sx q[3];
rz(0.21361175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13852791) q[0];
sx q[0];
rz(-0.55532885) q[0];
sx q[0];
rz(-1.7648765) q[0];
rz(1.6142913) q[1];
sx q[1];
rz(-1.9381899) q[1];
sx q[1];
rz(0.52070224) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4462858) q[0];
sx q[0];
rz(-2.5881672) q[0];
sx q[0];
rz(0.85352104) q[0];
x q[1];
rz(0.3600271) q[2];
sx q[2];
rz(-1.0485149) q[2];
sx q[2];
rz(1.479666) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5508073) q[1];
sx q[1];
rz(-1.1414764) q[1];
sx q[1];
rz(2.5270259) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52526955) q[3];
sx q[3];
rz(-0.97967463) q[3];
sx q[3];
rz(-1.0159462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.5248096) q[2];
sx q[2];
rz(-2.3526134) q[2];
sx q[2];
rz(1.3516124) q[2];
rz(2.1221519) q[3];
sx q[3];
rz(-0.56988684) q[3];
sx q[3];
rz(1.1714237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.549642) q[0];
sx q[0];
rz(-1.6787981) q[0];
sx q[0];
rz(0.098966448) q[0];
rz(2.4348266) q[1];
sx q[1];
rz(-0.87044972) q[1];
sx q[1];
rz(-0.4695355) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5564559) q[0];
sx q[0];
rz(-0.13988189) q[0];
sx q[0];
rz(-1.742247) q[0];
rz(0.97909285) q[2];
sx q[2];
rz(-0.62034494) q[2];
sx q[2];
rz(1.9949335) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.27891913) q[1];
sx q[1];
rz(-1.9544744) q[1];
sx q[1];
rz(2.5088599) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23271493) q[3];
sx q[3];
rz(-0.56833) q[3];
sx q[3];
rz(-1.1956904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.88064757) q[2];
sx q[2];
rz(-1.8069043) q[2];
sx q[2];
rz(1.3060695) q[2];
rz(0.4246873) q[3];
sx q[3];
rz(-1.3771907) q[3];
sx q[3];
rz(2.2556321) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0223618) q[0];
sx q[0];
rz(-0.62218085) q[0];
sx q[0];
rz(1.8192044) q[0];
rz(-2.0139096) q[1];
sx q[1];
rz(-1.6219982) q[1];
sx q[1];
rz(1.3816396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2974325) q[0];
sx q[0];
rz(-1.6571665) q[0];
sx q[0];
rz(-2.3303836) q[0];
x q[1];
rz(0.99389771) q[2];
sx q[2];
rz(-2.0303147) q[2];
sx q[2];
rz(-1.3764868) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4037212) q[1];
sx q[1];
rz(-0.779169) q[1];
sx q[1];
rz(2.8761151) q[1];
rz(-0.30140169) q[3];
sx q[3];
rz(-1.6414789) q[3];
sx q[3];
rz(1.7168728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7606925) q[2];
sx q[2];
rz(-1.6221294) q[2];
sx q[2];
rz(-0.48119989) q[2];
rz(2.6324658) q[3];
sx q[3];
rz(-0.23734084) q[3];
sx q[3];
rz(2.7595162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5612438) q[0];
sx q[0];
rz(-0.48155293) q[0];
sx q[0];
rz(0.089381889) q[0];
rz(-1.0786169) q[1];
sx q[1];
rz(-1.6856472) q[1];
sx q[1];
rz(-2.1481029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8427994) q[0];
sx q[0];
rz(-1.6090819) q[0];
sx q[0];
rz(0.78193112) q[0];
x q[1];
rz(-2.3553576) q[2];
sx q[2];
rz(-1.0362384) q[2];
sx q[2];
rz(1.7766118) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9844173) q[1];
rz(-pi) q[2];
rz(-2.4782789) q[3];
sx q[3];
rz(-0.74771008) q[3];
sx q[3];
rz(1.472689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.79954687) q[2];
sx q[2];
rz(-0.61508238) q[2];
sx q[2];
rz(-0.40204027) q[2];
rz(-1.0953995) q[3];
sx q[3];
rz(-1.2145372) q[3];
sx q[3];
rz(2.5636165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47356975) q[0];
sx q[0];
rz(-2.3325925) q[0];
sx q[0];
rz(1.3336257) q[0];
rz(0.85583055) q[1];
sx q[1];
rz(-1.4060833) q[1];
sx q[1];
rz(2.2241101) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8348351) q[0];
sx q[0];
rz(-1.4687755) q[0];
sx q[0];
rz(-1.7708885) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11882527) q[2];
sx q[2];
rz(-2.1911494) q[2];
sx q[2];
rz(2.9147749) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6010103) q[1];
sx q[1];
rz(-1.9153908) q[1];
sx q[1];
rz(0.2356727) q[1];
rz(0.64669426) q[3];
sx q[3];
rz(-1.1521253) q[3];
sx q[3];
rz(-0.44548098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5069919) q[2];
sx q[2];
rz(-1.736182) q[2];
sx q[2];
rz(1.290192) q[2];
rz(-0.90977943) q[3];
sx q[3];
rz(-1.4434283) q[3];
sx q[3];
rz(-3.1006052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23583394) q[0];
sx q[0];
rz(-1.9941149) q[0];
sx q[0];
rz(-2.8644417) q[0];
rz(-1.9174891) q[1];
sx q[1];
rz(-1.6033019) q[1];
sx q[1];
rz(1.7701497) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087589892) q[0];
sx q[0];
rz(-1.7642907) q[0];
sx q[0];
rz(0.67275472) q[0];
x q[1];
rz(0.54746898) q[2];
sx q[2];
rz(-1.1927422) q[2];
sx q[2];
rz(-2.7965656) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.692457) q[1];
sx q[1];
rz(-2.4583011) q[1];
sx q[1];
rz(-1.2399252) q[1];
rz(-pi) q[2];
rz(1.7684494) q[3];
sx q[3];
rz(-2.4098793) q[3];
sx q[3];
rz(-2.5885575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.203043) q[2];
sx q[2];
rz(-0.86157346) q[2];
sx q[2];
rz(0.68823632) q[2];
rz(-2.8271683) q[3];
sx q[3];
rz(-2.7721072) q[3];
sx q[3];
rz(-1.2615874) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37453434) q[0];
sx q[0];
rz(-2.2735167) q[0];
sx q[0];
rz(-0.57149291) q[0];
rz(-2.4608965) q[1];
sx q[1];
rz(-1.9711767) q[1];
sx q[1];
rz(0.64819711) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4179736) q[0];
sx q[0];
rz(-0.026935808) q[0];
sx q[0];
rz(1.0307781) q[0];
rz(0.83309116) q[2];
sx q[2];
rz(-1.0467741) q[2];
sx q[2];
rz(-1.3307216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.56420621) q[1];
sx q[1];
rz(-1.5485475) q[1];
sx q[1];
rz(0.69907) q[1];
rz(-pi) q[2];
rz(-1.7675507) q[3];
sx q[3];
rz(-1.0621539) q[3];
sx q[3];
rz(-3.0439049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.28006831) q[2];
sx q[2];
rz(-2.0678949) q[2];
sx q[2];
rz(1.6746707) q[2];
rz(-0.17659771) q[3];
sx q[3];
rz(-2.4956775) q[3];
sx q[3];
rz(-1.2944029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27547729) q[0];
sx q[0];
rz(-1.7452411) q[0];
sx q[0];
rz(-1.2757975) q[0];
rz(2.7571309) q[1];
sx q[1];
rz(-1.6356331) q[1];
sx q[1];
rz(2.5148139) q[1];
rz(-0.72970978) q[2];
sx q[2];
rz(-2.6860102) q[2];
sx q[2];
rz(-0.44853733) q[2];
rz(-1.0064784) q[3];
sx q[3];
rz(-1.1075533) q[3];
sx q[3];
rz(1.5816734) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
