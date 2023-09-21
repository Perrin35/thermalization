OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.835445) q[0];
sx q[0];
rz(2.4581576) q[0];
sx q[0];
rz(12.0876) q[0];
rz(0.03102826) q[1];
sx q[1];
rz(-1.1614769) q[1];
sx q[1];
rz(-2.5002313) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35533479) q[0];
sx q[0];
rz(-1.9541652) q[0];
sx q[0];
rz(-1.3134365) q[0];
x q[1];
rz(-1.9986721) q[2];
sx q[2];
rz(-1.9170205) q[2];
sx q[2];
rz(1.8510173) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3798843) q[1];
sx q[1];
rz(-2.0239081) q[1];
sx q[1];
rz(-0.74700991) q[1];
rz(2.6061329) q[3];
sx q[3];
rz(-0.70264953) q[3];
sx q[3];
rz(0.040576064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.550094) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(2.6634898) q[2];
rz(1.6889307) q[3];
sx q[3];
rz(-1.0457467) q[3];
sx q[3];
rz(2.5959192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1801382) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(1.0189198) q[0];
rz(1.4787176) q[1];
sx q[1];
rz(-0.61518413) q[1];
sx q[1];
rz(-2.5085124) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79257747) q[0];
sx q[0];
rz(-1.4903869) q[0];
sx q[0];
rz(2.3809459) q[0];
rz(-pi) q[1];
rz(0.98302977) q[2];
sx q[2];
rz(-0.9126185) q[2];
sx q[2];
rz(2.7764729) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.96453297) q[1];
sx q[1];
rz(-2.6033083) q[1];
sx q[1];
rz(2.2798377) q[1];
x q[2];
rz(-0.31140621) q[3];
sx q[3];
rz(-2.2283471) q[3];
sx q[3];
rz(0.94535512) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7818266) q[2];
sx q[2];
rz(-2.7414913) q[2];
sx q[2];
rz(0.11745545) q[2];
rz(-0.30101267) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(-1.7787748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85787073) q[0];
sx q[0];
rz(-2.4283333) q[0];
sx q[0];
rz(3.0531847) q[0];
rz(1.1075426) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(2.450768) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4118977) q[0];
sx q[0];
rz(-1.6459961) q[0];
sx q[0];
rz(-2.9883283) q[0];
x q[1];
rz(1.3356528) q[2];
sx q[2];
rz(-2.3217839) q[2];
sx q[2];
rz(0.90607925) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1580116) q[1];
sx q[1];
rz(-2.9576655) q[1];
sx q[1];
rz(1.4278825) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7615943) q[3];
sx q[3];
rz(-1.899154) q[3];
sx q[3];
rz(-2.2896555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2174125) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(1.4364093) q[3];
sx q[3];
rz(-0.36542106) q[3];
sx q[3];
rz(-1.6586554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0728264) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(0.52247125) q[0];
rz(-0.32896313) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(2.5879588) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41933665) q[0];
sx q[0];
rz(-1.1124848) q[0];
sx q[0];
rz(-0.76461794) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.372924) q[2];
sx q[2];
rz(-0.93066356) q[2];
sx q[2];
rz(-1.6878355) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1045038) q[1];
sx q[1];
rz(-1.4375086) q[1];
sx q[1];
rz(2.3549805) q[1];
x q[2];
rz(3.1133075) q[3];
sx q[3];
rz(-2.3445498) q[3];
sx q[3];
rz(2.459211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5840977) q[2];
sx q[2];
rz(-0.9157052) q[2];
sx q[2];
rz(0.49475691) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(1.2899227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6506127) q[0];
sx q[0];
rz(-0.74478331) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(-2.1814573) q[1];
sx q[1];
rz(-1.1791869) q[1];
sx q[1];
rz(-0.18403149) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9404011) q[0];
sx q[0];
rz(-1.2815676) q[0];
sx q[0];
rz(-1.4616696) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8475624) q[2];
sx q[2];
rz(-1.2185316) q[2];
sx q[2];
rz(0.89154746) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4689456) q[1];
sx q[1];
rz(-0.99860672) q[1];
sx q[1];
rz(-2.1395626) q[1];
rz(-pi) q[2];
rz(-1.3695413) q[3];
sx q[3];
rz(-2.2866837) q[3];
sx q[3];
rz(-1.9191238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2400143) q[2];
sx q[2];
rz(-2.7928536) q[2];
sx q[2];
rz(2.0007755) q[2];
rz(2.7187637) q[3];
sx q[3];
rz(-1.6641649) q[3];
sx q[3];
rz(-1.1269349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35247701) q[0];
sx q[0];
rz(-2.0334091) q[0];
sx q[0];
rz(-0.88678962) q[0];
rz(2.9011762) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(-2.9930847) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0200955) q[0];
sx q[0];
rz(-1.4098865) q[0];
sx q[0];
rz(0.73385977) q[0];
x q[1];
rz(2.912942) q[2];
sx q[2];
rz(-1.6778523) q[2];
sx q[2];
rz(0.46496898) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.442246) q[1];
sx q[1];
rz(-1.8298501) q[1];
sx q[1];
rz(0.11286347) q[1];
x q[2];
rz(-2.6753747) q[3];
sx q[3];
rz(-1.9293474) q[3];
sx q[3];
rz(2.8311604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.285816) q[2];
sx q[2];
rz(-0.4824051) q[2];
sx q[2];
rz(-0.4883858) q[2];
rz(-2.6521902) q[3];
sx q[3];
rz(-0.92178744) q[3];
sx q[3];
rz(1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4483036) q[0];
sx q[0];
rz(-1.8087837) q[0];
sx q[0];
rz(-0.81800246) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-0.39223448) q[1];
sx q[1];
rz(-1.12524) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.274652) q[0];
sx q[0];
rz(-0.38892239) q[0];
sx q[0];
rz(-1.2961943) q[0];
x q[1];
rz(-1.869404) q[2];
sx q[2];
rz(-1.719279) q[2];
sx q[2];
rz(-0.4543002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6553073) q[1];
sx q[1];
rz(-2.7318622) q[1];
sx q[1];
rz(2.3182931) q[1];
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
rz(-2.1606074) q[2];
sx q[2];
rz(0.64458624) q[2];
rz(-1.4792431) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(1.4054327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1709764) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(-1.53565) q[0];
rz(1.9203141) q[1];
sx q[1];
rz(-1.5131283) q[1];
sx q[1];
rz(-2.1309526) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3686417) q[0];
sx q[0];
rz(-2.3803582) q[0];
sx q[0];
rz(1.8038473) q[0];
x q[1];
rz(-1.0755195) q[2];
sx q[2];
rz(-1.8349378) q[2];
sx q[2];
rz(1.8758945) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9606564) q[1];
sx q[1];
rz(-1.3410543) q[1];
sx q[1];
rz(-0.22682637) q[1];
rz(-pi) q[2];
rz(-1.478785) q[3];
sx q[3];
rz(-1.1717456) q[3];
sx q[3];
rz(1.9012746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5236139) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(-0.69407216) q[2];
rz(-2.5726035) q[3];
sx q[3];
rz(-1.6945972) q[3];
sx q[3];
rz(-0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(-2.6468497) q[0];
rz(-0.61839473) q[1];
sx q[1];
rz(-1.4952375) q[1];
sx q[1];
rz(-0.075597413) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9915174) q[0];
sx q[0];
rz(-2.4872094) q[0];
sx q[0];
rz(2.0983216) q[0];
x q[1];
rz(-1.3952903) q[2];
sx q[2];
rz(-1.6298461) q[2];
sx q[2];
rz(-2.4621071) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.96502111) q[1];
sx q[1];
rz(-0.88639835) q[1];
sx q[1];
rz(-1.7040764) q[1];
x q[2];
rz(-2.5209386) q[3];
sx q[3];
rz(-1.422158) q[3];
sx q[3];
rz(-1.8295446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33346924) q[2];
sx q[2];
rz(-1.7227017) q[2];
sx q[2];
rz(1.3405651) q[2];
rz(2.8373485) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(2.5568967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2952404) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(-2.9647968) q[0];
rz(1.2416174) q[1];
sx q[1];
rz(-2.5071564) q[1];
sx q[1];
rz(-2.7005844) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1045751) q[0];
sx q[0];
rz(-1.772176) q[0];
sx q[0];
rz(-1.7778394) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15372865) q[2];
sx q[2];
rz(-1.6445465) q[2];
sx q[2];
rz(-1.7388625) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1383789) q[1];
sx q[1];
rz(-0.56204501) q[1];
sx q[1];
rz(-2.2467062) q[1];
rz(-2.7367758) q[3];
sx q[3];
rz(-2.4768618) q[3];
sx q[3];
rz(2.3354195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5796154) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-2.8397172) q[2];
rz(0.89312303) q[3];
sx q[3];
rz(-1.2877269) q[3];
sx q[3];
rz(2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4176035) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(3.1148615) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-0.46840657) q[2];
sx q[2];
rz(-0.9006587) q[2];
sx q[2];
rz(-2.547154) q[2];
rz(-2.5857153) q[3];
sx q[3];
rz(-2.0049958) q[3];
sx q[3];
rz(-0.47331664) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];