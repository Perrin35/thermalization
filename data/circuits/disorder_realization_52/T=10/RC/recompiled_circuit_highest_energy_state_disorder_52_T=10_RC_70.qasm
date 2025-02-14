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
rz(-2.4177457) q[0];
sx q[0];
rz(-1.9404193) q[0];
sx q[0];
rz(2.6475651) q[0];
rz(-1.5863034) q[1];
sx q[1];
rz(-2.2145693) q[1];
sx q[1];
rz(2.291099) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14334003) q[0];
sx q[0];
rz(-1.7368642) q[0];
sx q[0];
rz(-0.89675007) q[0];
rz(0.57588864) q[2];
sx q[2];
rz(-2.401899) q[2];
sx q[2];
rz(2.0832555) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8181818) q[1];
sx q[1];
rz(-1.2480624) q[1];
sx q[1];
rz(1.8052015) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36601074) q[3];
sx q[3];
rz(-1.5029534) q[3];
sx q[3];
rz(0.74857601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4832619) q[2];
sx q[2];
rz(-2.5265145) q[2];
sx q[2];
rz(-2.1966546) q[2];
rz(-3.1350709) q[3];
sx q[3];
rz(-2.3801453) q[3];
sx q[3];
rz(-2.884088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0354804) q[0];
sx q[0];
rz(-0.98863125) q[0];
sx q[0];
rz(-0.60580564) q[0];
rz(-2.2094191) q[1];
sx q[1];
rz(-1.4235539) q[1];
sx q[1];
rz(3.0359643) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52077928) q[0];
sx q[0];
rz(-2.438669) q[0];
sx q[0];
rz(-0.23971324) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1828853) q[2];
sx q[2];
rz(-2.5664133) q[2];
sx q[2];
rz(1.5782331) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3117506) q[1];
sx q[1];
rz(-1.1195604) q[1];
sx q[1];
rz(1.2430906) q[1];
rz(-1.1607882) q[3];
sx q[3];
rz(-1.9191529) q[3];
sx q[3];
rz(-0.94545555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.85670829) q[2];
sx q[2];
rz(-1.7785037) q[2];
sx q[2];
rz(-1.6178097) q[2];
rz(0.81426042) q[3];
sx q[3];
rz(-1.7819449) q[3];
sx q[3];
rz(2.2566569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.098323671) q[0];
sx q[0];
rz(-2.3461778) q[0];
sx q[0];
rz(1.5650308) q[0];
rz(-2.1463429) q[1];
sx q[1];
rz(-2.1627656) q[1];
sx q[1];
rz(-2.3547122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1724479) q[0];
sx q[0];
rz(-0.24227628) q[0];
sx q[0];
rz(0.97357018) q[0];
rz(1.4866531) q[2];
sx q[2];
rz(-1.9168789) q[2];
sx q[2];
rz(-0.10135191) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2128558) q[1];
sx q[1];
rz(-2.6784705) q[1];
sx q[1];
rz(-0.76174163) q[1];
x q[2];
rz(-1.5110817) q[3];
sx q[3];
rz(-1.51126) q[3];
sx q[3];
rz(1.2950031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3028822) q[2];
sx q[2];
rz(-1.6532712) q[2];
sx q[2];
rz(0.6380471) q[2];
rz(-0.4979411) q[3];
sx q[3];
rz(-1.0718071) q[3];
sx q[3];
rz(1.0897442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.2447253) q[0];
sx q[0];
rz(-1.3571955) q[0];
sx q[0];
rz(-0.063752739) q[0];
rz(1.4490734) q[1];
sx q[1];
rz(-1.8359102) q[1];
sx q[1];
rz(1.8316899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5777912) q[0];
sx q[0];
rz(-3.0534857) q[0];
sx q[0];
rz(2.7035575) q[0];
rz(-0.74964995) q[2];
sx q[2];
rz(-0.26826619) q[2];
sx q[2];
rz(-1.1531354) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79015649) q[1];
sx q[1];
rz(-2.3871941) q[1];
sx q[1];
rz(-2.960207) q[1];
rz(2.2579262) q[3];
sx q[3];
rz(-1.2901297) q[3];
sx q[3];
rz(1.3786045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.070907585) q[2];
sx q[2];
rz(-2.0347774) q[2];
sx q[2];
rz(-0.10565383) q[2];
rz(1.1834772) q[3];
sx q[3];
rz(-0.89619291) q[3];
sx q[3];
rz(2.4699874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.79528177) q[0];
sx q[0];
rz(-0.91049894) q[0];
sx q[0];
rz(1.3443391) q[0];
rz(-2.4686939) q[1];
sx q[1];
rz(-1.8310603) q[1];
sx q[1];
rz(-3.1281298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1707662) q[0];
sx q[0];
rz(-1.5970236) q[0];
sx q[0];
rz(0.88084014) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1273753) q[2];
sx q[2];
rz(-1.1566312) q[2];
sx q[2];
rz(-0.80314512) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9078482) q[1];
sx q[1];
rz(-1.9861167) q[1];
sx q[1];
rz(-1.5563346) q[1];
rz(-0.3810639) q[3];
sx q[3];
rz(-2.237202) q[3];
sx q[3];
rz(-1.8997418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5992735) q[2];
sx q[2];
rz(-2.4794674) q[2];
sx q[2];
rz(-1.8947961) q[2];
rz(0.072619297) q[3];
sx q[3];
rz(-2.9473372) q[3];
sx q[3];
rz(2.8323925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4319864) q[0];
sx q[0];
rz(-1.1259587) q[0];
sx q[0];
rz(0.11269888) q[0];
rz(0.20507774) q[1];
sx q[1];
rz(-1.3321184) q[1];
sx q[1];
rz(0.90788666) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31602177) q[0];
sx q[0];
rz(-0.96222034) q[0];
sx q[0];
rz(2.2951503) q[0];
x q[1];
rz(1.142318) q[2];
sx q[2];
rz(-1.0016164) q[2];
sx q[2];
rz(2.6250397) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7026313) q[1];
sx q[1];
rz(-1.4794032) q[1];
sx q[1];
rz(-0.63321873) q[1];
rz(-pi) q[2];
rz(-2.870918) q[3];
sx q[3];
rz(-2.4587016) q[3];
sx q[3];
rz(2.0840933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9997361) q[2];
sx q[2];
rz(-2.8047968) q[2];
sx q[2];
rz(-3.0736308) q[2];
rz(1.147602) q[3];
sx q[3];
rz(-1.8736898) q[3];
sx q[3];
rz(1.8806774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1934018) q[0];
sx q[0];
rz(-1.3113439) q[0];
sx q[0];
rz(-0.46352682) q[0];
rz(-2.9947128) q[1];
sx q[1];
rz(-1.905922) q[1];
sx q[1];
rz(0.99259496) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6382321) q[0];
sx q[0];
rz(-1.8192756) q[0];
sx q[0];
rz(2.4186224) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8568749) q[2];
sx q[2];
rz(-1.5062638) q[2];
sx q[2];
rz(2.4308824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0907572) q[1];
sx q[1];
rz(-2.2244284) q[1];
sx q[1];
rz(-0.98934503) q[1];
rz(2.4086508) q[3];
sx q[3];
rz(-1.2228726) q[3];
sx q[3];
rz(-2.3331353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2831882) q[2];
sx q[2];
rz(-1.7121199) q[2];
sx q[2];
rz(-2.6962213) q[2];
rz(1.3238268) q[3];
sx q[3];
rz(-2.0711074) q[3];
sx q[3];
rz(-1.0858735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0367947) q[0];
sx q[0];
rz(-3.1322271) q[0];
sx q[0];
rz(-0.15583663) q[0];
rz(1.2771295) q[1];
sx q[1];
rz(-1.3469478) q[1];
sx q[1];
rz(-3.1256622) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30980047) q[0];
sx q[0];
rz(-2.3224717) q[0];
sx q[0];
rz(-0.33162923) q[0];
rz(-pi) q[1];
rz(2.7819949) q[2];
sx q[2];
rz(-2.3162875) q[2];
sx q[2];
rz(-1.2957089) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.585938) q[1];
sx q[1];
rz(-0.92472142) q[1];
sx q[1];
rz(-2.487961) q[1];
rz(-2.5941284) q[3];
sx q[3];
rz(-2.4261835) q[3];
sx q[3];
rz(2.5304731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5126123) q[2];
sx q[2];
rz(-3.0292558) q[2];
sx q[2];
rz(-3.0158499) q[2];
rz(-0.9564774) q[3];
sx q[3];
rz(-1.7185017) q[3];
sx q[3];
rz(2.8023348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9823031) q[0];
sx q[0];
rz(-2.90137) q[0];
sx q[0];
rz(-0.45641986) q[0];
rz(-1.5727111) q[1];
sx q[1];
rz(-1.0036889) q[1];
sx q[1];
rz(0.68663866) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7266885) q[0];
sx q[0];
rz(-0.81302887) q[0];
sx q[0];
rz(0.70720478) q[0];
rz(-pi) q[1];
rz(-1.7438269) q[2];
sx q[2];
rz(-1.695172) q[2];
sx q[2];
rz(1.591452) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9514044) q[1];
sx q[1];
rz(-1.7177637) q[1];
sx q[1];
rz(1.4495871) q[1];
rz(-pi) q[2];
rz(-0.17466361) q[3];
sx q[3];
rz(-1.6399989) q[3];
sx q[3];
rz(-0.098002794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3829284) q[2];
sx q[2];
rz(-1.1229346) q[2];
sx q[2];
rz(-2.0056966) q[2];
rz(2.1618333) q[3];
sx q[3];
rz(-1.3330678) q[3];
sx q[3];
rz(-2.4433344) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45282388) q[0];
sx q[0];
rz(-1.8166421) q[0];
sx q[0];
rz(0.45595566) q[0];
rz(-3.0988354) q[1];
sx q[1];
rz(-1.9330934) q[1];
sx q[1];
rz(2.4483689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9915145) q[0];
sx q[0];
rz(-2.743609) q[0];
sx q[0];
rz(1.8282169) q[0];
x q[1];
rz(2.3568527) q[2];
sx q[2];
rz(-2.5973136) q[2];
sx q[2];
rz(-1.2960824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9849423) q[1];
sx q[1];
rz(-1.6688111) q[1];
sx q[1];
rz(-3.1344862) q[1];
x q[2];
rz(-2.5467123) q[3];
sx q[3];
rz(-0.52973807) q[3];
sx q[3];
rz(1.4929508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.91528714) q[2];
sx q[2];
rz(-1.8733571) q[2];
sx q[2];
rz(-2.477296) q[2];
rz(1.213446) q[3];
sx q[3];
rz(-2.4197141) q[3];
sx q[3];
rz(-2.2693995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6912457) q[0];
sx q[0];
rz(-2.3007614) q[0];
sx q[0];
rz(1.7578516) q[0];
rz(1.6830403) q[1];
sx q[1];
rz(-2.8589307) q[1];
sx q[1];
rz(0.9160441) q[1];
rz(-3.0677879) q[2];
sx q[2];
rz(-0.37141411) q[2];
sx q[2];
rz(0.85167533) q[2];
rz(-2.1324329) q[3];
sx q[3];
rz(-2.1957993) q[3];
sx q[3];
rz(-1.5002863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
