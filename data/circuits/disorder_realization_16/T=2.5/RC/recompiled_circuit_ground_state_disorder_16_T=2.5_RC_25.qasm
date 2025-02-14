OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.032441) q[0];
sx q[0];
rz(-1.1530131) q[0];
sx q[0];
rz(-0.074540019) q[0];
rz(1.6115161) q[1];
sx q[1];
rz(-3.0822152) q[1];
sx q[1];
rz(-0.52409726) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7439047) q[0];
sx q[0];
rz(-2.4389704) q[0];
sx q[0];
rz(-0.98708679) q[0];
rz(0.33676001) q[2];
sx q[2];
rz(-1.3408061) q[2];
sx q[2];
rz(-2.5652792) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.300507) q[1];
sx q[1];
rz(-0.99113722) q[1];
sx q[1];
rz(2.4546844) q[1];
x q[2];
rz(-0.040848537) q[3];
sx q[3];
rz(-1.4788879) q[3];
sx q[3];
rz(0.81113863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8147627) q[2];
sx q[2];
rz(-0.49850285) q[2];
sx q[2];
rz(0.89547431) q[2];
rz(2.879066) q[3];
sx q[3];
rz(-0.34590507) q[3];
sx q[3];
rz(3.053022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.217591) q[0];
sx q[0];
rz(-1.1955248) q[0];
sx q[0];
rz(-3.0254645) q[0];
rz(0.63236347) q[1];
sx q[1];
rz(-0.26591161) q[1];
sx q[1];
rz(1.6557453) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0362186) q[0];
sx q[0];
rz(-1.6375132) q[0];
sx q[0];
rz(1.5479888) q[0];
x q[1];
rz(-2.112442) q[2];
sx q[2];
rz(-1.7746762) q[2];
sx q[2];
rz(-1.980305) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9159269) q[1];
sx q[1];
rz(-1.6767595) q[1];
sx q[1];
rz(3.0132921) q[1];
x q[2];
rz(1.7255351) q[3];
sx q[3];
rz(-0.70634281) q[3];
sx q[3];
rz(1.4541436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.31132013) q[2];
sx q[2];
rz(-2.1612284) q[2];
sx q[2];
rz(-2.2118528) q[2];
rz(2.7021507) q[3];
sx q[3];
rz(-1.5403055) q[3];
sx q[3];
rz(2.3846227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.065652549) q[0];
sx q[0];
rz(-2.3891698) q[0];
sx q[0];
rz(2.2546076) q[0];
rz(-1.2607964) q[1];
sx q[1];
rz(-2.0340684) q[1];
sx q[1];
rz(-1.6541803) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5607779) q[0];
sx q[0];
rz(-3.0240623) q[0];
sx q[0];
rz(1.5752026) q[0];
rz(-pi) q[1];
rz(-2.7486408) q[2];
sx q[2];
rz(-1.4825511) q[2];
sx q[2];
rz(1.5618344) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9076242) q[1];
sx q[1];
rz(-2.8115936) q[1];
sx q[1];
rz(0.054363175) q[1];
rz(0.47580274) q[3];
sx q[3];
rz(-2.5047205) q[3];
sx q[3];
rz(2.3726557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6416574) q[2];
sx q[2];
rz(-0.96330088) q[2];
sx q[2];
rz(-0.43528834) q[2];
rz(0.38517243) q[3];
sx q[3];
rz(-1.2110854) q[3];
sx q[3];
rz(0.96863532) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73548549) q[0];
sx q[0];
rz(-0.084097363) q[0];
sx q[0];
rz(0.054585833) q[0];
rz(-1.6361884) q[1];
sx q[1];
rz(-1.9278229) q[1];
sx q[1];
rz(-0.41753599) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97440244) q[0];
sx q[0];
rz(-1.1343114) q[0];
sx q[0];
rz(-1.5065326) q[0];
rz(1.2415177) q[2];
sx q[2];
rz(-2.1770475) q[2];
sx q[2];
rz(2.1700567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3247973) q[1];
sx q[1];
rz(-2.9831605) q[1];
sx q[1];
rz(-0.043912391) q[1];
rz(3.0746704) q[3];
sx q[3];
rz(-1.8324781) q[3];
sx q[3];
rz(0.24623888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2893082) q[2];
sx q[2];
rz(-0.70312971) q[2];
sx q[2];
rz(-0.0011778041) q[2];
rz(-1.7581455) q[3];
sx q[3];
rz(-0.26418424) q[3];
sx q[3];
rz(2.0980825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99172878) q[0];
sx q[0];
rz(-2.9809451) q[0];
sx q[0];
rz(-2.7657261) q[0];
rz(2.1916892) q[1];
sx q[1];
rz(-1.4774731) q[1];
sx q[1];
rz(0.72521597) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8729018) q[0];
sx q[0];
rz(-0.16624545) q[0];
sx q[0];
rz(0.77456559) q[0];
rz(2.2967417) q[2];
sx q[2];
rz(-2.2965659) q[2];
sx q[2];
rz(1.8595075) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9053697) q[1];
sx q[1];
rz(-0.097644383) q[1];
sx q[1];
rz(-3.0821441) q[1];
rz(2.83738) q[3];
sx q[3];
rz(-0.72644573) q[3];
sx q[3];
rz(-3.1151047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3788562) q[2];
sx q[2];
rz(-1.4521658) q[2];
sx q[2];
rz(0.47259304) q[2];
rz(-3.0224814) q[3];
sx q[3];
rz(-0.62370682) q[3];
sx q[3];
rz(-0.81010336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8238207) q[0];
sx q[0];
rz(-2.9404984) q[0];
sx q[0];
rz(0.066545181) q[0];
rz(-0.19509527) q[1];
sx q[1];
rz(-2.0654443) q[1];
sx q[1];
rz(-1.1544352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6267363) q[0];
sx q[0];
rz(-1.8033228) q[0];
sx q[0];
rz(2.0621081) q[0];
x q[1];
rz(2.965314) q[2];
sx q[2];
rz(-0.46277324) q[2];
sx q[2];
rz(2.0001992) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5408386) q[1];
sx q[1];
rz(-0.37950613) q[1];
sx q[1];
rz(-1.2787766) q[1];
rz(-1.5510173) q[3];
sx q[3];
rz(-0.58073509) q[3];
sx q[3];
rz(1.6466717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2164312) q[2];
sx q[2];
rz(-1.5387646) q[2];
sx q[2];
rz(-0.72539854) q[2];
rz(-1.409449) q[3];
sx q[3];
rz(-1.1812482) q[3];
sx q[3];
rz(-0.60666549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0856536) q[0];
sx q[0];
rz(-0.5760718) q[0];
sx q[0];
rz(0.90580171) q[0];
rz(0.55465758) q[1];
sx q[1];
rz(-2.3034818) q[1];
sx q[1];
rz(2.99756) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5921123) q[0];
sx q[0];
rz(-1.0735816) q[0];
sx q[0];
rz(-0.29659941) q[0];
x q[1];
rz(1.7212825) q[2];
sx q[2];
rz(-0.51739365) q[2];
sx q[2];
rz(2.4491893) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1293948) q[1];
sx q[1];
rz(-1.8499814) q[1];
sx q[1];
rz(1.9703945) q[1];
rz(-pi) q[2];
rz(1.0353851) q[3];
sx q[3];
rz(-0.28942063) q[3];
sx q[3];
rz(-2.0339703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4732699) q[2];
sx q[2];
rz(-2.747135) q[2];
sx q[2];
rz(-2.1486166) q[2];
rz(-2.9571577) q[3];
sx q[3];
rz(-2.1858229) q[3];
sx q[3];
rz(-2.1485645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1297146) q[0];
sx q[0];
rz(-2.5373902) q[0];
sx q[0];
rz(-2.7469444) q[0];
rz(-0.53798419) q[1];
sx q[1];
rz(-0.6901651) q[1];
sx q[1];
rz(-1.949955) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2054128) q[0];
sx q[0];
rz(-1.8260341) q[0];
sx q[0];
rz(2.5827239) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.678474) q[2];
sx q[2];
rz(-2.7292433) q[2];
sx q[2];
rz(2.8099443) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0816457) q[1];
sx q[1];
rz(-2.4103202) q[1];
sx q[1];
rz(2.251723) q[1];
rz(-pi) q[2];
rz(-0.6362149) q[3];
sx q[3];
rz(-1.9466361) q[3];
sx q[3];
rz(-1.587589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.789232) q[2];
sx q[2];
rz(-1.612794) q[2];
sx q[2];
rz(-0.94190502) q[2];
rz(0.76147979) q[3];
sx q[3];
rz(-1.1250291) q[3];
sx q[3];
rz(-0.040508125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86904675) q[0];
sx q[0];
rz(-1.7750374) q[0];
sx q[0];
rz(-2.0696562) q[0];
rz(-2.5275224) q[1];
sx q[1];
rz(-2.2449988) q[1];
sx q[1];
rz(2.695172) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0907132) q[0];
sx q[0];
rz(-1.8514367) q[0];
sx q[0];
rz(-2.2801823) q[0];
rz(0.69009366) q[2];
sx q[2];
rz(-1.9999256) q[2];
sx q[2];
rz(-0.60455017) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0298456) q[1];
sx q[1];
rz(-0.66737755) q[1];
sx q[1];
rz(-2.9633629) q[1];
x q[2];
rz(2.6362801) q[3];
sx q[3];
rz(-1.4243505) q[3];
sx q[3];
rz(-3.05299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8322231) q[2];
sx q[2];
rz(-0.95755219) q[2];
sx q[2];
rz(1.8831801) q[2];
rz(1.9084515) q[3];
sx q[3];
rz(-0.62658739) q[3];
sx q[3];
rz(1.8241749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42245361) q[0];
sx q[0];
rz(-2.3261676) q[0];
sx q[0];
rz(-1.7183787) q[0];
rz(-0.24249679) q[1];
sx q[1];
rz(-0.47912326) q[1];
sx q[1];
rz(-1.2538145) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4187188) q[0];
sx q[0];
rz(-2.6986045) q[0];
sx q[0];
rz(-1.1305869) q[0];
rz(-pi) q[1];
rz(0.10842936) q[2];
sx q[2];
rz(-3.0745818) q[2];
sx q[2];
rz(-2.5052414) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.36404836) q[1];
sx q[1];
rz(-1.8692888) q[1];
sx q[1];
rz(1.0320028) q[1];
rz(-pi) q[2];
rz(-2.9217072) q[3];
sx q[3];
rz(-1.9337092) q[3];
sx q[3];
rz(0.21838926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8157876) q[2];
sx q[2];
rz(-2.3141404) q[2];
sx q[2];
rz(-3.117756) q[2];
rz(-1.8627953) q[3];
sx q[3];
rz(-0.33134225) q[3];
sx q[3];
rz(2.3648025) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9780289) q[0];
sx q[0];
rz(-1.704287) q[0];
sx q[0];
rz(-1.1821672) q[0];
rz(2.4143746) q[1];
sx q[1];
rz(-1.6738418) q[1];
sx q[1];
rz(2.3152836) q[1];
rz(1.6979065) q[2];
sx q[2];
rz(-0.15227507) q[2];
sx q[2];
rz(-1.1823552) q[2];
rz(-1.1990697) q[3];
sx q[3];
rz(-3.0109497) q[3];
sx q[3];
rz(2.0554832) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
