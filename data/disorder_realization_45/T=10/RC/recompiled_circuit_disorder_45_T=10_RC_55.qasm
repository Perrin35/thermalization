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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3135932) q[0];
sx q[0];
rz(-1.3324954) q[0];
sx q[0];
rz(0.39512623) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9986721) q[2];
sx q[2];
rz(-1.2245721) q[2];
sx q[2];
rz(-1.2905754) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76170834) q[1];
sx q[1];
rz(-2.0239081) q[1];
sx q[1];
rz(2.3945827) q[1];
rz(1.9786505) q[3];
sx q[3];
rz(-0.98148325) q[3];
sx q[3];
rz(0.70139635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59149867) q[2];
sx q[2];
rz(-1.2167565) q[2];
sx q[2];
rz(2.6634898) q[2];
rz(-1.452662) q[3];
sx q[3];
rz(-2.0958459) q[3];
sx q[3];
rz(12/(7*pi)) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96145445) q[0];
sx q[0];
rz(-2.366876) q[0];
sx q[0];
rz(-1.0189198) q[0];
rz(1.4787176) q[1];
sx q[1];
rz(-2.5264085) q[1];
sx q[1];
rz(2.5085124) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7019254) q[0];
sx q[0];
rz(-0.81322008) q[0];
sx q[0];
rz(-1.4600091) q[0];
rz(-pi) q[1];
rz(-2.3929246) q[2];
sx q[2];
rz(-1.1166755) q[2];
sx q[2];
rz(1.5926966) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9620348) q[1];
sx q[1];
rz(-1.9704559) q[1];
sx q[1];
rz(0.37079294) q[1];
x q[2];
rz(-2.2523316) q[3];
sx q[3];
rz(-1.325843) q[3];
sx q[3];
rz(-0.8196866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.35976609) q[2];
sx q[2];
rz(-0.40010139) q[2];
sx q[2];
rz(-0.11745545) q[2];
rz(2.84058) q[3];
sx q[3];
rz(-1.3344701) q[3];
sx q[3];
rz(1.3628179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(2.2837219) q[0];
sx q[0];
rz(-2.4283333) q[0];
sx q[0];
rz(-0.088407956) q[0];
rz(-1.1075426) q[1];
sx q[1];
rz(-2.2231367) q[1];
sx q[1];
rz(0.69082469) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.847825) q[0];
sx q[0];
rz(-2.9710037) q[0];
sx q[0];
rz(2.6831496) q[0];
x q[1];
rz(-1.3356528) q[2];
sx q[2];
rz(-0.81980875) q[2];
sx q[2];
rz(-2.2355134) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.01268) q[1];
sx q[1];
rz(-1.3887654) q[1];
sx q[1];
rz(0.026489594) q[1];
rz(-0.74294528) q[3];
sx q[3];
rz(-0.49693123) q[3];
sx q[3];
rz(-3.1019773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2174125) q[2];
sx q[2];
rz(-1.8841382) q[2];
sx q[2];
rz(0.10647354) q[2];
rz(1.7051833) q[3];
sx q[3];
rz(-2.7761716) q[3];
sx q[3];
rz(1.4829372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0728264) q[0];
sx q[0];
rz(-0.23385736) q[0];
sx q[0];
rz(2.6191214) q[0];
rz(0.32896313) q[1];
sx q[1];
rz(-1.6429699) q[1];
sx q[1];
rz(-2.5879588) q[1];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.81981084) q[2];
sx q[2];
rz(-2.1852583) q[2];
sx q[2];
rz(-2.7053506) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.33448996) q[1];
sx q[1];
rz(-2.346171) q[1];
sx q[1];
rz(-0.18717488) q[1];
rz(-0.028285154) q[3];
sx q[3];
rz(-2.3445498) q[3];
sx q[3];
rz(-0.68238168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.557495) q[2];
sx q[2];
rz(-2.2258874) q[2];
sx q[2];
rz(-2.6468357) q[2];
rz(2.2385712) q[3];
sx q[3];
rz(-1.5477864) q[3];
sx q[3];
rz(-1.8516699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49098) q[0];
sx q[0];
rz(-2.3968093) q[0];
sx q[0];
rz(-1.0850798) q[0];
rz(2.1814573) q[1];
sx q[1];
rz(-1.9624058) q[1];
sx q[1];
rz(2.9575612) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9404011) q[0];
sx q[0];
rz(-1.860025) q[0];
sx q[0];
rz(1.4616696) q[0];
rz(-pi) q[1];
rz(0.36501546) q[2];
sx q[2];
rz(-1.8301617) q[2];
sx q[2];
rz(-2.5600381) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3369493) q[1];
sx q[1];
rz(-0.78362432) q[1];
sx q[1];
rz(0.6964535) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4155951) q[3];
sx q[3];
rz(-1.4193924) q[3];
sx q[3];
rz(-2.9263673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.90157834) q[2];
sx q[2];
rz(-0.34873909) q[2];
sx q[2];
rz(-1.1408172) q[2];
rz(-2.7187637) q[3];
sx q[3];
rz(-1.4774277) q[3];
sx q[3];
rz(2.0146577) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7891156) q[0];
sx q[0];
rz(-1.1081835) q[0];
sx q[0];
rz(-2.254803) q[0];
rz(-0.24041644) q[1];
sx q[1];
rz(-2.1227032) q[1];
sx q[1];
rz(-2.9930847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37489031) q[0];
sx q[0];
rz(-2.3935211) q[0];
sx q[0];
rz(-2.9038249) q[0];
rz(-pi) q[1];
rz(-2.912942) q[2];
sx q[2];
rz(-1.6778523) q[2];
sx q[2];
rz(2.6766237) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6993466) q[1];
sx q[1];
rz(-1.3117426) q[1];
sx q[1];
rz(3.0287292) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46621795) q[3];
sx q[3];
rz(-1.9293474) q[3];
sx q[3];
rz(2.8311604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85577661) q[2];
sx q[2];
rz(-2.6591876) q[2];
sx q[2];
rz(2.6532069) q[2];
rz(-2.6521902) q[3];
sx q[3];
rz(-2.2198052) q[3];
sx q[3];
rz(-1.2667806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69328904) q[0];
sx q[0];
rz(-1.332809) q[0];
sx q[0];
rz(2.3235902) q[0];
rz(1.8891634) q[1];
sx q[1];
rz(-2.7493582) q[1];
sx q[1];
rz(1.12524) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1624807) q[0];
sx q[0];
rz(-1.9444124) q[0];
sx q[0];
rz(3.0309249) q[0];
x q[1];
rz(1.1004041) q[2];
sx q[2];
rz(-2.8090968) q[2];
sx q[2];
rz(-1.5770797) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.78912567) q[1];
sx q[1];
rz(-1.8450292) q[1];
sx q[1];
rz(1.2624361) q[1];
rz(-1.0057698) q[3];
sx q[3];
rz(-0.98494782) q[3];
sx q[3];
rz(2.1573531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0685048) q[2];
sx q[2];
rz(-2.1606074) q[2];
sx q[2];
rz(-0.64458624) q[2];
rz(1.4792431) q[3];
sx q[3];
rz(-0.95696604) q[3];
sx q[3];
rz(1.7361599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1709764) q[0];
sx q[0];
rz(-0.068844065) q[0];
sx q[0];
rz(1.53565) q[0];
rz(-1.2212785) q[1];
sx q[1];
rz(-1.6284643) q[1];
sx q[1];
rz(-1.01064) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7692524) q[0];
sx q[0];
rz(-1.7307889) q[0];
sx q[0];
rz(2.3183547) q[0];
x q[1];
rz(0.29823093) q[2];
sx q[2];
rz(-2.0474307) q[2];
sx q[2];
rz(-0.44520608) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6992221) q[1];
sx q[1];
rz(-1.3500299) q[1];
sx q[1];
rz(1.8063596) q[1];
x q[2];
rz(-2.9270494) q[3];
sx q[3];
rz(-2.7326267) q[3];
sx q[3];
rz(-1.4734801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6179787) q[2];
sx q[2];
rz(-2.8221059) q[2];
sx q[2];
rz(2.4475205) q[2];
rz(0.56898919) q[3];
sx q[3];
rz(-1.4469955) q[3];
sx q[3];
rz(0.93769658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086031832) q[0];
sx q[0];
rz(-1.9984351) q[0];
sx q[0];
rz(2.6468497) q[0];
rz(0.61839473) q[1];
sx q[1];
rz(-1.6463552) q[1];
sx q[1];
rz(3.0659952) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2878802) q[0];
sx q[0];
rz(-1.8822077) q[0];
sx q[0];
rz(0.98543723) q[0];
rz(1.2443301) q[2];
sx q[2];
rz(-2.9565161) q[2];
sx q[2];
rz(-1.9290123) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1739993) q[1];
sx q[1];
rz(-0.6951957) q[1];
sx q[1];
rz(-2.980152) q[1];
x q[2];
rz(-0.62065403) q[3];
sx q[3];
rz(-1.7194347) q[3];
sx q[3];
rz(1.312048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.33346924) q[2];
sx q[2];
rz(-1.418891) q[2];
sx q[2];
rz(1.3405651) q[2];
rz(2.8373485) q[3];
sx q[3];
rz(-0.86383581) q[3];
sx q[3];
rz(-0.58469599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2952404) q[0];
sx q[0];
rz(-1.8363991) q[0];
sx q[0];
rz(-2.9647968) q[0];
rz(1.8999752) q[1];
sx q[1];
rz(-0.63443628) q[1];
sx q[1];
rz(0.44100824) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2947185) q[0];
sx q[0];
rz(-0.28781578) q[0];
sx q[0];
rz(-0.78886445) q[0];
rz(-pi) q[1];
rz(0.15372865) q[2];
sx q[2];
rz(-1.6445465) q[2];
sx q[2];
rz(-1.7388625) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9778206) q[1];
sx q[1];
rz(-1.9107011) q[1];
sx q[1];
rz(1.1141067) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2714083) q[3];
sx q[3];
rz(-0.96794879) q[3];
sx q[3];
rz(2.833948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5796154) q[2];
sx q[2];
rz(-2.5728971) q[2];
sx q[2];
rz(-0.30187541) q[2];
rz(0.89312303) q[3];
sx q[3];
rz(-1.8538657) q[3];
sx q[3];
rz(-2.9368403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4176035) q[0];
sx q[0];
rz(-1.2932734) q[0];
sx q[0];
rz(-1.4751157) q[0];
rz(-0.026731116) q[1];
sx q[1];
rz(-1.6865128) q[1];
sx q[1];
rz(-1.7105688) q[1];
rz(-2.6731861) q[2];
sx q[2];
rz(-2.240934) q[2];
sx q[2];
rz(0.59443867) q[2];
rz(-0.7209575) q[3];
sx q[3];
rz(-0.69098916) q[3];
sx q[3];
rz(0.50222764) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
