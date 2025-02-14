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
rz(3.0348294) q[0];
sx q[0];
rz(-1.3314629) q[0];
sx q[0];
rz(1.3504299) q[0];
rz(0.30836937) q[1];
sx q[1];
rz(6.0344459) q[1];
sx q[1];
rz(11.587172) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0672507) q[0];
sx q[0];
rz(-2.3789211) q[0];
sx q[0];
rz(3.0347155) q[0];
x q[1];
rz(-2.0151695) q[2];
sx q[2];
rz(-2.4326486) q[2];
sx q[2];
rz(-0.63731784) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.82047082) q[1];
sx q[1];
rz(-2.781457) q[1];
sx q[1];
rz(-0.19345691) q[1];
rz(-pi) q[2];
rz(2.7704846) q[3];
sx q[3];
rz(-1.4083997) q[3];
sx q[3];
rz(-2.028167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2778492) q[2];
sx q[2];
rz(-1.2110445) q[2];
sx q[2];
rz(-0.63300526) q[2];
rz(-0.64166075) q[3];
sx q[3];
rz(-0.46052614) q[3];
sx q[3];
rz(-0.7707001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68384701) q[0];
sx q[0];
rz(-0.96413833) q[0];
sx q[0];
rz(0.82474166) q[0];
rz(1.9174346) q[1];
sx q[1];
rz(-0.85547248) q[1];
sx q[1];
rz(2.8214084) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63063568) q[0];
sx q[0];
rz(-0.63586006) q[0];
sx q[0];
rz(2.4102734) q[0];
rz(-1.0616706) q[2];
sx q[2];
rz(-1.3855644) q[2];
sx q[2];
rz(0.19134609) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6588179) q[1];
sx q[1];
rz(-1.59332) q[1];
sx q[1];
rz(3.060633) q[1];
rz(-pi) q[2];
rz(-1.792644) q[3];
sx q[3];
rz(-1.6288405) q[3];
sx q[3];
rz(2.0741418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1227293) q[2];
sx q[2];
rz(-1.9247232) q[2];
sx q[2];
rz(-1.2790206) q[2];
rz(-3.0458798) q[3];
sx q[3];
rz(-2.0666104) q[3];
sx q[3];
rz(1.3191282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9886446) q[0];
sx q[0];
rz(-0.64866346) q[0];
sx q[0];
rz(-1.0308107) q[0];
rz(0.55054322) q[1];
sx q[1];
rz(-2.2830453) q[1];
sx q[1];
rz(1.8969089) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36280945) q[0];
sx q[0];
rz(-0.16805695) q[0];
sx q[0];
rz(0.3248949) q[0];
rz(-pi) q[1];
rz(1.3366401) q[2];
sx q[2];
rz(-2.5114254) q[2];
sx q[2];
rz(0.62335287) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0708187) q[1];
sx q[1];
rz(-2.4386775) q[1];
sx q[1];
rz(1.8902814) q[1];
rz(1.1641509) q[3];
sx q[3];
rz(-0.79504025) q[3];
sx q[3];
rz(0.94889489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9024258) q[2];
sx q[2];
rz(-0.5639762) q[2];
sx q[2];
rz(-1.9695367) q[2];
rz(2.3160882) q[3];
sx q[3];
rz(-1.2647311) q[3];
sx q[3];
rz(2.4765305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0498407) q[0];
sx q[0];
rz(-1.6850152) q[0];
sx q[0];
rz(-0.41393429) q[0];
rz(-2.6992758) q[1];
sx q[1];
rz(-1.0041753) q[1];
sx q[1];
rz(2.5962459) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1236432) q[0];
sx q[0];
rz(-1.845903) q[0];
sx q[0];
rz(-1.6226044) q[0];
rz(0.77540437) q[2];
sx q[2];
rz(-2.0606962) q[2];
sx q[2];
rz(-2.6012914) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.28081799) q[1];
sx q[1];
rz(-1.0549666) q[1];
sx q[1];
rz(-0.3192374) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.993587) q[3];
sx q[3];
rz(-2.0750351) q[3];
sx q[3];
rz(-1.809169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.44117323) q[2];
sx q[2];
rz(-1.2607231) q[2];
sx q[2];
rz(-0.1612266) q[2];
rz(1.553933) q[3];
sx q[3];
rz(-0.59993887) q[3];
sx q[3];
rz(-2.5509295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0597543) q[0];
sx q[0];
rz(-1.2677544) q[0];
sx q[0];
rz(2.4181714) q[0];
rz(1.8266034) q[1];
sx q[1];
rz(-1.2774717) q[1];
sx q[1];
rz(2.1037197) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68907952) q[0];
sx q[0];
rz(-0.85001365) q[0];
sx q[0];
rz(1.6759273) q[0];
rz(2.1234799) q[2];
sx q[2];
rz(-0.92490126) q[2];
sx q[2];
rz(1.8551265) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2813869) q[1];
sx q[1];
rz(-1.1347927) q[1];
sx q[1];
rz(-1.866775) q[1];
rz(-pi) q[2];
x q[2];
rz(0.24017374) q[3];
sx q[3];
rz(-2.5824912) q[3];
sx q[3];
rz(2.9438275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7448685) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(-2.569516) q[2];
rz(2.5165596) q[3];
sx q[3];
rz(-1.0038989) q[3];
sx q[3];
rz(1.9151789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6201651) q[0];
sx q[0];
rz(-0.033997424) q[0];
sx q[0];
rz(-0.52325621) q[0];
rz(-1.322586) q[1];
sx q[1];
rz(-1.4679642) q[1];
sx q[1];
rz(0.62017131) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70411982) q[0];
sx q[0];
rz(-1.0123555) q[0];
sx q[0];
rz(1.062514) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.46862009) q[2];
sx q[2];
rz(-1.6362564) q[2];
sx q[2];
rz(-0.58131389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6002308) q[1];
sx q[1];
rz(-2.5087452) q[1];
sx q[1];
rz(-1.4221627) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67542757) q[3];
sx q[3];
rz(-2.8656061) q[3];
sx q[3];
rz(2.7819862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4124477) q[2];
sx q[2];
rz(-2.3690988) q[2];
sx q[2];
rz(1.283851) q[2];
rz(2.1733952) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(1.5865954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.17484434) q[0];
sx q[0];
rz(-2.0059678) q[0];
sx q[0];
rz(1.211776) q[0];
rz(-2.1719596) q[1];
sx q[1];
rz(-1.3215093) q[1];
sx q[1];
rz(1.8671573) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8961372) q[0];
sx q[0];
rz(-1.4198167) q[0];
sx q[0];
rz(2.2666988) q[0];
rz(-pi) q[1];
rz(-3.0465925) q[2];
sx q[2];
rz(-2.3534905) q[2];
sx q[2];
rz(1.011445) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74531534) q[1];
sx q[1];
rz(-2.4419129) q[1];
sx q[1];
rz(-2.0513644) q[1];
rz(-pi) q[2];
rz(1.2479374) q[3];
sx q[3];
rz(-2.0104694) q[3];
sx q[3];
rz(-2.6722801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4659287) q[2];
sx q[2];
rz(-1.482654) q[2];
sx q[2];
rz(2.8181804) q[2];
rz(3.1392426) q[3];
sx q[3];
rz(-1.4582783) q[3];
sx q[3];
rz(2.5200444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.49331409) q[0];
sx q[0];
rz(-1.4819773) q[0];
sx q[0];
rz(1.7907273) q[0];
rz(1.3510652) q[1];
sx q[1];
rz(-1.2850782) q[1];
sx q[1];
rz(-1.2877119) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0553594) q[0];
sx q[0];
rz(-1.4533193) q[0];
sx q[0];
rz(1.8553977) q[0];
rz(-pi) q[1];
rz(0.28952629) q[2];
sx q[2];
rz(-2.2643201) q[2];
sx q[2];
rz(-2.5159581) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2877174) q[1];
sx q[1];
rz(-2.7734904) q[1];
sx q[1];
rz(-2.6516857) q[1];
x q[2];
rz(-1.5519556) q[3];
sx q[3];
rz(-1.2219489) q[3];
sx q[3];
rz(-2.2486389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6172341) q[2];
sx q[2];
rz(-1.989216) q[2];
sx q[2];
rz(-2.1481245) q[2];
rz(1.0348381) q[3];
sx q[3];
rz(-0.60641685) q[3];
sx q[3];
rz(-2.9724227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8096823) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(1.6140953) q[0];
rz(0.88036674) q[1];
sx q[1];
rz(-2.7208734) q[1];
sx q[1];
rz(0.31164935) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0880222) q[0];
sx q[0];
rz(-2.7854475) q[0];
sx q[0];
rz(2.1689264) q[0];
rz(-3.1252091) q[2];
sx q[2];
rz(-1.3394622) q[2];
sx q[2];
rz(0.45778679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6039201) q[1];
sx q[1];
rz(-2.2419856) q[1];
sx q[1];
rz(-2.0449941) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1470471) q[3];
sx q[3];
rz(-0.69931385) q[3];
sx q[3];
rz(-0.55845234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.46939048) q[2];
sx q[2];
rz(-2.2180874) q[2];
sx q[2];
rz(1.7507318) q[2];
rz(-1.5270799) q[3];
sx q[3];
rz(-1.047784) q[3];
sx q[3];
rz(2.0670149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63956082) q[0];
sx q[0];
rz(-1.9968411) q[0];
sx q[0];
rz(0.61258739) q[0];
rz(0.9901498) q[1];
sx q[1];
rz(-1.9691111) q[1];
sx q[1];
rz(-1.6506857) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25967071) q[0];
sx q[0];
rz(-0.90357354) q[0];
sx q[0];
rz(1.0259969) q[0];
rz(-pi) q[1];
rz(1.2464929) q[2];
sx q[2];
rz(-0.68585472) q[2];
sx q[2];
rz(-1.7795479) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4348169) q[1];
sx q[1];
rz(-1.56456) q[1];
sx q[1];
rz(-1.732279) q[1];
x q[2];
rz(-1.8205216) q[3];
sx q[3];
rz(-2.2533544) q[3];
sx q[3];
rz(-0.54271316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.426173) q[2];
sx q[2];
rz(-0.82948589) q[2];
sx q[2];
rz(3.1066011) q[2];
rz(1.70111) q[3];
sx q[3];
rz(-0.8927497) q[3];
sx q[3];
rz(-0.54097241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69915019) q[0];
sx q[0];
rz(-2.2461666) q[0];
sx q[0];
rz(2.9641892) q[0];
rz(-1.198248) q[1];
sx q[1];
rz(-1.6234963) q[1];
sx q[1];
rz(-2.1877098) q[1];
rz(2.4631029) q[2];
sx q[2];
rz(-1.0889006) q[2];
sx q[2];
rz(2.0550645) q[2];
rz(-2.6260637) q[3];
sx q[3];
rz(-2.3524108) q[3];
sx q[3];
rz(-2.4084211) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
