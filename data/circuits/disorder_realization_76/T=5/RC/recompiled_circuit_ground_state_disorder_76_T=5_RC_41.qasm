OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4467093) q[0];
sx q[0];
rz(-2.0599685) q[0];
sx q[0];
rz(-2.4021436) q[0];
rz(2.3820355) q[1];
sx q[1];
rz(-1.3146725) q[1];
sx q[1];
rz(1.4256328) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48910558) q[0];
sx q[0];
rz(-1.28363) q[0];
sx q[0];
rz(2.1225568) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47503586) q[2];
sx q[2];
rz(-0.54752195) q[2];
sx q[2];
rz(2.5753367) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5162075) q[1];
sx q[1];
rz(-0.16730669) q[1];
sx q[1];
rz(-2.0545261) q[1];
x q[2];
rz(1.5208779) q[3];
sx q[3];
rz(-2.422159) q[3];
sx q[3];
rz(2.0572544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90929675) q[2];
sx q[2];
rz(-0.84315073) q[2];
sx q[2];
rz(2.9006531) q[2];
rz(-0.029189261) q[3];
sx q[3];
rz(-1.3386644) q[3];
sx q[3];
rz(1.1791112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9772684) q[0];
sx q[0];
rz(-1.203953) q[0];
sx q[0];
rz(2.6569195) q[0];
rz(-0.67972216) q[1];
sx q[1];
rz(-1.2815963) q[1];
sx q[1];
rz(-1.9814804) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.401448) q[0];
sx q[0];
rz(-2.8370246) q[0];
sx q[0];
rz(-0.061621678) q[0];
rz(-pi) q[1];
rz(-0.87748973) q[2];
sx q[2];
rz(-2.1602767) q[2];
sx q[2];
rz(-1.1699007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96148032) q[1];
sx q[1];
rz(-2.4870076) q[1];
sx q[1];
rz(-2.727319) q[1];
rz(-0.5228225) q[3];
sx q[3];
rz(-1.2827875) q[3];
sx q[3];
rz(-3.0436181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.29740563) q[2];
sx q[2];
rz(-0.30218267) q[2];
sx q[2];
rz(2.6837132) q[2];
rz(-2.0186021) q[3];
sx q[3];
rz(-2.1419958) q[3];
sx q[3];
rz(-0.56639731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26894012) q[0];
sx q[0];
rz(-2.432423) q[0];
sx q[0];
rz(0.031524468) q[0];
rz(2.8541376) q[1];
sx q[1];
rz(-0.87702409) q[1];
sx q[1];
rz(-1.9256176) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368204) q[0];
sx q[0];
rz(-2.2520814) q[0];
sx q[0];
rz(0.70685203) q[0];
rz(-pi) q[1];
rz(-2.7599665) q[2];
sx q[2];
rz(-2.6256621) q[2];
sx q[2];
rz(-2.6037773) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3423914) q[1];
sx q[1];
rz(-1.2974129) q[1];
sx q[1];
rz(2.8552516) q[1];
x q[2];
rz(-0.9312882) q[3];
sx q[3];
rz(-2.1289325) q[3];
sx q[3];
rz(-2.1344413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1027801) q[2];
sx q[2];
rz(-1.5422042) q[2];
sx q[2];
rz(0.83596027) q[2];
rz(1.3577667) q[3];
sx q[3];
rz(-2.416555) q[3];
sx q[3];
rz(2.3939705) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2826071) q[0];
sx q[0];
rz(-1.7361807) q[0];
sx q[0];
rz(-3.0614241) q[0];
rz(2.0591586) q[1];
sx q[1];
rz(-0.23324649) q[1];
sx q[1];
rz(-1.3287883) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69705582) q[0];
sx q[0];
rz(-1.2460684) q[0];
sx q[0];
rz(2.6022807) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3983509) q[2];
sx q[2];
rz(-1.8891462) q[2];
sx q[2];
rz(0.47296745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.13212407) q[1];
sx q[1];
rz(-1.4259725) q[1];
sx q[1];
rz(-1.138746) q[1];
rz(2.6065663) q[3];
sx q[3];
rz(-2.7561829) q[3];
sx q[3];
rz(-2.7148394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3749915) q[2];
sx q[2];
rz(-2.1583755) q[2];
sx q[2];
rz(-0.90744606) q[2];
rz(0.67029101) q[3];
sx q[3];
rz(-0.82796103) q[3];
sx q[3];
rz(-1.7214187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9012673) q[0];
sx q[0];
rz(-2.2690052) q[0];
sx q[0];
rz(-2.7101044) q[0];
rz(2.0101428) q[1];
sx q[1];
rz(-1.1791469) q[1];
sx q[1];
rz(0.44050899) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9517072) q[0];
sx q[0];
rz(-1.110297) q[0];
sx q[0];
rz(-3.0104464) q[0];
x q[1];
rz(-1.126664) q[2];
sx q[2];
rz(-2.6878549) q[2];
sx q[2];
rz(2.1778088) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4104328) q[1];
sx q[1];
rz(-1.9010547) q[1];
sx q[1];
rz(1.6124658) q[1];
rz(-pi) q[2];
rz(-2.2045361) q[3];
sx q[3];
rz(-1.0770633) q[3];
sx q[3];
rz(-2.1449094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0195007) q[2];
sx q[2];
rz(-0.92635218) q[2];
sx q[2];
rz(2.7276373) q[2];
rz(-1.3881989) q[3];
sx q[3];
rz(-1.5940758) q[3];
sx q[3];
rz(3.0164914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6269161) q[0];
sx q[0];
rz(-1.4056982) q[0];
sx q[0];
rz(-2.2077014) q[0];
rz(-0.67032188) q[1];
sx q[1];
rz(-1.2443845) q[1];
sx q[1];
rz(-1.8062887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40437296) q[0];
sx q[0];
rz(-1.6707194) q[0];
sx q[0];
rz(2.1735682) q[0];
rz(-pi) q[1];
rz(-1.5455294) q[2];
sx q[2];
rz(-1.7284365) q[2];
sx q[2];
rz(2.9767286) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.858243) q[1];
sx q[1];
rz(-1.5933541) q[1];
sx q[1];
rz(1.4580887) q[1];
x q[2];
rz(-2.504566) q[3];
sx q[3];
rz(-1.5777794) q[3];
sx q[3];
rz(1.03656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.89207092) q[2];
sx q[2];
rz(-2.8688909) q[2];
sx q[2];
rz(1.8708694) q[2];
rz(-1.3696085) q[3];
sx q[3];
rz(-1.9000051) q[3];
sx q[3];
rz(-0.52186596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18913604) q[0];
sx q[0];
rz(-0.58681762) q[0];
sx q[0];
rz(-2.9879046) q[0];
rz(-0.20248374) q[1];
sx q[1];
rz(-1.9053562) q[1];
sx q[1];
rz(2.1930146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0517387) q[0];
sx q[0];
rz(-1.7719381) q[0];
sx q[0];
rz(2.5178943) q[0];
x q[1];
rz(-1.9948629) q[2];
sx q[2];
rz(-2.4575007) q[2];
sx q[2];
rz(-2.7685194) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7708009) q[1];
sx q[1];
rz(-1.3389412) q[1];
sx q[1];
rz(1.3045909) q[1];
rz(-3.0074869) q[3];
sx q[3];
rz(-0.92536345) q[3];
sx q[3];
rz(-0.2303309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1367246) q[2];
sx q[2];
rz(-2.0437045) q[2];
sx q[2];
rz(-0.0014121545) q[2];
rz(-0.062156113) q[3];
sx q[3];
rz(-1.4096189) q[3];
sx q[3];
rz(2.1803161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38178) q[0];
sx q[0];
rz(-0.75730046) q[0];
sx q[0];
rz(-1.2148452) q[0];
rz(0.75792056) q[1];
sx q[1];
rz(-2.6140723) q[1];
sx q[1];
rz(-0.55799276) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3787631) q[0];
sx q[0];
rz(-0.42233322) q[0];
sx q[0];
rz(-0.93393737) q[0];
rz(-pi) q[1];
rz(-2.9079014) q[2];
sx q[2];
rz(-2.5776641) q[2];
sx q[2];
rz(0.83966161) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.8415738) q[1];
sx q[1];
rz(-1.8487367) q[1];
sx q[1];
rz(2.071408) q[1];
rz(3.1147546) q[3];
sx q[3];
rz(-1.8478113) q[3];
sx q[3];
rz(-2.7153496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5743635) q[2];
sx q[2];
rz(-1.1939253) q[2];
sx q[2];
rz(-2.8723259) q[2];
rz(1.1632129) q[3];
sx q[3];
rz(-2.5389157) q[3];
sx q[3];
rz(0.24063024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6729386) q[0];
sx q[0];
rz(-2.1042295) q[0];
sx q[0];
rz(1.0733676) q[0];
rz(1.1277554) q[1];
sx q[1];
rz(-0.22722166) q[1];
sx q[1];
rz(-2.5849297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2909235) q[0];
sx q[0];
rz(-1.8890831) q[0];
sx q[0];
rz(0.14936635) q[0];
rz(-pi) q[1];
rz(2.5684909) q[2];
sx q[2];
rz(-2.0656043) q[2];
sx q[2];
rz(3.0798517) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5334415) q[1];
sx q[1];
rz(-1.9284298) q[1];
sx q[1];
rz(0.74700345) q[1];
x q[2];
rz(-0.23454097) q[3];
sx q[3];
rz(-1.6681156) q[3];
sx q[3];
rz(3.0458801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2399981) q[2];
sx q[2];
rz(-1.0866714) q[2];
sx q[2];
rz(1.979801) q[2];
rz(-2.6084172) q[3];
sx q[3];
rz(-0.52459255) q[3];
sx q[3];
rz(0.1575135) q[3];
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
rz(-pi/2) q[3];
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
rz(-0.48851442) q[0];
sx q[0];
rz(-0.97452679) q[0];
sx q[0];
rz(2.5467806) q[0];
rz(-2.5626903) q[1];
sx q[1];
rz(-1.4756823) q[1];
sx q[1];
rz(-2.4822809) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66856495) q[0];
sx q[0];
rz(-2.1921232) q[0];
sx q[0];
rz(-1.7354119) q[0];
rz(-pi) q[1];
rz(-1.1289146) q[2];
sx q[2];
rz(-1.0143447) q[2];
sx q[2];
rz(2.93612) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6281471) q[1];
sx q[1];
rz(-1.4288578) q[1];
sx q[1];
rz(1.8177086) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8774977) q[3];
sx q[3];
rz(-1.3066829) q[3];
sx q[3];
rz(0.3981638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.79002964) q[2];
sx q[2];
rz(-1.943925) q[2];
sx q[2];
rz(0.053038049) q[2];
rz(1.7985581) q[3];
sx q[3];
rz(-1.8648632) q[3];
sx q[3];
rz(-3.067335) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85528436) q[0];
sx q[0];
rz(-1.7779779) q[0];
sx q[0];
rz(1.5386982) q[0];
rz(0.098943624) q[1];
sx q[1];
rz(-1.3700486) q[1];
sx q[1];
rz(0.055421967) q[1];
rz(-0.40545736) q[2];
sx q[2];
rz(-1.5787072) q[2];
sx q[2];
rz(-2.2070259) q[2];
rz(-2.4002038) q[3];
sx q[3];
rz(-1.6869873) q[3];
sx q[3];
rz(1.6968873) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
