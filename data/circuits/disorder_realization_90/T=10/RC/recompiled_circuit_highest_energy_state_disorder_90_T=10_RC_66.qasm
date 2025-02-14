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
rz(-2.443402) q[0];
sx q[0];
rz(-2.8289822) q[0];
sx q[0];
rz(-2.1460549) q[0];
rz(1.1065464) q[1];
sx q[1];
rz(-0.76835978) q[1];
sx q[1];
rz(2.8093657) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4058283) q[0];
sx q[0];
rz(-1.1678742) q[0];
sx q[0];
rz(-0.71056239) q[0];
rz(-pi) q[1];
x q[1];
rz(1.840749) q[2];
sx q[2];
rz(-2.6630262) q[2];
sx q[2];
rz(-1.1978483) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2549887) q[1];
sx q[1];
rz(-1.1082778) q[1];
sx q[1];
rz(-2.6673893) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99523441) q[3];
sx q[3];
rz(-2.0166631) q[3];
sx q[3];
rz(-0.80328548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60653162) q[2];
sx q[2];
rz(-2.1273095) q[2];
sx q[2];
rz(0.052058546) q[2];
rz(-1.082513) q[3];
sx q[3];
rz(-0.94729298) q[3];
sx q[3];
rz(-1.989971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5169736) q[0];
sx q[0];
rz(-0.032051429) q[0];
sx q[0];
rz(0.33729851) q[0];
rz(2.9889122) q[1];
sx q[1];
rz(-0.76714194) q[1];
sx q[1];
rz(-1.5229567) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47877676) q[0];
sx q[0];
rz(-2.1838476) q[0];
sx q[0];
rz(-1.7752035) q[0];
rz(-pi) q[1];
rz(0.27045336) q[2];
sx q[2];
rz(-0.65545299) q[2];
sx q[2];
rz(1.356911) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7584658) q[1];
sx q[1];
rz(-0.43769203) q[1];
sx q[1];
rz(-1.0785224) q[1];
rz(-pi) q[2];
rz(-0.2225575) q[3];
sx q[3];
rz(-0.61545505) q[3];
sx q[3];
rz(-0.30729242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0308257) q[2];
sx q[2];
rz(-2.775707) q[2];
sx q[2];
rz(-0.3863253) q[2];
rz(1.857916) q[3];
sx q[3];
rz(-1.9648896) q[3];
sx q[3];
rz(1.9234689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.06685) q[0];
sx q[0];
rz(-2.061494) q[0];
sx q[0];
rz(-2.1225488) q[0];
rz(-2.3661803) q[1];
sx q[1];
rz(-1.8653899) q[1];
sx q[1];
rz(-0.48616854) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3306823) q[0];
sx q[0];
rz(-0.031477246) q[0];
sx q[0];
rz(0.093491836) q[0];
rz(-pi) q[1];
rz(0.48392754) q[2];
sx q[2];
rz(-2.292715) q[2];
sx q[2];
rz(-2.1579822) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0224494) q[1];
sx q[1];
rz(-1.0876942) q[1];
sx q[1];
rz(-1.3234322) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73681946) q[3];
sx q[3];
rz(-0.42356682) q[3];
sx q[3];
rz(-1.119348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.39573085) q[2];
sx q[2];
rz(-1.9755325) q[2];
sx q[2];
rz(-0.80500785) q[2];
rz(-2.4896367) q[3];
sx q[3];
rz(-1.1019573) q[3];
sx q[3];
rz(-2.2074047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
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
rz(-2.4790799) q[0];
sx q[0];
rz(-1.3211687) q[0];
sx q[0];
rz(1.3385734) q[0];
rz(1.592912) q[1];
sx q[1];
rz(-2.4376696) q[1];
sx q[1];
rz(-2.7159363) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.073624728) q[0];
sx q[0];
rz(-1.5789766) q[0];
sx q[0];
rz(-3.0497396) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8786483) q[2];
sx q[2];
rz(-1.6037707) q[2];
sx q[2];
rz(-2.3127382) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50652177) q[1];
sx q[1];
rz(-2.2622411) q[1];
sx q[1];
rz(-2.1060645) q[1];
rz(-pi) q[2];
x q[2];
rz(0.71721806) q[3];
sx q[3];
rz(-1.3558431) q[3];
sx q[3];
rz(1.7058829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.272133) q[2];
sx q[2];
rz(-1.774186) q[2];
sx q[2];
rz(2.7222705) q[2];
rz(-1.1767293) q[3];
sx q[3];
rz(-0.71507016) q[3];
sx q[3];
rz(-0.56593219) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6856573) q[0];
sx q[0];
rz(-0.74860191) q[0];
sx q[0];
rz(-1.5022044) q[0];
rz(-1.7078851) q[1];
sx q[1];
rz(-1.2805484) q[1];
sx q[1];
rz(-2.7388403) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.086632816) q[0];
sx q[0];
rz(-1.3142337) q[0];
sx q[0];
rz(-0.94166468) q[0];
rz(0.5647955) q[2];
sx q[2];
rz(-1.5207054) q[2];
sx q[2];
rz(-1.7518856) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6943446) q[1];
sx q[1];
rz(-1.4741886) q[1];
sx q[1];
rz(3.0098947) q[1];
rz(-pi) q[2];
x q[2];
rz(0.14842965) q[3];
sx q[3];
rz(-2.469329) q[3];
sx q[3];
rz(-0.82151247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5002284) q[2];
sx q[2];
rz(-2.0945175) q[2];
sx q[2];
rz(-0.25263986) q[2];
rz(-2.3371475) q[3];
sx q[3];
rz(-2.6201456) q[3];
sx q[3];
rz(-2.6430602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8135524) q[0];
sx q[0];
rz(-3.033162) q[0];
sx q[0];
rz(0.88388467) q[0];
rz(-0.49993316) q[1];
sx q[1];
rz(-1.4686613) q[1];
sx q[1];
rz(0.51441851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.894569) q[0];
sx q[0];
rz(-1.9049562) q[0];
sx q[0];
rz(2.1804125) q[0];
rz(-pi) q[1];
rz(2.6096658) q[2];
sx q[2];
rz(-1.7830323) q[2];
sx q[2];
rz(-1.3440107) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0658256) q[1];
sx q[1];
rz(-1.7890686) q[1];
sx q[1];
rz(1.1595919) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3522369) q[3];
sx q[3];
rz(-1.5522) q[3];
sx q[3];
rz(-2.9531053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1540692) q[2];
sx q[2];
rz(-1.6776626) q[2];
sx q[2];
rz(3.0164111) q[2];
rz(-2.3212738) q[3];
sx q[3];
rz(-1.1743952) q[3];
sx q[3];
rz(0.34631795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5888551) q[0];
sx q[0];
rz(-0.6468361) q[0];
sx q[0];
rz(0.11189017) q[0];
rz(2.0018068) q[1];
sx q[1];
rz(-0.21509376) q[1];
sx q[1];
rz(1.2824167) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2578141) q[0];
sx q[0];
rz(-0.997966) q[0];
sx q[0];
rz(-2.8158067) q[0];
rz(-2.5865104) q[2];
sx q[2];
rz(-2.8709445) q[2];
sx q[2];
rz(2.3570983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0541898) q[1];
sx q[1];
rz(-2.3034181) q[1];
sx q[1];
rz(1.0637295) q[1];
rz(0.070738465) q[3];
sx q[3];
rz(-0.68545464) q[3];
sx q[3];
rz(0.56927337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.68975) q[2];
sx q[2];
rz(-2.8359783) q[2];
sx q[2];
rz(-1.3172733) q[2];
rz(-0.02056038) q[3];
sx q[3];
rz(-2.8097184) q[3];
sx q[3];
rz(2.1338972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2059712) q[0];
sx q[0];
rz(-0.99140778) q[0];
sx q[0];
rz(-2.8427065) q[0];
rz(1.0609974) q[1];
sx q[1];
rz(-1.7540878) q[1];
sx q[1];
rz(-1.908173) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73288146) q[0];
sx q[0];
rz(-1.2378843) q[0];
sx q[0];
rz(3.0391418) q[0];
x q[1];
rz(1.9154847) q[2];
sx q[2];
rz(-3.0200501) q[2];
sx q[2];
rz(-0.28608382) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3626676) q[1];
sx q[1];
rz(-2.2120419) q[1];
sx q[1];
rz(2.769313) q[1];
rz(2.4806907) q[3];
sx q[3];
rz(-1.1422302) q[3];
sx q[3];
rz(2.1235583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3491106) q[2];
sx q[2];
rz(-1.0976378) q[2];
sx q[2];
rz(0.20723542) q[2];
rz(-0.66050291) q[3];
sx q[3];
rz(-1.0005755) q[3];
sx q[3];
rz(-1.2787261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5527363) q[0];
sx q[0];
rz(-1.739946) q[0];
sx q[0];
rz(-1.2218342) q[0];
rz(-0.51512042) q[1];
sx q[1];
rz(-3.0312067) q[1];
sx q[1];
rz(-0.55330223) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.57896) q[0];
sx q[0];
rz(-1.4427358) q[0];
sx q[0];
rz(1.0109896) q[0];
rz(3.0764941) q[2];
sx q[2];
rz(-1.1506478) q[2];
sx q[2];
rz(-1.1769899) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.33722222) q[1];
sx q[1];
rz(-0.17717277) q[1];
sx q[1];
rz(1.5371662) q[1];
rz(-pi) q[2];
rz(-2.4724602) q[3];
sx q[3];
rz(-1.5520943) q[3];
sx q[3];
rz(-1.6220055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8052266) q[2];
sx q[2];
rz(-1.1104501) q[2];
sx q[2];
rz(0.52555788) q[2];
rz(-1.6622274) q[3];
sx q[3];
rz(-0.95366228) q[3];
sx q[3];
rz(2.3814538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0568327) q[0];
sx q[0];
rz(-2.3469717) q[0];
sx q[0];
rz(-0.42924616) q[0];
rz(1.5224573) q[1];
sx q[1];
rz(-1.2702962) q[1];
sx q[1];
rz(-2.2116908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36741716) q[0];
sx q[0];
rz(-0.40433592) q[0];
sx q[0];
rz(1.0003759) q[0];
rz(1.7903114) q[2];
sx q[2];
rz(-1.6782614) q[2];
sx q[2];
rz(-1.4023413) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7333755) q[1];
sx q[1];
rz(-1.853051) q[1];
sx q[1];
rz(-0.1677081) q[1];
rz(0.553073) q[3];
sx q[3];
rz(-1.6216842) q[3];
sx q[3];
rz(2.1257876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4055736) q[2];
sx q[2];
rz(-2.6655727) q[2];
sx q[2];
rz(2.948577) q[2];
rz(-3.0613101) q[3];
sx q[3];
rz(-2.2716227) q[3];
sx q[3];
rz(1.3840626) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3748462) q[0];
sx q[0];
rz(-1.9369047) q[0];
sx q[0];
rz(1.9932224) q[0];
rz(-2.171352) q[1];
sx q[1];
rz(-1.932825) q[1];
sx q[1];
rz(1.8995151) q[1];
rz(2.316723) q[2];
sx q[2];
rz(-1.3536659) q[2];
sx q[2];
rz(-2.547551) q[2];
rz(-2.8945782) q[3];
sx q[3];
rz(-0.09709662) q[3];
sx q[3];
rz(0.63963565) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
