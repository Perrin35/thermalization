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
rz(2.9958148) q[0];
sx q[0];
rz(-1.79359) q[0];
sx q[0];
rz(-0.3682799) q[0];
rz(-0.085973099) q[1];
sx q[1];
rz(-2.2987125) q[1];
sx q[1];
rz(-2.892363) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1220142) q[0];
sx q[0];
rz(-1.6824357) q[0];
sx q[0];
rz(0.48099244) q[0];
rz(-pi) q[1];
rz(2.8299242) q[2];
sx q[2];
rz(-1.0578007) q[2];
sx q[2];
rz(2.5530961) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.27654821) q[1];
sx q[1];
rz(-0.88403349) q[1];
sx q[1];
rz(-1.8774607) q[1];
x q[2];
rz(-1.6843819) q[3];
sx q[3];
rz(-1.2143597) q[3];
sx q[3];
rz(0.64768744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2117846) q[2];
sx q[2];
rz(-2.4346508) q[2];
sx q[2];
rz(2.9014273) q[2];
rz(2.5943878) q[3];
sx q[3];
rz(-1.6271084) q[3];
sx q[3];
rz(2.1521294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3419679) q[0];
sx q[0];
rz(-0.37536055) q[0];
sx q[0];
rz(-0.64390916) q[0];
rz(-1.0470692) q[1];
sx q[1];
rz(-1.964485) q[1];
sx q[1];
rz(0.26328304) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41516427) q[0];
sx q[0];
rz(-1.5405354) q[0];
sx q[0];
rz(-1.5983768) q[0];
rz(3.075454) q[2];
sx q[2];
rz(-1.5060584) q[2];
sx q[2];
rz(1.7867225) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4810679) q[1];
sx q[1];
rz(-0.98036375) q[1];
sx q[1];
rz(-0.51312889) q[1];
x q[2];
rz(-2.7277522) q[3];
sx q[3];
rz(-0.97025052) q[3];
sx q[3];
rz(-0.59970705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.457966) q[2];
sx q[2];
rz(-1.3764952) q[2];
sx q[2];
rz(2.4158884) q[2];
rz(-0.30722412) q[3];
sx q[3];
rz(-1.8351464) q[3];
sx q[3];
rz(2.0871625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3211408) q[0];
sx q[0];
rz(-1.4844002) q[0];
sx q[0];
rz(2.3725574) q[0];
rz(3.0342195) q[1];
sx q[1];
rz(-1.4391876) q[1];
sx q[1];
rz(0.69951397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61460068) q[0];
sx q[0];
rz(-1.8934947) q[0];
sx q[0];
rz(-1.931577) q[0];
x q[1];
rz(1.5753463) q[2];
sx q[2];
rz(-2.1855436) q[2];
sx q[2];
rz(0.67753032) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2809063) q[1];
sx q[1];
rz(-1.6526892) q[1];
sx q[1];
rz(2.0211092) q[1];
x q[2];
rz(-2.2992976) q[3];
sx q[3];
rz(-2.1811322) q[3];
sx q[3];
rz(-1.6323324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.37381441) q[2];
sx q[2];
rz(-0.57743293) q[2];
sx q[2];
rz(2.2334297) q[2];
rz(-2.5670037) q[3];
sx q[3];
rz(-1.2535973) q[3];
sx q[3];
rz(-0.38578924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3278367) q[0];
sx q[0];
rz(-1.8151374) q[0];
sx q[0];
rz(0.5089708) q[0];
rz(-0.058188997) q[1];
sx q[1];
rz(-1.6845082) q[1];
sx q[1];
rz(2.0740654) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71004407) q[0];
sx q[0];
rz(-0.87418927) q[0];
sx q[0];
rz(2.0569274) q[0];
rz(-1.9827438) q[2];
sx q[2];
rz(-1.5016404) q[2];
sx q[2];
rz(-2.1092887) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6696251) q[1];
sx q[1];
rz(-1.4635007) q[1];
sx q[1];
rz(-1.7194242) q[1];
rz(-0.068693585) q[3];
sx q[3];
rz(-1.8634081) q[3];
sx q[3];
rz(0.13394395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1099781) q[2];
sx q[2];
rz(-1.5441394) q[2];
sx q[2];
rz(-2.1515089) q[2];
rz(1.7456985) q[3];
sx q[3];
rz(-3.0315704) q[3];
sx q[3];
rz(1.8949738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9083967) q[0];
sx q[0];
rz(-1.3968503) q[0];
sx q[0];
rz(-1.0594077) q[0];
rz(1.1147095) q[1];
sx q[1];
rz(-1.6672986) q[1];
sx q[1];
rz(1.4310736) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7641033) q[0];
sx q[0];
rz(-1.7400578) q[0];
sx q[0];
rz(2.9891123) q[0];
rz(-pi) q[1];
rz(1.0658979) q[2];
sx q[2];
rz(-1.3422728) q[2];
sx q[2];
rz(1.8375719) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.614735) q[1];
sx q[1];
rz(-2.4522374) q[1];
sx q[1];
rz(-2.7562642) q[1];
rz(-pi) q[2];
rz(-1.4214755) q[3];
sx q[3];
rz(-1.7826555) q[3];
sx q[3];
rz(2.6184591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7913671) q[2];
sx q[2];
rz(-2.3855049) q[2];
sx q[2];
rz(1.4678601) q[2];
rz(1.0427467) q[3];
sx q[3];
rz(-2.0839033) q[3];
sx q[3];
rz(0.79152542) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70822155) q[0];
sx q[0];
rz(-1.8480166) q[0];
sx q[0];
rz(-0.21155393) q[0];
rz(-2.4987706) q[1];
sx q[1];
rz(-2.0170409) q[1];
sx q[1];
rz(-1.0483673) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12473561) q[0];
sx q[0];
rz(-0.87364679) q[0];
sx q[0];
rz(1.7378896) q[0];
rz(-2.8564212) q[2];
sx q[2];
rz(-0.18643555) q[2];
sx q[2];
rz(0.16229072) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0901186) q[1];
sx q[1];
rz(-1.5940545) q[1];
sx q[1];
rz(0.6384797) q[1];
rz(-0.048158483) q[3];
sx q[3];
rz(-0.95880303) q[3];
sx q[3];
rz(2.4186866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.50531498) q[2];
sx q[2];
rz(-1.9376398) q[2];
sx q[2];
rz(2.8875276) q[2];
rz(0.17357477) q[3];
sx q[3];
rz(-0.55145276) q[3];
sx q[3];
rz(-1.9614296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59488615) q[0];
sx q[0];
rz(-2.7017024) q[0];
sx q[0];
rz(-0.79676262) q[0];
rz(-0.88675371) q[1];
sx q[1];
rz(-1.3462857) q[1];
sx q[1];
rz(-0.60060445) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2830505) q[0];
sx q[0];
rz(-1.5158537) q[0];
sx q[0];
rz(-0.29015975) q[0];
x q[1];
rz(-1.4378635) q[2];
sx q[2];
rz(-0.80604751) q[2];
sx q[2];
rz(1.9423167) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.826304) q[1];
sx q[1];
rz(-2.4361594) q[1];
sx q[1];
rz(1.822173) q[1];
rz(-pi) q[2];
rz(1.0135456) q[3];
sx q[3];
rz(-0.60769121) q[3];
sx q[3];
rz(0.26085873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0299006) q[2];
sx q[2];
rz(-1.4375968) q[2];
sx q[2];
rz(-2.0491484) q[2];
rz(-1.7732636) q[3];
sx q[3];
rz(-2.5320801) q[3];
sx q[3];
rz(-0.4121367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71548897) q[0];
sx q[0];
rz(-2.0273209) q[0];
sx q[0];
rz(2.6901167) q[0];
rz(0.077797912) q[1];
sx q[1];
rz(-1.0323689) q[1];
sx q[1];
rz(2.480004) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5104354) q[0];
sx q[0];
rz(-1.9826188) q[0];
sx q[0];
rz(-2.9020082) q[0];
x q[1];
rz(-1.6986474) q[2];
sx q[2];
rz(-2.6397954) q[2];
sx q[2];
rz(1.0867829) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.80967951) q[1];
sx q[1];
rz(-0.47985489) q[1];
sx q[1];
rz(-1.859568) q[1];
rz(-pi) q[2];
rz(-2.745146) q[3];
sx q[3];
rz(-1.8155193) q[3];
sx q[3];
rz(0.94371599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0563858) q[2];
sx q[2];
rz(-1.2678009) q[2];
sx q[2];
rz(0.80643225) q[2];
rz(-1.8993529) q[3];
sx q[3];
rz(-0.16568383) q[3];
sx q[3];
rz(3.0012259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8388782) q[0];
sx q[0];
rz(-1.1352204) q[0];
sx q[0];
rz(0.43933991) q[0];
rz(1.9413403) q[1];
sx q[1];
rz(-2.3019583) q[1];
sx q[1];
rz(0.95019788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5419358) q[0];
sx q[0];
rz(-2.1996351) q[0];
sx q[0];
rz(0.46577439) q[0];
rz(-pi) q[1];
rz(0.15301159) q[2];
sx q[2];
rz(-1.8540314) q[2];
sx q[2];
rz(0.52545122) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9199782) q[1];
sx q[1];
rz(-1.6755381) q[1];
sx q[1];
rz(1.2135189) q[1];
rz(-pi) q[2];
rz(-2.1860649) q[3];
sx q[3];
rz(-2.5078109) q[3];
sx q[3];
rz(-1.2891226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.83303678) q[2];
sx q[2];
rz(-2.6335282) q[2];
sx q[2];
rz(1.1888095) q[2];
rz(-3.0360119) q[3];
sx q[3];
rz(-0.9413541) q[3];
sx q[3];
rz(2.1281435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.062716) q[0];
sx q[0];
rz(-2.833241) q[0];
sx q[0];
rz(2.4421332) q[0];
rz(1.4106916) q[1];
sx q[1];
rz(-2.3532093) q[1];
sx q[1];
rz(2.9404822) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5749436) q[0];
sx q[0];
rz(-2.0127949) q[0];
sx q[0];
rz(0.7676564) q[0];
x q[1];
rz(-0.91014782) q[2];
sx q[2];
rz(-2.4322699) q[2];
sx q[2];
rz(2.7924479) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8751018) q[1];
sx q[1];
rz(-1.3492341) q[1];
sx q[1];
rz(-2.2699577) q[1];
x q[2];
rz(2.6789959) q[3];
sx q[3];
rz(-0.2004846) q[3];
sx q[3];
rz(1.6311262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.19645277) q[2];
sx q[2];
rz(-1.8818734) q[2];
sx q[2];
rz(0.067616612) q[2];
rz(2.0130646) q[3];
sx q[3];
rz(-1.0849378) q[3];
sx q[3];
rz(1.6245925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6832798) q[0];
sx q[0];
rz(-1.0673609) q[0];
sx q[0];
rz(1.3491032) q[0];
rz(1.6364527) q[1];
sx q[1];
rz(-0.22947336) q[1];
sx q[1];
rz(-0.089692399) q[1];
rz(2.5386621) q[2];
sx q[2];
rz(-1.3984192) q[2];
sx q[2];
rz(2.0655823) q[2];
rz(-2.454461) q[3];
sx q[3];
rz(-2.7239068) q[3];
sx q[3];
rz(2.1129029) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
