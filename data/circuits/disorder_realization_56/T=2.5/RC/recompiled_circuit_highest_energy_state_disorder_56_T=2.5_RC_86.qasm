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
rz(4.4895953) q[0];
sx q[0];
rz(9.0564981) q[0];
rz(3.0556196) q[1];
sx q[1];
rz(-0.84288016) q[1];
sx q[1];
rz(2.892363) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0195784) q[0];
sx q[0];
rz(-1.459157) q[0];
sx q[0];
rz(0.48099244) q[0];
rz(-pi) q[1];
rz(2.8299242) q[2];
sx q[2];
rz(-2.083792) q[2];
sx q[2];
rz(-2.5530961) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.8650444) q[1];
sx q[1];
rz(-0.88403349) q[1];
sx q[1];
rz(1.264132) q[1];
rz(-pi) q[2];
rz(-0.29549142) q[3];
sx q[3];
rz(-0.37335941) q[3];
sx q[3];
rz(-0.96366027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.92980802) q[2];
sx q[2];
rz(-2.4346508) q[2];
sx q[2];
rz(-2.9014273) q[2];
rz(0.54720488) q[3];
sx q[3];
rz(-1.6271084) q[3];
sx q[3];
rz(0.98946324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3419679) q[0];
sx q[0];
rz(-0.37536055) q[0];
sx q[0];
rz(0.64390916) q[0];
rz(1.0470692) q[1];
sx q[1];
rz(-1.964485) q[1];
sx q[1];
rz(-0.26328304) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8174648) q[0];
sx q[0];
rz(-3.1006515) q[0];
sx q[0];
rz(-0.73887478) q[0];
rz(-3.075454) q[2];
sx q[2];
rz(-1.5060584) q[2];
sx q[2];
rz(-1.7867225) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66052478) q[1];
sx q[1];
rz(-2.1612289) q[1];
sx q[1];
rz(0.51312889) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1016779) q[3];
sx q[3];
rz(-0.71456075) q[3];
sx q[3];
rz(-1.8811864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6836267) q[2];
sx q[2];
rz(-1.7650975) q[2];
sx q[2];
rz(-2.4158884) q[2];
rz(-2.8343685) q[3];
sx q[3];
rz(-1.3064462) q[3];
sx q[3];
rz(-1.0544302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3211408) q[0];
sx q[0];
rz(-1.6571925) q[0];
sx q[0];
rz(-0.76903525) q[0];
rz(-3.0342195) q[1];
sx q[1];
rz(-1.4391876) q[1];
sx q[1];
rz(-0.69951397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8843667) q[0];
sx q[0];
rz(-0.47927893) q[0];
sx q[0];
rz(-2.3291161) q[0];
rz(-pi) q[1];
rz(3.1351481) q[2];
sx q[2];
rz(-2.5268308) q[2];
sx q[2];
rz(-2.4561735) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2809063) q[1];
sx q[1];
rz(-1.6526892) q[1];
sx q[1];
rz(2.0211092) q[1];
rz(0.75306692) q[3];
sx q[3];
rz(-0.99374607) q[3];
sx q[3];
rz(-2.6073539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37381441) q[2];
sx q[2];
rz(-0.57743293) q[2];
sx q[2];
rz(2.2334297) q[2];
rz(2.5670037) q[3];
sx q[3];
rz(-1.8879954) q[3];
sx q[3];
rz(-0.38578924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.813756) q[0];
sx q[0];
rz(-1.3264553) q[0];
sx q[0];
rz(-0.5089708) q[0];
rz(-0.058188997) q[1];
sx q[1];
rz(-1.4570844) q[1];
sx q[1];
rz(1.0675272) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4315486) q[0];
sx q[0];
rz(-2.2674034) q[0];
sx q[0];
rz(-2.0569274) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7420962) q[2];
sx q[2];
rz(-0.41738415) q[2];
sx q[2];
rz(0.69533747) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4219027) q[1];
sx q[1];
rz(-2.9585144) q[1];
sx q[1];
rz(0.94193926) q[1];
rz(-pi) q[2];
rz(-1.3467784) q[3];
sx q[3];
rz(-0.30034143) q[3];
sx q[3];
rz(-3.0413922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.031614583) q[2];
sx q[2];
rz(-1.5441394) q[2];
sx q[2];
rz(-2.1515089) q[2];
rz(-1.7456985) q[3];
sx q[3];
rz(-3.0315704) q[3];
sx q[3];
rz(-1.8949738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9083967) q[0];
sx q[0];
rz(-1.3968503) q[0];
sx q[0];
rz(-2.0821849) q[0];
rz(2.0268832) q[1];
sx q[1];
rz(-1.474294) q[1];
sx q[1];
rz(1.4310736) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9741668) q[0];
sx q[0];
rz(-1.7210809) q[0];
sx q[0];
rz(1.3995863) q[0];
rz(-2.8818508) q[2];
sx q[2];
rz(-1.0802104) q[2];
sx q[2];
rz(0.39133137) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0986745) q[1];
sx q[1];
rz(-0.94047767) q[1];
sx q[1];
rz(-1.2703671) q[1];
rz(-pi) q[2];
rz(2.5364872) q[3];
sx q[3];
rz(-0.25854585) q[3];
sx q[3];
rz(3.0437146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7913671) q[2];
sx q[2];
rz(-2.3855049) q[2];
sx q[2];
rz(-1.4678601) q[2];
rz(2.098846) q[3];
sx q[3];
rz(-2.0839033) q[3];
sx q[3];
rz(2.3500672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70822155) q[0];
sx q[0];
rz(-1.2935761) q[0];
sx q[0];
rz(0.21155393) q[0];
rz(0.64282203) q[1];
sx q[1];
rz(-1.1245518) q[1];
sx q[1];
rz(1.0483673) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7599568) q[0];
sx q[0];
rz(-0.71361976) q[0];
sx q[0];
rz(2.9455393) q[0];
rz(-pi) q[1];
rz(2.9625234) q[2];
sx q[2];
rz(-1.5186276) q[2];
sx q[2];
rz(1.6889926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6395289) q[1];
sx q[1];
rz(-2.2090753) q[1];
sx q[1];
rz(-1.5997574) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.95825742) q[3];
sx q[3];
rz(-1.6102092) q[3];
sx q[3];
rz(2.266021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.50531498) q[2];
sx q[2];
rz(-1.9376398) q[2];
sx q[2];
rz(0.25406507) q[2];
rz(-0.17357477) q[3];
sx q[3];
rz(-0.55145276) q[3];
sx q[3];
rz(-1.180163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59488615) q[0];
sx q[0];
rz(-0.43989023) q[0];
sx q[0];
rz(2.34483) q[0];
rz(0.88675371) q[1];
sx q[1];
rz(-1.795307) q[1];
sx q[1];
rz(-0.60060445) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69585872) q[0];
sx q[0];
rz(-1.2810871) q[0];
sx q[0];
rz(-1.5134619) q[0];
rz(-0.76917665) q[2];
sx q[2];
rz(-1.4750136) q[2];
sx q[2];
rz(-2.8623919) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4485885) q[1];
sx q[1];
rz(-1.4088165) q[1];
sx q[1];
rz(-0.88108351) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1280471) q[3];
sx q[3];
rz(-2.5339014) q[3];
sx q[3];
rz(-0.26085873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0299006) q[2];
sx q[2];
rz(-1.4375968) q[2];
sx q[2];
rz(1.0924443) q[2];
rz(1.7732636) q[3];
sx q[3];
rz(-0.60951257) q[3];
sx q[3];
rz(2.7294559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71548897) q[0];
sx q[0];
rz(-2.0273209) q[0];
sx q[0];
rz(-0.45147595) q[0];
rz(-3.0637947) q[1];
sx q[1];
rz(-2.1092238) q[1];
sx q[1];
rz(-2.480004) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2994227) q[0];
sx q[0];
rz(-1.3515858) q[0];
sx q[0];
rz(1.1482394) q[0];
rz(-0.069839283) q[2];
sx q[2];
rz(-1.0734715) q[2];
sx q[2];
rz(-2.2004009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6381423) q[1];
sx q[1];
rz(-1.7026445) q[1];
sx q[1];
rz(2.0335458) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5675232) q[3];
sx q[3];
rz(-0.46246734) q[3];
sx q[3];
rz(-1.9898349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0563858) q[2];
sx q[2];
rz(-1.2678009) q[2];
sx q[2];
rz(-2.3351604) q[2];
rz(-1.2422397) q[3];
sx q[3];
rz(-2.9759088) q[3];
sx q[3];
rz(3.0012259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8388782) q[0];
sx q[0];
rz(-1.1352204) q[0];
sx q[0];
rz(0.43933991) q[0];
rz(-1.2002523) q[1];
sx q[1];
rz(-2.3019583) q[1];
sx q[1];
rz(-2.1913948) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59965682) q[0];
sx q[0];
rz(-0.94195752) q[0];
sx q[0];
rz(2.6758183) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15301159) q[2];
sx q[2];
rz(-1.2875612) q[2];
sx q[2];
rz(-0.52545122) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2216144) q[1];
sx q[1];
rz(-1.6755381) q[1];
sx q[1];
rz(1.2135189) q[1];
rz(2.7404157) q[3];
sx q[3];
rz(-1.0660352) q[3];
sx q[3];
rz(1.1324319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3085559) q[2];
sx q[2];
rz(-2.6335282) q[2];
sx q[2];
rz(-1.1888095) q[2];
rz(-3.0360119) q[3];
sx q[3];
rz(-2.2002386) q[3];
sx q[3];
rz(-2.1281435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078876615) q[0];
sx q[0];
rz(-2.833241) q[0];
sx q[0];
rz(-2.4421332) q[0];
rz(-1.4106916) q[1];
sx q[1];
rz(-0.78838333) q[1];
sx q[1];
rz(-0.20111045) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7459262) q[0];
sx q[0];
rz(-2.2494083) q[0];
sx q[0];
rz(-0.9890438) q[0];
x q[1];
rz(2.1664326) q[2];
sx q[2];
rz(-1.1596408) q[2];
sx q[2];
rz(1.3871297) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.048677532) q[1];
sx q[1];
rz(-0.72775048) q[1];
sx q[1];
rz(-1.2341094) q[1];
rz(-pi) q[2];
rz(-2.6789959) q[3];
sx q[3];
rz(-0.2004846) q[3];
sx q[3];
rz(-1.6311262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9451399) q[2];
sx q[2];
rz(-1.8818734) q[2];
sx q[2];
rz(-0.067616612) q[2];
rz(-2.0130646) q[3];
sx q[3];
rz(-1.0849378) q[3];
sx q[3];
rz(-1.6245925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.6832798) q[0];
sx q[0];
rz(-2.0742317) q[0];
sx q[0];
rz(-1.7924894) q[0];
rz(-1.5051399) q[1];
sx q[1];
rz(-0.22947336) q[1];
sx q[1];
rz(-0.089692399) q[1];
rz(2.8436974) q[2];
sx q[2];
rz(-0.6241326) q[2];
sx q[2];
rz(0.2506102) q[2];
rz(-2.811089) q[3];
sx q[3];
rz(-1.8310343) q[3];
sx q[3];
rz(-1.9559947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
