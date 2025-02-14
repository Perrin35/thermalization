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
rz(3.0577793) q[0];
sx q[0];
rz(-0.35141355) q[0];
sx q[0];
rz(1.0454398) q[0];
rz(0.84402973) q[1];
sx q[1];
rz(-2.5764155) q[1];
sx q[1];
rz(-1.7796109) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0365252) q[0];
sx q[0];
rz(-1.0191259) q[0];
sx q[0];
rz(1.9776634) q[0];
rz(0.092081618) q[2];
sx q[2];
rz(-1.144334) q[2];
sx q[2];
rz(1.2464166) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2842362) q[1];
sx q[1];
rz(-1.7928837) q[1];
sx q[1];
rz(2.7276993) q[1];
rz(-pi) q[2];
rz(-2.8027439) q[3];
sx q[3];
rz(-2.04755) q[3];
sx q[3];
rz(-2.2173405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7094946) q[2];
sx q[2];
rz(-0.39958909) q[2];
sx q[2];
rz(2.438365) q[2];
rz(-2.382459) q[3];
sx q[3];
rz(-0.96733624) q[3];
sx q[3];
rz(2.4659992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264939) q[0];
sx q[0];
rz(-0.7248942) q[0];
sx q[0];
rz(0.79459992) q[0];
rz(-0.82851797) q[1];
sx q[1];
rz(-2.0683894) q[1];
sx q[1];
rz(0.4812831) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1409765) q[0];
sx q[0];
rz(-3.0191026) q[0];
sx q[0];
rz(0.23167892) q[0];
rz(-pi) q[1];
rz(2.6120466) q[2];
sx q[2];
rz(-0.82223071) q[2];
sx q[2];
rz(1.032743) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.1294205) q[1];
sx q[1];
rz(-0.24298619) q[1];
sx q[1];
rz(1.245973) q[1];
rz(-pi) q[2];
rz(2.4614905) q[3];
sx q[3];
rz(-0.77723336) q[3];
sx q[3];
rz(-0.2118563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1543697) q[2];
sx q[2];
rz(-0.96174806) q[2];
sx q[2];
rz(-2.6340384) q[2];
rz(0.65732035) q[3];
sx q[3];
rz(-2.1828914) q[3];
sx q[3];
rz(-1.9067732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2051314) q[0];
sx q[0];
rz(-1.6070123) q[0];
sx q[0];
rz(-0.24612799) q[0];
rz(-0.71324619) q[1];
sx q[1];
rz(-1.5245707) q[1];
sx q[1];
rz(-2.2770142) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.789398) q[0];
sx q[0];
rz(-1.3653132) q[0];
sx q[0];
rz(1.2557538) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6106493) q[2];
sx q[2];
rz(-2.1111672) q[2];
sx q[2];
rz(-0.8671538) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.1684226) q[1];
sx q[1];
rz(-1.4900786) q[1];
sx q[1];
rz(-1.1087085) q[1];
x q[2];
rz(0.15154408) q[3];
sx q[3];
rz(-2.3484548) q[3];
sx q[3];
rz(-2.0878937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41512179) q[2];
sx q[2];
rz(-0.39602009) q[2];
sx q[2];
rz(1.2064639) q[2];
rz(0.11780277) q[3];
sx q[3];
rz(-0.67540568) q[3];
sx q[3];
rz(-3.0142768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29663157) q[0];
sx q[0];
rz(-1.4609818) q[0];
sx q[0];
rz(0.46517459) q[0];
rz(2.2370715) q[1];
sx q[1];
rz(-2.1280839) q[1];
sx q[1];
rz(1.7101589) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97665107) q[0];
sx q[0];
rz(-1.8969844) q[0];
sx q[0];
rz(1.1378764) q[0];
x q[1];
rz(-1.2063556) q[2];
sx q[2];
rz(-2.6195457) q[2];
sx q[2];
rz(-2.4625157) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8697813) q[1];
sx q[1];
rz(-1.2220739) q[1];
sx q[1];
rz(0.70088245) q[1];
rz(-pi) q[2];
rz(1.9172098) q[3];
sx q[3];
rz(-0.87408057) q[3];
sx q[3];
rz(-0.55075607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.64100921) q[2];
sx q[2];
rz(-0.21037978) q[2];
sx q[2];
rz(1.6453936) q[2];
rz(3.1032041) q[3];
sx q[3];
rz(-1.3209141) q[3];
sx q[3];
rz(-0.7588318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1372304) q[0];
sx q[0];
rz(-0.69843233) q[0];
sx q[0];
rz(1.5414365) q[0];
rz(-1.6372708) q[1];
sx q[1];
rz(-1.5802822) q[1];
sx q[1];
rz(1.6391594) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3454895) q[0];
sx q[0];
rz(-2.0292583) q[0];
sx q[0];
rz(-1.0572516) q[0];
rz(-pi) q[1];
rz(-0.2404332) q[2];
sx q[2];
rz(-2.5747262) q[2];
sx q[2];
rz(0.64449233) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1790601) q[1];
sx q[1];
rz(-2.0174562) q[1];
sx q[1];
rz(1.0960177) q[1];
rz(-1.6385025) q[3];
sx q[3];
rz(-1.4622524) q[3];
sx q[3];
rz(-0.53311611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0043103546) q[2];
sx q[2];
rz(-2.8305125) q[2];
sx q[2];
rz(-0.38055554) q[2];
rz(2.6015094) q[3];
sx q[3];
rz(-1.8916847) q[3];
sx q[3];
rz(-1.165747) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49749097) q[0];
sx q[0];
rz(-0.044450132) q[0];
sx q[0];
rz(-0.38462001) q[0];
rz(-0.050845536) q[1];
sx q[1];
rz(-0.47639242) q[1];
sx q[1];
rz(-1.4643889) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58005262) q[0];
sx q[0];
rz(-2.0140352) q[0];
sx q[0];
rz(2.0361986) q[0];
rz(-pi) q[1];
rz(2.5821677) q[2];
sx q[2];
rz(-2.17387) q[2];
sx q[2];
rz(-0.81808264) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.44838833) q[1];
sx q[1];
rz(-0.12388661) q[1];
sx q[1];
rz(-2.3883567) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7767653) q[3];
sx q[3];
rz(-2.245083) q[3];
sx q[3];
rz(-0.29335653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9836318) q[2];
sx q[2];
rz(-2.3775358) q[2];
sx q[2];
rz(-2.2264886) q[2];
rz(-2.8369331) q[3];
sx q[3];
rz(-0.72158486) q[3];
sx q[3];
rz(1.6662395) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7810998) q[0];
sx q[0];
rz(-2.2486794) q[0];
sx q[0];
rz(-2.0258946) q[0];
rz(0.64104331) q[1];
sx q[1];
rz(-1.8843001) q[1];
sx q[1];
rz(1.3263652) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5005701) q[0];
sx q[0];
rz(-1.7020853) q[0];
sx q[0];
rz(2.6528051) q[0];
rz(0.4037598) q[2];
sx q[2];
rz(-0.40382622) q[2];
sx q[2];
rz(2.5346839) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7515727) q[1];
sx q[1];
rz(-2.2604001) q[1];
sx q[1];
rz(-0.71038664) q[1];
x q[2];
rz(0.85878813) q[3];
sx q[3];
rz(-0.988171) q[3];
sx q[3];
rz(1.7395541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.556584) q[2];
sx q[2];
rz(-0.98238397) q[2];
sx q[2];
rz(-1.2674241) q[2];
rz(3.0804539) q[3];
sx q[3];
rz(-1.2631402) q[3];
sx q[3];
rz(-2.256573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1214445) q[0];
sx q[0];
rz(-2.3763438) q[0];
sx q[0];
rz(2.6159317) q[0];
rz(-1.6540182) q[1];
sx q[1];
rz(-1.1143538) q[1];
sx q[1];
rz(0.18255998) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2935242) q[0];
sx q[0];
rz(-1.7901359) q[0];
sx q[0];
rz(0.42869403) q[0];
rz(-2.0261954) q[2];
sx q[2];
rz(-1.2898852) q[2];
sx q[2];
rz(1.8473491) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5044587) q[1];
sx q[1];
rz(-0.65385039) q[1];
sx q[1];
rz(1.8549396) q[1];
rz(-pi) q[2];
rz(1.0973516) q[3];
sx q[3];
rz(-1.3149188) q[3];
sx q[3];
rz(2.8476343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.93099) q[2];
sx q[2];
rz(-1.3347722) q[2];
sx q[2];
rz(-0.1087428) q[2];
rz(-0.10793081) q[3];
sx q[3];
rz(-0.35522541) q[3];
sx q[3];
rz(-1.7598553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0906319) q[0];
sx q[0];
rz(-1.9405631) q[0];
sx q[0];
rz(2.4884124) q[0];
rz(1.8494891) q[1];
sx q[1];
rz(-0.2457681) q[1];
sx q[1];
rz(-0.076315708) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.171577) q[0];
sx q[0];
rz(-1.2467915) q[0];
sx q[0];
rz(0.25361639) q[0];
rz(-pi) q[1];
rz(1.7387975) q[2];
sx q[2];
rz(-1.6036878) q[2];
sx q[2];
rz(2.5891182) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6725823) q[1];
sx q[1];
rz(-0.82090506) q[1];
sx q[1];
rz(-1.7609079) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1006484) q[3];
sx q[3];
rz(-0.68504928) q[3];
sx q[3];
rz(2.4026826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2448347) q[2];
sx q[2];
rz(-1.1847757) q[2];
sx q[2];
rz(0.99311382) q[2];
rz(-1.7874329) q[3];
sx q[3];
rz(-0.5391776) q[3];
sx q[3];
rz(1.7206934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3058474) q[0];
sx q[0];
rz(-1.638224) q[0];
sx q[0];
rz(0.29356965) q[0];
rz(-1.0319895) q[1];
sx q[1];
rz(-0.76621619) q[1];
sx q[1];
rz(2.8242677) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85451305) q[0];
sx q[0];
rz(-1.6995158) q[0];
sx q[0];
rz(-2.5230056) q[0];
x q[1];
rz(-1.3506838) q[2];
sx q[2];
rz(-2.1145472) q[2];
sx q[2];
rz(-1.3771283) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.36655246) q[1];
sx q[1];
rz(-1.9138819) q[1];
sx q[1];
rz(-0.45963414) q[1];
rz(1.6280319) q[3];
sx q[3];
rz(-1.8678819) q[3];
sx q[3];
rz(0.12606584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5130634) q[2];
sx q[2];
rz(-0.3231914) q[2];
sx q[2];
rz(2.826706) q[2];
rz(-3.1079187) q[3];
sx q[3];
rz(-2.0458872) q[3];
sx q[3];
rz(-1.7033738) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(2.2122129) q[0];
sx q[0];
rz(-2.0912981) q[0];
sx q[0];
rz(-0.97846497) q[0];
rz(2.9411511) q[1];
sx q[1];
rz(-1.0841752) q[1];
sx q[1];
rz(-0.60842327) q[1];
rz(-0.54666109) q[2];
sx q[2];
rz(-0.056686747) q[2];
sx q[2];
rz(-2.0591339) q[2];
rz(-1.325812) q[3];
sx q[3];
rz(-2.1143338) q[3];
sx q[3];
rz(-0.59296617) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
