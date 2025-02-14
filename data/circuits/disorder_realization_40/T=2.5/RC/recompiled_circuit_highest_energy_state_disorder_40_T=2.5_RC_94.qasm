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
rz(1.9873729) q[0];
sx q[0];
rz(4.6648751) q[0];
sx q[0];
rz(9.2821791) q[0];
rz(0.59347403) q[1];
sx q[1];
rz(3.7945336) q[1];
sx q[1];
rz(8.3795587) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77989376) q[0];
sx q[0];
rz(-2.4068659) q[0];
sx q[0];
rz(-0.72025062) q[0];
rz(0.93306834) q[2];
sx q[2];
rz(-1.5412113) q[2];
sx q[2];
rz(-2.3552908) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3821311) q[1];
sx q[1];
rz(-2.9540017) q[1];
sx q[1];
rz(1.0008903) q[1];
x q[2];
rz(1.7552283) q[3];
sx q[3];
rz(-0.57653713) q[3];
sx q[3];
rz(0.54631305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1656701) q[2];
sx q[2];
rz(-1.0079577) q[2];
sx q[2];
rz(2.9284076) q[2];
rz(-2.2255157) q[3];
sx q[3];
rz(-2.7824184) q[3];
sx q[3];
rz(0.38755125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857472) q[0];
sx q[0];
rz(-2.2540932) q[0];
sx q[0];
rz(0.47072738) q[0];
rz(-2.0384906) q[1];
sx q[1];
rz(-0.50869894) q[1];
sx q[1];
rz(-2.5618166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8033898) q[0];
sx q[0];
rz(-1.0472327) q[0];
sx q[0];
rz(-1.5016599) q[0];
rz(1.2389555) q[2];
sx q[2];
rz(-2.08925) q[2];
sx q[2];
rz(1.4084852) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6439542) q[1];
sx q[1];
rz(-2.4259287) q[1];
sx q[1];
rz(-2.1791451) q[1];
x q[2];
rz(0.26269368) q[3];
sx q[3];
rz(-1.212271) q[3];
sx q[3];
rz(2.3101447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9819928) q[2];
sx q[2];
rz(-0.85593587) q[2];
sx q[2];
rz(3.0465872) q[2];
rz(0.57560086) q[3];
sx q[3];
rz(-2.3592981) q[3];
sx q[3];
rz(-1.6949722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5461102) q[0];
sx q[0];
rz(-0.70850104) q[0];
sx q[0];
rz(0.27221671) q[0];
rz(-0.1046003) q[1];
sx q[1];
rz(-2.101254) q[1];
sx q[1];
rz(-1.6786172) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9538077) q[0];
sx q[0];
rz(-1.6208795) q[0];
sx q[0];
rz(-2.2950776) q[0];
x q[1];
rz(-0.78532312) q[2];
sx q[2];
rz(-1.5253272) q[2];
sx q[2];
rz(-2.3617705) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.453664) q[1];
sx q[1];
rz(-1.8519823) q[1];
sx q[1];
rz(-2.4635876) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6312863) q[3];
sx q[3];
rz(-1.20245) q[3];
sx q[3];
rz(2.9591536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0766803) q[2];
sx q[2];
rz(-1.7408337) q[2];
sx q[2];
rz(0.40985516) q[2];
rz(0.05154933) q[3];
sx q[3];
rz(-0.99286538) q[3];
sx q[3];
rz(-0.73369098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16875295) q[0];
sx q[0];
rz(-0.77744716) q[0];
sx q[0];
rz(0.69946104) q[0];
rz(0.93060023) q[1];
sx q[1];
rz(-1.7761296) q[1];
sx q[1];
rz(1.0994937) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4874026) q[0];
sx q[0];
rz(-2.497597) q[0];
sx q[0];
rz(1.0985159) q[0];
rz(2.1191988) q[2];
sx q[2];
rz(-1.363593) q[2];
sx q[2];
rz(2.0275379) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18106724) q[1];
sx q[1];
rz(-1.8904422) q[1];
sx q[1];
rz(-1.8962217) q[1];
rz(-pi) q[2];
x q[2];
rz(0.52249281) q[3];
sx q[3];
rz(-2.429649) q[3];
sx q[3];
rz(-0.57761907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.76116556) q[2];
sx q[2];
rz(-1.0720422) q[2];
sx q[2];
rz(2.2011444) q[2];
rz(0.56194168) q[3];
sx q[3];
rz(-0.95293957) q[3];
sx q[3];
rz(2.565062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1171653) q[0];
sx q[0];
rz(-3.0553525) q[0];
sx q[0];
rz(-1.0166136) q[0];
rz(-0.92574614) q[1];
sx q[1];
rz(-2.3912906) q[1];
sx q[1];
rz(0.076676682) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31875089) q[0];
sx q[0];
rz(-2.7116576) q[0];
sx q[0];
rz(-1.9775526) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3252649) q[2];
sx q[2];
rz(-1.7707728) q[2];
sx q[2];
rz(-0.94975805) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8342585) q[1];
sx q[1];
rz(-1.3594207) q[1];
sx q[1];
rz(1.0049051) q[1];
rz(-pi) q[2];
rz(1.2604976) q[3];
sx q[3];
rz(-2.3633133) q[3];
sx q[3];
rz(-0.33184127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3735247) q[2];
sx q[2];
rz(-0.7879476) q[2];
sx q[2];
rz(-0.12316556) q[2];
rz(-0.67433107) q[3];
sx q[3];
rz(-2.6682523) q[3];
sx q[3];
rz(0.82790747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3097565) q[0];
sx q[0];
rz(-0.18853822) q[0];
sx q[0];
rz(-2.6605666) q[0];
rz(0.24093974) q[1];
sx q[1];
rz(-0.53673458) q[1];
sx q[1];
rz(3.0066838) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9438926) q[0];
sx q[0];
rz(-1.1158082) q[0];
sx q[0];
rz(-0.83857341) q[0];
x q[1];
rz(-1.6740225) q[2];
sx q[2];
rz(-1.0380259) q[2];
sx q[2];
rz(-2.4599883) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72249216) q[1];
sx q[1];
rz(-2.3416391) q[1];
sx q[1];
rz(1.9582307) q[1];
rz(-0.44996275) q[3];
sx q[3];
rz(-1.6358161) q[3];
sx q[3];
rz(1.3715594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9350819) q[2];
sx q[2];
rz(-0.93588459) q[2];
sx q[2];
rz(1.8180465) q[2];
rz(0.27213085) q[3];
sx q[3];
rz(-0.85796732) q[3];
sx q[3];
rz(-2.9959397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018933522) q[0];
sx q[0];
rz(-0.7974565) q[0];
sx q[0];
rz(-2.4830699) q[0];
rz(1.5918484) q[1];
sx q[1];
rz(-2.1287983) q[1];
sx q[1];
rz(0.50643593) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50473266) q[0];
sx q[0];
rz(-1.2018645) q[0];
sx q[0];
rz(-1.4999267) q[0];
rz(-1.5120878) q[2];
sx q[2];
rz(-0.48248267) q[2];
sx q[2];
rz(-2.5168632) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4567007) q[1];
sx q[1];
rz(-1.646767) q[1];
sx q[1];
rz(-0.09346813) q[1];
rz(-0.5793367) q[3];
sx q[3];
rz(-1.5361353) q[3];
sx q[3];
rz(3.0140585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8062313) q[2];
sx q[2];
rz(-0.74593097) q[2];
sx q[2];
rz(1.8719505) q[2];
rz(-1.3658203) q[3];
sx q[3];
rz(-3.1277872) q[3];
sx q[3];
rz(0.55240101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7954471) q[0];
sx q[0];
rz(-0.2427635) q[0];
sx q[0];
rz(-0.41918293) q[0];
rz(-0.090713352) q[1];
sx q[1];
rz(-0.9181298) q[1];
sx q[1];
rz(-1.027164) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2614505) q[0];
sx q[0];
rz(-2.1322726) q[0];
sx q[0];
rz(-1.0358769) q[0];
x q[1];
rz(1.5628603) q[2];
sx q[2];
rz(-1.5936226) q[2];
sx q[2];
rz(-0.16032444) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.724714) q[1];
sx q[1];
rz(-1.6771349) q[1];
sx q[1];
rz(1.357253) q[1];
rz(-pi) q[2];
rz(0.054612463) q[3];
sx q[3];
rz(-0.88597466) q[3];
sx q[3];
rz(-2.7582061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1586228) q[2];
sx q[2];
rz(-2.136844) q[2];
sx q[2];
rz(2.1515004) q[2];
rz(-2.8236735) q[3];
sx q[3];
rz(-2.3293994) q[3];
sx q[3];
rz(3.1032491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7541499) q[0];
sx q[0];
rz(-0.70873547) q[0];
sx q[0];
rz(-2.5842174) q[0];
rz(2.7793461) q[1];
sx q[1];
rz(-1.7049512) q[1];
sx q[1];
rz(-0.16709669) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8734735) q[0];
sx q[0];
rz(-0.41129204) q[0];
sx q[0];
rz(-0.71843546) q[0];
x q[1];
rz(-2.9252824) q[2];
sx q[2];
rz(-0.71291332) q[2];
sx q[2];
rz(1.1110493) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5846338) q[1];
sx q[1];
rz(-1.2642197) q[1];
sx q[1];
rz(2.4584998) q[1];
x q[2];
rz(-1.6623508) q[3];
sx q[3];
rz(-2.4560438) q[3];
sx q[3];
rz(0.26622546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.96357137) q[2];
sx q[2];
rz(-0.19266291) q[2];
sx q[2];
rz(2.7131405) q[2];
rz(-2.4109449) q[3];
sx q[3];
rz(-2.0615536) q[3];
sx q[3];
rz(0.60034269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.44883248) q[0];
sx q[0];
rz(-1.6535783) q[0];
sx q[0];
rz(2.0220508) q[0];
rz(0.66850942) q[1];
sx q[1];
rz(-0.57066494) q[1];
sx q[1];
rz(2.8817435) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4762806) q[0];
sx q[0];
rz(-0.57054936) q[0];
sx q[0];
rz(-2.0670939) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2462559) q[2];
sx q[2];
rz(-2.5684212) q[2];
sx q[2];
rz(-1.2473388) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4230792) q[1];
sx q[1];
rz(-1.9155518) q[1];
sx q[1];
rz(-2.7026947) q[1];
x q[2];
rz(-0.88093984) q[3];
sx q[3];
rz(-1.4487113) q[3];
sx q[3];
rz(1.087904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4708289) q[2];
sx q[2];
rz(-1.908611) q[2];
sx q[2];
rz(-1.0864351) q[2];
rz(0.49232617) q[3];
sx q[3];
rz(-2.2958675) q[3];
sx q[3];
rz(2.4194748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.54871854) q[0];
sx q[0];
rz(-1.508779) q[0];
sx q[0];
rz(-1.4887703) q[0];
rz(1.2552352) q[1];
sx q[1];
rz(-1.8240758) q[1];
sx q[1];
rz(0.12771894) q[1];
rz(1.1716539) q[2];
sx q[2];
rz(-1.658434) q[2];
sx q[2];
rz(0.077736248) q[2];
rz(3.0782128) q[3];
sx q[3];
rz(-2.2324003) q[3];
sx q[3];
rz(1.5008055) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
