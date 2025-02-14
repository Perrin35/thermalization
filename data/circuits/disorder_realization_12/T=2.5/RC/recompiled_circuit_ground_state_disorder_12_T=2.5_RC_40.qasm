OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1801017) q[0];
sx q[0];
rz(-0.56668007) q[0];
sx q[0];
rz(2.3508747) q[0];
rz(2.2136731) q[1];
sx q[1];
rz(3.1766422) q[1];
sx q[1];
rz(9.1215134) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7105616) q[0];
sx q[0];
rz(-0.91791422) q[0];
sx q[0];
rz(1.8588572) q[0];
x q[1];
rz(-2.319807) q[2];
sx q[2];
rz(-2.6394834) q[2];
sx q[2];
rz(-0.443845) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.11285148) q[1];
sx q[1];
rz(-2.0926706) q[1];
sx q[1];
rz(-0.5571837) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.73533006) q[3];
sx q[3];
rz(-1.6323252) q[3];
sx q[3];
rz(1.3593591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.37320331) q[2];
sx q[2];
rz(-2.781227) q[2];
sx q[2];
rz(0.68161327) q[2];
rz(-0.56053376) q[3];
sx q[3];
rz(-2.3571641) q[3];
sx q[3];
rz(2.4973629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19175567) q[0];
sx q[0];
rz(-0.90832174) q[0];
sx q[0];
rz(-2.7857696) q[0];
rz(-2.8882354) q[1];
sx q[1];
rz(-0.8496049) q[1];
sx q[1];
rz(-1.8960948) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49911526) q[0];
sx q[0];
rz(-2.0428951) q[0];
sx q[0];
rz(-1.8224707) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0456845) q[2];
sx q[2];
rz(-1.0042937) q[2];
sx q[2];
rz(-2.1691466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.624234) q[1];
sx q[1];
rz(-1.1453724) q[1];
sx q[1];
rz(-2.5555771) q[1];
rz(-pi) q[2];
rz(0.83578556) q[3];
sx q[3];
rz(-2.4525149) q[3];
sx q[3];
rz(-0.101735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7443098) q[2];
sx q[2];
rz(-0.77703589) q[2];
sx q[2];
rz(-3.0139319) q[2];
rz(-2.4658261) q[3];
sx q[3];
rz(-2.6830169) q[3];
sx q[3];
rz(2.7202386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449988) q[0];
sx q[0];
rz(-2.8969722) q[0];
sx q[0];
rz(1.1340207) q[0];
rz(-2.2484312) q[1];
sx q[1];
rz(-0.90407073) q[1];
sx q[1];
rz(1.8697416) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29933481) q[0];
sx q[0];
rz(-1.0082196) q[0];
sx q[0];
rz(-3.1104053) q[0];
rz(-pi) q[1];
rz(-1.3273622) q[2];
sx q[2];
rz(-0.31844246) q[2];
sx q[2];
rz(1.3322347) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0772897) q[1];
sx q[1];
rz(-1.3659119) q[1];
sx q[1];
rz(-0.85475342) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9748114) q[3];
sx q[3];
rz(-1.7875009) q[3];
sx q[3];
rz(3.0971017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.21358718) q[2];
sx q[2];
rz(-0.86557937) q[2];
sx q[2];
rz(1.8013087) q[2];
rz(-1.2876997) q[3];
sx q[3];
rz(-1.6940593) q[3];
sx q[3];
rz(-2.652216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5260148) q[0];
sx q[0];
rz(-2.5647793) q[0];
sx q[0];
rz(-2.4082129) q[0];
rz(-2.3461657) q[1];
sx q[1];
rz(-0.19618244) q[1];
sx q[1];
rz(-1.9733852) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0709658) q[0];
sx q[0];
rz(-2.7160108) q[0];
sx q[0];
rz(-1.312567) q[0];
x q[1];
rz(-0.85251405) q[2];
sx q[2];
rz(-0.77272431) q[2];
sx q[2];
rz(2.1022405) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8454947) q[1];
sx q[1];
rz(-1.5278399) q[1];
sx q[1];
rz(1.5404433) q[1];
rz(-pi) q[2];
rz(-2.4766604) q[3];
sx q[3];
rz(-1.0626317) q[3];
sx q[3];
rz(2.6817961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1802133) q[2];
sx q[2];
rz(-0.81575477) q[2];
sx q[2];
rz(-0.95480314) q[2];
rz(1.0500326) q[3];
sx q[3];
rz(-0.78767109) q[3];
sx q[3];
rz(2.5900456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0410205) q[0];
sx q[0];
rz(-0.53871483) q[0];
sx q[0];
rz(0.69514298) q[0];
rz(1.2160542) q[1];
sx q[1];
rz(-1.0168409) q[1];
sx q[1];
rz(3.1207808) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1095525) q[0];
sx q[0];
rz(-1.4281245) q[0];
sx q[0];
rz(-0.90937664) q[0];
rz(2.8146652) q[2];
sx q[2];
rz(-0.86092868) q[2];
sx q[2];
rz(-2.8500789) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5095651) q[1];
sx q[1];
rz(-0.23938017) q[1];
sx q[1];
rz(-1.521895) q[1];
rz(-1.6324051) q[3];
sx q[3];
rz(-1.6635259) q[3];
sx q[3];
rz(-1.9237067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1702599) q[2];
sx q[2];
rz(-0.33997619) q[2];
sx q[2];
rz(2.5374832) q[2];
rz(0.64770925) q[3];
sx q[3];
rz(-0.52102399) q[3];
sx q[3];
rz(-0.46085301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0172417) q[0];
sx q[0];
rz(-2.0483973) q[0];
sx q[0];
rz(2.7830615) q[0];
rz(2.5784946) q[1];
sx q[1];
rz(-0.9022572) q[1];
sx q[1];
rz(-2.8360352) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76680032) q[0];
sx q[0];
rz(-1.9309485) q[0];
sx q[0];
rz(2.4712579) q[0];
rz(-pi) q[1];
x q[1];
rz(0.18230171) q[2];
sx q[2];
rz(-1.7392225) q[2];
sx q[2];
rz(1.4252848) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.36781947) q[1];
sx q[1];
rz(-1.0350643) q[1];
sx q[1];
rz(2.1526835) q[1];
rz(-0.73718585) q[3];
sx q[3];
rz(-3.0372826) q[3];
sx q[3];
rz(0.3462458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.87865889) q[2];
sx q[2];
rz(-2.1808251) q[2];
sx q[2];
rz(-3.089454) q[2];
rz(-2.7122998) q[3];
sx q[3];
rz(-1.7310127) q[3];
sx q[3];
rz(-0.2521635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6843863) q[0];
sx q[0];
rz(-0.52427137) q[0];
sx q[0];
rz(0.23505178) q[0];
rz(-0.60335195) q[1];
sx q[1];
rz(-0.82913202) q[1];
sx q[1];
rz(0.57650173) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.910945) q[0];
sx q[0];
rz(-2.0752788) q[0];
sx q[0];
rz(2.3672558) q[0];
rz(0.34502657) q[2];
sx q[2];
rz(-2.3357311) q[2];
sx q[2];
rz(2.8553183) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.602142) q[1];
sx q[1];
rz(-1.8906525) q[1];
sx q[1];
rz(1.1509239) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1436809) q[3];
sx q[3];
rz(-2.0485123) q[3];
sx q[3];
rz(2.9766072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.06360402) q[2];
sx q[2];
rz(-0.84090191) q[2];
sx q[2];
rz(1.3464751) q[2];
rz(-0.50955647) q[3];
sx q[3];
rz(-1.8815123) q[3];
sx q[3];
rz(2.8349561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52617306) q[0];
sx q[0];
rz(-1.0831447) q[0];
sx q[0];
rz(-2.7857067) q[0];
rz(-2.4399759) q[1];
sx q[1];
rz(-2.9569148) q[1];
sx q[1];
rz(0.96431771) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50011694) q[0];
sx q[0];
rz(-1.7044412) q[0];
sx q[0];
rz(1.1125314) q[0];
rz(-pi) q[1];
rz(1.6379328) q[2];
sx q[2];
rz(-0.29342857) q[2];
sx q[2];
rz(2.1591883) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0810534) q[1];
sx q[1];
rz(-2.6463228) q[1];
sx q[1];
rz(-0.041306382) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6526493) q[3];
sx q[3];
rz(-1.2602031) q[3];
sx q[3];
rz(-1.6037343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.648916) q[2];
sx q[2];
rz(-2.4339088) q[2];
sx q[2];
rz(-0.81563449) q[2];
rz(0.91551578) q[3];
sx q[3];
rz(-2.2736277) q[3];
sx q[3];
rz(0.22667949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30064073) q[0];
sx q[0];
rz(-2.9673567) q[0];
sx q[0];
rz(-1.9006282) q[0];
rz(-3.0366483) q[1];
sx q[1];
rz(-2.8009156) q[1];
sx q[1];
rz(-0.45936432) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7442087) q[0];
sx q[0];
rz(-0.43699118) q[0];
sx q[0];
rz(2.1676012) q[0];
rz(-1.2392339) q[2];
sx q[2];
rz(-1.3156097) q[2];
sx q[2];
rz(0.58411264) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33973673) q[1];
sx q[1];
rz(-1.1366399) q[1];
sx q[1];
rz(-1.804002) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4308484) q[3];
sx q[3];
rz(-1.0010747) q[3];
sx q[3];
rz(2.8960814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0910054) q[2];
sx q[2];
rz(-0.688474) q[2];
sx q[2];
rz(0.090593226) q[2];
rz(-0.76720864) q[3];
sx q[3];
rz(-1.6588666) q[3];
sx q[3];
rz(1.4513133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7200658) q[0];
sx q[0];
rz(-0.92362112) q[0];
sx q[0];
rz(-1.5371171) q[0];
rz(-2.1859258) q[1];
sx q[1];
rz(-1.4612528) q[1];
sx q[1];
rz(2.59424) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.065009251) q[0];
sx q[0];
rz(-0.07993789) q[0];
sx q[0];
rz(-1.4968833) q[0];
x q[1];
rz(-0.53031594) q[2];
sx q[2];
rz(-1.2387453) q[2];
sx q[2];
rz(0.9584934) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.66710752) q[1];
sx q[1];
rz(-0.82208668) q[1];
sx q[1];
rz(-0.32331031) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0738952) q[3];
sx q[3];
rz(-2.2184531) q[3];
sx q[3];
rz(1.5941509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1110111) q[2];
sx q[2];
rz(-2.9647398) q[2];
sx q[2];
rz(-0.64876968) q[2];
rz(-0.23160058) q[3];
sx q[3];
rz(-0.88445556) q[3];
sx q[3];
rz(-3.054936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34407525) q[0];
sx q[0];
rz(-1.2114914) q[0];
sx q[0];
rz(-1.399566) q[0];
rz(2.4823785) q[1];
sx q[1];
rz(-1.8623687) q[1];
sx q[1];
rz(1.6946793) q[1];
rz(-1.7991015) q[2];
sx q[2];
rz(-1.5226428) q[2];
sx q[2];
rz(-1.5626283) q[2];
rz(1.9271196) q[3];
sx q[3];
rz(-0.15542843) q[3];
sx q[3];
rz(-0.036416362) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
