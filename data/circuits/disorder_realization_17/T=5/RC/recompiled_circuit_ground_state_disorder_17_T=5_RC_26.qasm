OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.021304) q[0];
sx q[0];
rz(-0.53417438) q[0];
sx q[0];
rz(1.0502653) q[0];
rz(1.8771111) q[1];
sx q[1];
rz(-0.85327947) q[1];
sx q[1];
rz(-2.5929911) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3162606) q[0];
sx q[0];
rz(-1.863593) q[0];
sx q[0];
rz(-1.2626075) q[0];
rz(-pi) q[1];
rz(3.110496) q[2];
sx q[2];
rz(-0.59760909) q[2];
sx q[2];
rz(-0.55127599) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42014112) q[1];
sx q[1];
rz(-0.68183091) q[1];
sx q[1];
rz(-0.89116606) q[1];
rz(1.886142) q[3];
sx q[3];
rz(-1.7996801) q[3];
sx q[3];
rz(1.6381519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41142622) q[2];
sx q[2];
rz(-2.7226518) q[2];
sx q[2];
rz(2.1298998) q[2];
rz(-2.2569979) q[3];
sx q[3];
rz(-1.9832289) q[3];
sx q[3];
rz(-1.2456892) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0678134) q[0];
sx q[0];
rz(-2.1429006) q[0];
sx q[0];
rz(-0.78773898) q[0];
rz(0.9785606) q[1];
sx q[1];
rz(-1.1509044) q[1];
sx q[1];
rz(-1.2501134) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12560168) q[0];
sx q[0];
rz(-2.0569394) q[0];
sx q[0];
rz(-0.12195964) q[0];
x q[1];
rz(2.1627126) q[2];
sx q[2];
rz(-2.667281) q[2];
sx q[2];
rz(-1.7537376) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3722222) q[1];
sx q[1];
rz(-1.5466855) q[1];
sx q[1];
rz(0.80609821) q[1];
rz(-pi) q[2];
rz(1.3504418) q[3];
sx q[3];
rz(-1.7856132) q[3];
sx q[3];
rz(-0.99220105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1193739) q[2];
sx q[2];
rz(-1.4639857) q[2];
sx q[2];
rz(3.019943) q[2];
rz(2.4308448) q[3];
sx q[3];
rz(-2.9080279) q[3];
sx q[3];
rz(2.5392883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0372593) q[0];
sx q[0];
rz(-2.4132044) q[0];
sx q[0];
rz(-0.65993586) q[0];
rz(0.51689369) q[1];
sx q[1];
rz(-0.70534244) q[1];
sx q[1];
rz(1.0924115) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4683974) q[0];
sx q[0];
rz(-1.930113) q[0];
sx q[0];
rz(2.8617489) q[0];
rz(-pi) q[1];
rz(-1.925611) q[2];
sx q[2];
rz(-0.72941226) q[2];
sx q[2];
rz(1.5340005) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9398492) q[1];
sx q[1];
rz(-2.1437316) q[1];
sx q[1];
rz(-2.4948289) q[1];
rz(-2.4928983) q[3];
sx q[3];
rz(-1.3645384) q[3];
sx q[3];
rz(1.0252531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2163781) q[2];
sx q[2];
rz(-0.34319147) q[2];
sx q[2];
rz(1.7197616) q[2];
rz(-2.6575994) q[3];
sx q[3];
rz(-1.4839987) q[3];
sx q[3];
rz(-1.4181731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5487109) q[0];
sx q[0];
rz(-3.0573248) q[0];
sx q[0];
rz(-1.3053869) q[0];
rz(3.1343754) q[1];
sx q[1];
rz(-2.9199298) q[1];
sx q[1];
rz(-0.73297393) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3696267) q[0];
sx q[0];
rz(-1.4384171) q[0];
sx q[0];
rz(0.11970971) q[0];
rz(-pi) q[1];
rz(0.29455955) q[2];
sx q[2];
rz(-1.4617051) q[2];
sx q[2];
rz(2.8677169) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.32506714) q[1];
sx q[1];
rz(-0.96615929) q[1];
sx q[1];
rz(0.52765347) q[1];
rz(-1.4603314) q[3];
sx q[3];
rz(-0.85064954) q[3];
sx q[3];
rz(-0.11321774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2923773) q[2];
sx q[2];
rz(-1.4039618) q[2];
sx q[2];
rz(-0.89548573) q[2];
rz(2.605947) q[3];
sx q[3];
rz(-1.5629385) q[3];
sx q[3];
rz(-3.079788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.8931005) q[0];
sx q[0];
rz(-2.94815) q[0];
sx q[0];
rz(-0.50791159) q[0];
rz(-1.6912564) q[1];
sx q[1];
rz(-2.129887) q[1];
sx q[1];
rz(1.3571665) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0471418) q[0];
sx q[0];
rz(-3.019382) q[0];
sx q[0];
rz(-2.727174) q[0];
rz(2.4139348) q[2];
sx q[2];
rz(-0.59077677) q[2];
sx q[2];
rz(-2.1075005) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0101945) q[1];
sx q[1];
rz(-1.6784188) q[1];
sx q[1];
rz(-1.1788998) q[1];
rz(-1.0413707) q[3];
sx q[3];
rz(-0.66877194) q[3];
sx q[3];
rz(-3.1094299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42673972) q[2];
sx q[2];
rz(-1.4740976) q[2];
sx q[2];
rz(1.1023785) q[2];
rz(-2.0057996) q[3];
sx q[3];
rz(-1.2165242) q[3];
sx q[3];
rz(-0.77146012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.256846) q[0];
sx q[0];
rz(-2.1717635) q[0];
sx q[0];
rz(-0.92887512) q[0];
rz(0.38201395) q[1];
sx q[1];
rz(-2.1451352) q[1];
sx q[1];
rz(1.6928203) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0576768) q[0];
sx q[0];
rz(-1.8734212) q[0];
sx q[0];
rz(2.641119) q[0];
rz(-1.5097627) q[2];
sx q[2];
rz(-1.0866797) q[2];
sx q[2];
rz(2.8060437) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9517437) q[1];
sx q[1];
rz(-2.7286224) q[1];
sx q[1];
rz(-1.8443405) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9881416) q[3];
sx q[3];
rz(-2.7355237) q[3];
sx q[3];
rz(-1.4514597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6021619) q[2];
sx q[2];
rz(-0.346589) q[2];
sx q[2];
rz(-0.76751417) q[2];
rz(1.2933939) q[3];
sx q[3];
rz(-0.69197217) q[3];
sx q[3];
rz(-0.20492157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9354189) q[0];
sx q[0];
rz(-0.17429166) q[0];
sx q[0];
rz(0.25099227) q[0];
rz(0.17768606) q[1];
sx q[1];
rz(-1.6975941) q[1];
sx q[1];
rz(-2.4846855) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0278695) q[0];
sx q[0];
rz(-1.8060469) q[0];
sx q[0];
rz(0.21224888) q[0];
rz(1.7288293) q[2];
sx q[2];
rz(-0.92301805) q[2];
sx q[2];
rz(-0.72731804) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.73093534) q[1];
sx q[1];
rz(-1.1852131) q[1];
sx q[1];
rz(1.247184) q[1];
x q[2];
rz(-1.1360938) q[3];
sx q[3];
rz(-2.1699749) q[3];
sx q[3];
rz(-0.51927943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3333007) q[2];
sx q[2];
rz(-2.5623645) q[2];
sx q[2];
rz(-0.91326886) q[2];
rz(-2.463786) q[3];
sx q[3];
rz(-1.5014239) q[3];
sx q[3];
rz(-1.0031797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
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
rz(1.0388357) q[0];
sx q[0];
rz(-1.1433733) q[0];
sx q[0];
rz(0.58297771) q[0];
rz(-1.4684756) q[1];
sx q[1];
rz(-2.1491094) q[1];
sx q[1];
rz(-0.34117064) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7852029) q[0];
sx q[0];
rz(-1.5283397) q[0];
sx q[0];
rz(1.3713942) q[0];
rz(-pi) q[1];
rz(-2.274373) q[2];
sx q[2];
rz(-0.3647621) q[2];
sx q[2];
rz(-2.4976969) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9205089) q[1];
sx q[1];
rz(-0.53295164) q[1];
sx q[1];
rz(-2.3531662) q[1];
rz(2.2489585) q[3];
sx q[3];
rz(-1.2481108) q[3];
sx q[3];
rz(1.762515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.94656452) q[2];
sx q[2];
rz(-0.82411689) q[2];
sx q[2];
rz(0.80671802) q[2];
rz(-2.6840456) q[3];
sx q[3];
rz(-2.1541607) q[3];
sx q[3];
rz(-2.257982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6954527) q[0];
sx q[0];
rz(-2.0841632) q[0];
sx q[0];
rz(-0.17644185) q[0];
rz(0.31013075) q[1];
sx q[1];
rz(-1.6519494) q[1];
sx q[1];
rz(2.2055221) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1131206) q[0];
sx q[0];
rz(-0.33361379) q[0];
sx q[0];
rz(-0.46974515) q[0];
rz(-pi) q[1];
rz(1.1175198) q[2];
sx q[2];
rz(-1.0220811) q[2];
sx q[2];
rz(-1.1191776) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8011978) q[1];
sx q[1];
rz(-2.7615393) q[1];
sx q[1];
rz(-1.317522) q[1];
x q[2];
rz(-0.12054875) q[3];
sx q[3];
rz(-1.0998187) q[3];
sx q[3];
rz(0.69191832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.04756847) q[2];
sx q[2];
rz(-1.3022283) q[2];
sx q[2];
rz(3.1330718) q[2];
rz(2.408037) q[3];
sx q[3];
rz(-2.3365884) q[3];
sx q[3];
rz(0.12695299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-1.9200639) q[0];
sx q[0];
rz(-0.41291741) q[0];
sx q[0];
rz(2.371696) q[0];
rz(0.13433111) q[1];
sx q[1];
rz(-2.3513992) q[1];
sx q[1];
rz(1.92164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8304571) q[0];
sx q[0];
rz(-2.756098) q[0];
sx q[0];
rz(2.4804082) q[0];
rz(-pi) q[1];
rz(-2.7292615) q[2];
sx q[2];
rz(-2.4930232) q[2];
sx q[2];
rz(2.1438053) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.091948) q[1];
sx q[1];
rz(-1.94749) q[1];
sx q[1];
rz(2.6833181) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0073529) q[3];
sx q[3];
rz(-2.9749492) q[3];
sx q[3];
rz(1.8145916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8613345) q[2];
sx q[2];
rz(-0.65797776) q[2];
sx q[2];
rz(-0.78557837) q[2];
rz(-0.33513364) q[3];
sx q[3];
rz(-2.1527055) q[3];
sx q[3];
rz(-2.3194763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32421865) q[0];
sx q[0];
rz(-2.2786409) q[0];
sx q[0];
rz(1.8550158) q[0];
rz(-2.4556976) q[1];
sx q[1];
rz(-0.58882014) q[1];
sx q[1];
rz(-1.3480766) q[1];
rz(1.5790719) q[2];
sx q[2];
rz(-0.58115056) q[2];
sx q[2];
rz(-1.5528408) q[2];
rz(-3.104628) q[3];
sx q[3];
rz(-1.2782469) q[3];
sx q[3];
rz(-0.20991355) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
