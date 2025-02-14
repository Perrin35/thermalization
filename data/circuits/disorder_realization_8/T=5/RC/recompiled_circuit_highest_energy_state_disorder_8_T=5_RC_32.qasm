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
rz(0.60628211) q[0];
sx q[0];
rz(4.192306) q[0];
sx q[0];
rz(11.769073) q[0];
rz(0.94766131) q[1];
sx q[1];
rz(-1.892765) q[1];
sx q[1];
rz(0.35511425) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.170318) q[0];
sx q[0];
rz(-1.5040893) q[0];
sx q[0];
rz(2.8913777) q[0];
rz(-1.5325512) q[2];
sx q[2];
rz(-2.6638263) q[2];
sx q[2];
rz(-2.2528799) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0715636) q[1];
sx q[1];
rz(-1.3243081) q[1];
sx q[1];
rz(0.60004289) q[1];
rz(-pi) q[2];
x q[2];
rz(2.303351) q[3];
sx q[3];
rz(-2.5154193) q[3];
sx q[3];
rz(-0.29919344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1261824) q[2];
sx q[2];
rz(-2.5842032) q[2];
sx q[2];
rz(-1.1861447) q[2];
rz(1.0960389) q[3];
sx q[3];
rz(-0.78117424) q[3];
sx q[3];
rz(-1.621915) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1315411) q[0];
sx q[0];
rz(-3.0942656) q[0];
sx q[0];
rz(-1.9897687) q[0];
rz(2.8861619) q[1];
sx q[1];
rz(-1.3574182) q[1];
sx q[1];
rz(0.40886042) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6050312) q[0];
sx q[0];
rz(-0.6694109) q[0];
sx q[0];
rz(2.1971357) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4299655) q[2];
sx q[2];
rz(-2.4438876) q[2];
sx q[2];
rz(-2.066246) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.042603167) q[1];
sx q[1];
rz(-0.87638646) q[1];
sx q[1];
rz(-2.8907772) q[1];
x q[2];
rz(-2.3005083) q[3];
sx q[3];
rz(-0.99811038) q[3];
sx q[3];
rz(-3.0149079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.46201181) q[2];
sx q[2];
rz(-2.2363594) q[2];
sx q[2];
rz(0.26642624) q[2];
rz(2.5816495) q[3];
sx q[3];
rz(-0.67548776) q[3];
sx q[3];
rz(2.717836) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68488455) q[0];
sx q[0];
rz(-2.279778) q[0];
sx q[0];
rz(1.3927214) q[0];
rz(-2.4640153) q[1];
sx q[1];
rz(-1.1303439) q[1];
sx q[1];
rz(2.5052501) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0669495) q[0];
sx q[0];
rz(-2.136629) q[0];
sx q[0];
rz(2.0412579) q[0];
x q[1];
rz(-2.8767831) q[2];
sx q[2];
rz(-2.2356459) q[2];
sx q[2];
rz(-2.2547371) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.174068) q[1];
sx q[1];
rz(-1.2651769) q[1];
sx q[1];
rz(-0.6966865) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78181617) q[3];
sx q[3];
rz(-1.2788427) q[3];
sx q[3];
rz(-1.7192732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59844083) q[2];
sx q[2];
rz(-1.2739015) q[2];
sx q[2];
rz(0.72371036) q[2];
rz(1.2094469) q[3];
sx q[3];
rz(-0.8330605) q[3];
sx q[3];
rz(2.3256653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4010312) q[0];
sx q[0];
rz(-2.7771692) q[0];
sx q[0];
rz(-2.5357699) q[0];
rz(1.1653028) q[1];
sx q[1];
rz(-1.106768) q[1];
sx q[1];
rz(0.63444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.070557) q[0];
sx q[0];
rz(-1.0939176) q[0];
sx q[0];
rz(-1.6183075) q[0];
rz(-pi) q[1];
rz(2.6003378) q[2];
sx q[2];
rz(-1.6879904) q[2];
sx q[2];
rz(-1.8126496) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0418718) q[1];
sx q[1];
rz(-2.1812154) q[1];
sx q[1];
rz(2.1919151) q[1];
x q[2];
rz(3.0511749) q[3];
sx q[3];
rz(-1.3141229) q[3];
sx q[3];
rz(2.924447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4623462) q[2];
sx q[2];
rz(-0.99455589) q[2];
sx q[2];
rz(1.7328978) q[2];
rz(-0.6399706) q[3];
sx q[3];
rz(-2.7100115) q[3];
sx q[3];
rz(-1.7190759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6022279) q[0];
sx q[0];
rz(-0.27502763) q[0];
sx q[0];
rz(-1.4825024) q[0];
rz(-1.1461343) q[1];
sx q[1];
rz(-2.4303552) q[1];
sx q[1];
rz(-1.3292638) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6791554) q[0];
sx q[0];
rz(-0.94171038) q[0];
sx q[0];
rz(-0.85926677) q[0];
rz(-0.48945697) q[2];
sx q[2];
rz(-0.93867368) q[2];
sx q[2];
rz(2.7165627) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.98745495) q[1];
sx q[1];
rz(-1.7510268) q[1];
sx q[1];
rz(1.3280416) q[1];
rz(-pi) q[2];
x q[2];
rz(0.063956694) q[3];
sx q[3];
rz(-1.8956479) q[3];
sx q[3];
rz(-2.9717546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9596404) q[2];
sx q[2];
rz(-2.1722309) q[2];
sx q[2];
rz(2.1592965) q[2];
rz(-2.9082409) q[3];
sx q[3];
rz(-1.2956649) q[3];
sx q[3];
rz(-0.54096925) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0503466) q[0];
sx q[0];
rz(-2.8627658) q[0];
sx q[0];
rz(2.4009551) q[0];
rz(-2.8431173) q[1];
sx q[1];
rz(-1.1497755) q[1];
sx q[1];
rz(2.1508335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6753526) q[0];
sx q[0];
rz(-1.9431433) q[0];
sx q[0];
rz(-1.5574934) q[0];
rz(-pi) q[1];
rz(-1.2392943) q[2];
sx q[2];
rz(-1.1109034) q[2];
sx q[2];
rz(-2.1761328) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37153944) q[1];
sx q[1];
rz(-0.75870514) q[1];
sx q[1];
rz(-0.9758238) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3902682) q[3];
sx q[3];
rz(-1.2653143) q[3];
sx q[3];
rz(0.014367291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6159181) q[2];
sx q[2];
rz(-2.2423223) q[2];
sx q[2];
rz(-0.57596469) q[2];
rz(1.2032262) q[3];
sx q[3];
rz(-1.4148182) q[3];
sx q[3];
rz(-1.5244213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7791157) q[0];
sx q[0];
rz(-2.8562107) q[0];
sx q[0];
rz(1.9051636) q[0];
rz(-1.3800887) q[1];
sx q[1];
rz(-0.80685902) q[1];
sx q[1];
rz(-2.7469452) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1427396) q[0];
sx q[0];
rz(-1.79302) q[0];
sx q[0];
rz(3.1059899) q[0];
rz(1.790843) q[2];
sx q[2];
rz(-1.4119822) q[2];
sx q[2];
rz(-1.0823298) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.27543008) q[1];
sx q[1];
rz(-2.1808561) q[1];
sx q[1];
rz(-0.21924217) q[1];
x q[2];
rz(2.5684082) q[3];
sx q[3];
rz(-0.34332192) q[3];
sx q[3];
rz(-1.1487701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3820162) q[2];
sx q[2];
rz(-0.76618177) q[2];
sx q[2];
rz(-1.7570599) q[2];
rz(-1.1197439) q[3];
sx q[3];
rz(-1.1625959) q[3];
sx q[3];
rz(-2.3118481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8501594) q[0];
sx q[0];
rz(-2.9740574) q[0];
sx q[0];
rz(-0.39655381) q[0];
rz(1.5810168) q[1];
sx q[1];
rz(-2.719559) q[1];
sx q[1];
rz(-2.5975554) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85910857) q[0];
sx q[0];
rz(-0.73900925) q[0];
sx q[0];
rz(-3.0942731) q[0];
x q[1];
rz(1.3048473) q[2];
sx q[2];
rz(-1.4575934) q[2];
sx q[2];
rz(1.2987657) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99847356) q[1];
sx q[1];
rz(-1.0556738) q[1];
sx q[1];
rz(-1.7410941) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3633779) q[3];
sx q[3];
rz(-0.44379297) q[3];
sx q[3];
rz(-1.3794514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3364526) q[2];
sx q[2];
rz(-0.84638941) q[2];
sx q[2];
rz(1.7479755) q[2];
rz(-2.9584068) q[3];
sx q[3];
rz(-1.8905247) q[3];
sx q[3];
rz(0.92060554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2038323) q[0];
sx q[0];
rz(-2.5354009) q[0];
sx q[0];
rz(-3.0699741) q[0];
rz(2.3764745) q[1];
sx q[1];
rz(-1.7444997) q[1];
sx q[1];
rz(-0.58806932) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7808666) q[0];
sx q[0];
rz(-0.67029291) q[0];
sx q[0];
rz(1.0295342) q[0];
rz(-pi) q[1];
rz(-1.4977354) q[2];
sx q[2];
rz(-0.68314766) q[2];
sx q[2];
rz(-1.1047995) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7797554) q[1];
sx q[1];
rz(-1.4762403) q[1];
sx q[1];
rz(-0.66715099) q[1];
rz(-pi) q[2];
rz(2.189013) q[3];
sx q[3];
rz(-1.8164571) q[3];
sx q[3];
rz(0.21354833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.4582943) q[2];
sx q[2];
rz(-0.27687645) q[2];
sx q[2];
rz(-1.5923306) q[2];
rz(1.27093) q[3];
sx q[3];
rz(-1.9586321) q[3];
sx q[3];
rz(0.32196885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1147605) q[0];
sx q[0];
rz(-2.8883567) q[0];
sx q[0];
rz(-0.34227553) q[0];
rz(0.5510785) q[1];
sx q[1];
rz(-1.9218048) q[1];
sx q[1];
rz(-2.6384027) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4069243) q[0];
sx q[0];
rz(-0.97678629) q[0];
sx q[0];
rz(1.9918562) q[0];
rz(-pi) q[1];
rz(3.0936623) q[2];
sx q[2];
rz(-0.88719545) q[2];
sx q[2];
rz(-1.1178818) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3437475) q[1];
sx q[1];
rz(-0.78727951) q[1];
sx q[1];
rz(-0.78664732) q[1];
rz(-pi) q[2];
rz(-0.14519174) q[3];
sx q[3];
rz(-0.54779966) q[3];
sx q[3];
rz(-2.5992268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8503348) q[2];
sx q[2];
rz(-1.0185654) q[2];
sx q[2];
rz(0.84602082) q[2];
rz(0.033898354) q[3];
sx q[3];
rz(-0.97550052) q[3];
sx q[3];
rz(-0.94023824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7593183) q[0];
sx q[0];
rz(-1.2713534) q[0];
sx q[0];
rz(-0.8937339) q[0];
rz(-1.8779124) q[1];
sx q[1];
rz(-2.3448941) q[1];
sx q[1];
rz(5/(14*pi)) q[1];
rz(-1.0882174) q[2];
sx q[2];
rz(-0.84090424) q[2];
sx q[2];
rz(3.1089208) q[2];
rz(1.9119156) q[3];
sx q[3];
rz(-2.1003953) q[3];
sx q[3];
rz(0.25672124) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
