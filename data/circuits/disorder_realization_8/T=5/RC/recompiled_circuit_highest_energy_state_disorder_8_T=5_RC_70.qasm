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
rz(-2.1939313) q[1];
sx q[1];
rz(-1.2488276) q[1];
sx q[1];
rz(-0.35511425) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14534403) q[0];
sx q[0];
rz(-0.25877413) q[0];
sx q[0];
rz(-2.8780659) q[0];
rz(-pi) q[1];
x q[1];
rz(0.019794982) q[2];
sx q[2];
rz(-1.0934085) q[2];
sx q[2];
rz(-2.2098178) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0700291) q[1];
sx q[1];
rz(-1.3243081) q[1];
sx q[1];
rz(0.60004289) q[1];
x q[2];
rz(2.6910681) q[3];
sx q[3];
rz(-2.0216216) q[3];
sx q[3];
rz(-2.6032347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.015410272) q[2];
sx q[2];
rz(-0.55738944) q[2];
sx q[2];
rz(-1.9554479) q[2];
rz(1.0960389) q[3];
sx q[3];
rz(-0.78117424) q[3];
sx q[3];
rz(-1.621915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.010051589) q[0];
sx q[0];
rz(-0.047327071) q[0];
sx q[0];
rz(-1.151824) q[0];
rz(-0.25543073) q[1];
sx q[1];
rz(-1.7841745) q[1];
sx q[1];
rz(2.7327322) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6050312) q[0];
sx q[0];
rz(-2.4721818) q[0];
sx q[0];
rz(0.94445697) q[0];
rz(-pi) q[1];
rz(-1.7116272) q[2];
sx q[2];
rz(-0.69770506) q[2];
sx q[2];
rz(1.0753466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0989895) q[1];
sx q[1];
rz(-2.2652062) q[1];
sx q[1];
rz(2.8907772) q[1];
rz(2.4284372) q[3];
sx q[3];
rz(-0.9760614) q[3];
sx q[3];
rz(-2.1487873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46201181) q[2];
sx q[2];
rz(-2.2363594) q[2];
sx q[2];
rz(0.26642624) q[2];
rz(0.5599432) q[3];
sx q[3];
rz(-0.67548776) q[3];
sx q[3];
rz(0.42375666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-0.68488455) q[0];
sx q[0];
rz(-2.279778) q[0];
sx q[0];
rz(-1.7488712) q[0];
rz(2.4640153) q[1];
sx q[1];
rz(-2.0112487) q[1];
sx q[1];
rz(-0.63634253) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9116117) q[0];
sx q[0];
rz(-1.9634569) q[0];
sx q[0];
rz(2.5225181) q[0];
rz(-pi) q[1];
rz(1.8930334) q[2];
sx q[2];
rz(-0.70813484) q[2];
sx q[2];
rz(1.3009877) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8835523) q[1];
sx q[1];
rz(-0.75037724) q[1];
sx q[1];
rz(0.45697327) q[1];
rz(-pi) q[2];
rz(-0.40319188) q[3];
sx q[3];
rz(-2.3180214) q[3];
sx q[3];
rz(0.43063569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59844083) q[2];
sx q[2];
rz(-1.8676912) q[2];
sx q[2];
rz(-0.72371036) q[2];
rz(1.9321457) q[3];
sx q[3];
rz(-0.8330605) q[3];
sx q[3];
rz(0.81592733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4010312) q[0];
sx q[0];
rz(-2.7771692) q[0];
sx q[0];
rz(2.5357699) q[0];
rz(-1.9762899) q[1];
sx q[1];
rz(-1.106768) q[1];
sx q[1];
rz(0.63444) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96781743) q[0];
sx q[0];
rz(-0.47905827) q[0];
sx q[0];
rz(3.0499247) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4342821) q[2];
sx q[2];
rz(-2.1079328) q[2];
sx q[2];
rz(2.8295663) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.099720864) q[1];
sx q[1];
rz(-0.96037727) q[1];
sx q[1];
rz(2.1919151) q[1];
rz(-pi) q[2];
rz(-1.9021413) q[3];
sx q[3];
rz(-2.8697909) q[3];
sx q[3];
rz(0.1258752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.67924649) q[2];
sx q[2];
rz(-0.99455589) q[2];
sx q[2];
rz(-1.7328978) q[2];
rz(2.5016221) q[3];
sx q[3];
rz(-2.7100115) q[3];
sx q[3];
rz(-1.7190759) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6022279) q[0];
sx q[0];
rz(-0.27502763) q[0];
sx q[0];
rz(-1.4825024) q[0];
rz(1.1461343) q[1];
sx q[1];
rz(-2.4303552) q[1];
sx q[1];
rz(1.3292638) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6791554) q[0];
sx q[0];
rz(-2.1998823) q[0];
sx q[0];
rz(2.2823259) q[0];
x q[1];
rz(0.87814754) q[2];
sx q[2];
rz(-1.1817538) q[2];
sx q[2];
rz(1.6908975) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.043221323) q[1];
sx q[1];
rz(-0.30128208) q[1];
sx q[1];
rz(-2.2194018) q[1];
rz(-pi) q[2];
rz(3.077636) q[3];
sx q[3];
rz(-1.2459447) q[3];
sx q[3];
rz(0.16983804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1819522) q[2];
sx q[2];
rz(-0.96936172) q[2];
sx q[2];
rz(0.98229614) q[2];
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
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.091246009) q[0];
sx q[0];
rz(-2.8627658) q[0];
sx q[0];
rz(-2.4009551) q[0];
rz(-2.8431173) q[1];
sx q[1];
rz(-1.9918171) q[1];
sx q[1];
rz(-2.1508335) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0418763) q[0];
sx q[0];
rz(-1.5831876) q[0];
sx q[0];
rz(-0.37237694) q[0];
rz(-0.58133881) q[2];
sx q[2];
rz(-2.5816988) q[2];
sx q[2];
rz(2.8357504) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.37153944) q[1];
sx q[1];
rz(-0.75870514) q[1];
sx q[1];
rz(-0.9758238) q[1];
rz(-2.7088183) q[3];
sx q[3];
rz(-2.3419341) q[3];
sx q[3];
rz(-1.8679269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.52567452) q[2];
sx q[2];
rz(-0.89927036) q[2];
sx q[2];
rz(-2.565628) q[2];
rz(-1.9383664) q[3];
sx q[3];
rz(-1.4148182) q[3];
sx q[3];
rz(1.6171713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3624769) q[0];
sx q[0];
rz(-0.285382) q[0];
sx q[0];
rz(-1.2364291) q[0];
rz(1.3800887) q[1];
sx q[1];
rz(-2.3347336) q[1];
sx q[1];
rz(0.39464748) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1427396) q[0];
sx q[0];
rz(-1.79302) q[0];
sx q[0];
rz(3.1059899) q[0];
rz(-pi) q[1];
rz(-0.93776607) q[2];
sx q[2];
rz(-2.870976) q[2];
sx q[2];
rz(-1.1039162) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1683936) q[1];
sx q[1];
rz(-1.3915807) q[1];
sx q[1];
rz(2.1922796) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29177417) q[3];
sx q[3];
rz(-1.387216) q[3];
sx q[3];
rz(-2.1734299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75957647) q[2];
sx q[2];
rz(-0.76618177) q[2];
sx q[2];
rz(-1.3845328) q[2];
rz(2.0218487) q[3];
sx q[3];
rz(-1.1625959) q[3];
sx q[3];
rz(-2.3118481) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8501594) q[0];
sx q[0];
rz(-2.9740574) q[0];
sx q[0];
rz(0.39655381) q[0];
rz(1.5810168) q[1];
sx q[1];
rz(-0.4220337) q[1];
sx q[1];
rz(2.5975554) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67670031) q[0];
sx q[0];
rz(-1.5389305) q[0];
sx q[0];
rz(-2.403141) q[0];
rz(-pi) q[1];
rz(1.3048473) q[2];
sx q[2];
rz(-1.6839993) q[2];
sx q[2];
rz(1.8428269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.48780826) q[1];
sx q[1];
rz(-1.422773) q[1];
sx q[1];
rz(2.6201999) q[1];
x q[2];
rz(-2.0062449) q[3];
sx q[3];
rz(-1.6593336) q[3];
sx q[3];
rz(0.003525894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3364526) q[2];
sx q[2];
rz(-2.2952032) q[2];
sx q[2];
rz(1.7479755) q[2];
rz(0.18318583) q[3];
sx q[3];
rz(-1.2510679) q[3];
sx q[3];
rz(-0.92060554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2038323) q[0];
sx q[0];
rz(-2.5354009) q[0];
sx q[0];
rz(-3.0699741) q[0];
rz(0.76511818) q[1];
sx q[1];
rz(-1.3970929) q[1];
sx q[1];
rz(-0.58806932) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.360726) q[0];
sx q[0];
rz(-2.4712997) q[0];
sx q[0];
rz(2.1120584) q[0];
x q[1];
rz(2.2526365) q[2];
sx q[2];
rz(-1.6168904) q[2];
sx q[2];
rz(-0.40929133) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.32826281) q[1];
sx q[1];
rz(-2.4687903) q[1];
sx q[1];
rz(-2.9895014) q[1];
x q[2];
rz(-0.29847904) q[3];
sx q[3];
rz(-2.1678181) q[3];
sx q[3];
rz(-1.1859758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6832983) q[2];
sx q[2];
rz(-2.8647162) q[2];
sx q[2];
rz(1.5923306) q[2];
rz(-1.8706627) q[3];
sx q[3];
rz(-1.9586321) q[3];
sx q[3];
rz(-2.8196238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1147605) q[0];
sx q[0];
rz(-0.25323594) q[0];
sx q[0];
rz(-2.7993171) q[0];
rz(2.5905142) q[1];
sx q[1];
rz(-1.9218048) q[1];
sx q[1];
rz(2.6384027) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4095112) q[0];
sx q[0];
rz(-2.4284673) q[0];
sx q[0];
rz(-2.5973707) q[0];
rz(3.0936623) q[2];
sx q[2];
rz(-0.88719545) q[2];
sx q[2];
rz(2.0237109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3863661) q[1];
sx q[1];
rz(-2.0947573) q[1];
sx q[1];
rz(-0.95295277) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9964009) q[3];
sx q[3];
rz(-2.593793) q[3];
sx q[3];
rz(-0.54236587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2912579) q[2];
sx q[2];
rz(-2.1230272) q[2];
sx q[2];
rz(-0.84602082) q[2];
rz(3.1076943) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(2.2013544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3822744) q[0];
sx q[0];
rz(-1.8702392) q[0];
sx q[0];
rz(2.2478588) q[0];
rz(1.8779124) q[1];
sx q[1];
rz(-0.79669851) q[1];
sx q[1];
rz(-3.0279108) q[1];
rz(-0.47847099) q[2];
sx q[2];
rz(-0.84979117) q[2];
sx q[2];
rz(-0.69862943) q[2];
rz(2.6223948) q[3];
sx q[3];
rz(-2.52057) q[3];
sx q[3];
rz(-2.2723335) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
