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
rz(-2.5353105) q[0];
sx q[0];
rz(-1.0507133) q[0];
sx q[0];
rz(-2.3442955) q[0];
rz(0.94766131) q[1];
sx q[1];
rz(-1.892765) q[1];
sx q[1];
rz(0.35511425) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41751244) q[0];
sx q[0];
rz(-1.3211498) q[0];
sx q[0];
rz(1.5019519) q[0];
x q[1];
rz(2.0482641) q[2];
sx q[2];
rz(-1.588378) q[2];
sx q[2];
rz(2.4934752) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0700291) q[1];
sx q[1];
rz(-1.8172846) q[1];
sx q[1];
rz(0.60004289) q[1];
rz(2.0641692) q[3];
sx q[3];
rz(-1.1680654) q[3];
sx q[3];
rz(-1.9014408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.015410272) q[2];
sx q[2];
rz(-2.5842032) q[2];
sx q[2];
rz(-1.1861447) q[2];
rz(-1.0960389) q[3];
sx q[3];
rz(-0.78117424) q[3];
sx q[3];
rz(-1.5196777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010051589) q[0];
sx q[0];
rz(-3.0942656) q[0];
sx q[0];
rz(1.9897687) q[0];
rz(0.25543073) q[1];
sx q[1];
rz(-1.7841745) q[1];
sx q[1];
rz(0.40886042) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20862234) q[0];
sx q[0];
rz(-1.044036) q[0];
sx q[0];
rz(0.43430682) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11714012) q[2];
sx q[2];
rz(-0.88136281) q[2];
sx q[2];
rz(1.8833138) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3656816) q[1];
sx q[1];
rz(-1.3789007) q[1];
sx q[1];
rz(-0.86071162) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84108431) q[3];
sx q[3];
rz(-2.1434823) q[3];
sx q[3];
rz(0.12668474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6795808) q[2];
sx q[2];
rz(-2.2363594) q[2];
sx q[2];
rz(-0.26642624) q[2];
rz(2.5816495) q[3];
sx q[3];
rz(-0.67548776) q[3];
sx q[3];
rz(-0.42375666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68488455) q[0];
sx q[0];
rz(-2.279778) q[0];
sx q[0];
rz(1.7488712) q[0];
rz(-2.4640153) q[1];
sx q[1];
rz(-1.1303439) q[1];
sx q[1];
rz(-0.63634253) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22998097) q[0];
sx q[0];
rz(-1.1781357) q[0];
sx q[0];
rz(-0.61907452) q[0];
x q[1];
rz(0.88865033) q[2];
sx q[2];
rz(-1.7782619) q[2];
sx q[2];
rz(2.6234174) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84987931) q[1];
sx q[1];
rz(-0.91227898) q[1];
sx q[1];
rz(1.9610599) q[1];
rz(-pi) q[2];
rz(1.1701835) q[3];
sx q[3];
rz(-0.83016443) q[3];
sx q[3];
rz(-0.12987353) q[3];
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
rz(1.9321457) q[3];
sx q[3];
rz(-2.3085322) q[3];
sx q[3];
rz(-0.81592733) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4010312) q[0];
sx q[0];
rz(-2.7771692) q[0];
sx q[0];
rz(0.60582274) q[0];
rz(-1.1653028) q[1];
sx q[1];
rz(-2.0348246) q[1];
sx q[1];
rz(-2.5071526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.620011) q[0];
sx q[0];
rz(-1.6130035) q[0];
sx q[0];
rz(-0.47733929) q[0];
rz(-pi) q[1];
rz(1.4342821) q[2];
sx q[2];
rz(-1.0336598) q[2];
sx q[2];
rz(-0.31202635) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.099720864) q[1];
sx q[1];
rz(-2.1812154) q[1];
sx q[1];
rz(-0.94967752) q[1];
rz(-0.090417763) q[3];
sx q[3];
rz(-1.3141229) q[3];
sx q[3];
rz(2.924447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67924649) q[2];
sx q[2];
rz(-2.1470368) q[2];
sx q[2];
rz(-1.4086949) q[2];
rz(2.5016221) q[3];
sx q[3];
rz(-2.7100115) q[3];
sx q[3];
rz(1.4225167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6022279) q[0];
sx q[0];
rz(-2.866565) q[0];
sx q[0];
rz(-1.6590903) q[0];
rz(1.1461343) q[1];
sx q[1];
rz(-2.4303552) q[1];
sx q[1];
rz(-1.8123288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70721165) q[0];
sx q[0];
rz(-0.9117306) q[0];
sx q[0];
rz(-0.7313318) q[0];
rz(2.2634451) q[2];
sx q[2];
rz(-1.9598388) q[2];
sx q[2];
rz(-1.4506952) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.043221323) q[1];
sx q[1];
rz(-2.8403106) q[1];
sx q[1];
rz(-2.2194018) q[1];
rz(-pi) q[2];
rz(-1.8962675) q[3];
sx q[3];
rz(-1.6314037) q[3];
sx q[3];
rz(1.7610723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1819522) q[2];
sx q[2];
rz(-0.96936172) q[2];
sx q[2];
rz(-2.1592965) q[2];
rz(0.2333518) q[3];
sx q[3];
rz(-1.2956649) q[3];
sx q[3];
rz(2.6006234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0503466) q[0];
sx q[0];
rz(-2.8627658) q[0];
sx q[0];
rz(-2.4009551) q[0];
rz(0.29847538) q[1];
sx q[1];
rz(-1.1497755) q[1];
sx q[1];
rz(-0.99075913) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6388004) q[0];
sx q[0];
rz(-0.37257344) q[0];
sx q[0];
rz(-0.034046455) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.58133881) q[2];
sx q[2];
rz(-2.5816988) q[2];
sx q[2];
rz(2.8357504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7700532) q[1];
sx q[1];
rz(-2.3828875) q[1];
sx q[1];
rz(0.9758238) q[1];
x q[2];
rz(-1.9781817) q[3];
sx q[3];
rz(-2.2797027) q[3];
sx q[3];
rz(-1.2825625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6159181) q[2];
sx q[2];
rz(-0.89927036) q[2];
sx q[2];
rz(2.565628) q[2];
rz(1.2032262) q[3];
sx q[3];
rz(-1.4148182) q[3];
sx q[3];
rz(1.6171713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7791157) q[0];
sx q[0];
rz(-0.285382) q[0];
sx q[0];
rz(-1.9051636) q[0];
rz(1.7615039) q[1];
sx q[1];
rz(-2.3347336) q[1];
sx q[1];
rz(-0.39464748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7056859) q[0];
sx q[0];
rz(-1.5360695) q[0];
sx q[0];
rz(-1.7931563) q[0];
x q[1];
rz(0.93776607) q[2];
sx q[2];
rz(-2.870976) q[2];
sx q[2];
rz(1.1039162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9731991) q[1];
sx q[1];
rz(-1.3915807) q[1];
sx q[1];
rz(-2.1922796) q[1];
rz(-1.7622855) q[3];
sx q[3];
rz(-1.8575274) q[3];
sx q[3];
rz(-0.54785997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.75957647) q[2];
sx q[2];
rz(-0.76618177) q[2];
sx q[2];
rz(1.7570599) q[2];
rz(2.0218487) q[3];
sx q[3];
rz(-1.9789968) q[3];
sx q[3];
rz(-0.82974452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(1.8501594) q[0];
sx q[0];
rz(-2.9740574) q[0];
sx q[0];
rz(2.7450388) q[0];
rz(1.5605759) q[1];
sx q[1];
rz(-0.4220337) q[1];
sx q[1];
rz(0.54403725) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4648923) q[0];
sx q[0];
rz(-1.6026622) q[0];
sx q[0];
rz(-0.7384517) q[0];
x q[1];
rz(1.3048473) q[2];
sx q[2];
rz(-1.6839993) q[2];
sx q[2];
rz(-1.2987657) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48780826) q[1];
sx q[1];
rz(-1.7188196) q[1];
sx q[1];
rz(-2.6201999) q[1];
rz(-0.097594768) q[3];
sx q[3];
rz(-1.1371693) q[3];
sx q[3];
rz(-1.6083838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3364526) q[2];
sx q[2];
rz(-0.84638941) q[2];
sx q[2];
rz(1.7479755) q[2];
rz(-2.9584068) q[3];
sx q[3];
rz(-1.2510679) q[3];
sx q[3];
rz(2.2209871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377604) q[0];
sx q[0];
rz(-2.5354009) q[0];
sx q[0];
rz(3.0699741) q[0];
rz(0.76511818) q[1];
sx q[1];
rz(-1.3970929) q[1];
sx q[1];
rz(2.5535233) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.360726) q[0];
sx q[0];
rz(-2.4712997) q[0];
sx q[0];
rz(1.0295342) q[0];
rz(-pi) q[1];
rz(-2.2526365) q[2];
sx q[2];
rz(-1.6168904) q[2];
sx q[2];
rz(0.40929133) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7797554) q[1];
sx q[1];
rz(-1.6653524) q[1];
sx q[1];
rz(-2.4744417) q[1];
x q[2];
rz(-1.9790827) q[3];
sx q[3];
rz(-0.65925778) q[3];
sx q[3];
rz(1.4548649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6832983) q[2];
sx q[2];
rz(-0.27687645) q[2];
sx q[2];
rz(1.5923306) q[2];
rz(-1.8706627) q[3];
sx q[3];
rz(-1.1829605) q[3];
sx q[3];
rz(2.8196238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1147605) q[0];
sx q[0];
rz(-0.25323594) q[0];
sx q[0];
rz(-0.34227553) q[0];
rz(2.5905142) q[1];
sx q[1];
rz(-1.2197878) q[1];
sx q[1];
rz(-2.6384027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4095112) q[0];
sx q[0];
rz(-2.4284673) q[0];
sx q[0];
rz(0.54422191) q[0];
x q[1];
rz(0.047930389) q[2];
sx q[2];
rz(-0.88719545) q[2];
sx q[2];
rz(-2.0237109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79784517) q[1];
sx q[1];
rz(-0.78727951) q[1];
sx q[1];
rz(2.3549453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6588361) q[3];
sx q[3];
rz(-2.1121893) q[3];
sx q[3];
rz(-0.71200395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2912579) q[2];
sx q[2];
rz(-2.1230272) q[2];
sx q[2];
rz(-2.2955718) q[2];
rz(-0.033898354) q[3];
sx q[3];
rz(-0.97550052) q[3];
sx q[3];
rz(-2.2013544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7593183) q[0];
sx q[0];
rz(-1.8702392) q[0];
sx q[0];
rz(2.2478588) q[0];
rz(-1.2636802) q[1];
sx q[1];
rz(-0.79669851) q[1];
sx q[1];
rz(-3.0279108) q[1];
rz(2.6631217) q[2];
sx q[2];
rz(-0.84979117) q[2];
sx q[2];
rz(-0.69862943) q[2];
rz(0.51919785) q[3];
sx q[3];
rz(-0.62102269) q[3];
sx q[3];
rz(0.86925918) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
