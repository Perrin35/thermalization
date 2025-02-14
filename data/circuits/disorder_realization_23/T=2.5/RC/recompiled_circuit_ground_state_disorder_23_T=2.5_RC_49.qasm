OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.58148122) q[0];
sx q[0];
rz(-1.886263) q[0];
sx q[0];
rz(-0.49129593) q[0];
rz(2.8331941) q[1];
sx q[1];
rz(-1.3003132) q[1];
sx q[1];
rz(0.058606776) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8347625) q[0];
sx q[0];
rz(-2.1170724) q[0];
sx q[0];
rz(-0.52301126) q[0];
rz(0.7841447) q[2];
sx q[2];
rz(-1.9040739) q[2];
sx q[2];
rz(-3.0212457) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3112009) q[1];
sx q[1];
rz(-1.7313384) q[1];
sx q[1];
rz(2.83648) q[1];
rz(-pi) q[2];
rz(-1.3335854) q[3];
sx q[3];
rz(-2.8136642) q[3];
sx q[3];
rz(-1.6746317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.10597524) q[2];
sx q[2];
rz(-0.79088894) q[2];
sx q[2];
rz(2.7604575) q[2];
rz(3.1086339) q[3];
sx q[3];
rz(-1.6842753) q[3];
sx q[3];
rz(-1.0491252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1586665) q[0];
sx q[0];
rz(-2.6148112) q[0];
sx q[0];
rz(0.43989936) q[0];
rz(-0.061821763) q[1];
sx q[1];
rz(-1.5301751) q[1];
sx q[1];
rz(0.27438146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8500811) q[0];
sx q[0];
rz(-2.3604104) q[0];
sx q[0];
rz(2.2456332) q[0];
x q[1];
rz(2.0437061) q[2];
sx q[2];
rz(-0.47240546) q[2];
sx q[2];
rz(2.1771569) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8148796) q[1];
sx q[1];
rz(-0.45456262) q[1];
sx q[1];
rz(1.7437032) q[1];
rz(-pi) q[2];
rz(0.37526826) q[3];
sx q[3];
rz(-0.57884848) q[3];
sx q[3];
rz(-1.3001833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8476094) q[2];
sx q[2];
rz(-0.01394883) q[2];
sx q[2];
rz(0.070276109) q[2];
rz(-0.62486068) q[3];
sx q[3];
rz(-1.3288386) q[3];
sx q[3];
rz(-0.074987324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6845508) q[0];
sx q[0];
rz(-1.5032737) q[0];
sx q[0];
rz(0.3918089) q[0];
rz(-2.1458972) q[1];
sx q[1];
rz(-2.3052146) q[1];
sx q[1];
rz(0.044274274) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5747815) q[0];
sx q[0];
rz(-2.1315082) q[0];
sx q[0];
rz(-0.36578481) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4344014) q[2];
sx q[2];
rz(-0.91333616) q[2];
sx q[2];
rz(0.75188504) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4391195) q[1];
sx q[1];
rz(-0.28077341) q[1];
sx q[1];
rz(-2.4149815) q[1];
rz(1.6674394) q[3];
sx q[3];
rz(-2.1657888) q[3];
sx q[3];
rz(-0.71221029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6280881) q[2];
sx q[2];
rz(-2.2254483) q[2];
sx q[2];
rz(1.8872895) q[2];
rz(-0.11145505) q[3];
sx q[3];
rz(-0.97193757) q[3];
sx q[3];
rz(-2.286262) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3067538) q[0];
sx q[0];
rz(-0.67973891) q[0];
sx q[0];
rz(0.87598959) q[0];
rz(2.7442979) q[1];
sx q[1];
rz(-1.5700424) q[1];
sx q[1];
rz(1.809459) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0007325) q[0];
sx q[0];
rz(-2.6291641) q[0];
sx q[0];
rz(-0.39649244) q[0];
rz(-pi) q[1];
rz(0.45355637) q[2];
sx q[2];
rz(-0.95090129) q[2];
sx q[2];
rz(-0.93728055) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.5046103) q[1];
sx q[1];
rz(-0.17269293) q[1];
sx q[1];
rz(2.6673698) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0334084) q[3];
sx q[3];
rz(-1.7741446) q[3];
sx q[3];
rz(2.5739079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9925655) q[2];
sx q[2];
rz(-1.6529447) q[2];
sx q[2];
rz(1.7472902) q[2];
rz(0.99327883) q[3];
sx q[3];
rz(-1.9781878) q[3];
sx q[3];
rz(-1.2629898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8472327) q[0];
sx q[0];
rz(-2.9672186) q[0];
sx q[0];
rz(0.85884035) q[0];
rz(0.42293388) q[1];
sx q[1];
rz(-2.6507288) q[1];
sx q[1];
rz(-1.4609969) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57371712) q[0];
sx q[0];
rz(-1.6934614) q[0];
sx q[0];
rz(0.50267641) q[0];
rz(-pi) q[1];
rz(2.016045) q[2];
sx q[2];
rz(-1.5684874) q[2];
sx q[2];
rz(-0.18773676) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7778421) q[1];
sx q[1];
rz(-1.3652752) q[1];
sx q[1];
rz(1.7670034) q[1];
x q[2];
rz(2.9988838) q[3];
sx q[3];
rz(-2.0721216) q[3];
sx q[3];
rz(-0.10979788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.4451663) q[2];
sx q[2];
rz(-0.40881613) q[2];
sx q[2];
rz(-1.6430631) q[2];
rz(-1.1284575) q[3];
sx q[3];
rz(-2.0339637) q[3];
sx q[3];
rz(2.4988153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.4437101) q[0];
sx q[0];
rz(-1.3297798) q[0];
sx q[0];
rz(-0.71769303) q[0];
rz(-2.7806661) q[1];
sx q[1];
rz(-2.2310427) q[1];
sx q[1];
rz(-2.8275729) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12355328) q[0];
sx q[0];
rz(-1.6121056) q[0];
sx q[0];
rz(2.6459751) q[0];
x q[1];
rz(1.5271565) q[2];
sx q[2];
rz(-0.72481643) q[2];
sx q[2];
rz(0.61312719) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5962299) q[1];
sx q[1];
rz(-2.0343421) q[1];
sx q[1];
rz(0.35407663) q[1];
x q[2];
rz(-0.76992294) q[3];
sx q[3];
rz(-1.7020546) q[3];
sx q[3];
rz(2.8974887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.47739425) q[2];
sx q[2];
rz(-1.7919284) q[2];
sx q[2];
rz(2.9273709) q[2];
rz(2.2231806) q[3];
sx q[3];
rz(-0.94614202) q[3];
sx q[3];
rz(1.4001747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65900954) q[0];
sx q[0];
rz(-0.56203401) q[0];
sx q[0];
rz(-2.5157978) q[0];
rz(1.2311426) q[1];
sx q[1];
rz(-1.8742671) q[1];
sx q[1];
rz(2.321718) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5974665) q[0];
sx q[0];
rz(-1.0080999) q[0];
sx q[0];
rz(1.1981127) q[0];
x q[1];
rz(-0.49687044) q[2];
sx q[2];
rz(-2.3665603) q[2];
sx q[2];
rz(-2.1642223) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.35988369) q[1];
sx q[1];
rz(-1.4952342) q[1];
sx q[1];
rz(1.9025926) q[1];
x q[2];
rz(-1.5890454) q[3];
sx q[3];
rz(-1.4762763) q[3];
sx q[3];
rz(2.8121473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1890373) q[2];
sx q[2];
rz(-1.2574715) q[2];
sx q[2];
rz(0.38743585) q[2];
rz(2.6881325) q[3];
sx q[3];
rz(-0.87415868) q[3];
sx q[3];
rz(0.15997729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3979935) q[0];
sx q[0];
rz(-2.0026119) q[0];
sx q[0];
rz(1.1238264) q[0];
rz(0.75346142) q[1];
sx q[1];
rz(-1.1898142) q[1];
sx q[1];
rz(1.9728647) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3982464) q[0];
sx q[0];
rz(-0.69730824) q[0];
sx q[0];
rz(1.7137609) q[0];
rz(0.0059466023) q[2];
sx q[2];
rz(-1.1134673) q[2];
sx q[2];
rz(-1.3471239) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.67186165) q[1];
sx q[1];
rz(-2.1657092) q[1];
sx q[1];
rz(2.7487526) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99645741) q[3];
sx q[3];
rz(-2.399181) q[3];
sx q[3];
rz(1.3848828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0231102) q[2];
sx q[2];
rz(-2.1239943) q[2];
sx q[2];
rz(-2.6712096) q[2];
rz(-1.6810301) q[3];
sx q[3];
rz(-1.8813671) q[3];
sx q[3];
rz(2.1605261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39101446) q[0];
sx q[0];
rz(-1.7711201) q[0];
sx q[0];
rz(1.4960666) q[0];
rz(1.1356614) q[1];
sx q[1];
rz(-1.2874425) q[1];
sx q[1];
rz(-0.48887238) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35621168) q[0];
sx q[0];
rz(-1.4309034) q[0];
sx q[0];
rz(3.1091305) q[0];
rz(-pi) q[1];
rz(1.4004009) q[2];
sx q[2];
rz(-2.6238447) q[2];
sx q[2];
rz(-2.2108159) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.006039) q[1];
sx q[1];
rz(-0.054243739) q[1];
sx q[1];
rz(2.7179621) q[1];
rz(-pi) q[2];
rz(3.0472894) q[3];
sx q[3];
rz(-0.62178388) q[3];
sx q[3];
rz(-1.0434542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6966072) q[2];
sx q[2];
rz(-0.61034909) q[2];
sx q[2];
rz(1.102591) q[2];
rz(1.5251112) q[3];
sx q[3];
rz(-2.3307255) q[3];
sx q[3];
rz(1.471224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1333604) q[0];
sx q[0];
rz(-0.52953774) q[0];
sx q[0];
rz(1.3326921) q[0];
rz(-2.7545199) q[1];
sx q[1];
rz(-1.8862855) q[1];
sx q[1];
rz(2.0852087) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7232976) q[0];
sx q[0];
rz(-1.3271062) q[0];
sx q[0];
rz(2.9979628) q[0];
rz(1.6956639) q[2];
sx q[2];
rz(-1.1556632) q[2];
sx q[2];
rz(-1.268569) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38195189) q[1];
sx q[1];
rz(-2.3084062) q[1];
sx q[1];
rz(-2.6159389) q[1];
rz(-pi) q[2];
rz(-0.71223082) q[3];
sx q[3];
rz(-2.3281186) q[3];
sx q[3];
rz(1.2016736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.010783823) q[2];
sx q[2];
rz(-0.84785145) q[2];
sx q[2];
rz(-1.7424142) q[2];
rz(-1.0731267) q[3];
sx q[3];
rz(-1.224204) q[3];
sx q[3];
rz(2.6789902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.241093) q[0];
sx q[0];
rz(-1.4971965) q[0];
sx q[0];
rz(1.5370488) q[0];
rz(0.72147876) q[1];
sx q[1];
rz(-0.57054467) q[1];
sx q[1];
rz(1.9346938) q[1];
rz(-0.73312326) q[2];
sx q[2];
rz(-0.55126581) q[2];
sx q[2];
rz(0.049964213) q[2];
rz(-1.6695475) q[3];
sx q[3];
rz(-0.80174123) q[3];
sx q[3];
rz(1.5855736) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
