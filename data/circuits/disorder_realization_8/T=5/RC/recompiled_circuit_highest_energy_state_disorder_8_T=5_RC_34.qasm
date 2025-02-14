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
rz(-2.0908794) q[0];
sx q[0];
rz(2.3442955) q[0];
rz(-2.1939313) q[1];
sx q[1];
rz(5.0343577) q[1];
sx q[1];
rz(9.0696637) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9962486) q[0];
sx q[0];
rz(-0.25877413) q[0];
sx q[0];
rz(0.26352675) q[0];
x q[1];
rz(3.1217977) q[2];
sx q[2];
rz(-2.0481841) q[2];
sx q[2];
rz(-2.2098178) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.15820643) q[1];
sx q[1];
rz(-0.64288989) q[1];
sx q[1];
rz(-2.7224273) q[1];
rz(-pi) q[2];
x q[2];
rz(2.303351) q[3];
sx q[3];
rz(-2.5154193) q[3];
sx q[3];
rz(2.8423992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1261824) q[2];
sx q[2];
rz(-0.55738944) q[2];
sx q[2];
rz(-1.9554479) q[2];
rz(-2.0455537) q[3];
sx q[3];
rz(-2.3604184) q[3];
sx q[3];
rz(1.621915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-3.1315411) q[0];
sx q[0];
rz(-3.0942656) q[0];
sx q[0];
rz(-1.151824) q[0];
rz(-2.8861619) q[1];
sx q[1];
rz(-1.3574182) q[1];
sx q[1];
rz(-0.40886042) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20862234) q[0];
sx q[0];
rz(-1.044036) q[0];
sx q[0];
rz(-0.43430682) q[0];
rz(-0.87798543) q[2];
sx q[2];
rz(-1.6610985) q[2];
sx q[2];
rz(2.754359) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0989895) q[1];
sx q[1];
rz(-2.2652062) q[1];
sx q[1];
rz(-0.25081544) q[1];
rz(-pi) q[2];
rz(2.4284372) q[3];
sx q[3];
rz(-2.1655313) q[3];
sx q[3];
rz(-0.99280533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.46201181) q[2];
sx q[2];
rz(-0.90523326) q[2];
sx q[2];
rz(0.26642624) q[2];
rz(2.5816495) q[3];
sx q[3];
rz(-2.4661049) q[3];
sx q[3];
rz(-2.717836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4567081) q[0];
sx q[0];
rz(-0.86181462) q[0];
sx q[0];
rz(1.3927214) q[0];
rz(0.67757732) q[1];
sx q[1];
rz(-1.1303439) q[1];
sx q[1];
rz(2.5052501) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9116117) q[0];
sx q[0];
rz(-1.9634569) q[0];
sx q[0];
rz(-0.61907452) q[0];
rz(-1.2485593) q[2];
sx q[2];
rz(-0.70813484) q[2];
sx q[2];
rz(1.3009877) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2917133) q[1];
sx q[1];
rz(-0.91227898) q[1];
sx q[1];
rz(1.1805328) q[1];
rz(-2.3597765) q[3];
sx q[3];
rz(-1.2788427) q[3];
sx q[3];
rz(-1.4223195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59844083) q[2];
sx q[2];
rz(-1.8676912) q[2];
sx q[2];
rz(2.4178823) q[2];
rz(-1.2094469) q[3];
sx q[3];
rz(-0.8330605) q[3];
sx q[3];
rz(0.81592733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
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
rz(1.1653028) q[1];
sx q[1];
rz(-1.106768) q[1];
sx q[1];
rz(-2.5071526) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1737752) q[0];
sx q[0];
rz(-0.47905827) q[0];
sx q[0];
rz(-3.0499247) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7073105) q[2];
sx q[2];
rz(-1.0336598) q[2];
sx q[2];
rz(2.8295663) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.099720864) q[1];
sx q[1];
rz(-2.1812154) q[1];
sx q[1];
rz(-2.1919151) q[1];
rz(-pi) q[2];
rz(-1.2394514) q[3];
sx q[3];
rz(-2.8697909) q[3];
sx q[3];
rz(3.0157175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4623462) q[2];
sx q[2];
rz(-0.99455589) q[2];
sx q[2];
rz(-1.4086949) q[2];
rz(0.6399706) q[3];
sx q[3];
rz(-2.7100115) q[3];
sx q[3];
rz(-1.4225167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53936476) q[0];
sx q[0];
rz(-2.866565) q[0];
sx q[0];
rz(-1.4825024) q[0];
rz(1.1461343) q[1];
sx q[1];
rz(-2.4303552) q[1];
sx q[1];
rz(-1.8123288) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6791554) q[0];
sx q[0];
rz(-0.94171038) q[0];
sx q[0];
rz(0.85926677) q[0];
rz(-2.1414928) q[2];
sx q[2];
rz(-2.3632103) q[2];
sx q[2];
rz(0.308643) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5138897) q[1];
sx q[1];
rz(-1.3320508) q[1];
sx q[1];
rz(0.18555218) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.063956694) q[3];
sx q[3];
rz(-1.2459447) q[3];
sx q[3];
rz(0.16983804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9596404) q[2];
sx q[2];
rz(-2.1722309) q[2];
sx q[2];
rz(-0.98229614) q[2];
rz(0.2333518) q[3];
sx q[3];
rz(-1.8459277) q[3];
sx q[3];
rz(0.54096925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.091246009) q[0];
sx q[0];
rz(-0.27882689) q[0];
sx q[0];
rz(-2.4009551) q[0];
rz(2.8431173) q[1];
sx q[1];
rz(-1.1497755) q[1];
sx q[1];
rz(-2.1508335) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0418763) q[0];
sx q[0];
rz(-1.558405) q[0];
sx q[0];
rz(0.37237694) q[0];
x q[1];
rz(-2.6590587) q[2];
sx q[2];
rz(-1.8667456) q[2];
sx q[2];
rz(-0.75693989) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0195398) q[1];
sx q[1];
rz(-0.96457997) q[1];
sx q[1];
rz(2.6531924) q[1];
rz(-pi) q[2];
rz(2.7088183) q[3];
sx q[3];
rz(-0.7996586) q[3];
sx q[3];
rz(-1.8679269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52567452) q[2];
sx q[2];
rz(-2.2423223) q[2];
sx q[2];
rz(-2.565628) q[2];
rz(-1.9383664) q[3];
sx q[3];
rz(-1.7267745) q[3];
sx q[3];
rz(-1.6171713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7791157) q[0];
sx q[0];
rz(-2.8562107) q[0];
sx q[0];
rz(1.9051636) q[0];
rz(-1.7615039) q[1];
sx q[1];
rz(-0.80685902) q[1];
sx q[1];
rz(2.7469452) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1590736) q[0];
sx q[0];
rz(-0.22501105) q[0];
sx q[0];
rz(1.7270442) q[0];
rz(0.93776607) q[2];
sx q[2];
rz(-2.870976) q[2];
sx q[2];
rz(1.1039162) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4952364) q[1];
sx q[1];
rz(-2.4980821) q[1];
sx q[1];
rz(-1.8724426) q[1];
x q[2];
rz(0.5731845) q[3];
sx q[3];
rz(-0.34332192) q[3];
sx q[3];
rz(1.1487701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3820162) q[2];
sx q[2];
rz(-2.3754109) q[2];
sx q[2];
rz(-1.7570599) q[2];
rz(-2.0218487) q[3];
sx q[3];
rz(-1.9789968) q[3];
sx q[3];
rz(0.82974452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2914332) q[0];
sx q[0];
rz(-2.9740574) q[0];
sx q[0];
rz(-0.39655381) q[0];
rz(1.5605759) q[1];
sx q[1];
rz(-0.4220337) q[1];
sx q[1];
rz(-2.5975554) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85910857) q[0];
sx q[0];
rz(-0.73900925) q[0];
sx q[0];
rz(0.047319558) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9790573) q[2];
sx q[2];
rz(-2.8530793) q[2];
sx q[2];
rz(2.4764594) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.99847356) q[1];
sx q[1];
rz(-1.0556738) q[1];
sx q[1];
rz(-1.7410941) q[1];
rz(-pi) q[2];
rz(-0.097594768) q[3];
sx q[3];
rz(-2.0044234) q[3];
sx q[3];
rz(1.6083838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3364526) q[2];
sx q[2];
rz(-0.84638941) q[2];
sx q[2];
rz(-1.3936172) q[2];
rz(-0.18318583) q[3];
sx q[3];
rz(-1.2510679) q[3];
sx q[3];
rz(-2.2209871) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9377604) q[0];
sx q[0];
rz(-0.60619175) q[0];
sx q[0];
rz(3.0699741) q[0];
rz(-0.76511818) q[1];
sx q[1];
rz(-1.7444997) q[1];
sx q[1];
rz(2.5535233) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65031728) q[0];
sx q[0];
rz(-1.2450019) q[0];
sx q[0];
rz(-0.97401818) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0822524) q[2];
sx q[2];
rz(-2.2517747) q[2];
sx q[2];
rz(-1.1988893) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0068609) q[1];
sx q[1];
rz(-2.2344338) q[1];
sx q[1];
rz(-1.4506542) q[1];
rz(-pi) q[2];
rz(-0.95257967) q[3];
sx q[3];
rz(-1.8164571) q[3];
sx q[3];
rz(0.21354833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4582943) q[2];
sx q[2];
rz(-2.8647162) q[2];
sx q[2];
rz(1.5923306) q[2];
rz(-1.8706627) q[3];
sx q[3];
rz(-1.9586321) q[3];
sx q[3];
rz(0.32196885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0268321) q[0];
sx q[0];
rz(-0.25323594) q[0];
sx q[0];
rz(-0.34227553) q[0];
rz(-2.5905142) q[1];
sx q[1];
rz(-1.2197878) q[1];
sx q[1];
rz(2.6384027) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7321271) q[0];
sx q[0];
rz(-1.9163462) q[0];
sx q[0];
rz(2.5045128) q[0];
rz(-pi) q[1];
rz(1.6295428) q[2];
sx q[2];
rz(-0.68500942) q[2];
sx q[2];
rz(2.0995121) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.15718225) q[1];
sx q[1];
rz(-1.0453893) q[1];
sx q[1];
rz(-2.5249283) q[1];
x q[2];
rz(-0.54310829) q[3];
sx q[3];
rz(-1.4953729) q[3];
sx q[3];
rz(2.2373445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2912579) q[2];
sx q[2];
rz(-1.0185654) q[2];
sx q[2];
rz(-0.84602082) q[2];
rz(-3.1076943) q[3];
sx q[3];
rz(-0.97550052) q[3];
sx q[3];
rz(2.2013544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.7593183) q[0];
sx q[0];
rz(-1.8702392) q[0];
sx q[0];
rz(2.2478588) q[0];
rz(1.2636802) q[1];
sx q[1];
rz(-2.3448941) q[1];
sx q[1];
rz(5/(14*pi)) q[1];
rz(2.0533753) q[2];
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
