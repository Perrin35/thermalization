OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5034135) q[0];
sx q[0];
rz(-1.2794275) q[0];
sx q[0];
rz(-2.3655565) q[0];
rz(0.70932055) q[1];
sx q[1];
rz(-1.6242937) q[1];
sx q[1];
rz(-0.59901839) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4871019) q[0];
sx q[0];
rz(-1.4485437) q[0];
sx q[0];
rz(-2.046665) q[0];
rz(-1.232595) q[2];
sx q[2];
rz(-0.71824348) q[2];
sx q[2];
rz(2.3423549) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.30556074) q[1];
sx q[1];
rz(-1.2183883) q[1];
sx q[1];
rz(0.52959092) q[1];
rz(0.78035155) q[3];
sx q[3];
rz(-1.8858741) q[3];
sx q[3];
rz(-1.9596069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0029995) q[2];
sx q[2];
rz(-2.004576) q[2];
sx q[2];
rz(2.9453759) q[2];
rz(1.044322) q[3];
sx q[3];
rz(-2.315867) q[3];
sx q[3];
rz(0.31061068) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1021378) q[0];
sx q[0];
rz(-2.4882443) q[0];
sx q[0];
rz(0.85773221) q[0];
rz(2.6013382) q[1];
sx q[1];
rz(-1.0539571) q[1];
sx q[1];
rz(0.78261715) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40561134) q[0];
sx q[0];
rz(-2.1850039) q[0];
sx q[0];
rz(-0.10615291) q[0];
rz(-1.247632) q[2];
sx q[2];
rz(-0.64715451) q[2];
sx q[2];
rz(1.1101071) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2007717) q[1];
sx q[1];
rz(-1.7634974) q[1];
sx q[1];
rz(-0.60639834) q[1];
rz(2.2223118) q[3];
sx q[3];
rz(-0.79943919) q[3];
sx q[3];
rz(-2.8795867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.1648078) q[2];
sx q[2];
rz(-0.79443496) q[2];
sx q[2];
rz(2.5505193) q[2];
rz(-1.0278541) q[3];
sx q[3];
rz(-2.5058392) q[3];
sx q[3];
rz(3.0052321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22981055) q[0];
sx q[0];
rz(-0.52017838) q[0];
sx q[0];
rz(-2.0853364) q[0];
rz(-0.99110574) q[1];
sx q[1];
rz(-1.7678363) q[1];
sx q[1];
rz(2.3371005) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5493889) q[0];
sx q[0];
rz(-2.443586) q[0];
sx q[0];
rz(-2.7207359) q[0];
rz(-1.7011004) q[2];
sx q[2];
rz(-2.5087803) q[2];
sx q[2];
rz(-1.9595343) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.75315216) q[1];
sx q[1];
rz(-2.8049251) q[1];
sx q[1];
rz(0.80516385) q[1];
rz(-pi) q[2];
rz(-1.3785754) q[3];
sx q[3];
rz(-0.9323191) q[3];
sx q[3];
rz(2.9137127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.450401) q[2];
sx q[2];
rz(-0.85019008) q[2];
sx q[2];
rz(0.56785339) q[2];
rz(1.9495226) q[3];
sx q[3];
rz(-0.53693938) q[3];
sx q[3];
rz(0.15538628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(1.6464226) q[0];
sx q[0];
rz(-1.7429054) q[0];
sx q[0];
rz(1.1871185) q[0];
rz(2.8861956) q[1];
sx q[1];
rz(-1.4030158) q[1];
sx q[1];
rz(0.38937169) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0856649) q[0];
sx q[0];
rz(-1.7048258) q[0];
sx q[0];
rz(-0.22179752) q[0];
rz(-pi) q[1];
x q[1];
rz(0.84367532) q[2];
sx q[2];
rz(-1.5932343) q[2];
sx q[2];
rz(-2.0539935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.90044777) q[1];
sx q[1];
rz(-2.4551848) q[1];
sx q[1];
rz(1.871984) q[1];
rz(-3.0087972) q[3];
sx q[3];
rz(-1.7657874) q[3];
sx q[3];
rz(-2.8094069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4046459) q[2];
sx q[2];
rz(-0.32361042) q[2];
sx q[2];
rz(1.7741989) q[2];
rz(0.5395475) q[3];
sx q[3];
rz(-1.5893785) q[3];
sx q[3];
rz(-1.8168943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.0747727) q[0];
sx q[0];
rz(-2.9434151) q[0];
sx q[0];
rz(2.5812126) q[0];
rz(2.9226774) q[1];
sx q[1];
rz(-1.0812662) q[1];
sx q[1];
rz(-0.17094368) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31256277) q[0];
sx q[0];
rz(-2.0382904) q[0];
sx q[0];
rz(-0.74689052) q[0];
rz(-0.60336171) q[2];
sx q[2];
rz(-2.8953279) q[2];
sx q[2];
rz(-0.45349301) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7485908) q[1];
sx q[1];
rz(-2.1868863) q[1];
sx q[1];
rz(2.2115117) q[1];
x q[2];
rz(0.070329026) q[3];
sx q[3];
rz(-2.189895) q[3];
sx q[3];
rz(-0.17507565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.35974744) q[2];
sx q[2];
rz(-1.1283987) q[2];
sx q[2];
rz(-2.183059) q[2];
rz(2.9790699) q[3];
sx q[3];
rz(-0.84238094) q[3];
sx q[3];
rz(1.5490279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.573134) q[0];
sx q[0];
rz(-0.064366654) q[0];
sx q[0];
rz(-0.74813133) q[0];
rz(-2.023078) q[1];
sx q[1];
rz(-2.7225814) q[1];
sx q[1];
rz(2.8245139) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1823163) q[0];
sx q[0];
rz(-1.8166564) q[0];
sx q[0];
rz(-1.6410752) q[0];
x q[1];
rz(0.24131624) q[2];
sx q[2];
rz(-0.20212999) q[2];
sx q[2];
rz(1.2723107) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2851781) q[1];
sx q[1];
rz(-0.76907255) q[1];
sx q[1];
rz(-1.4682707) q[1];
x q[2];
rz(-1.5311144) q[3];
sx q[3];
rz(-0.79426605) q[3];
sx q[3];
rz(-0.021707699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3743484) q[2];
sx q[2];
rz(-0.92864645) q[2];
sx q[2];
rz(1.413215) q[2];
rz(-0.080538571) q[3];
sx q[3];
rz(-1.7929411) q[3];
sx q[3];
rz(0.22182375) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0291075) q[0];
sx q[0];
rz(-1.2487829) q[0];
sx q[0];
rz(1.8674194) q[0];
rz(1.3451276) q[1];
sx q[1];
rz(-2.3990264) q[1];
sx q[1];
rz(1.3412195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2182512) q[0];
sx q[0];
rz(-0.78235432) q[0];
sx q[0];
rz(-2.4104491) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6244782) q[2];
sx q[2];
rz(-0.54143751) q[2];
sx q[2];
rz(2.0064029) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0772895) q[1];
sx q[1];
rz(-1.7619362) q[1];
sx q[1];
rz(-0.15644507) q[1];
rz(-pi) q[2];
rz(-1.4169832) q[3];
sx q[3];
rz(-2.4857368) q[3];
sx q[3];
rz(-0.7747618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.93366569) q[2];
sx q[2];
rz(-1.9542481) q[2];
sx q[2];
rz(-0.46736091) q[2];
rz(2.3815239) q[3];
sx q[3];
rz(-1.4773388) q[3];
sx q[3];
rz(-1.2490341) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46045983) q[0];
sx q[0];
rz(-1.0204027) q[0];
sx q[0];
rz(-1.6946174) q[0];
rz(-2.815822) q[1];
sx q[1];
rz(-0.21369801) q[1];
sx q[1];
rz(-1.6206585) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64031271) q[0];
sx q[0];
rz(-0.58459832) q[0];
sx q[0];
rz(2.6452097) q[0];
x q[1];
rz(0.09163945) q[2];
sx q[2];
rz(-1.5704201) q[2];
sx q[2];
rz(-2.3918652) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.28315) q[1];
sx q[1];
rz(-1.2179188) q[1];
sx q[1];
rz(-2.9550885) q[1];
rz(-1.3255903) q[3];
sx q[3];
rz(-0.71906656) q[3];
sx q[3];
rz(1.9089501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1763566) q[2];
sx q[2];
rz(-0.74702817) q[2];
sx q[2];
rz(0.61526862) q[2];
rz(1.8687013) q[3];
sx q[3];
rz(-0.26510173) q[3];
sx q[3];
rz(0.42241514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1212921) q[0];
sx q[0];
rz(-1.2739807) q[0];
sx q[0];
rz(2.2913388) q[0];
rz(-2.0203159) q[1];
sx q[1];
rz(-1.774615) q[1];
sx q[1];
rz(-2.3698295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61311713) q[0];
sx q[0];
rz(-2.4474047) q[0];
sx q[0];
rz(-1.608169) q[0];
rz(1.6408368) q[2];
sx q[2];
rz(-1.3904499) q[2];
sx q[2];
rz(1.548962) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38852019) q[1];
sx q[1];
rz(-0.24986146) q[1];
sx q[1];
rz(1.7953963) q[1];
rz(-0.27098591) q[3];
sx q[3];
rz(-1.8838722) q[3];
sx q[3];
rz(-2.3199758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0261592) q[2];
sx q[2];
rz(-1.5449646) q[2];
sx q[2];
rz(-2.3010632) q[2];
rz(0.080951512) q[3];
sx q[3];
rz(-0.5778802) q[3];
sx q[3];
rz(-0.24278434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7131272) q[0];
sx q[0];
rz(-2.9436538) q[0];
sx q[0];
rz(-1.9848829) q[0];
rz(2.1139862) q[1];
sx q[1];
rz(-1.3778069) q[1];
sx q[1];
rz(-2.3311232) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9543957) q[0];
sx q[0];
rz(-2.3141197) q[0];
sx q[0];
rz(-1.275497) q[0];
x q[1];
rz(-0.54525156) q[2];
sx q[2];
rz(-1.8601189) q[2];
sx q[2];
rz(-0.74301737) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.6582074) q[1];
sx q[1];
rz(-1.864363) q[1];
sx q[1];
rz(-2.844643) q[1];
rz(-pi) q[2];
rz(0.23328554) q[3];
sx q[3];
rz(-0.57396997) q[3];
sx q[3];
rz(-0.026653224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.345574) q[2];
sx q[2];
rz(-1.2556262) q[2];
sx q[2];
rz(1.5239117) q[2];
rz(-1.1003305) q[3];
sx q[3];
rz(-1.9610145) q[3];
sx q[3];
rz(-0.17980096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1179467) q[0];
sx q[0];
rz(-1.7272341) q[0];
sx q[0];
rz(-0.9247307) q[0];
rz(0.042451518) q[1];
sx q[1];
rz(-1.1141384) q[1];
sx q[1];
rz(1.3420807) q[1];
rz(-3.0855865) q[2];
sx q[2];
rz(-1.1650892) q[2];
sx q[2];
rz(2.7964609) q[2];
rz(-1.5347979) q[3];
sx q[3];
rz(-1.1136342) q[3];
sx q[3];
rz(-1.3200091) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
