OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(-1.9999737) q[1];
sx q[1];
rz(3.5715754) q[1];
sx q[1];
rz(6.9663098) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0177512) q[0];
sx q[0];
rz(-0.089028247) q[0];
sx q[0];
rz(2.428399) q[0];
rz(-1.2201266) q[2];
sx q[2];
rz(-1.3436683) q[2];
sx q[2];
rz(0.35958689) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5176757) q[1];
sx q[1];
rz(-0.78849925) q[1];
sx q[1];
rz(-0.6448402) q[1];
x q[2];
rz(-2.1360374) q[3];
sx q[3];
rz(-1.2473277) q[3];
sx q[3];
rz(-2.2068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1258939) q[2];
sx q[2];
rz(-1.761972) q[2];
sx q[2];
rz(1.0985628) q[2];
rz(2.0627608) q[3];
sx q[3];
rz(-0.94509411) q[3];
sx q[3];
rz(-1.1014972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.9803479) q[0];
sx q[0];
rz(-2.8256567) q[0];
sx q[0];
rz(-2.9336477) q[0];
rz(2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.4651441) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11008308) q[0];
sx q[0];
rz(-0.84658909) q[0];
sx q[0];
rz(-1.5741041) q[0];
rz(-pi) q[1];
rz(1.1128088) q[2];
sx q[2];
rz(-2.6070242) q[2];
sx q[2];
rz(2.2528258) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0456902) q[1];
sx q[1];
rz(-2.0369139) q[1];
sx q[1];
rz(2.5065266) q[1];
rz(-pi) q[2];
x q[2];
rz(1.915669) q[3];
sx q[3];
rz(-0.85540918) q[3];
sx q[3];
rz(-1.8638924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1229822) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(-3.0318276) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-2.770335) q[3];
sx q[3];
rz(1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-2.3410472) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1313989) q[0];
sx q[0];
rz(-2.7079765) q[0];
sx q[0];
rz(1.9171159) q[0];
rz(-pi) q[1];
rz(-2.3163047) q[2];
sx q[2];
rz(-0.85916677) q[2];
sx q[2];
rz(-0.7691783) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2785981) q[1];
sx q[1];
rz(-1.9837712) q[1];
sx q[1];
rz(-1.5763271) q[1];
x q[2];
rz(-2.839746) q[3];
sx q[3];
rz(-2.2491124) q[3];
sx q[3];
rz(0.18920004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.67733726) q[2];
sx q[2];
rz(-0.3192454) q[2];
sx q[2];
rz(-1.3595954) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.8680957) q[3];
sx q[3];
rz(1.4250071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1490705) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(2.676679) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(-2.0565313) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5269055) q[0];
sx q[0];
rz(-2.5589716) q[0];
sx q[0];
rz(-2.7132062) q[0];
rz(-1.2041353) q[2];
sx q[2];
rz(-1.8030093) q[2];
sx q[2];
rz(1.0166849) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.76772056) q[1];
sx q[1];
rz(-1.671119) q[1];
sx q[1];
rz(-2.9994681) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.61120175) q[3];
sx q[3];
rz(-1.9107606) q[3];
sx q[3];
rz(-1.1838278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-0.51182169) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27914771) q[0];
sx q[0];
rz(-2.6036766) q[0];
sx q[0];
rz(1.7329247) q[0];
rz(2.7092343) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(-2.1599105) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0504793) q[0];
sx q[0];
rz(-1.0266773) q[0];
sx q[0];
rz(-1.0247466) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9972516) q[2];
sx q[2];
rz(-0.28595668) q[2];
sx q[2];
rz(0.85948247) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7179467) q[1];
sx q[1];
rz(-1.4232993) q[1];
sx q[1];
rz(0.46848483) q[1];
rz(-pi) q[2];
rz(3.0675689) q[3];
sx q[3];
rz(-0.97462666) q[3];
sx q[3];
rz(-2.1285469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0354054) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(-0.48941082) q[2];
rz(-2.126746) q[3];
sx q[3];
rz(-0.5265407) q[3];
sx q[3];
rz(1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.6569825) q[0];
sx q[0];
rz(-0.15743142) q[0];
sx q[0];
rz(-0.37242517) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(2.9763124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0342456) q[0];
sx q[0];
rz(-3.0421071) q[0];
sx q[0];
rz(-2.8427567) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6406306) q[2];
sx q[2];
rz(-0.88431057) q[2];
sx q[2];
rz(-1.7718466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.11152553) q[1];
sx q[1];
rz(-2.1038052) q[1];
sx q[1];
rz(-1.0982606) q[1];
rz(-2.9061756) q[3];
sx q[3];
rz(-1.132292) q[3];
sx q[3];
rz(0.32270839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.91339397) q[2];
sx q[2];
rz(-2.1134351) q[2];
sx q[2];
rz(-1.4432663) q[2];
rz(-0.36744395) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(-2.929556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3570324) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(-0.34926397) q[0];
rz(2.3941984) q[1];
sx q[1];
rz(-0.29577297) q[1];
sx q[1];
rz(-0.73648891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9518785) q[0];
sx q[0];
rz(-1.5218381) q[0];
sx q[0];
rz(3.1244856) q[0];
rz(-pi) q[1];
rz(2.1963504) q[2];
sx q[2];
rz(-2.7316299) q[2];
sx q[2];
rz(1.9117102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.10662096) q[1];
sx q[1];
rz(-1.7502516) q[1];
sx q[1];
rz(-2.3751395) q[1];
rz(-pi) q[2];
rz(-0.50076671) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(0.70481833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13654576) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(0.53156701) q[2];
rz(0.68938869) q[3];
sx q[3];
rz(-1.4644943) q[3];
sx q[3];
rz(-1.846107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078995973) q[0];
sx q[0];
rz(-1.9364708) q[0];
sx q[0];
rz(-1.2980365) q[0];
rz(-0.80728665) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(-0.92179006) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73496504) q[0];
sx q[0];
rz(-1.4359183) q[0];
sx q[0];
rz(0.69358967) q[0];
rz(-1.8838896) q[2];
sx q[2];
rz(-0.76258341) q[2];
sx q[2];
rz(1.2459754) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8242952) q[1];
sx q[1];
rz(-1.3166787) q[1];
sx q[1];
rz(-1.7829249) q[1];
rz(-pi) q[2];
rz(-0.17955762) q[3];
sx q[3];
rz(-2.2710544) q[3];
sx q[3];
rz(-1.2683271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.90199295) q[2];
sx q[2];
rz(-0.56100503) q[2];
sx q[2];
rz(1.9699338) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.6849018) q[3];
sx q[3];
rz(-1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(0.25306025) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(-1.4755479) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(0.67970651) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98201671) q[0];
sx q[0];
rz(-1.7734818) q[0];
sx q[0];
rz(-2.6152339) q[0];
x q[1];
rz(-1.3557415) q[2];
sx q[2];
rz(-0.70677033) q[2];
sx q[2];
rz(3.0422473) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.45021536) q[1];
sx q[1];
rz(-2.2274349) q[1];
sx q[1];
rz(-2.0037829) q[1];
x q[2];
rz(0.67919517) q[3];
sx q[3];
rz(-1.067357) q[3];
sx q[3];
rz(-0.60590832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0150962) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(2.2311907) q[2];
rz(2.4662468) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-2.2909686) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0891721) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(2.4269379) q[0];
rz(0.71406281) q[1];
sx q[1];
rz(-2.1866182) q[1];
sx q[1];
rz(1.1766599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0929716) q[0];
sx q[0];
rz(-1.7341359) q[0];
sx q[0];
rz(-3.0924348) q[0];
rz(-pi) q[1];
rz(-0.54729692) q[2];
sx q[2];
rz(-2.1314869) q[2];
sx q[2];
rz(1.0123569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1365876) q[1];
sx q[1];
rz(-2.1274381) q[1];
sx q[1];
rz(3.0535166) q[1];
x q[2];
rz(0.15575274) q[3];
sx q[3];
rz(-2.6253346) q[3];
sx q[3];
rz(1.932365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.24361336) q[2];
sx q[2];
rz(-1.4695797) q[2];
sx q[2];
rz(2.0142377) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(-1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7174299) q[0];
sx q[0];
rz(-1.2107727) q[0];
sx q[0];
rz(-2.6834224) q[0];
rz(0.39623109) q[1];
sx q[1];
rz(-0.025066499) q[1];
sx q[1];
rz(0.33096663) q[1];
rz(-1.2926119) q[2];
sx q[2];
rz(-0.85360151) q[2];
sx q[2];
rz(-3.1384946) q[2];
rz(-3.1191961) q[3];
sx q[3];
rz(-0.35084421) q[3];
sx q[3];
rz(2.4435333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
