OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5873592) q[0];
sx q[0];
rz(-0.98266196) q[0];
sx q[0];
rz(2.3798556) q[0];
rz(0.76454437) q[1];
sx q[1];
rz(4.2188797) q[1];
sx q[1];
rz(10.168434) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24510358) q[0];
sx q[0];
rz(-1.6365543) q[0];
sx q[0];
rz(-1.8114281) q[0];
x q[1];
rz(-2.488963) q[2];
sx q[2];
rz(-1.4748117) q[2];
sx q[2];
rz(1.3053615) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8929157) q[1];
sx q[1];
rz(-0.49202737) q[1];
sx q[1];
rz(-2.1881275) q[1];
rz(-pi) q[2];
rz(-1.3841964) q[3];
sx q[3];
rz(-0.36986923) q[3];
sx q[3];
rz(-2.6922525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7068229) q[2];
sx q[2];
rz(-2.7396024) q[2];
sx q[2];
rz(3.0337231) q[2];
rz(-2.9989631) q[3];
sx q[3];
rz(-1.4167891) q[3];
sx q[3];
rz(-0.67255783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07664872) q[0];
sx q[0];
rz(-0.80524421) q[0];
sx q[0];
rz(-2.9108677) q[0];
rz(-1.8143066) q[1];
sx q[1];
rz(-2.4704411) q[1];
sx q[1];
rz(0.040963106) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7728535) q[0];
sx q[0];
rz(-2.8452747) q[0];
sx q[0];
rz(-2.1570737) q[0];
rz(1.2626921) q[2];
sx q[2];
rz(-1.3034046) q[2];
sx q[2];
rz(-1.845713) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.37521024) q[1];
sx q[1];
rz(-2.3751405) q[1];
sx q[1];
rz(-2.4541928) q[1];
rz(-pi) q[2];
rz(-0.81539865) q[3];
sx q[3];
rz(-2.6147463) q[3];
sx q[3];
rz(0.53838733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9791947) q[2];
sx q[2];
rz(-1.4731864) q[2];
sx q[2];
rz(-0.56817788) q[2];
rz(-2.6290821) q[3];
sx q[3];
rz(-0.5957225) q[3];
sx q[3];
rz(-0.56186831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7115962) q[0];
sx q[0];
rz(-1.5189518) q[0];
sx q[0];
rz(2.6913753) q[0];
rz(-1.846116) q[1];
sx q[1];
rz(-1.9772915) q[1];
sx q[1];
rz(-2.4643262) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2579502) q[0];
sx q[0];
rz(-1.6227239) q[0];
sx q[0];
rz(0.16846637) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0063653221) q[2];
sx q[2];
rz(-1.5696988) q[2];
sx q[2];
rz(1.7249677) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.87264204) q[1];
sx q[1];
rz(-2.9193455) q[1];
sx q[1];
rz(2.2058669) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.47795313) q[3];
sx q[3];
rz(-1.9216929) q[3];
sx q[3];
rz(-1.4357476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9230817) q[2];
sx q[2];
rz(-1.9873025) q[2];
sx q[2];
rz(-3.0947321) q[2];
rz(0.81165195) q[3];
sx q[3];
rz(-0.67957687) q[3];
sx q[3];
rz(0.17015447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1470404) q[0];
sx q[0];
rz(-0.11163286) q[0];
sx q[0];
rz(-3.0901093) q[0];
rz(2.6507846) q[1];
sx q[1];
rz(-2.1665116) q[1];
sx q[1];
rz(-1.1725918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8797982) q[0];
sx q[0];
rz(-2.2248631) q[0];
sx q[0];
rz(-1.7394702) q[0];
rz(-1.555221) q[2];
sx q[2];
rz(-2.1767463) q[2];
sx q[2];
rz(-1.7473999) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.51845156) q[1];
sx q[1];
rz(-1.8703096) q[1];
sx q[1];
rz(-1.9407942) q[1];
x q[2];
rz(-3.067279) q[3];
sx q[3];
rz(-1.1399817) q[3];
sx q[3];
rz(-0.58882344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.88901687) q[2];
sx q[2];
rz(-0.93095195) q[2];
sx q[2];
rz(-2.6546997) q[2];
rz(1.1934818) q[3];
sx q[3];
rz(-2.8119757) q[3];
sx q[3];
rz(-0.025432767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30381969) q[0];
sx q[0];
rz(-1.0205512) q[0];
sx q[0];
rz(-1.3254962) q[0];
rz(-2.5505113) q[1];
sx q[1];
rz(-2.4609844) q[1];
sx q[1];
rz(-2.4868734) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39010534) q[0];
sx q[0];
rz(-2.562398) q[0];
sx q[0];
rz(-2.4672227) q[0];
rz(-pi) q[1];
rz(-0.40277092) q[2];
sx q[2];
rz(-2.6908814) q[2];
sx q[2];
rz(2.7223301) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5419546) q[1];
sx q[1];
rz(-1.9497739) q[1];
sx q[1];
rz(0.71211262) q[1];
rz(-pi) q[2];
rz(-1.6642903) q[3];
sx q[3];
rz(-1.1870541) q[3];
sx q[3];
rz(-1.5941217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6902265) q[2];
sx q[2];
rz(-1.9407242) q[2];
sx q[2];
rz(-2.5704685) q[2];
rz(-0.58978224) q[3];
sx q[3];
rz(-2.6815806) q[3];
sx q[3];
rz(-3.0736198) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94474435) q[0];
sx q[0];
rz(-1.9510883) q[0];
sx q[0];
rz(2.8425472) q[0];
rz(1.8213182) q[1];
sx q[1];
rz(-0.25779217) q[1];
sx q[1];
rz(-1.4978131) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59693903) q[0];
sx q[0];
rz(-2.4785846) q[0];
sx q[0];
rz(1.6556428) q[0];
x q[1];
rz(2.1645855) q[2];
sx q[2];
rz(-2.0679681) q[2];
sx q[2];
rz(-1.1360816) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.3469997) q[1];
sx q[1];
rz(-2.0302677) q[1];
sx q[1];
rz(2.8192892) q[1];
x q[2];
rz(-2.4160556) q[3];
sx q[3];
rz(-2.7751338) q[3];
sx q[3];
rz(-0.08034245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.51612878) q[2];
sx q[2];
rz(-1.4279782) q[2];
sx q[2];
rz(-3.1141172) q[2];
rz(-2.6190858) q[3];
sx q[3];
rz(-2.3404739) q[3];
sx q[3];
rz(0.81645042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41247535) q[0];
sx q[0];
rz(-2.0232047) q[0];
sx q[0];
rz(0.11418848) q[0];
rz(-0.97822899) q[1];
sx q[1];
rz(-2.5949635) q[1];
sx q[1];
rz(2.3506929) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40490155) q[0];
sx q[0];
rz(-1.6653403) q[0];
sx q[0];
rz(2.1928284) q[0];
x q[1];
rz(-2.3979264) q[2];
sx q[2];
rz(-1.8218092) q[2];
sx q[2];
rz(0.92168346) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0930209) q[1];
sx q[1];
rz(-1.2067817) q[1];
sx q[1];
rz(3.0572592) q[1];
rz(-pi) q[2];
rz(-2.7780276) q[3];
sx q[3];
rz(-2.225038) q[3];
sx q[3];
rz(3.0594861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.87166446) q[2];
sx q[2];
rz(-1.1130788) q[2];
sx q[2];
rz(2.2765735) q[2];
rz(-0.67251742) q[3];
sx q[3];
rz(-1.1285684) q[3];
sx q[3];
rz(-3.045936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6487811) q[0];
sx q[0];
rz(-0.5195986) q[0];
sx q[0];
rz(0.46736026) q[0];
rz(-0.53721792) q[1];
sx q[1];
rz(-0.98058128) q[1];
sx q[1];
rz(-0.25407243) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2124436) q[0];
sx q[0];
rz(-1.375276) q[0];
sx q[0];
rz(-1.6385727) q[0];
rz(-pi) q[1];
rz(1.5723096) q[2];
sx q[2];
rz(-1.5655893) q[2];
sx q[2];
rz(-1.7092012) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2368187) q[1];
sx q[1];
rz(-2.8042951) q[1];
sx q[1];
rz(2.7493613) q[1];
rz(0.59928008) q[3];
sx q[3];
rz(-1.5044971) q[3];
sx q[3];
rz(1.9950206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6415928) q[2];
sx q[2];
rz(-2.3708512) q[2];
sx q[2];
rz(-2.8059778) q[2];
rz(-2.8149758) q[3];
sx q[3];
rz(-2.2461522) q[3];
sx q[3];
rz(1.4183104) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0582054) q[0];
sx q[0];
rz(-0.21033062) q[0];
sx q[0];
rz(-2.0876419) q[0];
rz(-0.15696934) q[1];
sx q[1];
rz(-1.418768) q[1];
sx q[1];
rz(0.98186791) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4440585) q[0];
sx q[0];
rz(-0.42675787) q[0];
sx q[0];
rz(-0.27192893) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0332546) q[2];
sx q[2];
rz(-1.0045369) q[2];
sx q[2];
rz(-2.1926751) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5280142) q[1];
sx q[1];
rz(-0.19436969) q[1];
sx q[1];
rz(-1.8730875) q[1];
rz(2.8544159) q[3];
sx q[3];
rz(-0.73162006) q[3];
sx q[3];
rz(-1.8457796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2173569) q[2];
sx q[2];
rz(-2.0841667) q[2];
sx q[2];
rz(2.1441933) q[2];
rz(-0.50619566) q[3];
sx q[3];
rz(-2.180407) q[3];
sx q[3];
rz(3.1072646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.85703325) q[0];
sx q[0];
rz(-0.71145809) q[0];
sx q[0];
rz(-0.51666623) q[0];
rz(2.754028) q[1];
sx q[1];
rz(-1.060408) q[1];
sx q[1];
rz(-0.26836747) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70290138) q[0];
sx q[0];
rz(-2.1500906) q[0];
sx q[0];
rz(-0.86433522) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1860113) q[2];
sx q[2];
rz(-2.4590543) q[2];
sx q[2];
rz(-0.62894422) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6827991) q[1];
sx q[1];
rz(-1.8384215) q[1];
sx q[1];
rz(2.5330179) q[1];
x q[2];
rz(0.99276944) q[3];
sx q[3];
rz(-1.7939995) q[3];
sx q[3];
rz(0.85152599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9364075) q[2];
sx q[2];
rz(-0.78339094) q[2];
sx q[2];
rz(-2.5777585) q[2];
rz(1.9514203) q[3];
sx q[3];
rz(-2.1824013) q[3];
sx q[3];
rz(0.64175516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.664809) q[0];
sx q[0];
rz(-1.2510779) q[0];
sx q[0];
rz(-1.0673987) q[0];
rz(-1.339636) q[1];
sx q[1];
rz(-1.7032774) q[1];
sx q[1];
rz(1.3443321) q[1];
rz(2.6396991) q[2];
sx q[2];
rz(-0.5770275) q[2];
sx q[2];
rz(1.8736476) q[2];
rz(2.5793544) q[3];
sx q[3];
rz(-0.36459618) q[3];
sx q[3];
rz(0.96095745) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];