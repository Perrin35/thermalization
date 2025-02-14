OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.073590241) q[0];
sx q[0];
rz(-1.6728787) q[0];
sx q[0];
rz(-0.52748632) q[0];
rz(-2.5975851) q[1];
sx q[1];
rz(-0.44096947) q[1];
sx q[1];
rz(-3.1402332) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6327335) q[0];
sx q[0];
rz(-1.5616882) q[0];
sx q[0];
rz(1.5004116) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5293526) q[2];
sx q[2];
rz(-1.7667102) q[2];
sx q[2];
rz(1.6454091) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7391752) q[1];
sx q[1];
rz(-1.1858984) q[1];
sx q[1];
rz(-0.56300122) q[1];
rz(2.9980761) q[3];
sx q[3];
rz(-2.9401192) q[3];
sx q[3];
rz(-0.62769393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.15343836) q[2];
sx q[2];
rz(-0.64186382) q[2];
sx q[2];
rz(1.7786857) q[2];
rz(3.0617132) q[3];
sx q[3];
rz(-0.86013836) q[3];
sx q[3];
rz(-0.0063272198) q[3];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7774571) q[0];
sx q[0];
rz(-1.9123001) q[0];
sx q[0];
rz(2.6312713) q[0];
rz(1.3279042) q[1];
sx q[1];
rz(-2.4102305) q[1];
sx q[1];
rz(-0.43689835) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53117673) q[0];
sx q[0];
rz(-1.0836837) q[0];
sx q[0];
rz(-0.25419828) q[0];
rz(-pi) q[1];
rz(-2.6969621) q[2];
sx q[2];
rz(-1.7845588) q[2];
sx q[2];
rz(-1.7997326) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4480312) q[1];
sx q[1];
rz(-0.80151075) q[1];
sx q[1];
rz(-1.526725) q[1];
x q[2];
rz(1.1673959) q[3];
sx q[3];
rz(-1.6296436) q[3];
sx q[3];
rz(-1.4368536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3029311) q[2];
sx q[2];
rz(-2.624056) q[2];
sx q[2];
rz(-2.9663864) q[2];
rz(-3.0025499) q[3];
sx q[3];
rz(-1.4929205) q[3];
sx q[3];
rz(2.2429332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43056968) q[0];
sx q[0];
rz(-2.2570606) q[0];
sx q[0];
rz(2.7810466) q[0];
rz(-1.4973466) q[1];
sx q[1];
rz(-1.9162354) q[1];
sx q[1];
rz(-2.5490733) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1688878) q[0];
sx q[0];
rz(-2.1967432) q[0];
sx q[0];
rz(-0.79373534) q[0];
rz(-1.4176263) q[2];
sx q[2];
rz(-0.94050759) q[2];
sx q[2];
rz(3.1068124) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.46099037) q[1];
sx q[1];
rz(-2.7719738) q[1];
sx q[1];
rz(1.9113944) q[1];
x q[2];
rz(-2.6733562) q[3];
sx q[3];
rz(-1.9561202) q[3];
sx q[3];
rz(-1.4778397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3561919) q[2];
sx q[2];
rz(-2.5742026) q[2];
sx q[2];
rz(-2.1991335) q[2];
rz(2.5475907) q[3];
sx q[3];
rz(-2.3022251) q[3];
sx q[3];
rz(-3.0187637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9451611) q[0];
sx q[0];
rz(-2.5158947) q[0];
sx q[0];
rz(0.64836597) q[0];
rz(-2.2658589) q[1];
sx q[1];
rz(-1.6301194) q[1];
sx q[1];
rz(-1.3513563) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195456) q[0];
sx q[0];
rz(-0.81455671) q[0];
sx q[0];
rz(0.053653663) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9586323) q[2];
sx q[2];
rz(-0.81987971) q[2];
sx q[2];
rz(-0.41868107) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0160603) q[1];
sx q[1];
rz(-2.2439306) q[1];
sx q[1];
rz(-0.036733997) q[1];
x q[2];
rz(1.625678) q[3];
sx q[3];
rz(-1.1883231) q[3];
sx q[3];
rz(0.11688133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.15328345) q[2];
sx q[2];
rz(-0.51258665) q[2];
sx q[2];
rz(1.2549866) q[2];
rz(-0.60162383) q[3];
sx q[3];
rz(-0.87698495) q[3];
sx q[3];
rz(1.6833444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0772142) q[0];
sx q[0];
rz(-2.0266396) q[0];
sx q[0];
rz(0.13378046) q[0];
rz(-1.6556219) q[1];
sx q[1];
rz(-0.32951117) q[1];
sx q[1];
rz(-2.9697184) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3456589) q[0];
sx q[0];
rz(-2.5194114) q[0];
sx q[0];
rz(2.2337929) q[0];
rz(1.1823229) q[2];
sx q[2];
rz(-2.3459593) q[2];
sx q[2];
rz(2.1284136) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.682595) q[1];
sx q[1];
rz(-0.38308278) q[1];
sx q[1];
rz(-1.9869845) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.71545728) q[3];
sx q[3];
rz(-0.60327941) q[3];
sx q[3];
rz(2.6048599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4922501) q[2];
sx q[2];
rz(-2.616373) q[2];
sx q[2];
rz(-1.4225175) q[2];
rz(2.9366142) q[3];
sx q[3];
rz(-2.180438) q[3];
sx q[3];
rz(-0.80406308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.499046) q[0];
sx q[0];
rz(-2.5048984) q[0];
sx q[0];
rz(-2.9472886) q[0];
rz(2.4957472) q[1];
sx q[1];
rz(-1.1767118) q[1];
sx q[1];
rz(-2.7913854) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80794835) q[0];
sx q[0];
rz(-2.5910113) q[0];
sx q[0];
rz(-0.26170611) q[0];
rz(-2.5123213) q[2];
sx q[2];
rz(-2.2221178) q[2];
sx q[2];
rz(2.4547612) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3694075) q[1];
sx q[1];
rz(-0.52474743) q[1];
sx q[1];
rz(0.16903846) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4713481) q[3];
sx q[3];
rz(-1.6247162) q[3];
sx q[3];
rz(2.7516493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0575867) q[2];
sx q[2];
rz(-2.7454594) q[2];
sx q[2];
rz(0.32220379) q[2];
rz(2.3816439) q[3];
sx q[3];
rz(-0.51637572) q[3];
sx q[3];
rz(0.050839067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8754804) q[0];
sx q[0];
rz(-0.66829824) q[0];
sx q[0];
rz(-0.71568263) q[0];
rz(-0.066787668) q[1];
sx q[1];
rz(-2.5804434) q[1];
sx q[1];
rz(-2.7776048) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9371295) q[0];
sx q[0];
rz(-1.9162906) q[0];
sx q[0];
rz(-2.1895467) q[0];
rz(-pi) q[1];
rz(0.90708676) q[2];
sx q[2];
rz(-1.0561475) q[2];
sx q[2];
rz(-2.3163925) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5426483) q[1];
sx q[1];
rz(-2.3348044) q[1];
sx q[1];
rz(-1.5933977) q[1];
x q[2];
rz(-0.19192275) q[3];
sx q[3];
rz(-2.4867184) q[3];
sx q[3];
rz(-2.690243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.62465847) q[2];
sx q[2];
rz(-2.2397569) q[2];
sx q[2];
rz(-1.9067524) q[2];
rz(0.45461795) q[3];
sx q[3];
rz(-1.2841299) q[3];
sx q[3];
rz(0.07479085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55998498) q[0];
sx q[0];
rz(-3.0261664) q[0];
sx q[0];
rz(2.5009632) q[0];
rz(-0.36264125) q[1];
sx q[1];
rz(-0.33932149) q[1];
sx q[1];
rz(1.5193264) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4962728) q[0];
sx q[0];
rz(-2.4970166) q[0];
sx q[0];
rz(-0.98981895) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73070406) q[2];
sx q[2];
rz(-1.3514697) q[2];
sx q[2];
rz(0.64576572) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1040428) q[1];
sx q[1];
rz(-0.68244237) q[1];
sx q[1];
rz(-1.4242092) q[1];
rz(1.0985116) q[3];
sx q[3];
rz(-1.660822) q[3];
sx q[3];
rz(1.93879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2210239) q[2];
sx q[2];
rz(-1.2732882) q[2];
sx q[2];
rz(-1.1576687) q[2];
rz(-2.8122592) q[3];
sx q[3];
rz(-2.9521827) q[3];
sx q[3];
rz(-2.2277189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90410239) q[0];
sx q[0];
rz(-0.79039031) q[0];
sx q[0];
rz(-0.41123408) q[0];
rz(2.5143738) q[1];
sx q[1];
rz(-0.5694446) q[1];
sx q[1];
rz(1.8597182) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86253765) q[0];
sx q[0];
rz(-1.2077966) q[0];
sx q[0];
rz(-1.8022861) q[0];
rz(-pi) q[1];
rz(0.87780805) q[2];
sx q[2];
rz(-2.3690732) q[2];
sx q[2];
rz(2.7038717) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8258649) q[1];
sx q[1];
rz(-2.4701354) q[1];
sx q[1];
rz(-1.3888098) q[1];
rz(-pi) q[2];
rz(-0.19880812) q[3];
sx q[3];
rz(-1.947177) q[3];
sx q[3];
rz(-2.876296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7810479) q[2];
sx q[2];
rz(-1.5881528) q[2];
sx q[2];
rz(-0.9786728) q[2];
rz(2.9141038) q[3];
sx q[3];
rz(-2.3922908) q[3];
sx q[3];
rz(-0.45660517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.6351673) q[0];
sx q[0];
rz(-0.058110617) q[0];
sx q[0];
rz(-0.79570049) q[0];
rz(-2.6129163) q[1];
sx q[1];
rz(-2.1541336) q[1];
sx q[1];
rz(0.19435571) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.07819019) q[0];
sx q[0];
rz(-1.1991797) q[0];
sx q[0];
rz(-0.097784575) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8647381) q[2];
sx q[2];
rz(-1.1938098) q[2];
sx q[2];
rz(1.8238719) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4116507) q[1];
sx q[1];
rz(-0.46272181) q[1];
sx q[1];
rz(0.28280854) q[1];
x q[2];
rz(-2.0967257) q[3];
sx q[3];
rz(-1.1218438) q[3];
sx q[3];
rz(2.8766905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.84870321) q[2];
sx q[2];
rz(-1.6699426) q[2];
sx q[2];
rz(-0.29937747) q[2];
rz(0.40144604) q[3];
sx q[3];
rz(-0.58599389) q[3];
sx q[3];
rz(2.6011023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9837579) q[0];
sx q[0];
rz(-0.83136375) q[0];
sx q[0];
rz(-1.2663483) q[0];
rz(3.0345295) q[1];
sx q[1];
rz(-2.3657847) q[1];
sx q[1];
rz(1.8484144) q[1];
rz(-2.5994095) q[2];
sx q[2];
rz(-1.5246921) q[2];
sx q[2];
rz(0.17490457) q[2];
rz(0.23672531) q[3];
sx q[3];
rz(-1.4631347) q[3];
sx q[3];
rz(2.2377228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
