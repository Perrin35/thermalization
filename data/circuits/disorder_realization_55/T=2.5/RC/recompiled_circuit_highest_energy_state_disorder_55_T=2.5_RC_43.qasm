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
rz(-1.363938) q[0];
sx q[0];
rz(3.5600297) q[0];
sx q[0];
rz(10.726396) q[0];
rz(-2.9786181) q[1];
sx q[1];
rz(-1.6956704) q[1];
sx q[1];
rz(3.1248098) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7977475) q[0];
sx q[0];
rz(-0.36917403) q[0];
sx q[0];
rz(-1.7789715) q[0];
rz(-pi) q[1];
rz(1.4799825) q[2];
sx q[2];
rz(-0.1802643) q[2];
sx q[2];
rz(-1.6423051) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7960733) q[1];
sx q[1];
rz(-1.5658448) q[1];
sx q[1];
rz(-1.5657182) q[1];
x q[2];
rz(-2.9892606) q[3];
sx q[3];
rz(-1.6248584) q[3];
sx q[3];
rz(-3.1076589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5590543) q[2];
sx q[2];
rz(-0.6607008) q[2];
sx q[2];
rz(-1.5793229) q[2];
rz(2.9301379) q[3];
sx q[3];
rz(-0.00051694218) q[3];
sx q[3];
rz(2.9805984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.529539) q[0];
sx q[0];
rz(-2.8654629) q[0];
sx q[0];
rz(-1.3268693) q[0];
rz(2.5621085) q[1];
sx q[1];
rz(-0.0038298413) q[1];
sx q[1];
rz(-0.63900596) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9352183) q[0];
sx q[0];
rz(-0.91418302) q[0];
sx q[0];
rz(-0.73120631) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5875568) q[2];
sx q[2];
rz(-1.4498561) q[2];
sx q[2];
rz(0.018509381) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7182448) q[1];
sx q[1];
rz(-3.1210174) q[1];
sx q[1];
rz(0.59845509) q[1];
rz(-0.78044807) q[3];
sx q[3];
rz(-1.6326346) q[3];
sx q[3];
rz(-2.0831747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7063893) q[2];
sx q[2];
rz(-0.13627626) q[2];
sx q[2];
rz(1.5380247) q[2];
rz(1.5806574) q[3];
sx q[3];
rz(-3.1272562) q[3];
sx q[3];
rz(-3.1106136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8921709) q[0];
sx q[0];
rz(-2.6264661) q[0];
sx q[0];
rz(0.38145915) q[0];
rz(2.4341266) q[1];
sx q[1];
rz(-3.1222157) q[1];
sx q[1];
rz(2.0170508) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8302001) q[0];
sx q[0];
rz(-1.8233577) q[0];
sx q[0];
rz(0.24858944) q[0];
rz(-pi) q[1];
rz(-1.3493645) q[2];
sx q[2];
rz(-3.0237758) q[2];
sx q[2];
rz(3.0750781) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2644653) q[1];
sx q[1];
rz(-1.6337602) q[1];
sx q[1];
rz(-1.5829045) q[1];
rz(-2.0865404) q[3];
sx q[3];
rz(-2.0655144) q[3];
sx q[3];
rz(-0.62913857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4418929) q[2];
sx q[2];
rz(-0.012233891) q[2];
sx q[2];
rz(3.0921248) q[2];
rz(-0.60702819) q[3];
sx q[3];
rz(-3.1403465) q[3];
sx q[3];
rz(-1.9342669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98168755) q[0];
sx q[0];
rz(-0.1683546) q[0];
sx q[0];
rz(-3.1244151) q[0];
rz(-2.848564) q[1];
sx q[1];
rz(-0.79048645) q[1];
sx q[1];
rz(-1.5944098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7538504) q[0];
sx q[0];
rz(-2.1374636) q[0];
sx q[0];
rz(-0.64105861) q[0];
x q[1];
rz(-0.36357673) q[2];
sx q[2];
rz(-0.57206094) q[2];
sx q[2];
rz(-0.043302082) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.80611011) q[1];
sx q[1];
rz(-1.4467738) q[1];
sx q[1];
rz(1.5107442) q[1];
x q[2];
rz(-0.64597102) q[3];
sx q[3];
rz(-1.5785909) q[3];
sx q[3];
rz(-0.92168671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25345099) q[2];
sx q[2];
rz(-0.47051045) q[2];
sx q[2];
rz(-2.4593501) q[2];
rz(0.055179723) q[3];
sx q[3];
rz(-3.1339055) q[3];
sx q[3];
rz(-1.8365708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76160112) q[0];
sx q[0];
rz(-0.43552265) q[0];
sx q[0];
rz(-0.4250266) q[0];
rz(-1.5399326) q[1];
sx q[1];
rz(-0.48307499) q[1];
sx q[1];
rz(-0.81533283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9352933) q[0];
sx q[0];
rz(-1.5816798) q[0];
sx q[0];
rz(0.061200415) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.710395) q[2];
sx q[2];
rz(-0.023471467) q[2];
sx q[2];
rz(-1.0271629) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.40774959) q[1];
sx q[1];
rz(-1.4516648) q[1];
sx q[1];
rz(1.4834189) q[1];
rz(0.32834239) q[3];
sx q[3];
rz(-0.36188618) q[3];
sx q[3];
rz(-0.34235024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1347947) q[2];
sx q[2];
rz(-0.012601348) q[2];
sx q[2];
rz(-1.6689782) q[2];
rz(2.2115479) q[3];
sx q[3];
rz(-3.127122) q[3];
sx q[3];
rz(0.84021935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8188266) q[0];
sx q[0];
rz(-3.1184986) q[0];
sx q[0];
rz(-1.7148788) q[0];
rz(2.4115883) q[1];
sx q[1];
rz(-2.5529824) q[1];
sx q[1];
rz(1.1013365) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0520605) q[0];
sx q[0];
rz(-0.87498481) q[0];
sx q[0];
rz(-0.3541644) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8011212) q[2];
sx q[2];
rz(-1.4061951) q[2];
sx q[2];
rz(2.4245461) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.28336542) q[1];
sx q[1];
rz(-1.4472597) q[1];
sx q[1];
rz(3.0459411) q[1];
rz(-pi) q[2];
rz(1.8227449) q[3];
sx q[3];
rz(-0.12853208) q[3];
sx q[3];
rz(2.4737308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.71658984) q[2];
sx q[2];
rz(-0.060881946) q[2];
sx q[2];
rz(1.3072183) q[2];
rz(2.763125) q[3];
sx q[3];
rz(-3.1186447) q[3];
sx q[3];
rz(2.4881261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
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
rz(-1.3017479) q[0];
sx q[0];
rz(-1.8740338) q[0];
sx q[0];
rz(0.92754716) q[0];
rz(1.357366) q[1];
sx q[1];
rz(-0.83186847) q[1];
sx q[1];
rz(-1.5313139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6283693) q[0];
sx q[0];
rz(-1.5931061) q[0];
sx q[0];
rz(1.6091899) q[0];
x q[1];
rz(-1.8157829) q[2];
sx q[2];
rz(-1.9568482) q[2];
sx q[2];
rz(2.1419142) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5591874) q[1];
sx q[1];
rz(-1.4417164) q[1];
sx q[1];
rz(-3.1388793) q[1];
rz(0.31503265) q[3];
sx q[3];
rz(-2.0612392) q[3];
sx q[3];
rz(-0.88446188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7808468) q[2];
sx q[2];
rz(-3.1372034) q[2];
sx q[2];
rz(1.2537664) q[2];
rz(0.69416657) q[3];
sx q[3];
rz(-0.73918754) q[3];
sx q[3];
rz(-2.9085801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53741443) q[0];
sx q[0];
rz(-2.1359213) q[0];
sx q[0];
rz(-1.0373212) q[0];
rz(1.6088156) q[1];
sx q[1];
rz(-0.2205801) q[1];
sx q[1];
rz(1.467009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5418422) q[0];
sx q[0];
rz(-1.6381725) q[0];
sx q[0];
rz(-1.3753433) q[0];
rz(-3.0612512) q[2];
sx q[2];
rz(-0.27896491) q[2];
sx q[2];
rz(0.34948784) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.92865366) q[1];
sx q[1];
rz(-0.0012081971) q[1];
sx q[1];
rz(-1.1986284) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9503533) q[3];
sx q[3];
rz(-0.46267022) q[3];
sx q[3];
rz(-1.3706051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.1303225) q[2];
sx q[2];
rz(-2.932817) q[2];
sx q[2];
rz(-0.043896349) q[2];
rz(-0.52055001) q[3];
sx q[3];
rz(-3.136941) q[3];
sx q[3];
rz(-1.6827778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540102) q[0];
sx q[0];
rz(-0.0024604877) q[0];
sx q[0];
rz(1.822923) q[0];
rz(1.7240546) q[1];
sx q[1];
rz(-2.8520165) q[1];
sx q[1];
rz(1.5444548) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1653571) q[0];
sx q[0];
rz(-2.6476323) q[0];
sx q[0];
rz(-0.24458347) q[0];
x q[1];
rz(2.4840646) q[2];
sx q[2];
rz(-1.3597466) q[2];
sx q[2];
rz(-1.576265) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3785951) q[1];
sx q[1];
rz(-1.7471658) q[1];
sx q[1];
rz(1.2552983) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61267743) q[3];
sx q[3];
rz(-0.78734382) q[3];
sx q[3];
rz(-1.2932216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3370207) q[2];
sx q[2];
rz(-1.8420409) q[2];
sx q[2];
rz(-2.9447832) q[2];
rz(1.192441) q[3];
sx q[3];
rz(-2.9350023) q[3];
sx q[3];
rz(2.9406252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30109677) q[0];
sx q[0];
rz(-1.7844642) q[0];
sx q[0];
rz(-1.949973) q[0];
rz(-1.6169029) q[1];
sx q[1];
rz(-0.646851) q[1];
sx q[1];
rz(-1.5651388) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3707664) q[0];
sx q[0];
rz(-1.8642555) q[0];
sx q[0];
rz(1.3956885) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0166111) q[2];
sx q[2];
rz(-2.6458394) q[2];
sx q[2];
rz(0.021136016) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8305739) q[1];
sx q[1];
rz(-1.5717447) q[1];
sx q[1];
rz(1.5702973) q[1];
x q[2];
rz(-3.0281046) q[3];
sx q[3];
rz(-1.4295414) q[3];
sx q[3];
rz(-0.065441386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9578751) q[2];
sx q[2];
rz(-0.58791939) q[2];
sx q[2];
rz(1.4584165) q[2];
rz(0.030473907) q[3];
sx q[3];
rz(-0.009549791) q[3];
sx q[3];
rz(-0.20235801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6644345) q[0];
sx q[0];
rz(-1.3344593) q[0];
sx q[0];
rz(1.6819171) q[0];
rz(1.5674113) q[1];
sx q[1];
rz(-1.3290783) q[1];
sx q[1];
rz(-3.0507416) q[1];
rz(-1.5051928) q[2];
sx q[2];
rz(-0.075724307) q[2];
sx q[2];
rz(-2.9157467) q[2];
rz(2.6371757) q[3];
sx q[3];
rz(-2.258959) q[3];
sx q[3];
rz(-2.7213982) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
