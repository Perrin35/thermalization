OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(5.3263721) q[0];
sx q[0];
rz(10.858067) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8960436) q[0];
sx q[0];
rz(-2.890812) q[0];
sx q[0];
rz(1.8401237) q[0];
rz(-1.990591) q[2];
sx q[2];
rz(-1.8325873) q[2];
sx q[2];
rz(-0.37444886) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6822328) q[1];
sx q[1];
rz(-1.1182922) q[1];
sx q[1];
rz(-3.1205936) q[1];
x q[2];
rz(1.8688698) q[3];
sx q[3];
rz(-1.0469336) q[3];
sx q[3];
rz(0.03447547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0248489) q[2];
sx q[2];
rz(-1.8086834) q[2];
sx q[2];
rz(-1.2343181) q[2];
rz(-1.0553137) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(-1.9887259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18594436) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(2.3117075) q[0];
rz(-0.25289598) q[1];
sx q[1];
rz(-2.3650832) q[1];
sx q[1];
rz(-0.84709644) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0704437) q[0];
sx q[0];
rz(-1.1304454) q[0];
sx q[0];
rz(2.0106993) q[0];
rz(2.3624624) q[2];
sx q[2];
rz(-1.7841633) q[2];
sx q[2];
rz(1.2029635) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4699128) q[1];
sx q[1];
rz(-2.5922853) q[1];
sx q[1];
rz(0.13110199) q[1];
rz(-pi) q[2];
x q[2];
rz(2.794572) q[3];
sx q[3];
rz(-1.1565398) q[3];
sx q[3];
rz(1.8300717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0236686) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(0.71933293) q[2];
rz(1.2891278) q[3];
sx q[3];
rz(-0.23580655) q[3];
sx q[3];
rz(-1.4891362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711202) q[0];
sx q[0];
rz(-1.2501165) q[0];
sx q[0];
rz(-2.5701994) q[0];
rz(2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(-2.6142696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8674832) q[0];
sx q[0];
rz(-1.880078) q[0];
sx q[0];
rz(1.5779737) q[0];
rz(-pi) q[1];
rz(-0.6730404) q[2];
sx q[2];
rz(-1.3125784) q[2];
sx q[2];
rz(-0.25707993) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5273683) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(2.7781092) q[1];
rz(-2.4089291) q[3];
sx q[3];
rz(-1.3244197) q[3];
sx q[3];
rz(-0.25783595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30806914) q[2];
sx q[2];
rz(-1.6922733) q[2];
sx q[2];
rz(2.4148338) q[2];
rz(-0.99749342) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(-1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8961287) q[0];
sx q[0];
rz(-1.4828869) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(2.3153268) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4429312) q[0];
sx q[0];
rz(-2.1974265) q[0];
sx q[0];
rz(0.71039623) q[0];
rz(-1.516953) q[2];
sx q[2];
rz(-0.87140897) q[2];
sx q[2];
rz(2.2211423) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5845881) q[1];
sx q[1];
rz(-1.4950206) q[1];
sx q[1];
rz(1.9083379) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6310001) q[3];
sx q[3];
rz(-0.69805745) q[3];
sx q[3];
rz(-1.5654246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.56759175) q[2];
sx q[2];
rz(-1.8969798) q[2];
sx q[2];
rz(-2.2272002) q[2];
rz(-3.0363723) q[3];
sx q[3];
rz(-1.5172232) q[3];
sx q[3];
rz(2.5523394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88702622) q[0];
sx q[0];
rz(-1.5019324) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(0.0056313593) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(-2.5240135) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3546346) q[0];
sx q[0];
rz(-1.600453) q[0];
sx q[0];
rz(-1.6164854) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1281625) q[2];
sx q[2];
rz(-0.97634041) q[2];
sx q[2];
rz(0.15304676) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4434112) q[1];
sx q[1];
rz(-0.65084208) q[1];
sx q[1];
rz(2.0425668) q[1];
x q[2];
rz(0.072237416) q[3];
sx q[3];
rz(-0.17599711) q[3];
sx q[3];
rz(-0.76398677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(-2.9079672) q[2];
rz(-0.99308333) q[3];
sx q[3];
rz(-2.1916094) q[3];
sx q[3];
rz(0.63794199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96930209) q[0];
sx q[0];
rz(-3.0811221) q[0];
sx q[0];
rz(-3.0517975) q[0];
rz(-0.88109294) q[1];
sx q[1];
rz(-2.6125364) q[1];
sx q[1];
rz(-2.9615013) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0463379) q[0];
sx q[0];
rz(-2.1349847) q[0];
sx q[0];
rz(0.5395704) q[0];
rz(2.9526268) q[2];
sx q[2];
rz(-0.97550387) q[2];
sx q[2];
rz(-1.8408066) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.99299586) q[1];
sx q[1];
rz(-2.8088048) q[1];
sx q[1];
rz(-2.6246043) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0780219) q[3];
sx q[3];
rz(-1.5544976) q[3];
sx q[3];
rz(-1.8340045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5774612) q[2];
sx q[2];
rz(1.1759261) q[2];
rz(0.073143395) q[3];
sx q[3];
rz(-0.72435838) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(-1.1649959) q[1];
sx q[1];
rz(-0.55924758) q[1];
sx q[1];
rz(-0.040239008) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2049853) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(-1.3315014) q[0];
rz(1.6778498) q[2];
sx q[2];
rz(-2.7333626) q[2];
sx q[2];
rz(1.5232616) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7525967) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(-0.1715626) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3761446) q[3];
sx q[3];
rz(-0.42740373) q[3];
sx q[3];
rz(-3.1020853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-2.3562727) q[2];
sx q[2];
rz(0.17383943) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(-0.85913908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180856) q[0];
sx q[0];
rz(-2.9243587) q[0];
sx q[0];
rz(-1.7370976) q[0];
rz(1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(1.5054024) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0818427) q[0];
sx q[0];
rz(-1.8847701) q[0];
sx q[0];
rz(2.7008675) q[0];
rz(2.3178188) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(1.6549695) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1464403) q[1];
sx q[1];
rz(-0.69680981) q[1];
sx q[1];
rz(3.0909096) q[1];
x q[2];
rz(-0.96945854) q[3];
sx q[3];
rz(-1.2974032) q[3];
sx q[3];
rz(2.6858342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3796842) q[2];
sx q[2];
rz(-2.6613993) q[2];
sx q[2];
rz(2.9910679) q[2];
rz(1.5395509) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(0.62301821) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6983011) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(-0.64569965) q[0];
rz(0.49939108) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(-1.8766778) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7495959) q[0];
sx q[0];
rz(-2.2889334) q[0];
sx q[0];
rz(-1.0438265) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0090946) q[2];
sx q[2];
rz(-1.0602789) q[2];
sx q[2];
rz(0.27062182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1688655) q[1];
sx q[1];
rz(-1.7164413) q[1];
sx q[1];
rz(1.5345854) q[1];
rz(-pi) q[2];
rz(2.7154269) q[3];
sx q[3];
rz(-0.45939547) q[3];
sx q[3];
rz(0.029749425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.50679961) q[2];
sx q[2];
rz(-2.1719077) q[2];
sx q[2];
rz(2.6169422) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.5126712) q[3];
sx q[3];
rz(-2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5234579) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(2.8961704) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(-1.6773178) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2357764) q[0];
sx q[0];
rz(-0.24484867) q[0];
sx q[0];
rz(2.66923) q[0];
rz(-pi) q[1];
rz(0.22428959) q[2];
sx q[2];
rz(-0.99594342) q[2];
sx q[2];
rz(-1.526051) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6357395) q[1];
sx q[1];
rz(-2.3931599) q[1];
sx q[1];
rz(1.9745419) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7467935) q[3];
sx q[3];
rz(-1.4868075) q[3];
sx q[3];
rz(-1.2916816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2422553) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(2.396092) q[2];
rz(-1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(1.1283114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2765008) q[0];
sx q[0];
rz(-1.6456974) q[0];
sx q[0];
rz(0.67281848) q[0];
rz(-1.8854234) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(2.0586661) q[2];
sx q[2];
rz(-0.91960533) q[2];
sx q[2];
rz(-1.4056924) q[2];
rz(1.1938548) q[3];
sx q[3];
rz(-2.1115163) q[3];
sx q[3];
rz(1.9797309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
