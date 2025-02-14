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
rz(3.0348294) q[0];
sx q[0];
rz(-1.3314629) q[0];
sx q[0];
rz(1.3504299) q[0];
rz(-2.8332233) q[1];
sx q[1];
rz(-2.8928533) q[1];
sx q[1];
rz(0.97919908) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.072974394) q[0];
sx q[0];
rz(-0.81356293) q[0];
sx q[0];
rz(1.6723775) q[0];
x q[1];
rz(2.2297284) q[2];
sx q[2];
rz(-1.2871337) q[2];
sx q[2];
rz(0.58667573) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1147656) q[1];
sx q[1];
rz(-1.2176759) q[1];
sx q[1];
rz(-1.4985282) q[1];
rz(-1.7448252) q[3];
sx q[3];
rz(-1.2048033) q[3];
sx q[3];
rz(-2.7470565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2778492) q[2];
sx q[2];
rz(-1.2110445) q[2];
sx q[2];
rz(0.63300526) q[2];
rz(-2.4999319) q[3];
sx q[3];
rz(-2.6810665) q[3];
sx q[3];
rz(2.3708926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68384701) q[0];
sx q[0];
rz(-0.96413833) q[0];
sx q[0];
rz(0.82474166) q[0];
rz(-1.224158) q[1];
sx q[1];
rz(-2.2861202) q[1];
sx q[1];
rz(0.32018426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.510957) q[0];
sx q[0];
rz(-0.63586006) q[0];
sx q[0];
rz(0.73131928) q[0];
rz(-pi) q[1];
rz(-1.2037781) q[2];
sx q[2];
rz(-0.53897714) q[2];
sx q[2];
rz(2.0808329) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9588543) q[1];
sx q[1];
rz(-3.0575648) q[1];
sx q[1];
rz(-2.8699204) q[1];
x q[2];
rz(-1.8289964) q[3];
sx q[3];
rz(-0.22919433) q[3];
sx q[3];
rz(2.3865139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0188633) q[2];
sx q[2];
rz(-1.9247232) q[2];
sx q[2];
rz(-1.2790206) q[2];
rz(-0.095712885) q[3];
sx q[3];
rz(-1.0749823) q[3];
sx q[3];
rz(-1.8224645) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1529481) q[0];
sx q[0];
rz(-2.4929292) q[0];
sx q[0];
rz(1.0308107) q[0];
rz(-2.5910494) q[1];
sx q[1];
rz(-2.2830453) q[1];
sx q[1];
rz(1.8969089) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7787832) q[0];
sx q[0];
rz(-0.16805695) q[0];
sx q[0];
rz(-0.3248949) q[0];
x q[1];
rz(2.9739506) q[2];
sx q[2];
rz(-2.1812005) q[2];
sx q[2];
rz(2.8053225) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0708187) q[1];
sx q[1];
rz(-0.70291513) q[1];
sx q[1];
rz(1.2513112) q[1];
rz(-pi) q[2];
rz(-2.3233285) q[3];
sx q[3];
rz(-1.8570559) q[3];
sx q[3];
rz(0.32901007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9024258) q[2];
sx q[2];
rz(-2.5776165) q[2];
sx q[2];
rz(1.9695367) q[2];
rz(-2.3160882) q[3];
sx q[3];
rz(-1.8768616) q[3];
sx q[3];
rz(-0.66506213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.091752) q[0];
sx q[0];
rz(-1.6850152) q[0];
sx q[0];
rz(0.41393429) q[0];
rz(2.6992758) q[1];
sx q[1];
rz(-2.1374173) q[1];
sx q[1];
rz(-0.54534674) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53876153) q[0];
sx q[0];
rz(-1.6206546) q[0];
sx q[0];
rz(-0.2754579) q[0];
rz(-pi) q[1];
x q[1];
rz(0.92939922) q[2];
sx q[2];
rz(-0.9050194) q[2];
sx q[2];
rz(-2.5432472) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0132349) q[1];
sx q[1];
rz(-1.2942794) q[1];
sx q[1];
rz(2.1091631) q[1];
rz(-pi) q[2];
rz(-1.8319301) q[3];
sx q[3];
rz(-0.52370893) q[3];
sx q[3];
rz(-1.6317489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7004194) q[2];
sx q[2];
rz(-1.8808695) q[2];
sx q[2];
rz(2.9803661) q[2];
rz(-1.553933) q[3];
sx q[3];
rz(-0.59993887) q[3];
sx q[3];
rz(2.5509295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0597543) q[0];
sx q[0];
rz(-1.8738382) q[0];
sx q[0];
rz(-0.72342122) q[0];
rz(1.8266034) q[1];
sx q[1];
rz(-1.2774717) q[1];
sx q[1];
rz(2.1037197) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4525131) q[0];
sx q[0];
rz(-2.291579) q[0];
sx q[0];
rz(1.4656653) q[0];
rz(-pi) q[1];
rz(-0.60836253) q[2];
sx q[2];
rz(-2.3181097) q[2];
sx q[2];
rz(-1.057511) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.41748504) q[1];
sx q[1];
rz(-1.3032262) q[1];
sx q[1];
rz(-0.45326434) q[1];
rz(-pi) q[2];
rz(-2.9014189) q[3];
sx q[3];
rz(-2.5824912) q[3];
sx q[3];
rz(-0.19776519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7448685) q[2];
sx q[2];
rz(-2.0266271) q[2];
sx q[2];
rz(0.57207668) q[2];
rz(2.5165596) q[3];
sx q[3];
rz(-2.1376938) q[3];
sx q[3];
rz(-1.9151789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52142757) q[0];
sx q[0];
rz(-3.1075952) q[0];
sx q[0];
rz(0.52325621) q[0];
rz(1.322586) q[1];
sx q[1];
rz(-1.4679642) q[1];
sx q[1];
rz(-0.62017131) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70411982) q[0];
sx q[0];
rz(-1.0123555) q[0];
sx q[0];
rz(-2.0790786) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6441395) q[2];
sx q[2];
rz(-2.0383325) q[2];
sx q[2];
rz(-2.1852124) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0508681) q[1];
sx q[1];
rz(-1.4830989) q[1];
sx q[1];
rz(2.1983653) q[1];
rz(-pi) q[2];
rz(0.2175339) q[3];
sx q[3];
rz(-1.399588) q[3];
sx q[3];
rz(1.273716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.72914499) q[2];
sx q[2];
rz(-2.3690988) q[2];
sx q[2];
rz(-1.8577417) q[2];
rz(-2.1733952) q[3];
sx q[3];
rz(-1.7627534) q[3];
sx q[3];
rz(-1.5865954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17484434) q[0];
sx q[0];
rz(-2.0059678) q[0];
sx q[0];
rz(1.9298166) q[0];
rz(-0.9696331) q[1];
sx q[1];
rz(-1.8200834) q[1];
sx q[1];
rz(1.8671573) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2454555) q[0];
sx q[0];
rz(-1.4198167) q[0];
sx q[0];
rz(-2.2666988) q[0];
rz(-pi) q[1];
rz(-1.6658804) q[2];
sx q[2];
rz(-0.78721744) q[2];
sx q[2];
rz(1.9958391) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.74531534) q[1];
sx q[1];
rz(-2.4419129) q[1];
sx q[1];
rz(-1.0902283) q[1];
x q[2];
rz(-0.59341615) q[3];
sx q[3];
rz(-0.53916603) q[3];
sx q[3];
rz(2.9447458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4659287) q[2];
sx q[2];
rz(-1.482654) q[2];
sx q[2];
rz(-0.3234123) q[2];
rz(-0.0023500738) q[3];
sx q[3];
rz(-1.6833143) q[3];
sx q[3];
rz(0.6215483) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49331409) q[0];
sx q[0];
rz(-1.6596153) q[0];
sx q[0];
rz(1.7907273) q[0];
rz(-1.3510652) q[1];
sx q[1];
rz(-1.2850782) q[1];
sx q[1];
rz(-1.8538808) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13437102) q[0];
sx q[0];
rz(-2.8343081) q[0];
sx q[0];
rz(1.968712) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2852934) q[2];
sx q[2];
rz(-1.3494455) q[2];
sx q[2];
rz(0.75698392) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.334633) q[1];
sx q[1];
rz(-1.8939085) q[1];
sx q[1];
rz(-1.3912702) q[1];
rz(-3.0898422) q[3];
sx q[3];
rz(-2.7922575) q[3];
sx q[3];
rz(-2.3037095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6172341) q[2];
sx q[2];
rz(-1.1523767) q[2];
sx q[2];
rz(-0.99346811) q[2];
rz(1.0348381) q[3];
sx q[3];
rz(-2.5351758) q[3];
sx q[3];
rz(-0.16916999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33191037) q[0];
sx q[0];
rz(-1.9669635) q[0];
sx q[0];
rz(1.5274973) q[0];
rz(2.2612259) q[1];
sx q[1];
rz(-0.4207193) q[1];
sx q[1];
rz(-2.8299433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051119251) q[0];
sx q[0];
rz(-1.3731806) q[0];
sx q[0];
rz(1.86905) q[0];
x q[1];
rz(0.016383532) q[2];
sx q[2];
rz(-1.8021305) q[2];
sx q[2];
rz(2.6838059) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.27585718) q[1];
sx q[1];
rz(-1.9364662) q[1];
sx q[1];
rz(-2.4128298) q[1];
rz(1.9945456) q[3];
sx q[3];
rz(-0.69931385) q[3];
sx q[3];
rz(2.5831403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6722022) q[2];
sx q[2];
rz(-0.92350525) q[2];
sx q[2];
rz(-1.7507318) q[2];
rz(1.6145128) q[3];
sx q[3];
rz(-1.047784) q[3];
sx q[3];
rz(-1.0745777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-0.63956082) q[0];
sx q[0];
rz(-1.1447516) q[0];
sx q[0];
rz(-0.61258739) q[0];
rz(2.1514429) q[1];
sx q[1];
rz(-1.1724816) q[1];
sx q[1];
rz(1.490907) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25967071) q[0];
sx q[0];
rz(-2.2380191) q[0];
sx q[0];
rz(1.0259969) q[0];
rz(-pi) q[1];
rz(2.886495) q[2];
sx q[2];
rz(-0.92689415) q[2];
sx q[2];
rz(-1.7718499) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.1742434) q[1];
sx q[1];
rz(-2.9799906) q[1];
sx q[1];
rz(-1.6095649) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8464523) q[3];
sx q[3];
rz(-0.71986976) q[3];
sx q[3];
rz(0.9269549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.426173) q[2];
sx q[2];
rz(-0.82948589) q[2];
sx q[2];
rz(-0.034991525) q[2];
rz(1.4404826) q[3];
sx q[3];
rz(-2.248843) q[3];
sx q[3];
rz(-0.54097241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4424425) q[0];
sx q[0];
rz(-0.89542605) q[0];
sx q[0];
rz(-0.17740346) q[0];
rz(-1.9433446) q[1];
sx q[1];
rz(-1.5180963) q[1];
sx q[1];
rz(0.95388283) q[1];
rz(2.162355) q[2];
sx q[2];
rz(-0.981058) q[2];
sx q[2];
rz(-2.2997591) q[2];
rz(-1.1097601) q[3];
sx q[3];
rz(-0.90519917) q[3];
sx q[3];
rz(-1.7310033) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
