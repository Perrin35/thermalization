OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2296978) q[0];
sx q[0];
rz(4.3877572) q[0];
sx q[0];
rz(10.945209) q[0];
rz(-2.7819832) q[1];
sx q[1];
rz(-0.26878992) q[1];
sx q[1];
rz(0.51528817) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7069076) q[0];
sx q[0];
rz(-1.3680172) q[0];
sx q[0];
rz(2.7903656) q[0];
rz(-1.2688387) q[2];
sx q[2];
rz(-1.0450604) q[2];
sx q[2];
rz(2.5765004) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1396619) q[1];
sx q[1];
rz(-1.7822305) q[1];
sx q[1];
rz(-2.7359208) q[1];
rz(-pi) q[2];
rz(2.4042729) q[3];
sx q[3];
rz(-2.2529369) q[3];
sx q[3];
rz(-1.2350262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.23060736) q[2];
sx q[2];
rz(-1.5585941) q[2];
sx q[2];
rz(-0.67414635) q[2];
rz(2.9437183) q[3];
sx q[3];
rz(-1.2400235) q[3];
sx q[3];
rz(-1.5052634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3182217) q[0];
sx q[0];
rz(-0.30105337) q[0];
sx q[0];
rz(-0.2163042) q[0];
rz(-3.0220616) q[1];
sx q[1];
rz(-1.1263589) q[1];
sx q[1];
rz(1.4215887) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5622373) q[0];
sx q[0];
rz(-2.6174712) q[0];
sx q[0];
rz(-1.7299132) q[0];
rz(-pi) q[1];
rz(0.40070285) q[2];
sx q[2];
rz(-1.7661912) q[2];
sx q[2];
rz(0.71824441) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6232963) q[1];
sx q[1];
rz(-0.5591679) q[1];
sx q[1];
rz(-0.95018799) q[1];
x q[2];
rz(-0.98540636) q[3];
sx q[3];
rz(-0.93275654) q[3];
sx q[3];
rz(0.93321484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8372832) q[2];
sx q[2];
rz(-1.195636) q[2];
sx q[2];
rz(0.99011123) q[2];
rz(-0.75634161) q[3];
sx q[3];
rz(-2.4939311) q[3];
sx q[3];
rz(-2.9653449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
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
rz(-2.9847617) q[0];
sx q[0];
rz(-0.75060833) q[0];
sx q[0];
rz(-1.1358776) q[0];
rz(2.1319977) q[1];
sx q[1];
rz(-1.1547836) q[1];
sx q[1];
rz(-0.44581595) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2899781) q[0];
sx q[0];
rz(-1.6400684) q[0];
sx q[0];
rz(-1.4659856) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4381329) q[2];
sx q[2];
rz(-1.5477053) q[2];
sx q[2];
rz(2.9863043) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.793632) q[1];
sx q[1];
rz(-0.98274454) q[1];
sx q[1];
rz(2.6403202) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0296311) q[3];
sx q[3];
rz(-2.6918133) q[3];
sx q[3];
rz(-1.8729608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3710215) q[2];
sx q[2];
rz(-1.0022481) q[2];
sx q[2];
rz(-0.64024964) q[2];
rz(-0.69784969) q[3];
sx q[3];
rz(-1.6777439) q[3];
sx q[3];
rz(-2.8016134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8401538) q[0];
sx q[0];
rz(-0.88742632) q[0];
sx q[0];
rz(-1.3288757) q[0];
rz(-1.9245194) q[1];
sx q[1];
rz(-0.71958676) q[1];
sx q[1];
rz(-0.77635366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2503795) q[0];
sx q[0];
rz(-0.73085143) q[0];
sx q[0];
rz(2.5688085) q[0];
rz(2.4716581) q[2];
sx q[2];
rz(-2.4342854) q[2];
sx q[2];
rz(-2.7279127) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.66527077) q[1];
sx q[1];
rz(-1.5821536) q[1];
sx q[1];
rz(-1.2922056) q[1];
rz(1.5723002) q[3];
sx q[3];
rz(-2.2553911) q[3];
sx q[3];
rz(0.94438121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2886469) q[2];
sx q[2];
rz(-0.32175803) q[2];
sx q[2];
rz(-1.9256437) q[2];
rz(-2.0208042) q[3];
sx q[3];
rz(-1.8782764) q[3];
sx q[3];
rz(2.2585675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9057366) q[0];
sx q[0];
rz(-2.164542) q[0];
sx q[0];
rz(-0.46347722) q[0];
rz(-1.0889168) q[1];
sx q[1];
rz(-1.2605647) q[1];
sx q[1];
rz(-1.3190528) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70196296) q[0];
sx q[0];
rz(-1.3101398) q[0];
sx q[0];
rz(-1.8219276) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3438101) q[2];
sx q[2];
rz(-1.451865) q[2];
sx q[2];
rz(0.2695131) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8221093) q[1];
sx q[1];
rz(-2.5527375) q[1];
sx q[1];
rz(1.4134779) q[1];
rz(-1.9897668) q[3];
sx q[3];
rz(-1.4052466) q[3];
sx q[3];
rz(-2.3922684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17135458) q[2];
sx q[2];
rz(-2.3562608) q[2];
sx q[2];
rz(0.63924092) q[2];
rz(-2.3302737) q[3];
sx q[3];
rz(-2.6526484) q[3];
sx q[3];
rz(1.607224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92523471) q[0];
sx q[0];
rz(-0.41330591) q[0];
sx q[0];
rz(-2.9226724) q[0];
rz(1.9256598) q[1];
sx q[1];
rz(-2.6775807) q[1];
sx q[1];
rz(1.4930412) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1136883) q[0];
sx q[0];
rz(-2.0311715) q[0];
sx q[0];
rz(1.7229863) q[0];
x q[1];
rz(2.7240407) q[2];
sx q[2];
rz(-1.5735911) q[2];
sx q[2];
rz(-0.33139834) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5801489) q[1];
sx q[1];
rz(-2.3030048) q[1];
sx q[1];
rz(1.5931908) q[1];
rz(-0.031367112) q[3];
sx q[3];
rz(-3.0423954) q[3];
sx q[3];
rz(1.8435696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0367937) q[2];
sx q[2];
rz(-1.7918189) q[2];
sx q[2];
rz(1.8297423) q[2];
rz(-1.5051684) q[3];
sx q[3];
rz(-1.2994095) q[3];
sx q[3];
rz(-2.6817491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63603193) q[0];
sx q[0];
rz(-2.6513031) q[0];
sx q[0];
rz(-1.7946515) q[0];
rz(-1.3268283) q[1];
sx q[1];
rz(-2.3074) q[1];
sx q[1];
rz(2.7562275) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0085394) q[0];
sx q[0];
rz(-2.7811249) q[0];
sx q[0];
rz(-3.051877) q[0];
rz(-2.6527638) q[2];
sx q[2];
rz(-1.0932473) q[2];
sx q[2];
rz(-1.0350641) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24033156) q[1];
sx q[1];
rz(-1.5442532) q[1];
sx q[1];
rz(2.6955797) q[1];
rz(2.6680384) q[3];
sx q[3];
rz(-1.1692759) q[3];
sx q[3];
rz(0.36916379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4555326) q[2];
sx q[2];
rz(-0.81642381) q[2];
sx q[2];
rz(-1.653999) q[2];
rz(0.16573302) q[3];
sx q[3];
rz(-2.3159852) q[3];
sx q[3];
rz(2.9674271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97891775) q[0];
sx q[0];
rz(-2.5161777) q[0];
sx q[0];
rz(2.9130574) q[0];
rz(2.8288815) q[1];
sx q[1];
rz(-0.87930185) q[1];
sx q[1];
rz(-1.762134) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7373567) q[0];
sx q[0];
rz(-1.1690869) q[0];
sx q[0];
rz(1.6549003) q[0];
rz(0.32416043) q[2];
sx q[2];
rz(-1.7334786) q[2];
sx q[2];
rz(-1.206347) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5818565) q[1];
sx q[1];
rz(-2.7908771) q[1];
sx q[1];
rz(-1.4697671) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58016915) q[3];
sx q[3];
rz(-2.3456367) q[3];
sx q[3];
rz(-2.7027276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10869965) q[2];
sx q[2];
rz(-2.097082) q[2];
sx q[2];
rz(-1.0408939) q[2];
rz(0.4852455) q[3];
sx q[3];
rz(-1.3121366) q[3];
sx q[3];
rz(-2.7237039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8682206) q[0];
sx q[0];
rz(-2.1594248) q[0];
sx q[0];
rz(2.3305273) q[0];
rz(-1.9048994) q[1];
sx q[1];
rz(-0.86634723) q[1];
sx q[1];
rz(-2.9467357) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9962353) q[0];
sx q[0];
rz(-0.35439098) q[0];
sx q[0];
rz(-2.3217391) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5783875) q[2];
sx q[2];
rz(-2.1426472) q[2];
sx q[2];
rz(-0.88746136) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.81842717) q[1];
sx q[1];
rz(-2.0173397) q[1];
sx q[1];
rz(0.38231229) q[1];
x q[2];
rz(-0.55826346) q[3];
sx q[3];
rz(-0.9953863) q[3];
sx q[3];
rz(0.093274506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1607232) q[2];
sx q[2];
rz(-1.8941433) q[2];
sx q[2];
rz(-0.54171872) q[2];
rz(1.8187652) q[3];
sx q[3];
rz(-1.3397237) q[3];
sx q[3];
rz(-2.8790348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59170359) q[0];
sx q[0];
rz(-1.6240969) q[0];
sx q[0];
rz(2.8926335) q[0];
rz(1.9242363) q[1];
sx q[1];
rz(-1.8739506) q[1];
sx q[1];
rz(0.41044661) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0276709) q[0];
sx q[0];
rz(-0.57909896) q[0];
sx q[0];
rz(0.76580745) q[0];
rz(-pi) q[1];
rz(-3.1308799) q[2];
sx q[2];
rz(-1.6617643) q[2];
sx q[2];
rz(-2.9988924) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.527566) q[1];
sx q[1];
rz(-2.0260794) q[1];
sx q[1];
rz(-2.8705995) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66907042) q[3];
sx q[3];
rz(-2.4134372) q[3];
sx q[3];
rz(1.6403584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2534788) q[2];
sx q[2];
rz(-1.1682744) q[2];
sx q[2];
rz(-2.0237563) q[2];
rz(0.11876336) q[3];
sx q[3];
rz(-0.6898841) q[3];
sx q[3];
rz(1.9411055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73040199) q[0];
sx q[0];
rz(-1.406519) q[0];
sx q[0];
rz(2.7997959) q[0];
rz(-3.0939915) q[1];
sx q[1];
rz(-1.2536512) q[1];
sx q[1];
rz(1.8258078) q[1];
rz(-0.45817057) q[2];
sx q[2];
rz(-0.86877634) q[2];
sx q[2];
rz(1.7553117) q[2];
rz(0.66524617) q[3];
sx q[3];
rz(-0.86626296) q[3];
sx q[3];
rz(-2.5588425) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
