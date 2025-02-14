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
rz(0.76437104) q[0];
sx q[0];
rz(-1.3410913) q[0];
sx q[0];
rz(2.2361225) q[0];
rz(1.5068997) q[1];
sx q[1];
rz(-2.1807179) q[1];
sx q[1];
rz(-1.6406055) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.577271) q[0];
sx q[0];
rz(-0.73350302) q[0];
sx q[0];
rz(-1.8689687) q[0];
rz(2.9653984) q[2];
sx q[2];
rz(-1.1374047) q[2];
sx q[2];
rz(-0.36461634) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6480744) q[1];
sx q[1];
rz(-1.5941125) q[1];
sx q[1];
rz(-3.0991395) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.34637536) q[3];
sx q[3];
rz(-1.0358255) q[3];
sx q[3];
rz(0.76081027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.91244873) q[2];
sx q[2];
rz(-1.392776) q[2];
sx q[2];
rz(0.32192117) q[2];
rz(2.1866482) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(1.9979075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8973273) q[0];
sx q[0];
rz(-2.0877512) q[0];
sx q[0];
rz(0.51189297) q[0];
rz(1.5533252) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(-1.1221251) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0034135) q[0];
sx q[0];
rz(-1.8818568) q[0];
sx q[0];
rz(1.6775234) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4830515) q[2];
sx q[2];
rz(-0.43102396) q[2];
sx q[2];
rz(-1.4813678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0341238) q[1];
sx q[1];
rz(-1.9975047) q[1];
sx q[1];
rz(-0.36860768) q[1];
x q[2];
rz(2.4347794) q[3];
sx q[3];
rz(-1.9403096) q[3];
sx q[3];
rz(-0.57716864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3731709) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(-0.46305099) q[2];
rz(0.57524663) q[3];
sx q[3];
rz(-1.7762215) q[3];
sx q[3];
rz(-3.0423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4082044) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(-0.43011618) q[0];
rz(-0.46562132) q[1];
sx q[1];
rz(-0.72048134) q[1];
sx q[1];
rz(2.5045085) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4642536) q[0];
sx q[0];
rz(-1.6116983) q[0];
sx q[0];
rz(-3.1275355) q[0];
rz(-pi) q[1];
rz(1.5970206) q[2];
sx q[2];
rz(-2.3560212) q[2];
sx q[2];
rz(-2.3674813) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7635258) q[1];
sx q[1];
rz(-0.77516205) q[1];
sx q[1];
rz(0.35825348) q[1];
rz(-2.5175321) q[3];
sx q[3];
rz(-0.61316031) q[3];
sx q[3];
rz(1.1414736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35885262) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(0.72506881) q[2];
rz(1.3310165) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(-0.63841188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8722039) q[0];
sx q[0];
rz(-1.4545472) q[0];
sx q[0];
rz(1.697502) q[0];
rz(-0.048642453) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(-0.28894249) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4699997) q[0];
sx q[0];
rz(-2.6942188) q[0];
sx q[0];
rz(2.1567221) q[0];
rz(-1.7554531) q[2];
sx q[2];
rz(-0.15310213) q[2];
sx q[2];
rz(-0.36400041) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1565981) q[1];
sx q[1];
rz(-0.20594507) q[1];
sx q[1];
rz(1.0099645) q[1];
rz(-1.5163317) q[3];
sx q[3];
rz(-0.90622444) q[3];
sx q[3];
rz(1.8502473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.85663969) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(2.7613769) q[2];
rz(-2.5034261) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(2.0130472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3087092) q[0];
sx q[0];
rz(-2.5888011) q[0];
sx q[0];
rz(2.047245) q[0];
rz(-0.87567466) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(2.0424776) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9834866) q[0];
sx q[0];
rz(-1.9560342) q[0];
sx q[0];
rz(1.7190821) q[0];
rz(-pi) q[1];
rz(1.9149295) q[2];
sx q[2];
rz(-1.8948855) q[2];
sx q[2];
rz(-1.4365172) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2420038) q[1];
sx q[1];
rz(-2.5676641) q[1];
sx q[1];
rz(1.8651047) q[1];
rz(-pi) q[2];
rz(-1.7775675) q[3];
sx q[3];
rz(-2.4768157) q[3];
sx q[3];
rz(2.1867276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8984453) q[2];
sx q[2];
rz(-1.3546319) q[2];
sx q[2];
rz(-0.97664991) q[2];
rz(-2.6099033) q[3];
sx q[3];
rz(-2.6952126) q[3];
sx q[3];
rz(0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831182) q[0];
sx q[0];
rz(-2.0396905) q[0];
sx q[0];
rz(-1.8699159) q[0];
rz(2.7187128) q[1];
sx q[1];
rz(-1.6115178) q[1];
sx q[1];
rz(1.5333102) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0853839) q[0];
sx q[0];
rz(-1.6374491) q[0];
sx q[0];
rz(0.012270836) q[0];
rz(-pi) q[1];
rz(3.1321636) q[2];
sx q[2];
rz(-1.6688188) q[2];
sx q[2];
rz(0.67154166) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.530513) q[1];
sx q[1];
rz(-2.8332498) q[1];
sx q[1];
rz(0.18402305) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9098784) q[3];
sx q[3];
rz(-2.3149256) q[3];
sx q[3];
rz(0.8207013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.61234683) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(0.38522729) q[2];
rz(-3.0569844) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(2.2658074) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92152921) q[0];
sx q[0];
rz(-2.0863057) q[0];
sx q[0];
rz(0.39749843) q[0];
rz(2.6037604) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(-0.0040815512) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86212457) q[0];
sx q[0];
rz(-1.6982268) q[0];
sx q[0];
rz(-2.3412933) q[0];
rz(2.7360544) q[2];
sx q[2];
rz(-0.92541646) q[2];
sx q[2];
rz(-2.1359512) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4392802) q[1];
sx q[1];
rz(-1.1902307) q[1];
sx q[1];
rz(0.61813942) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8936196) q[3];
sx q[3];
rz(-2.1359518) q[3];
sx q[3];
rz(0.86097417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5092545) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(0.21305591) q[2];
rz(0.80900711) q[3];
sx q[3];
rz(-0.041497858) q[3];
sx q[3];
rz(1.2808778) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1236561) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(-1.9785471) q[0];
rz(0.2991547) q[1];
sx q[1];
rz(-1.2048293) q[1];
sx q[1];
rz(2.1655653) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9508787) q[0];
sx q[0];
rz(-1.6813206) q[0];
sx q[0];
rz(2.8948363) q[0];
x q[1];
rz(1.7485745) q[2];
sx q[2];
rz(-1.7209098) q[2];
sx q[2];
rz(0.77265384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.4055172) q[1];
sx q[1];
rz(-1.3323405) q[1];
sx q[1];
rz(-1.9333436) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6292442) q[3];
sx q[3];
rz(-0.47790018) q[3];
sx q[3];
rz(1.1262788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.30190793) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(0.1304661) q[2];
rz(1.8179551) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2015304) q[0];
sx q[0];
rz(-2.764743) q[0];
sx q[0];
rz(1.4991722) q[0];
rz(2.6014853) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(1.3444208) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3657065) q[0];
sx q[0];
rz(-1.2822064) q[0];
sx q[0];
rz(-1.9802753) q[0];
x q[1];
rz(-0.70602472) q[2];
sx q[2];
rz(-0.75504061) q[2];
sx q[2];
rz(1.6598827) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0186179) q[1];
sx q[1];
rz(-2.5135165) q[1];
sx q[1];
rz(2.0043892) q[1];
rz(-pi) q[2];
rz(-0.62219859) q[3];
sx q[3];
rz(-2.1902764) q[3];
sx q[3];
rz(-0.066368997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53238791) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(-1.6877635) q[2];
rz(0.42803556) q[3];
sx q[3];
rz(-1.5902218) q[3];
sx q[3];
rz(-1.9868896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353272) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(-1.6499299) q[0];
rz(2.5904169) q[1];
sx q[1];
rz(-2.1291514) q[1];
sx q[1];
rz(2.5490882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036464036) q[0];
sx q[0];
rz(-1.8499287) q[0];
sx q[0];
rz(1.817784) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.76982381) q[2];
sx q[2];
rz(-0.69902674) q[2];
sx q[2];
rz(1.5848643) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4514102) q[1];
sx q[1];
rz(-1.8333149) q[1];
sx q[1];
rz(2.7464944) q[1];
rz(2.6423934) q[3];
sx q[3];
rz(-1.6147524) q[3];
sx q[3];
rz(-1.4756965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.42980117) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(-3.0832624) q[2];
rz(0.37060261) q[3];
sx q[3];
rz(-0.27675089) q[3];
sx q[3];
rz(0.0037732865) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8534828) q[0];
sx q[0];
rz(-1.3567038) q[0];
sx q[0];
rz(0.10467341) q[0];
rz(-1.7755605) q[1];
sx q[1];
rz(-2.3401101) q[1];
sx q[1];
rz(0.95536864) q[1];
rz(1.8283394) q[2];
sx q[2];
rz(-1.0124442) q[2];
sx q[2];
rz(-0.53895216) q[2];
rz(0.52363734) q[3];
sx q[3];
rz(-1.0900396) q[3];
sx q[3];
rz(-0.47245792) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
