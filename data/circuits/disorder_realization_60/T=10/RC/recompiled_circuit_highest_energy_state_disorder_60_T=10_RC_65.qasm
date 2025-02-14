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
rz(4.6484923) q[1];
sx q[1];
rz(5.3223106) q[1];
sx q[1];
rz(7.9237908) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.577271) q[0];
sx q[0];
rz(-2.4080896) q[0];
sx q[0];
rz(1.2726239) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1314279) q[2];
sx q[2];
rz(-1.7305534) q[2];
sx q[2];
rz(1.8607832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0653049) q[1];
sx q[1];
rz(-1.5283547) q[1];
sx q[1];
rz(-1.5941335) q[1];
rz(-pi) q[2];
rz(-2.7952173) q[3];
sx q[3];
rz(-2.1057671) q[3];
sx q[3];
rz(0.76081027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91244873) q[2];
sx q[2];
rz(-1.392776) q[2];
sx q[2];
rz(-2.8196715) q[2];
rz(-0.95494444) q[3];
sx q[3];
rz(-2.3117282) q[3];
sx q[3];
rz(1.9979075) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2442653) q[0];
sx q[0];
rz(-2.0877512) q[0];
sx q[0];
rz(0.51189297) q[0];
rz(1.5533252) q[1];
sx q[1];
rz(-2.6542108) q[1];
sx q[1];
rz(1.1221251) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47488362) q[0];
sx q[0];
rz(-0.32829744) q[0];
sx q[0];
rz(2.8216381) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1412337) q[2];
sx q[2];
rz(-1.5341752) q[2];
sx q[2];
rz(-0.0096732339) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1074688) q[1];
sx q[1];
rz(-1.144088) q[1];
sx q[1];
rz(0.36860768) q[1];
rz(-pi) q[2];
rz(0.70681326) q[3];
sx q[3];
rz(-1.2012831) q[3];
sx q[3];
rz(-0.57716864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.7684218) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(-0.46305099) q[2];
rz(-2.566346) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(-0.099253207) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4082044) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(0.43011618) q[0];
rz(2.6759713) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(0.63708416) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7953932) q[0];
sx q[0];
rz(-3.0983438) q[0];
sx q[0];
rz(-1.9016483) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3561959) q[2];
sx q[2];
rz(-1.5893418) q[2];
sx q[2];
rz(0.77814276) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0729941) q[1];
sx q[1];
rz(-1.322876) q[1];
sx q[1];
rz(2.3991755) q[1];
x q[2];
rz(-2.5175321) q[3];
sx q[3];
rz(-0.61316031) q[3];
sx q[3];
rz(1.1414736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.78274) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(0.72506881) q[2];
rz(-1.8105761) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(-0.63841188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8722039) q[0];
sx q[0];
rz(-1.4545472) q[0];
sx q[0];
rz(-1.697502) q[0];
rz(3.0929502) q[1];
sx q[1];
rz(-1.8154058) q[1];
sx q[1];
rz(0.28894249) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83544448) q[0];
sx q[0];
rz(-1.939538) q[0];
sx q[0];
rz(-0.25935092) q[0];
rz(1.3861395) q[2];
sx q[2];
rz(-2.9884905) q[2];
sx q[2];
rz(0.36400041) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7270941) q[1];
sx q[1];
rz(-1.3967522) q[1];
sx q[1];
rz(-3.0309309) q[1];
x q[2];
rz(-1.5163317) q[3];
sx q[3];
rz(-0.90622444) q[3];
sx q[3];
rz(-1.2913454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.85663969) q[2];
sx q[2];
rz(-2.608947) q[2];
sx q[2];
rz(0.3802158) q[2];
rz(0.63816655) q[3];
sx q[3];
rz(-1.9117982) q[3];
sx q[3];
rz(1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83288348) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(-2.047245) q[0];
rz(-2.265918) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(-2.0424776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7849784) q[0];
sx q[0];
rz(-1.4334502) q[0];
sx q[0];
rz(0.38909586) q[0];
rz(-pi) q[1];
rz(1.2266632) q[2];
sx q[2];
rz(-1.2467071) q[2];
sx q[2];
rz(-1.4365172) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2420038) q[1];
sx q[1];
rz(-2.5676641) q[1];
sx q[1];
rz(1.276488) q[1];
x q[2];
rz(-2.9820479) q[3];
sx q[3];
rz(-2.2189848) q[3];
sx q[3];
rz(1.9262579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8984453) q[2];
sx q[2];
rz(-1.3546319) q[2];
sx q[2];
rz(0.97664991) q[2];
rz(-2.6099033) q[3];
sx q[3];
rz(-2.6952126) q[3];
sx q[3];
rz(0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584745) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(1.2716768) q[0];
rz(0.42287982) q[1];
sx q[1];
rz(-1.5300749) q[1];
sx q[1];
rz(1.5333102) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0853839) q[0];
sx q[0];
rz(-1.6374491) q[0];
sx q[0];
rz(0.012270836) q[0];
rz(-pi) q[1];
rz(-3.1321636) q[2];
sx q[2];
rz(-1.4727739) q[2];
sx q[2];
rz(0.67154166) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3376006) q[1];
sx q[1];
rz(-1.873766) q[1];
sx q[1];
rz(1.5125809) q[1];
rz(-pi) q[2];
rz(1.815237) q[3];
sx q[3];
rz(-2.3688753) q[3];
sx q[3];
rz(-2.6561007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5292458) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(2.7563654) q[2];
rz(-3.0569844) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(-0.87578526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2200634) q[0];
sx q[0];
rz(-2.0863057) q[0];
sx q[0];
rz(-2.7440942) q[0];
rz(2.6037604) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(-0.0040815512) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83140552) q[0];
sx q[0];
rz(-0.80813974) q[0];
sx q[0];
rz(-0.17669295) q[0];
rz(-pi) q[1];
rz(0.88433684) q[2];
sx q[2];
rz(-1.2501653) q[2];
sx q[2];
rz(0.31242958) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7917408) q[1];
sx q[1];
rz(-2.428973) q[1];
sx q[1];
rz(0.60421677) q[1];
rz(2.3621906) q[3];
sx q[3];
rz(-0.85242295) q[3];
sx q[3];
rz(-1.8441594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6323382) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(-0.21305591) q[2];
rz(-0.80900711) q[3];
sx q[3];
rz(-0.041497858) q[3];
sx q[3];
rz(-1.2808778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1236561) q[0];
sx q[0];
rz(-1.2024711) q[0];
sx q[0];
rz(-1.9785471) q[0];
rz(0.2991547) q[1];
sx q[1];
rz(-1.9367633) q[1];
sx q[1];
rz(0.97602731) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7892864) q[0];
sx q[0];
rz(-1.3255766) q[0];
sx q[0];
rz(1.6847436) q[0];
rz(-pi) q[1];
x q[1];
rz(0.86324228) q[2];
sx q[2];
rz(-0.23216557) q[2];
sx q[2];
rz(-1.4922929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.27777003) q[1];
sx q[1];
rz(-0.43102095) q[1];
sx q[1];
rz(-0.96993877) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51234849) q[3];
sx q[3];
rz(-2.6636925) q[3];
sx q[3];
rz(1.1262788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30190793) q[2];
sx q[2];
rz(-0.32543108) q[2];
sx q[2];
rz(0.1304661) q[2];
rz(1.8179551) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(2.8470993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2015304) q[0];
sx q[0];
rz(-2.764743) q[0];
sx q[0];
rz(1.4991722) q[0];
rz(-0.54010737) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(-1.7971719) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8136052) q[0];
sx q[0];
rz(-1.9623956) q[0];
sx q[0];
rz(2.8285976) q[0];
rz(2.4355679) q[2];
sx q[2];
rz(-0.75504061) q[2];
sx q[2];
rz(1.6598827) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.089382305) q[1];
sx q[1];
rz(-1.8202413) q[1];
sx q[1];
rz(0.98813842) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62219859) q[3];
sx q[3];
rz(-0.95131627) q[3];
sx q[3];
rz(3.0752237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6092047) q[2];
sx q[2];
rz(-0.33668533) q[2];
sx q[2];
rz(-1.6877635) q[2];
rz(0.42803556) q[3];
sx q[3];
rz(-1.5513709) q[3];
sx q[3];
rz(1.9868896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353272) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(1.4916627) q[0];
rz(2.5904169) q[1];
sx q[1];
rz(-1.0124413) q[1];
sx q[1];
rz(-2.5490882) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1051286) q[0];
sx q[0];
rz(-1.291664) q[0];
sx q[0];
rz(-1.3238086) q[0];
x q[1];
rz(-0.76982381) q[2];
sx q[2];
rz(-0.69902674) q[2];
sx q[2];
rz(1.5848643) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3687835) q[1];
sx q[1];
rz(-1.9516489) q[1];
sx q[1];
rz(1.2874777) q[1];
rz(-pi) q[2];
rz(-1.5207401) q[3];
sx q[3];
rz(-2.069469) q[3];
sx q[3];
rz(-3.0704481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7117915) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(-3.0832624) q[2];
rz(0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(3.1378194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8534828) q[0];
sx q[0];
rz(-1.3567038) q[0];
sx q[0];
rz(0.10467341) q[0];
rz(-1.3660322) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(-2.568062) q[2];
sx q[2];
rz(-1.7885359) q[2];
sx q[2];
rz(0.89319695) q[2];
rz(1.0287063) q[3];
sx q[3];
rz(-1.1114612) q[3];
sx q[3];
rz(-2.3041861) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
