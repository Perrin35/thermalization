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
rz(2.3693585) q[0];
sx q[0];
rz(-1.9263664) q[0];
sx q[0];
rz(1.1507432) q[0];
rz(2.8590705) q[1];
sx q[1];
rz(-1.1257659) q[1];
sx q[1];
rz(2.0292625) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5476569) q[0];
sx q[0];
rz(-0.2875244) q[0];
sx q[0];
rz(-0.42371807) q[0];
x q[1];
rz(-0.1126099) q[2];
sx q[2];
rz(-0.99163429) q[2];
sx q[2];
rz(-0.65161588) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2124651) q[1];
sx q[1];
rz(-2.0501137) q[1];
sx q[1];
rz(-2.1838837) q[1];
rz(-pi) q[2];
rz(-0.86701534) q[3];
sx q[3];
rz(-2.2261282) q[3];
sx q[3];
rz(2.1239779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.34899601) q[2];
sx q[2];
rz(-2.0670321) q[2];
sx q[2];
rz(-1.3953155) q[2];
rz(3.0830749) q[3];
sx q[3];
rz(-1.4165712) q[3];
sx q[3];
rz(1.831656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7175452) q[0];
sx q[0];
rz(-0.40322867) q[0];
sx q[0];
rz(2.6918217) q[0];
rz(-0.46478081) q[1];
sx q[1];
rz(-1.3100781) q[1];
sx q[1];
rz(-0.25258499) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16134444) q[0];
sx q[0];
rz(-0.35156116) q[0];
sx q[0];
rz(1.5241966) q[0];
rz(-1.1285601) q[2];
sx q[2];
rz(-2.1827843) q[2];
sx q[2];
rz(2.649518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.8014268) q[1];
sx q[1];
rz(-2.8599842) q[1];
sx q[1];
rz(1.5810709) q[1];
rz(2.3303495) q[3];
sx q[3];
rz(-2.0245123) q[3];
sx q[3];
rz(-2.8123901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2437336) q[2];
sx q[2];
rz(-0.89577883) q[2];
sx q[2];
rz(-0.016156999) q[2];
rz(-2.2195418) q[3];
sx q[3];
rz(-0.99370304) q[3];
sx q[3];
rz(0.13185681) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.09318) q[0];
sx q[0];
rz(-1.6572297) q[0];
sx q[0];
rz(0.68159252) q[0];
rz(0.59773481) q[1];
sx q[1];
rz(-1.455541) q[1];
sx q[1];
rz(0.32424232) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5094788) q[0];
sx q[0];
rz(-1.1369708) q[0];
sx q[0];
rz(-1.9471517) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6347972) q[2];
sx q[2];
rz(-1.629771) q[2];
sx q[2];
rz(2.5304194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8475839) q[1];
sx q[1];
rz(-1.4022744) q[1];
sx q[1];
rz(-0.042673341) q[1];
rz(-pi) q[2];
rz(-0.16894582) q[3];
sx q[3];
rz(-2.1732531) q[3];
sx q[3];
rz(1.4698616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9165667) q[2];
sx q[2];
rz(-0.39174199) q[2];
sx q[2];
rz(1.0507091) q[2];
rz(2.4010036) q[3];
sx q[3];
rz(-1.8480443) q[3];
sx q[3];
rz(3.0802687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447613) q[0];
sx q[0];
rz(-0.83503857) q[0];
sx q[0];
rz(2.187425) q[0];
rz(-1.1369811) q[1];
sx q[1];
rz(-1.515772) q[1];
sx q[1];
rz(2.383393) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9228464) q[0];
sx q[0];
rz(-2.8675666) q[0];
sx q[0];
rz(1.42204) q[0];
rz(-pi) q[1];
rz(1.3761282) q[2];
sx q[2];
rz(-0.83651453) q[2];
sx q[2];
rz(1.0502953) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0727194) q[1];
sx q[1];
rz(-2.2617635) q[1];
sx q[1];
rz(0.12703883) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.011298374) q[3];
sx q[3];
rz(-0.73888981) q[3];
sx q[3];
rz(2.9113462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.053146426) q[2];
sx q[2];
rz(-2.5909178) q[2];
sx q[2];
rz(1.3147563) q[2];
rz(-2.8407319) q[3];
sx q[3];
rz(-1.8010537) q[3];
sx q[3];
rz(0.68527591) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.143173) q[0];
sx q[0];
rz(-1.5851861) q[0];
sx q[0];
rz(-1.5013303) q[0];
rz(0.38652626) q[1];
sx q[1];
rz(-1.0685579) q[1];
sx q[1];
rz(1.5949257) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37171897) q[0];
sx q[0];
rz(-1.3981016) q[0];
sx q[0];
rz(2.4358389) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32471409) q[2];
sx q[2];
rz(-0.98718671) q[2];
sx q[2];
rz(1.8167855) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7666801) q[1];
sx q[1];
rz(-1.1385575) q[1];
sx q[1];
rz(-1.6633227) q[1];
rz(3.0647072) q[3];
sx q[3];
rz(-2.5664483) q[3];
sx q[3];
rz(-0.98371668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.74653643) q[2];
sx q[2];
rz(-2.655513) q[2];
sx q[2];
rz(2.0716095) q[2];
rz(-0.88548958) q[3];
sx q[3];
rz(-2.077379) q[3];
sx q[3];
rz(-1.989367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2162061) q[0];
sx q[0];
rz(-0.33846551) q[0];
sx q[0];
rz(2.0881407) q[0];
rz(-2.9523051) q[1];
sx q[1];
rz(-0.82866755) q[1];
sx q[1];
rz(-1.4922356) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1150619) q[0];
sx q[0];
rz(-1.4326864) q[0];
sx q[0];
rz(1.2006635) q[0];
rz(-1.8771421) q[2];
sx q[2];
rz(-1.3679081) q[2];
sx q[2];
rz(2.5032147) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.88613207) q[1];
sx q[1];
rz(-1.0696054) q[1];
sx q[1];
rz(-1.4574128) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1023472) q[3];
sx q[3];
rz(-2.3330894) q[3];
sx q[3];
rz(-1.8448585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.80387992) q[2];
sx q[2];
rz(-1.0151851) q[2];
sx q[2];
rz(1.052915) q[2];
rz(3.0247011) q[3];
sx q[3];
rz(-2.0785073) q[3];
sx q[3];
rz(1.4611999) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18148024) q[0];
sx q[0];
rz(-2.5616665) q[0];
sx q[0];
rz(2.5970698) q[0];
rz(0.25360423) q[1];
sx q[1];
rz(-0.2519775) q[1];
sx q[1];
rz(-0.25467083) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.666018) q[0];
sx q[0];
rz(-3.0694002) q[0];
sx q[0];
rz(-1.484314) q[0];
rz(-pi) q[1];
rz(-2.9355132) q[2];
sx q[2];
rz(-2.7408618) q[2];
sx q[2];
rz(-1.3361267) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0510301) q[1];
sx q[1];
rz(-1.6335618) q[1];
sx q[1];
rz(0.81845567) q[1];
rz(2.3319844) q[3];
sx q[3];
rz(-2.0008464) q[3];
sx q[3];
rz(2.3187217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2250681) q[2];
sx q[2];
rz(-0.93776408) q[2];
sx q[2];
rz(0.41898215) q[2];
rz(-0.032912832) q[3];
sx q[3];
rz(-1.907828) q[3];
sx q[3];
rz(-2.029443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58881775) q[0];
sx q[0];
rz(-0.7426312) q[0];
sx q[0];
rz(-1.8137929) q[0];
rz(2.1090419) q[1];
sx q[1];
rz(-1.3464709) q[1];
sx q[1];
rz(-0.58951497) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1731271) q[0];
sx q[0];
rz(-1.0778946) q[0];
sx q[0];
rz(-0.94287927) q[0];
rz(-2.8500798) q[2];
sx q[2];
rz(-1.8290724) q[2];
sx q[2];
rz(-1.8669389) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8679598) q[1];
sx q[1];
rz(-2.3936749) q[1];
sx q[1];
rz(-2.2272417) q[1];
rz(-0.39535268) q[3];
sx q[3];
rz(-2.0369923) q[3];
sx q[3];
rz(2.6564997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0699658) q[2];
sx q[2];
rz(-2.4667141) q[2];
sx q[2];
rz(-2.2281632) q[2];
rz(2.6895798) q[3];
sx q[3];
rz(-1.2830696) q[3];
sx q[3];
rz(-1.7970596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2695059) q[0];
sx q[0];
rz(-1.0259314) q[0];
sx q[0];
rz(-1.5103229) q[0];
rz(-1.3144846) q[1];
sx q[1];
rz(-1.038082) q[1];
sx q[1];
rz(-0.41183919) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.340419) q[0];
sx q[0];
rz(-2.1864751) q[0];
sx q[0];
rz(2.8007617) q[0];
x q[1];
rz(-2.9963295) q[2];
sx q[2];
rz(-0.32215873) q[2];
sx q[2];
rz(1.8069428) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5631521) q[1];
sx q[1];
rz(-0.65453374) q[1];
sx q[1];
rz(-2.5451979) q[1];
rz(-pi) q[2];
rz(1.8985073) q[3];
sx q[3];
rz(-1.2043038) q[3];
sx q[3];
rz(1.2372563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1371548) q[2];
sx q[2];
rz(-1.4919446) q[2];
sx q[2];
rz(0.23492661) q[2];
rz(-2.2233985) q[3];
sx q[3];
rz(-0.72489679) q[3];
sx q[3];
rz(2.0535859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
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
rz(-3.039577) q[0];
sx q[0];
rz(-0.35440847) q[0];
sx q[0];
rz(1.8812195) q[0];
rz(1.4541939) q[1];
sx q[1];
rz(-1.6834384) q[1];
sx q[1];
rz(1.6967324) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4871976) q[0];
sx q[0];
rz(-0.62614589) q[0];
sx q[0];
rz(2.3930413) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2606662) q[2];
sx q[2];
rz(-1.147715) q[2];
sx q[2];
rz(-1.0151052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9493172) q[1];
sx q[1];
rz(-2.2819464) q[1];
sx q[1];
rz(-1.7356731) q[1];
rz(-2.4244196) q[3];
sx q[3];
rz(-0.48472084) q[3];
sx q[3];
rz(-2.6074099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5626278) q[2];
sx q[2];
rz(-0.29890385) q[2];
sx q[2];
rz(-3.01801) q[2];
rz(0.28892162) q[3];
sx q[3];
rz(-1.8653423) q[3];
sx q[3];
rz(0.45946521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73744437) q[0];
sx q[0];
rz(-1.7973719) q[0];
sx q[0];
rz(1.6112882) q[0];
rz(0.96724802) q[1];
sx q[1];
rz(-2.1687242) q[1];
sx q[1];
rz(0.8868934) q[1];
rz(-3.0616765) q[2];
sx q[2];
rz(-2.5408404) q[2];
sx q[2];
rz(0.26206721) q[2];
rz(2.5206) q[3];
sx q[3];
rz(-2.1541193) q[3];
sx q[3];
rz(-0.79574058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
