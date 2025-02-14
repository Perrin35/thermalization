OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.90280688) q[0];
sx q[0];
rz(-1.1385318) q[0];
sx q[0];
rz(1.0983374) q[0];
rz(1.1605473) q[1];
sx q[1];
rz(-1.1359954) q[1];
sx q[1];
rz(2.5392037) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0011604) q[0];
sx q[0];
rz(-0.70297697) q[0];
sx q[0];
rz(2.4329937) q[0];
rz(0.7135993) q[2];
sx q[2];
rz(-1.5290135) q[2];
sx q[2];
rz(3.0609727) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1225886) q[1];
sx q[1];
rz(-1.0961696) q[1];
sx q[1];
rz(-3.0130375) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1423774) q[3];
sx q[3];
rz(-1.5366293) q[3];
sx q[3];
rz(2.7673801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.41523734) q[2];
sx q[2];
rz(-2.0550315) q[2];
sx q[2];
rz(-1.497867) q[2];
rz(-0.1114397) q[3];
sx q[3];
rz(-1.8569511) q[3];
sx q[3];
rz(2.1845412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9621256) q[0];
sx q[0];
rz(-0.92573708) q[0];
sx q[0];
rz(2.977648) q[0];
rz(-2.1666849) q[1];
sx q[1];
rz(-2.4539852) q[1];
sx q[1];
rz(-1.6418537) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3079118) q[0];
sx q[0];
rz(-1.320086) q[0];
sx q[0];
rz(-2.3019019) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2664505) q[2];
sx q[2];
rz(-0.89584699) q[2];
sx q[2];
rz(0.067719134) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.837599) q[1];
sx q[1];
rz(-1.6483867) q[1];
sx q[1];
rz(-0.71240387) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1466931) q[3];
sx q[3];
rz(-2.7585976) q[3];
sx q[3];
rz(-2.8378311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.46453309) q[2];
sx q[2];
rz(-2.634282) q[2];
sx q[2];
rz(2.9110009) q[2];
rz(-2.4584037) q[3];
sx q[3];
rz(-0.69470996) q[3];
sx q[3];
rz(-2.3901239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36000073) q[0];
sx q[0];
rz(-1.8311904) q[0];
sx q[0];
rz(-0.42136425) q[0];
rz(1.5158117) q[1];
sx q[1];
rz(-0.49585626) q[1];
sx q[1];
rz(0.96744195) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28855265) q[0];
sx q[0];
rz(-1.3167736) q[0];
sx q[0];
rz(-1.2328531) q[0];
rz(-pi) q[1];
rz(-2.0386364) q[2];
sx q[2];
rz(-1.7262344) q[2];
sx q[2];
rz(1.1303678) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.9989717) q[1];
sx q[1];
rz(-1.2983067) q[1];
sx q[1];
rz(2.9914844) q[1];
rz(-pi) q[2];
rz(2.580242) q[3];
sx q[3];
rz(-1.2863942) q[3];
sx q[3];
rz(1.4903414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7385873) q[2];
sx q[2];
rz(-2.1295197) q[2];
sx q[2];
rz(2.9834413) q[2];
rz(0.75913298) q[3];
sx q[3];
rz(-1.682155) q[3];
sx q[3];
rz(-0.37091836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89449969) q[0];
sx q[0];
rz(-0.58421725) q[0];
sx q[0];
rz(1.9547113) q[0];
rz(3.0184556) q[1];
sx q[1];
rz(-1.5703399) q[1];
sx q[1];
rz(-1.8298967) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9378915) q[0];
sx q[0];
rz(-2.9726056) q[0];
sx q[0];
rz(1.0615055) q[0];
rz(1.6475296) q[2];
sx q[2];
rz(-1.6729681) q[2];
sx q[2];
rz(-2.7213784) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0614491) q[1];
sx q[1];
rz(-1.9654566) q[1];
sx q[1];
rz(0.45545642) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.00057083227) q[3];
sx q[3];
rz(-1.5667241) q[3];
sx q[3];
rz(2.3696871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19795236) q[2];
sx q[2];
rz(-1.3200878) q[2];
sx q[2];
rz(-0.8521592) q[2];
rz(0.85396829) q[3];
sx q[3];
rz(-1.2205418) q[3];
sx q[3];
rz(-2.061969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78344408) q[0];
sx q[0];
rz(-1.7572948) q[0];
sx q[0];
rz(-3.070991) q[0];
rz(-2.103503) q[1];
sx q[1];
rz(-1.1425273) q[1];
sx q[1];
rz(0.80931726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.071595777) q[0];
sx q[0];
rz(-2.5629099) q[0];
sx q[0];
rz(-3.0093638) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4368183) q[2];
sx q[2];
rz(-0.88181978) q[2];
sx q[2];
rz(-2.5344684) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.63974158) q[1];
sx q[1];
rz(-2.1110404) q[1];
sx q[1];
rz(1.8748724) q[1];
rz(-pi) q[2];
rz(1.7820939) q[3];
sx q[3];
rz(-0.70970067) q[3];
sx q[3];
rz(3.078408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.91850963) q[2];
sx q[2];
rz(-0.56879908) q[2];
sx q[2];
rz(-1.1893547) q[2];
rz(0.14063028) q[3];
sx q[3];
rz(-0.46969241) q[3];
sx q[3];
rz(-0.3895337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4024432) q[0];
sx q[0];
rz(-2.2279255) q[0];
sx q[0];
rz(1.6074578) q[0];
rz(0.88273478) q[1];
sx q[1];
rz(-0.84461707) q[1];
sx q[1];
rz(-0.61558634) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7587335) q[0];
sx q[0];
rz(-0.40399562) q[0];
sx q[0];
rz(-2.9055779) q[0];
rz(-0.83258038) q[2];
sx q[2];
rz(-0.84917779) q[2];
sx q[2];
rz(2.941212) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1751354) q[1];
sx q[1];
rz(-1.7086204) q[1];
sx q[1];
rz(-1.3786999) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.011544051) q[3];
sx q[3];
rz(-2.1911616) q[3];
sx q[3];
rz(-0.64988713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9756056) q[2];
sx q[2];
rz(-0.96419445) q[2];
sx q[2];
rz(1.1791505) q[2];
rz(2.9278582) q[3];
sx q[3];
rz(-1.390018) q[3];
sx q[3];
rz(-0.24553044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24446503) q[0];
sx q[0];
rz(-1.1470969) q[0];
sx q[0];
rz(3.097528) q[0];
rz(2.7524718) q[1];
sx q[1];
rz(-0.48946425) q[1];
sx q[1];
rz(-2.3997831) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0866068) q[0];
sx q[0];
rz(-1.9339203) q[0];
sx q[0];
rz(-1.2753049) q[0];
rz(-pi) q[1];
rz(2.7905383) q[2];
sx q[2];
rz(-2.2864955) q[2];
sx q[2];
rz(-2.9911016) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1434113) q[1];
sx q[1];
rz(-0.3340958) q[1];
sx q[1];
rz(-0.27004529) q[1];
rz(-pi) q[2];
rz(2.4398285) q[3];
sx q[3];
rz(-1.7265111) q[3];
sx q[3];
rz(0.053618535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6645708) q[2];
sx q[2];
rz(-1.7844113) q[2];
sx q[2];
rz(0.50194293) q[2];
rz(1.4255514) q[3];
sx q[3];
rz(-0.52949667) q[3];
sx q[3];
rz(2.3486923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.75941706) q[0];
sx q[0];
rz(-1.734153) q[0];
sx q[0];
rz(2.7446246) q[0];
rz(1.601903) q[1];
sx q[1];
rz(-2.868728) q[1];
sx q[1];
rz(-0.67169619) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7714193) q[0];
sx q[0];
rz(-2.015593) q[0];
sx q[0];
rz(0.14011393) q[0];
x q[1];
rz(-0.36988389) q[2];
sx q[2];
rz(-1.1544168) q[2];
sx q[2];
rz(2.2781792) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5956143) q[1];
sx q[1];
rz(-0.71450662) q[1];
sx q[1];
rz(0.48518588) q[1];
x q[2];
rz(-0.58010871) q[3];
sx q[3];
rz(-0.91843361) q[3];
sx q[3];
rz(-2.7325055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.73528618) q[2];
sx q[2];
rz(-2.6267509) q[2];
sx q[2];
rz(1.7740645) q[2];
rz(-0.96274084) q[3];
sx q[3];
rz(-1.5509501) q[3];
sx q[3];
rz(-0.65403691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.830736) q[0];
sx q[0];
rz(-1.5850569) q[0];
sx q[0];
rz(-2.8209525) q[0];
rz(2.2036208) q[1];
sx q[1];
rz(-1.9357977) q[1];
sx q[1];
rz(-1.4042312) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7787301) q[0];
sx q[0];
rz(-1.6741279) q[0];
sx q[0];
rz(-0.10424239) q[0];
rz(-pi) q[1];
rz(-0.088555468) q[2];
sx q[2];
rz(-0.27805304) q[2];
sx q[2];
rz(-1.1403569) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.93130805) q[1];
sx q[1];
rz(-3.0979994) q[1];
sx q[1];
rz(1.0924073) q[1];
rz(-pi) q[2];
rz(-1.8102989) q[3];
sx q[3];
rz(-0.82297269) q[3];
sx q[3];
rz(-2.8209958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9414901) q[2];
sx q[2];
rz(-0.76277554) q[2];
sx q[2];
rz(-2.4007559) q[2];
rz(1.0624933) q[3];
sx q[3];
rz(-2.0334358) q[3];
sx q[3];
rz(-1.9443996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.2049388) q[0];
sx q[0];
rz(-0.42977253) q[0];
sx q[0];
rz(-2.3858261) q[0];
rz(-1.8273805) q[1];
sx q[1];
rz(-1.5145489) q[1];
sx q[1];
rz(0.72602138) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7814815) q[0];
sx q[0];
rz(-1.6051634) q[0];
sx q[0];
rz(0.78256677) q[0];
rz(-3.0376932) q[2];
sx q[2];
rz(-0.88822059) q[2];
sx q[2];
rz(2.7524591) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5948702) q[1];
sx q[1];
rz(-0.38179643) q[1];
sx q[1];
rz(-1.9893655) q[1];
rz(-1.8733379) q[3];
sx q[3];
rz(-0.46756755) q[3];
sx q[3];
rz(-2.9085161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8119592) q[2];
sx q[2];
rz(-2.6128431) q[2];
sx q[2];
rz(2.192396) q[2];
rz(2.9694929) q[3];
sx q[3];
rz(-0.97730079) q[3];
sx q[3];
rz(2.7026091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1127472) q[0];
sx q[0];
rz(-1.5929359) q[0];
sx q[0];
rz(1.5449217) q[0];
rz(1.9817837) q[1];
sx q[1];
rz(-1.3216959) q[1];
sx q[1];
rz(-1.6285223) q[1];
rz(1.7261577) q[2];
sx q[2];
rz(-1.0280357) q[2];
sx q[2];
rz(0.24161777) q[2];
rz(3.0712745) q[3];
sx q[3];
rz(-1.0828154) q[3];
sx q[3];
rz(1.2630496) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
