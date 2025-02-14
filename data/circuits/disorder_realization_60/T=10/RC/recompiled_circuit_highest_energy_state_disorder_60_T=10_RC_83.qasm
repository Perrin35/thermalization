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
rz(-1.634693) q[1];
sx q[1];
rz(-0.96087471) q[1];
sx q[1];
rz(-1.5009872) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.577271) q[0];
sx q[0];
rz(-2.4080896) q[0];
sx q[0];
rz(-1.8689687) q[0];
x q[1];
rz(2.9653984) q[2];
sx q[2];
rz(-1.1374047) q[2];
sx q[2];
rz(2.7769763) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.562362) q[1];
sx q[1];
rz(-0.048431245) q[1];
sx q[1];
rz(-0.50244759) q[1];
rz(-pi) q[2];
rz(-2.1330058) q[3];
sx q[3];
rz(-1.8671745) q[3];
sx q[3];
rz(2.5135771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2291439) q[2];
sx q[2];
rz(-1.392776) q[2];
sx q[2];
rz(-0.32192117) q[2];
rz(2.1866482) q[3];
sx q[3];
rz(-0.82986444) q[3];
sx q[3];
rz(1.1436852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.2442653) q[0];
sx q[0];
rz(-2.0877512) q[0];
sx q[0];
rz(2.6296997) q[0];
rz(-1.5533252) q[1];
sx q[1];
rz(-2.6542108) q[1];
sx q[1];
rz(2.0194676) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0034135) q[0];
sx q[0];
rz(-1.2597359) q[0];
sx q[0];
rz(-1.4640693) q[0];
rz(1.4830515) q[2];
sx q[2];
rz(-0.43102396) q[2];
sx q[2];
rz(1.4813678) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.304804) q[1];
sx q[1];
rz(-1.2366017) q[1];
sx q[1];
rz(-2.0242974) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.53776017) q[3];
sx q[3];
rz(-0.7825635) q[3];
sx q[3];
rz(1.7478706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.7684218) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(2.6785417) q[2];
rz(2.566346) q[3];
sx q[3];
rz(-1.7762215) q[3];
sx q[3];
rz(3.0423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7333882) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(-0.43011618) q[0];
rz(-2.6759713) q[1];
sx q[1];
rz(-2.4211113) q[1];
sx q[1];
rz(-0.63708416) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7953932) q[0];
sx q[0];
rz(-0.043248873) q[0];
sx q[0];
rz(1.9016483) q[0];
rz(-3.1153684) q[2];
sx q[2];
rz(-0.78556873) q[2];
sx q[2];
rz(0.81119591) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2807389) q[1];
sx q[1];
rz(-0.85608427) q[1];
sx q[1];
rz(1.2398941) q[1];
rz(-pi) q[2];
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
rz(0.35885262) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(-0.72506881) q[2];
rz(-1.8105761) q[3];
sx q[3];
rz(-1.015181) q[3];
sx q[3];
rz(2.5031808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8722039) q[0];
sx q[0];
rz(-1.4545472) q[0];
sx q[0];
rz(-1.697502) q[0];
rz(-3.0929502) q[1];
sx q[1];
rz(-1.3261869) q[1];
sx q[1];
rz(0.28894249) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.671593) q[0];
sx q[0];
rz(-2.6942188) q[0];
sx q[0];
rz(-0.98487052) q[0];
x q[1];
rz(1.7554531) q[2];
sx q[2];
rz(-0.15310213) q[2];
sx q[2];
rz(-2.7775922) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1565981) q[1];
sx q[1];
rz(-2.9356476) q[1];
sx q[1];
rz(1.0099645) q[1];
rz(0.069372481) q[3];
sx q[3];
rz(-0.66646229) q[3];
sx q[3];
rz(1.3795167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.85663969) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(2.7613769) q[2];
rz(2.5034261) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(-2.0130472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83288348) q[0];
sx q[0];
rz(-2.5888011) q[0];
sx q[0];
rz(2.047245) q[0];
rz(-2.265918) q[1];
sx q[1];
rz(-2.5119669) q[1];
sx q[1];
rz(1.099115) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1581061) q[0];
sx q[0];
rz(-1.9560342) q[0];
sx q[0];
rz(-1.7190821) q[0];
rz(-0.78754707) q[2];
sx q[2];
rz(-0.468245) q[2];
sx q[2];
rz(0.59216532) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.079541072) q[1];
sx q[1];
rz(-1.7289484) q[1];
sx q[1];
rz(-1.016721) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9820479) q[3];
sx q[3];
rz(-2.2189848) q[3];
sx q[3];
rz(1.2153347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8984453) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(-0.97664991) q[2];
rz(2.6099033) q[3];
sx q[3];
rz(-2.6952126) q[3];
sx q[3];
rz(2.4272052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0584745) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(1.2716768) q[0];
rz(2.7187128) q[1];
sx q[1];
rz(-1.5300749) q[1];
sx q[1];
rz(-1.5333102) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0853839) q[0];
sx q[0];
rz(-1.5041435) q[0];
sx q[0];
rz(-3.1293218) q[0];
x q[1];
rz(1.6663867) q[2];
sx q[2];
rz(-0.098473452) q[2];
sx q[2];
rz(2.5661039) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6110797) q[1];
sx q[1];
rz(-0.30834282) q[1];
sx q[1];
rz(0.18402305) q[1];
x q[2];
rz(1.3263557) q[3];
sx q[3];
rz(-0.77271739) q[3];
sx q[3];
rz(-2.6561007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61234683) q[2];
sx q[2];
rz(-1.8491448) q[2];
sx q[2];
rz(-2.7563654) q[2];
rz(-3.0569844) q[3];
sx q[3];
rz(-0.43360964) q[3];
sx q[3];
rz(0.87578526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2200634) q[0];
sx q[0];
rz(-2.0863057) q[0];
sx q[0];
rz(2.7440942) q[0];
rz(2.6037604) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(-0.0040815512) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794681) q[0];
sx q[0];
rz(-1.6982268) q[0];
sx q[0];
rz(2.3412933) q[0];
rz(2.053431) q[2];
sx q[2];
rz(-0.74649278) q[2];
sx q[2];
rz(-1.5160402) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38975484) q[1];
sx q[1];
rz(-1.0026649) q[1];
sx q[1];
rz(-2.0271432) q[1];
rz(-0.8936196) q[3];
sx q[3];
rz(-1.0056408) q[3];
sx q[3];
rz(2.2806185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5092545) q[2];
sx q[2];
rz(-1.1804487) q[2];
sx q[2];
rz(-2.9285367) q[2];
rz(-0.80900711) q[3];
sx q[3];
rz(-0.041497858) q[3];
sx q[3];
rz(1.8607148) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0179366) q[0];
sx q[0];
rz(-1.9391215) q[0];
sx q[0];
rz(-1.1630455) q[0];
rz(0.2991547) q[1];
sx q[1];
rz(-1.2048293) q[1];
sx q[1];
rz(2.1655653) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7892864) q[0];
sx q[0];
rz(-1.8160161) q[0];
sx q[0];
rz(-1.6847436) q[0];
rz(-pi) q[1];
rz(0.15248044) q[2];
sx q[2];
rz(-1.3950384) q[2];
sx q[2];
rz(-0.7712785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8638226) q[1];
sx q[1];
rz(-2.7105717) q[1];
sx q[1];
rz(0.96993877) q[1];
rz(1.3221413) q[3];
sx q[3];
rz(-1.1583405) q[3];
sx q[3];
rz(1.4506884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8396847) q[2];
sx q[2];
rz(-2.8161616) q[2];
sx q[2];
rz(3.0111266) q[2];
rz(-1.8179551) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(2.8470993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2015304) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(-1.4991722) q[0];
rz(-0.54010737) q[1];
sx q[1];
rz(-2.3441548) q[1];
sx q[1];
rz(-1.3444208) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3657065) q[0];
sx q[0];
rz(-1.2822064) q[0];
sx q[0];
rz(-1.9802753) q[0];
x q[1];
rz(0.62144582) q[2];
sx q[2];
rz(-1.1100195) q[2];
sx q[2];
rz(2.4968392) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6426797) q[1];
sx q[1];
rz(-2.1331926) q[1];
sx q[1];
rz(2.845473) q[1];
rz(-pi) q[2];
rz(-0.88560652) q[3];
sx q[3];
rz(-2.2937498) q[3];
sx q[3];
rz(2.1849439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.53238791) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(1.6877635) q[2];
rz(0.42803556) q[3];
sx q[3];
rz(-1.5513709) q[3];
sx q[3];
rz(1.9868896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55353272) q[0];
sx q[0];
rz(-2.689671) q[0];
sx q[0];
rz(-1.4916627) q[0];
rz(2.5904169) q[1];
sx q[1];
rz(-1.0124413) q[1];
sx q[1];
rz(0.59250441) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.036464036) q[0];
sx q[0];
rz(-1.291664) q[0];
sx q[0];
rz(-1.817784) q[0];
rz(-pi) q[1];
rz(-2.5985322) q[2];
sx q[2];
rz(-1.1064227) q[2];
sx q[2];
rz(2.4891702) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4372651) q[1];
sx q[1];
rz(-2.6710837) q[1];
sx q[1];
rz(0.60948845) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.09162147) q[3];
sx q[3];
rz(-2.6406248) q[3];
sx q[3];
rz(-2.9660781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.42980117) q[2];
sx q[2];
rz(-2.1000803) q[2];
sx q[2];
rz(-3.0832624) q[2];
rz(-2.77099) q[3];
sx q[3];
rz(-0.27675089) q[3];
sx q[3];
rz(-3.1378194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2881099) q[0];
sx q[0];
rz(-1.7848889) q[0];
sx q[0];
rz(-3.0369192) q[0];
rz(-1.7755605) q[1];
sx q[1];
rz(-2.3401101) q[1];
sx q[1];
rz(0.95536864) q[1];
rz(0.38717196) q[2];
sx q[2];
rz(-2.5324814) q[2];
sx q[2];
rz(-1.0003288) q[2];
rz(2.3351135) q[3];
sx q[3];
rz(-2.4462593) q[3];
sx q[3];
rz(-1.367955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
