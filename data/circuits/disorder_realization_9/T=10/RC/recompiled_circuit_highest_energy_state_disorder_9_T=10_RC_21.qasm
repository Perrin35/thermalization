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
rz(0.23586805) q[0];
sx q[0];
rz(-1.1455102) q[0];
sx q[0];
rz(0.048576485) q[0];
rz(-1.6254758) q[1];
sx q[1];
rz(-1.0928417) q[1];
sx q[1];
rz(-2.4258274) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5950105) q[0];
sx q[0];
rz(-1.5413741) q[0];
sx q[0];
rz(-1.741623) q[0];
rz(-pi) q[1];
rz(1.8507929) q[2];
sx q[2];
rz(-1.3716619) q[2];
sx q[2];
rz(-3.0906349) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.29808334) q[1];
sx q[1];
rz(-0.88551003) q[1];
sx q[1];
rz(-2.09649) q[1];
x q[2];
rz(1.1668901) q[3];
sx q[3];
rz(-2.0527186) q[3];
sx q[3];
rz(1.9306077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48892659) q[2];
sx q[2];
rz(-2.581614) q[2];
sx q[2];
rz(2.0860591) q[2];
rz(0.40600285) q[3];
sx q[3];
rz(-1.3136274) q[3];
sx q[3];
rz(0.82160151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15776289) q[0];
sx q[0];
rz(-1.0394179) q[0];
sx q[0];
rz(-2.5854172) q[0];
rz(1.3293386) q[1];
sx q[1];
rz(-0.61873299) q[1];
sx q[1];
rz(-0.083273085) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8502959) q[0];
sx q[0];
rz(-1.2477739) q[0];
sx q[0];
rz(2.4730541) q[0];
rz(-2.1883972) q[2];
sx q[2];
rz(-1.0824656) q[2];
sx q[2];
rz(0.40541475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1494807) q[1];
sx q[1];
rz(-0.50051033) q[1];
sx q[1];
rz(0.47187279) q[1];
rz(-pi) q[2];
rz(2.137763) q[3];
sx q[3];
rz(-1.1440082) q[3];
sx q[3];
rz(-2.4231644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9437774) q[2];
sx q[2];
rz(-1.752172) q[2];
sx q[2];
rz(-2.6488292) q[2];
rz(2.1150151) q[3];
sx q[3];
rz(-2.8387098) q[3];
sx q[3];
rz(-1.4125642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76348412) q[0];
sx q[0];
rz(-0.56270993) q[0];
sx q[0];
rz(0.43877959) q[0];
rz(-1.4525061) q[1];
sx q[1];
rz(-2.8197598) q[1];
sx q[1];
rz(2.7486393) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9071436) q[0];
sx q[0];
rz(-2.6496072) q[0];
sx q[0];
rz(-1.6299295) q[0];
rz(1.0424445) q[2];
sx q[2];
rz(-1.5990886) q[2];
sx q[2];
rz(1.9118146) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.71686983) q[1];
sx q[1];
rz(-2.7174414) q[1];
sx q[1];
rz(1.39759) q[1];
rz(-pi) q[2];
rz(0.64928253) q[3];
sx q[3];
rz(-1.6341795) q[3];
sx q[3];
rz(1.6398304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4517333) q[2];
sx q[2];
rz(-0.84443513) q[2];
sx q[2];
rz(0.27352697) q[2];
rz(0.61255974) q[3];
sx q[3];
rz(-1.4140244) q[3];
sx q[3];
rz(2.2851022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81371152) q[0];
sx q[0];
rz(-0.59309816) q[0];
sx q[0];
rz(-0.80075049) q[0];
rz(0.71209359) q[1];
sx q[1];
rz(-2.021603) q[1];
sx q[1];
rz(0.48316479) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8022435) q[0];
sx q[0];
rz(-0.74269811) q[0];
sx q[0];
rz(-0.00089641103) q[0];
x q[1];
rz(1.0739378) q[2];
sx q[2];
rz(-1.0845079) q[2];
sx q[2];
rz(-2.4118347) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2103211) q[1];
sx q[1];
rz(-0.64669791) q[1];
sx q[1];
rz(-3.0937322) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7399446) q[3];
sx q[3];
rz(-0.73026087) q[3];
sx q[3];
rz(0.83056565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.157865) q[2];
sx q[2];
rz(-0.71061504) q[2];
sx q[2];
rz(1.2845117) q[2];
rz(2.6537248) q[3];
sx q[3];
rz(-1.7043461) q[3];
sx q[3];
rz(1.454486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(2.3294285) q[0];
sx q[0];
rz(-2.282113) q[0];
sx q[0];
rz(-2.676945) q[0];
rz(0.94841415) q[1];
sx q[1];
rz(-2.3034838) q[1];
sx q[1];
rz(-1.6882187) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68935716) q[0];
sx q[0];
rz(-0.043137155) q[0];
sx q[0];
rz(2.1998911) q[0];
rz(-pi) q[1];
rz(0.74044944) q[2];
sx q[2];
rz(-1.7621104) q[2];
sx q[2];
rz(0.96093169) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.075045295) q[1];
sx q[1];
rz(-0.61457026) q[1];
sx q[1];
rz(3.0250508) q[1];
rz(2.5713163) q[3];
sx q[3];
rz(-2.3304984) q[3];
sx q[3];
rz(2.9549884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1571618) q[2];
sx q[2];
rz(-2.0941907) q[2];
sx q[2];
rz(-0.39443991) q[2];
rz(1.9858342) q[3];
sx q[3];
rz(-2.211536) q[3];
sx q[3];
rz(1.3084779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30094576) q[0];
sx q[0];
rz(-2.7257958) q[0];
sx q[0];
rz(-1.0619324) q[0];
rz(2.1901219) q[1];
sx q[1];
rz(-0.87363344) q[1];
sx q[1];
rz(-1.5207312) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5439107) q[0];
sx q[0];
rz(-1.8476227) q[0];
sx q[0];
rz(3.0863161) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1459424) q[2];
sx q[2];
rz(-1.0077927) q[2];
sx q[2];
rz(-0.90331385) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.86506063) q[1];
sx q[1];
rz(-0.78462696) q[1];
sx q[1];
rz(0.23179455) q[1];
x q[2];
rz(-0.49455182) q[3];
sx q[3];
rz(-2.7003717) q[3];
sx q[3];
rz(-2.0989204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.025853) q[2];
sx q[2];
rz(-1.226475) q[2];
sx q[2];
rz(2.8958877) q[2];
rz(-1.687441) q[3];
sx q[3];
rz(-1.5116296) q[3];
sx q[3];
rz(-0.83034849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-0.65569735) q[0];
sx q[0];
rz(-0.73807722) q[0];
sx q[0];
rz(0.6947211) q[0];
rz(3.0409536) q[1];
sx q[1];
rz(-2.1058319) q[1];
sx q[1];
rz(2.2763841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4072587) q[0];
sx q[0];
rz(-1.5293122) q[0];
sx q[0];
rz(2.7232552) q[0];
rz(-0.91036441) q[2];
sx q[2];
rz(-1.5857134) q[2];
sx q[2];
rz(-0.93176022) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.45564902) q[1];
sx q[1];
rz(-0.63243094) q[1];
sx q[1];
rz(-1.8258105) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1221625) q[3];
sx q[3];
rz(-0.97276455) q[3];
sx q[3];
rz(0.74265146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9575017) q[2];
sx q[2];
rz(-2.3881193) q[2];
sx q[2];
rz(2.2313879) q[2];
rz(1.4522067) q[3];
sx q[3];
rz(-1.2319177) q[3];
sx q[3];
rz(0.46428251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24377395) q[0];
sx q[0];
rz(-2.0938566) q[0];
sx q[0];
rz(-2.5183103) q[0];
rz(-2.9858164) q[1];
sx q[1];
rz(-2.364295) q[1];
sx q[1];
rz(-2.8782841) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0513228) q[0];
sx q[0];
rz(-1.3557845) q[0];
sx q[0];
rz(-0.010943576) q[0];
x q[1];
rz(-1.6313577) q[2];
sx q[2];
rz(-1.6471631) q[2];
sx q[2];
rz(1.5983943) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.840082) q[1];
sx q[1];
rz(-2.5599944) q[1];
sx q[1];
rz(0.14029293) q[1];
x q[2];
rz(-2.5978388) q[3];
sx q[3];
rz(-0.86287806) q[3];
sx q[3];
rz(-1.5073564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.89173633) q[2];
sx q[2];
rz(-0.24954924) q[2];
sx q[2];
rz(-0.88805324) q[2];
rz(2.213721) q[3];
sx q[3];
rz(-2.0136621) q[3];
sx q[3];
rz(-0.97308648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4230098) q[0];
sx q[0];
rz(-2.6319478) q[0];
sx q[0];
rz(2.6408559) q[0];
rz(-1.1205193) q[1];
sx q[1];
rz(-0.78103939) q[1];
sx q[1];
rz(0.80148554) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00068110739) q[0];
sx q[0];
rz(-1.6417703) q[0];
sx q[0];
rz(1.5727726) q[0];
rz(-2.5762465) q[2];
sx q[2];
rz(-1.1156811) q[2];
sx q[2];
rz(0.75789497) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7940143) q[1];
sx q[1];
rz(-1.1433351) q[1];
sx q[1];
rz(0.68750896) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3128223) q[3];
sx q[3];
rz(-1.7593062) q[3];
sx q[3];
rz(2.9303355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4215309) q[2];
sx q[2];
rz(-1.820182) q[2];
sx q[2];
rz(2.9113972) q[2];
rz(-0.0097533334) q[3];
sx q[3];
rz(-1.6246656) q[3];
sx q[3];
rz(-0.15596381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6941876) q[0];
sx q[0];
rz(-1.367584) q[0];
sx q[0];
rz(-0.76960027) q[0];
rz(0.96254483) q[1];
sx q[1];
rz(-1.8237518) q[1];
sx q[1];
rz(2.1544971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0113163) q[0];
sx q[0];
rz(-0.53142953) q[0];
sx q[0];
rz(1.7087144) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15322282) q[2];
sx q[2];
rz(-1.8109057) q[2];
sx q[2];
rz(2.3599412) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7915114) q[1];
sx q[1];
rz(-1.6120076) q[1];
sx q[1];
rz(3.1084395) q[1];
rz(-pi) q[2];
x q[2];
rz(2.065583) q[3];
sx q[3];
rz(-0.85799137) q[3];
sx q[3];
rz(3.0312169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4276245) q[2];
sx q[2];
rz(-0.93728137) q[2];
sx q[2];
rz(0.06812185) q[2];
rz(0.44975975) q[3];
sx q[3];
rz(-0.41523784) q[3];
sx q[3];
rz(-1.7033887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9925256) q[0];
sx q[0];
rz(-1.608792) q[0];
sx q[0];
rz(1.7444862) q[0];
rz(2.8515011) q[1];
sx q[1];
rz(-2.7128704) q[1];
sx q[1];
rz(-2.6147978) q[1];
rz(-1.9214966) q[2];
sx q[2];
rz(-2.0439707) q[2];
sx q[2];
rz(1.1743152) q[2];
rz(2.7431106) q[3];
sx q[3];
rz(-1.553411) q[3];
sx q[3];
rz(-0.51391272) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
