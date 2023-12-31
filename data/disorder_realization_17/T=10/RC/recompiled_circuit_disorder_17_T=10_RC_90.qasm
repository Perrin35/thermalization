OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2330115) q[0];
sx q[0];
rz(-1.1865948) q[0];
sx q[0];
rz(-1.9957805) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(-0.29247984) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8669517) q[0];
sx q[0];
rz(-0.39817087) q[0];
sx q[0];
rz(-1.4273594) q[0];
x q[1];
rz(-0.05188611) q[2];
sx q[2];
rz(-2.3615026) q[2];
sx q[2];
rz(2.8442596) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18947345) q[1];
sx q[1];
rz(-1.5249624) q[1];
sx q[1];
rz(-2.1650251) q[1];
x q[2];
rz(-3.1396477) q[3];
sx q[3];
rz(-1.0969775) q[3];
sx q[3];
rz(2.4634944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4108489) q[2];
sx q[2];
rz(-1.0935254) q[2];
sx q[2];
rz(0.4494108) q[2];
rz(2.4959026) q[3];
sx q[3];
rz(-0.68104762) q[3];
sx q[3];
rz(-1.8507563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9280424) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(-2.399562) q[0];
rz(1.4713326) q[1];
sx q[1];
rz(-2.4447618) q[1];
sx q[1];
rz(1.3630294) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8105158) q[0];
sx q[0];
rz(-1.9609945) q[0];
sx q[0];
rz(2.4742728) q[0];
rz(-pi) q[1];
rz(-1.7861373) q[2];
sx q[2];
rz(-1.1110348) q[2];
sx q[2];
rz(0.41972566) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0568309) q[1];
sx q[1];
rz(-2.8014604) q[1];
sx q[1];
rz(3.0840193) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8369254) q[3];
sx q[3];
rz(-1.7566655) q[3];
sx q[3];
rz(-0.64418018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.30119511) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(-1.2724686) q[2];
rz(2.8043591) q[3];
sx q[3];
rz(-1.1703346) q[3];
sx q[3];
rz(-0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.067588016) q[0];
sx q[0];
rz(-0.33960605) q[0];
sx q[0];
rz(0.2628251) q[0];
rz(2.7858531) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(1.9062818) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4848029) q[0];
sx q[0];
rz(-3.1233388) q[0];
sx q[0];
rz(2.546893) q[0];
rz(-pi) q[1];
rz(1.6791061) q[2];
sx q[2];
rz(-0.24178594) q[2];
sx q[2];
rz(-1.0847278) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1097818) q[1];
sx q[1];
rz(-2.8019252) q[1];
sx q[1];
rz(-0.32082816) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1615109) q[3];
sx q[3];
rz(-2.403879) q[3];
sx q[3];
rz(-2.0091332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.7586781) q[2];
sx q[2];
rz(0.18033218) q[2];
rz(-2.5056433) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(0.37160555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3746049) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(0.71907991) q[0];
rz(0.51302296) q[1];
sx q[1];
rz(-1.2027556) q[1];
sx q[1];
rz(1.0346574) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35095222) q[0];
sx q[0];
rz(-1.6296903) q[0];
sx q[0];
rz(-1.590593) q[0];
rz(1.9070508) q[2];
sx q[2];
rz(-0.14306919) q[2];
sx q[2];
rz(0.87749315) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0532916) q[1];
sx q[1];
rz(-2.3846855) q[1];
sx q[1];
rz(-0.41815586) q[1];
x q[2];
rz(-1.3961117) q[3];
sx q[3];
rz(-2.0377199) q[3];
sx q[3];
rz(-0.13388453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2085312) q[2];
sx q[2];
rz(-0.26991093) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(-1.158372) q[3];
sx q[3];
rz(-1.002243) q[3];
sx q[3];
rz(2.2896144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-0.029595705) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(2.239256) q[0];
rz(2.2205655) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(2.5193118) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8392004) q[0];
sx q[0];
rz(-1.3883739) q[0];
sx q[0];
rz(-1.8586041) q[0];
rz(-pi) q[1];
rz(3.0201966) q[2];
sx q[2];
rz(-0.65181323) q[2];
sx q[2];
rz(-0.35169841) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29187782) q[1];
sx q[1];
rz(-1.4018702) q[1];
sx q[1];
rz(-0.090976322) q[1];
rz(-pi) q[2];
rz(-2.4683687) q[3];
sx q[3];
rz(-2.554318) q[3];
sx q[3];
rz(-0.67929635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2174012) q[2];
sx q[2];
rz(-1.8175586) q[2];
sx q[2];
rz(0.021281555) q[2];
rz(1.6993258) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4845881) q[0];
sx q[0];
rz(-2.5452884) q[0];
sx q[0];
rz(-0.044629991) q[0];
rz(-2.1394829) q[1];
sx q[1];
rz(-2.1410746) q[1];
sx q[1];
rz(-3.1071641) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7444367) q[0];
sx q[0];
rz(-1.7560871) q[0];
sx q[0];
rz(3.1277083) q[0];
rz(-pi) q[1];
rz(1.0677412) q[2];
sx q[2];
rz(-1.6776553) q[2];
sx q[2];
rz(1.6018794) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.62528246) q[1];
sx q[1];
rz(-1.0370266) q[1];
sx q[1];
rz(-0.37955243) q[1];
rz(-2.4363082) q[3];
sx q[3];
rz(-1.4449638) q[3];
sx q[3];
rz(2.3182486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3520711) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(0.96674031) q[2];
rz(0.012185193) q[3];
sx q[3];
rz(-0.92958486) q[3];
sx q[3];
rz(-0.10908443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063342) q[0];
sx q[0];
rz(-0.2671347) q[0];
sx q[0];
rz(-0.99408856) q[0];
rz(3.0864691) q[1];
sx q[1];
rz(-1.2722641) q[1];
sx q[1];
rz(0.73928839) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4261739) q[0];
sx q[0];
rz(-0.16695515) q[0];
sx q[0];
rz(-2.1782317) q[0];
rz(-pi) q[1];
rz(2.5653966) q[2];
sx q[2];
rz(-1.9857166) q[2];
sx q[2];
rz(-0.31838271) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6608097) q[1];
sx q[1];
rz(-2.0936831) q[1];
sx q[1];
rz(-1.3516264) q[1];
x q[2];
rz(0.39412734) q[3];
sx q[3];
rz(-0.65791241) q[3];
sx q[3];
rz(0.037362785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58549515) q[2];
sx q[2];
rz(-1.7332417) q[2];
sx q[2];
rz(-2.5778256) q[2];
rz(2.901315) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(1.3747922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69650841) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(-0.58404303) q[0];
rz(2.7208327) q[1];
sx q[1];
rz(-1.0412443) q[1];
sx q[1];
rz(0.27841321) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1906492) q[0];
sx q[0];
rz(-1.4454495) q[0];
sx q[0];
rz(-0.9992674) q[0];
x q[1];
rz(-0.39324795) q[2];
sx q[2];
rz(-1.3967782) q[2];
sx q[2];
rz(-0.72088036) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.486057) q[1];
sx q[1];
rz(-1.2206519) q[1];
sx q[1];
rz(2.9085367) q[1];
rz(0.24448963) q[3];
sx q[3];
rz(-1.036662) q[3];
sx q[3];
rz(0.29230803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2723096) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(-0.38796866) q[2];
rz(-1.3502454) q[3];
sx q[3];
rz(-2.6781121) q[3];
sx q[3];
rz(2.8665682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(-0.7318837) q[0];
rz(-2.9751119) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(0.073721185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4972992) q[0];
sx q[0];
rz(-1.6719581) q[0];
sx q[0];
rz(2.0072719) q[0];
rz(-pi) q[1];
rz(1.45103) q[2];
sx q[2];
rz(-1.4904067) q[2];
sx q[2];
rz(0.62270852) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8718865) q[1];
sx q[1];
rz(-1.6082113) q[1];
sx q[1];
rz(-0.58362959) q[1];
rz(-1.681626) q[3];
sx q[3];
rz(-2.14114) q[3];
sx q[3];
rz(-0.092872083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59447294) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(-0.42868844) q[2];
rz(-1.8244913) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(-1.2111506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73228943) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(-2.0987341) q[0];
rz(1.5111142) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(1.3483378) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52085984) q[0];
sx q[0];
rz(-3.0254786) q[0];
sx q[0];
rz(1.7945047) q[0];
rz(0.83421631) q[2];
sx q[2];
rz(-1.6305106) q[2];
sx q[2];
rz(0.63876736) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9546982) q[1];
sx q[1];
rz(-0.3317301) q[1];
sx q[1];
rz(-1.7897254) q[1];
x q[2];
rz(-0.24095778) q[3];
sx q[3];
rz(-1.3853419) q[3];
sx q[3];
rz(-1.0518187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.80031359) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(-0.094816118) q[2];
rz(1.9684277) q[3];
sx q[3];
rz(-1.4011819) q[3];
sx q[3];
rz(-2.9828984) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.512758) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(1.5402773) q[1];
sx q[1];
rz(-1.5581144) q[1];
sx q[1];
rz(1.5099572) q[1];
rz(-1.7839292) q[2];
sx q[2];
rz(-1.6518946) q[2];
sx q[2];
rz(1.2988731) q[2];
rz(-2.1223162) q[3];
sx q[3];
rz(-0.84001361) q[3];
sx q[3];
rz(2.7915814) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
