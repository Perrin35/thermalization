OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(-1.4533071) q[0];
sx q[0];
rz(0.31153554) q[0];
rz(-0.43752924) q[1];
sx q[1];
rz(-1.8234) q[1];
sx q[1];
rz(0.55895609) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6214949) q[0];
sx q[0];
rz(-1.1557475) q[0];
sx q[0];
rz(-2.9893304) q[0];
rz(-pi) q[1];
rz(0.020157651) q[2];
sx q[2];
rz(-2.5841789) q[2];
sx q[2];
rz(-1.3486947) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54139582) q[1];
sx q[1];
rz(-2.2943498) q[1];
sx q[1];
rz(-2.3266836) q[1];
rz(-pi) q[2];
rz(2.7028014) q[3];
sx q[3];
rz(-1.3109129) q[3];
sx q[3];
rz(-2.4364542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3922334) q[2];
sx q[2];
rz(-1.2831251) q[2];
sx q[2];
rz(0.63670811) q[2];
rz(0.84896815) q[3];
sx q[3];
rz(-2.520112) q[3];
sx q[3];
rz(2.9076715) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7648776) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(2.9887181) q[0];
rz(-2.3846467) q[1];
sx q[1];
rz(-1.5544954) q[1];
sx q[1];
rz(-2.1551932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5312885) q[0];
sx q[0];
rz(-0.072040759) q[0];
sx q[0];
rz(-2.5487367) q[0];
rz(0.54665357) q[2];
sx q[2];
rz(-0.84178998) q[2];
sx q[2];
rz(2.5559049) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0609329) q[1];
sx q[1];
rz(-0.45127171) q[1];
sx q[1];
rz(-2.5732451) q[1];
x q[2];
rz(1.9513449) q[3];
sx q[3];
rz(-0.35223397) q[3];
sx q[3];
rz(-0.16570839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5622921) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(0.78312773) q[2];
rz(-0.018571818) q[3];
sx q[3];
rz(-1.5037856) q[3];
sx q[3];
rz(-0.40772453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0884393) q[0];
sx q[0];
rz(-0.29215559) q[0];
sx q[0];
rz(2.1799178) q[0];
rz(-2.7812474) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(-3.0128984) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6260687) q[0];
sx q[0];
rz(-1.4883853) q[0];
sx q[0];
rz(0.047973085) q[0];
rz(-0.80758904) q[2];
sx q[2];
rz(-1.1620887) q[2];
sx q[2];
rz(-1.7847716) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.80850959) q[1];
sx q[1];
rz(-1.2915478) q[1];
sx q[1];
rz(-2.6781494) q[1];
x q[2];
rz(0.80273654) q[3];
sx q[3];
rz(-1.667913) q[3];
sx q[3];
rz(0.10955284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.06015691) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(-1.8998247) q[2];
rz(0.5870108) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(-0.97755066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63242763) q[0];
sx q[0];
rz(-0.88212633) q[0];
sx q[0];
rz(-2.0571016) q[0];
rz(-1.4831316) q[1];
sx q[1];
rz(-2.5741534) q[1];
sx q[1];
rz(-0.09253563) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4400892) q[0];
sx q[0];
rz(-1.2883696) q[0];
sx q[0];
rz(0.28042067) q[0];
rz(1.091435) q[2];
sx q[2];
rz(-1.879868) q[2];
sx q[2];
rz(-2.022559) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9738237) q[1];
sx q[1];
rz(-1.9134221) q[1];
sx q[1];
rz(1.6325566) q[1];
rz(-0.16102287) q[3];
sx q[3];
rz(-2.2194214) q[3];
sx q[3];
rz(0.30740689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3080421) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(-0.33205024) q[2];
rz(-1.0559233) q[3];
sx q[3];
rz(-2.8639586) q[3];
sx q[3];
rz(-0.61029303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8191391) q[0];
sx q[0];
rz(-1.2912913) q[0];
sx q[0];
rz(2.9300368) q[0];
rz(1.8353204) q[1];
sx q[1];
rz(-1.897656) q[1];
sx q[1];
rz(-0.64770118) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0928597) q[0];
sx q[0];
rz(-1.6824241) q[0];
sx q[0];
rz(-0.1603006) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7251882) q[2];
sx q[2];
rz(-1.2263745) q[2];
sx q[2];
rz(-1.7248578) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.16072893) q[1];
sx q[1];
rz(-2.594922) q[1];
sx q[1];
rz(-1.9561808) q[1];
rz(-pi) q[2];
rz(0.019142814) q[3];
sx q[3];
rz(-1.585841) q[3];
sx q[3];
rz(-2.6263833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.56090474) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(2.4482751) q[2];
rz(-0.66926113) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(-0.0049237331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3174021) q[0];
sx q[0];
rz(-2.9142002) q[0];
sx q[0];
rz(-1.2325226) q[0];
rz(-1.0725853) q[1];
sx q[1];
rz(-1.0718081) q[1];
sx q[1];
rz(-0.17428621) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3586853) q[0];
sx q[0];
rz(-1.8795965) q[0];
sx q[0];
rz(-0.82682825) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2767378) q[2];
sx q[2];
rz(-2.1055429) q[2];
sx q[2];
rz(-0.45644444) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41387687) q[1];
sx q[1];
rz(-1.1912279) q[1];
sx q[1];
rz(0.80645251) q[1];
rz(-1.918119) q[3];
sx q[3];
rz(-2.8109549) q[3];
sx q[3];
rz(-2.2280681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.8217414) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(-2.9439587) q[2];
rz(0.28891426) q[3];
sx q[3];
rz(-0.87696004) q[3];
sx q[3];
rz(1.7355841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.063868) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(0.20859627) q[0];
rz(2.1754307) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(-1.6360412) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87058545) q[0];
sx q[0];
rz(-0.85575543) q[0];
sx q[0];
rz(-2.4164651) q[0];
rz(-pi) q[1];
rz(-3.0622919) q[2];
sx q[2];
rz(-0.43159017) q[2];
sx q[2];
rz(2.0183795) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.58563102) q[1];
sx q[1];
rz(-1.0502083) q[1];
sx q[1];
rz(-2.6160137) q[1];
rz(-pi) q[2];
rz(-2.8553477) q[3];
sx q[3];
rz(-1.6702594) q[3];
sx q[3];
rz(0.092982987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.85764) q[2];
sx q[2];
rz(-2.3985034) q[2];
sx q[2];
rz(3.0440142) q[2];
rz(1.7476667) q[3];
sx q[3];
rz(-1.3503617) q[3];
sx q[3];
rz(-0.71684366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-2.7512648) q[0];
sx q[0];
rz(-1.3289691) q[0];
sx q[0];
rz(-0.51399291) q[0];
rz(3.0184074) q[1];
sx q[1];
rz(-2.8942278) q[1];
sx q[1];
rz(-2.2095912) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1876353) q[0];
sx q[0];
rz(-1.4353416) q[0];
sx q[0];
rz(-2.016469) q[0];
x q[1];
rz(-1.2134238) q[2];
sx q[2];
rz(-2.4744611) q[2];
sx q[2];
rz(-0.92598976) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.92703687) q[1];
sx q[1];
rz(-0.80074691) q[1];
sx q[1];
rz(0.75896778) q[1];
x q[2];
rz(-0.85429116) q[3];
sx q[3];
rz(-1.4201418) q[3];
sx q[3];
rz(2.6137969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1121858) q[2];
sx q[2];
rz(-1.1108578) q[2];
sx q[2];
rz(2.712148) q[2];
rz(-1.9321692) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(0.69303524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76686239) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(1.3388348) q[0];
rz(-0.70612899) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(-2.1910117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50981748) q[0];
sx q[0];
rz(-1.063856) q[0];
sx q[0];
rz(0.99960534) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9865932) q[2];
sx q[2];
rz(-1.6396513) q[2];
sx q[2];
rz(-1.6378251) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1817158) q[1];
sx q[1];
rz(-2.1023589) q[1];
sx q[1];
rz(-2.9243484) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0334942) q[3];
sx q[3];
rz(-2.6712382) q[3];
sx q[3];
rz(0.62894097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.4799698) q[2];
sx q[2];
rz(-1.491549) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(2.2144923) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1442239) q[0];
sx q[0];
rz(-0.08865083) q[0];
sx q[0];
rz(0.22928672) q[0];
rz(-0.43481049) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(0.71892175) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8573498) q[0];
sx q[0];
rz(-1.6650668) q[0];
sx q[0];
rz(2.1988792) q[0];
rz(0.11833338) q[2];
sx q[2];
rz(-2.1552857) q[2];
sx q[2];
rz(0.10533939) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1246008) q[1];
sx q[1];
rz(-1.7182087) q[1];
sx q[1];
rz(0.37313811) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2273916) q[3];
sx q[3];
rz(-2.4061892) q[3];
sx q[3];
rz(-0.38287336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7583313) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-0.74679217) q[2];
rz(-2.2693999) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(0.94223589) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.223021) q[0];
sx q[0];
rz(-1.7445607) q[0];
sx q[0];
rz(1.8013409) q[0];
rz(0.37721286) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(-2.626426) q[2];
sx q[2];
rz(-2.5611521) q[2];
sx q[2];
rz(-0.44708154) q[2];
rz(-1.4217581) q[3];
sx q[3];
rz(-2.2825713) q[3];
sx q[3];
rz(1.5406516) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];