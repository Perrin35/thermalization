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
rz(1.963653) q[0];
sx q[0];
rz(-0.093129245) q[0];
sx q[0];
rz(7.0153305) q[0];
rz(-1.5692476) q[1];
sx q[1];
rz(-0.88580004) q[1];
sx q[1];
rz(-2.7028309) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2193091) q[0];
sx q[0];
rz(-1.5153335) q[0];
sx q[0];
rz(-1.903141) q[0];
rz(-pi) q[1];
x q[1];
rz(2.311539) q[2];
sx q[2];
rz(-2.7657653) q[2];
sx q[2];
rz(-0.91700441) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6631515) q[1];
sx q[1];
rz(-1.8290797) q[1];
sx q[1];
rz(1.5228935) q[1];
rz(-pi) q[2];
rz(3.0436727) q[3];
sx q[3];
rz(-0.46028462) q[3];
sx q[3];
rz(-1.8079881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0438805) q[2];
sx q[2];
rz(-0.040412929) q[2];
sx q[2];
rz(2.2229693) q[2];
rz(-1.0527323) q[3];
sx q[3];
rz(-2.2248) q[3];
sx q[3];
rz(0.91744939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87616462) q[0];
sx q[0];
rz(-2.3227203) q[0];
sx q[0];
rz(-2.2692666) q[0];
rz(-2.1288952) q[1];
sx q[1];
rz(-1.6302949) q[1];
sx q[1];
rz(-0.33908078) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6144086) q[0];
sx q[0];
rz(-1.1361343) q[0];
sx q[0];
rz(-1.4285157) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0929957) q[2];
sx q[2];
rz(-1.1638008) q[2];
sx q[2];
rz(-2.5841449) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2660632) q[1];
sx q[1];
rz(-1.5782981) q[1];
sx q[1];
rz(-0.041635978) q[1];
rz(1.2148772) q[3];
sx q[3];
rz(-2.7971075) q[3];
sx q[3];
rz(-0.075862715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.33112153) q[2];
sx q[2];
rz(-1.9942185) q[2];
sx q[2];
rz(-2.5362711) q[2];
rz(-0.33453861) q[3];
sx q[3];
rz(-0.3722705) q[3];
sx q[3];
rz(-0.076586671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89638585) q[0];
sx q[0];
rz(-0.67436445) q[0];
sx q[0];
rz(1.6687923) q[0];
rz(2.8145166) q[1];
sx q[1];
rz(-1.3027124) q[1];
sx q[1];
rz(2.8715141) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.208858) q[0];
sx q[0];
rz(-1.7678102) q[0];
sx q[0];
rz(-0.52893649) q[0];
rz(-pi) q[1];
x q[1];
rz(0.075614838) q[2];
sx q[2];
rz(-1.9984198) q[2];
sx q[2];
rz(-2.2958034) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1535608) q[1];
sx q[1];
rz(-2.2940293) q[1];
sx q[1];
rz(-2.5392697) q[1];
rz(-pi) q[2];
rz(1.5669426) q[3];
sx q[3];
rz(-1.650963) q[3];
sx q[3];
rz(-1.7146669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.104287) q[2];
sx q[2];
rz(-1.9992) q[2];
sx q[2];
rz(-1.2874862) q[2];
rz(1.2152952) q[3];
sx q[3];
rz(-1.3853962) q[3];
sx q[3];
rz(2.9638885) q[3];
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
rz(0.52008587) q[0];
sx q[0];
rz(-0.44317133) q[0];
sx q[0];
rz(-0.31975123) q[0];
rz(-2.7301835) q[1];
sx q[1];
rz(-1.5450059) q[1];
sx q[1];
rz(-3.06126) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92341766) q[0];
sx q[0];
rz(-1.8785155) q[0];
sx q[0];
rz(0.20826343) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80647237) q[2];
sx q[2];
rz(-1.6739168) q[2];
sx q[2];
rz(2.0224624) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2633966) q[1];
sx q[1];
rz(-0.85021633) q[1];
sx q[1];
rz(2.7769068) q[1];
rz(-pi) q[2];
rz(2.5307054) q[3];
sx q[3];
rz(-1.7912994) q[3];
sx q[3];
rz(1.0569416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8722998) q[2];
sx q[2];
rz(-0.23739693) q[2];
sx q[2];
rz(0.044895127) q[2];
rz(-0.53523713) q[3];
sx q[3];
rz(-1.5072631) q[3];
sx q[3];
rz(3.0285192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.991268) q[0];
sx q[0];
rz(-0.70882216) q[0];
sx q[0];
rz(2.2659361) q[0];
rz(-2.9278897) q[1];
sx q[1];
rz(-2.6468266) q[1];
sx q[1];
rz(1.5692086) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.039883651) q[0];
sx q[0];
rz(-1.3613762) q[0];
sx q[0];
rz(-1.8906192) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46342586) q[2];
sx q[2];
rz(-0.48115402) q[2];
sx q[2];
rz(1.1347186) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8991291) q[1];
sx q[1];
rz(-1.4664259) q[1];
sx q[1];
rz(-2.4422798) q[1];
rz(-1.476579) q[3];
sx q[3];
rz(-1.251128) q[3];
sx q[3];
rz(1.6188542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.956942) q[2];
sx q[2];
rz(-2.9968379) q[2];
sx q[2];
rz(-1.0682028) q[2];
rz(-0.60643658) q[3];
sx q[3];
rz(-1.2111827) q[3];
sx q[3];
rz(0.010312168) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40474263) q[0];
sx q[0];
rz(-0.24979845) q[0];
sx q[0];
rz(-1.3909719) q[0];
rz(-0.5237611) q[1];
sx q[1];
rz(-0.57558376) q[1];
sx q[1];
rz(1.9578804) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9212657) q[0];
sx q[0];
rz(-1.8588716) q[0];
sx q[0];
rz(-0.057805268) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.95518388) q[2];
sx q[2];
rz(-3.0490626) q[2];
sx q[2];
rz(1.1902155) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8642446) q[1];
sx q[1];
rz(-1.7226761) q[1];
sx q[1];
rz(2.1870696) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88504412) q[3];
sx q[3];
rz(-2.9687299) q[3];
sx q[3];
rz(-2.1751753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4411053) q[2];
sx q[2];
rz(-0.59868559) q[2];
sx q[2];
rz(2.109745) q[2];
rz(1.4619689) q[3];
sx q[3];
rz(-2.4470191) q[3];
sx q[3];
rz(1.7803378) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4782891) q[0];
sx q[0];
rz(-0.42600584) q[0];
sx q[0];
rz(1.7167094) q[0];
rz(2.5583963) q[1];
sx q[1];
rz(-1.9164663) q[1];
sx q[1];
rz(-3.084175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9641477) q[0];
sx q[0];
rz(-2.3741407) q[0];
sx q[0];
rz(-1.5758118) q[0];
rz(-1.0591828) q[2];
sx q[2];
rz(-1.0213189) q[2];
sx q[2];
rz(-1.8715931) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8805931) q[1];
sx q[1];
rz(-2.2003502) q[1];
sx q[1];
rz(-1.7884607) q[1];
x q[2];
rz(0.52773169) q[3];
sx q[3];
rz(-0.70042983) q[3];
sx q[3];
rz(2.79984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8251557) q[2];
sx q[2];
rz(-1.2978483) q[2];
sx q[2];
rz(-1.8334897) q[2];
rz(1.3206652) q[3];
sx q[3];
rz(-2.2949009) q[3];
sx q[3];
rz(-2.2804885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2945781) q[0];
sx q[0];
rz(-2.4569643) q[0];
sx q[0];
rz(-1.8735877) q[0];
rz(2.0199203) q[1];
sx q[1];
rz(-1.8662165) q[1];
sx q[1];
rz(2.4915288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2156159) q[0];
sx q[0];
rz(-2.4575666) q[0];
sx q[0];
rz(-2.2355272) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6261601) q[2];
sx q[2];
rz(-2.353172) q[2];
sx q[2];
rz(1.5962708) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.5019116) q[1];
sx q[1];
rz(-0.6957275) q[1];
sx q[1];
rz(-0.63668107) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1459747) q[3];
sx q[3];
rz(-1.1523968) q[3];
sx q[3];
rz(-2.0783238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0012297) q[2];
sx q[2];
rz(-1.4823806) q[2];
sx q[2];
rz(-0.39230997) q[2];
rz(-2.9366734) q[3];
sx q[3];
rz(-0.50060087) q[3];
sx q[3];
rz(2.4216381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0023163) q[0];
sx q[0];
rz(-1.5895695) q[0];
sx q[0];
rz(-1.4578777) q[0];
rz(-0.32683364) q[1];
sx q[1];
rz(-1.146233) q[1];
sx q[1];
rz(-2.0976417) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9245968) q[0];
sx q[0];
rz(-1.9925653) q[0];
sx q[0];
rz(2.934458) q[0];
rz(-1.7108261) q[2];
sx q[2];
rz(-1.1159889) q[2];
sx q[2];
rz(0.32250139) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3957246) q[1];
sx q[1];
rz(-1.1765132) q[1];
sx q[1];
rz(3.0529725) q[1];
rz(-pi) q[2];
rz(-0.92071988) q[3];
sx q[3];
rz(-1.8333922) q[3];
sx q[3];
rz(0.85185862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8331208) q[2];
sx q[2];
rz(-1.0288419) q[2];
sx q[2];
rz(-0.92588818) q[2];
rz(-1.067767) q[3];
sx q[3];
rz(-1.1238778) q[3];
sx q[3];
rz(-1.3656176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2866216) q[0];
sx q[0];
rz(-1.4737031) q[0];
sx q[0];
rz(-0.58706748) q[0];
rz(-2.558737) q[1];
sx q[1];
rz(-0.28935495) q[1];
sx q[1];
rz(2.6606182) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60763393) q[0];
sx q[0];
rz(-1.9780428) q[0];
sx q[0];
rz(-1.3429705) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9513861) q[2];
sx q[2];
rz(-1.9252732) q[2];
sx q[2];
rz(0.5917393) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.15156432) q[1];
sx q[1];
rz(-2.1089601) q[1];
sx q[1];
rz(-2.0694222) q[1];
x q[2];
rz(2.6943227) q[3];
sx q[3];
rz(-0.4382689) q[3];
sx q[3];
rz(0.13188383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.19758548) q[2];
sx q[2];
rz(-2.6052167) q[2];
sx q[2];
rz(-1.8120922) q[2];
rz(2.3594989) q[3];
sx q[3];
rz(-2.8160281) q[3];
sx q[3];
rz(1.9935002) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44163497) q[0];
sx q[0];
rz(-1.4186207) q[0];
sx q[0];
rz(1.6765208) q[0];
rz(-1.5326473) q[1];
sx q[1];
rz(-1.1679222) q[1];
sx q[1];
rz(2.4293778) q[1];
rz(2.4075422) q[2];
sx q[2];
rz(-1.0577591) q[2];
sx q[2];
rz(-1.9815097) q[2];
rz(1.731012) q[3];
sx q[3];
rz(-1.0777149) q[3];
sx q[3];
rz(3.1330681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
