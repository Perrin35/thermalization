OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(-1.7927875) q[1];
sx q[1];
rz(-0.92372149) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1169352) q[0];
sx q[0];
rz(-2.5664461) q[0];
sx q[0];
rz(2.2206109) q[0];
rz(-1.853763) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(1.1793009) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.87286283) q[1];
sx q[1];
rz(-1.1464719) q[1];
sx q[1];
rz(-0.55117589) q[1];
rz(-pi) q[2];
rz(0.49433319) q[3];
sx q[3];
rz(-1.2327854) q[3];
sx q[3];
rz(2.5288343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2279921) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(-1.9447928) q[0];
rz(-0.67990047) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.4555567) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50549492) q[0];
sx q[0];
rz(-2.539145) q[0];
sx q[0];
rz(-1.2732182) q[0];
x q[1];
rz(0.37462072) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(-0.74707109) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1728954) q[1];
sx q[1];
rz(-2.6186133) q[1];
sx q[1];
rz(1.5978659) q[1];
rz(0.060828408) q[3];
sx q[3];
rz(-0.41234499) q[3];
sx q[3];
rz(-2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882554) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(-2.341111) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.9690537) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5867509) q[0];
sx q[0];
rz(-1.2625853) q[0];
sx q[0];
rz(0.59535938) q[0];
rz(-pi) q[1];
rz(2.3861888) q[2];
sx q[2];
rz(-1.3284151) q[2];
sx q[2];
rz(0.6005477) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8923924) q[1];
sx q[1];
rz(-1.3274267) q[1];
sx q[1];
rz(-2.1057486) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3451194) q[3];
sx q[3];
rz(-0.25697069) q[3];
sx q[3];
rz(-2.8185609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(2.2303936) q[2];
rz(-2.1905812) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(-2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(3.0134841) q[0];
rz(-0.076106636) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(-2.6180843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9273705) q[0];
sx q[0];
rz(-2.4282051) q[0];
sx q[0];
rz(-2.5582696) q[0];
x q[1];
rz(2.8115494) q[2];
sx q[2];
rz(-1.651262) q[2];
sx q[2];
rz(2.6835494) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.88672968) q[1];
sx q[1];
rz(-2.2543395) q[1];
sx q[1];
rz(2.8121594) q[1];
rz(-pi) q[2];
rz(-1.773049) q[3];
sx q[3];
rz(-0.60086717) q[3];
sx q[3];
rz(-0.7522538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6161502) q[2];
sx q[2];
rz(-1.5972861) q[2];
sx q[2];
rz(-2.5775487) q[2];
rz(2.8530252) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(-2.2619757) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(-2.1496444) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0569699) q[0];
sx q[0];
rz(-1.3966494) q[0];
sx q[0];
rz(-1.4625545) q[0];
x q[1];
rz(1.5422103) q[2];
sx q[2];
rz(-2.8386142) q[2];
sx q[2];
rz(2.6945393) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.76584133) q[1];
sx q[1];
rz(-1.6110794) q[1];
sx q[1];
rz(0.30893107) q[1];
rz(-0.12985142) q[3];
sx q[3];
rz(-0.81749812) q[3];
sx q[3];
rz(-2.0714456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(2.8707855) q[2];
rz(2.9233542) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4145684) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(1.4250925) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093229175) q[0];
sx q[0];
rz(-1.3472124) q[0];
sx q[0];
rz(0.50888749) q[0];
x q[1];
rz(-1.8259949) q[2];
sx q[2];
rz(-1.4824502) q[2];
sx q[2];
rz(1.3510072) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7879494) q[1];
sx q[1];
rz(-2.0197581) q[1];
sx q[1];
rz(-0.073972703) q[1];
rz(-0.53374966) q[3];
sx q[3];
rz(-2.1042049) q[3];
sx q[3];
rz(2.651754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1340593) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(-2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5086223) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(-2.3197876) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9200631) q[0];
sx q[0];
rz(-2.1863345) q[0];
sx q[0];
rz(-2.2277742) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1173382) q[2];
sx q[2];
rz(-0.80596906) q[2];
sx q[2];
rz(2.1794127) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6268839) q[1];
sx q[1];
rz(-2.2885972) q[1];
sx q[1];
rz(-1.2674598) q[1];
rz(-pi) q[2];
rz(2.6131367) q[3];
sx q[3];
rz(-2.4848357) q[3];
sx q[3];
rz(0.17880759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(-0.24469963) q[2];
rz(-3.0120567) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(-1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(1.3051422) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1730723) q[0];
sx q[0];
rz(-2.1497796) q[0];
sx q[0];
rz(1.1770244) q[0];
rz(-pi) q[1];
rz(-2.4370286) q[2];
sx q[2];
rz(-1.2837871) q[2];
sx q[2];
rz(-1.2003843) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9040363) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(-2.0380286) q[1];
rz(1.6373789) q[3];
sx q[3];
rz(-0.3399907) q[3];
sx q[3];
rz(2.7293527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1398853) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(-1.4902327) q[2];
rz(-2.0643318) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-0.70294356) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84530572) q[0];
sx q[0];
rz(-1.1830813) q[0];
sx q[0];
rz(1.9592497) q[0];
x q[1];
rz(-1.0614971) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(-2.408037) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1785537) q[1];
sx q[1];
rz(-0.56561618) q[1];
sx q[1];
rz(-1.2600684) q[1];
rz(-1.1144936) q[3];
sx q[3];
rz(-0.66702402) q[3];
sx q[3];
rz(1.4240571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(0.30570269) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383485) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(-1.9158069) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(-2.6729565) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0271284) q[0];
sx q[0];
rz(-0.44136029) q[0];
sx q[0];
rz(-2.6370185) q[0];
rz(-pi) q[1];
rz(2.8154545) q[2];
sx q[2];
rz(-1.9378127) q[2];
sx q[2];
rz(2.0740167) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74433078) q[1];
sx q[1];
rz(-1.3052193) q[1];
sx q[1];
rz(-3.0708405) q[1];
rz(0.43379421) q[3];
sx q[3];
rz(-1.1489831) q[3];
sx q[3];
rz(1.0292366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6580711) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(-1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(-0.62190965) q[1];
sx q[1];
rz(-1.6786631) q[1];
sx q[1];
rz(2.8181029) q[1];
rz(2.3399578) q[2];
sx q[2];
rz(-2.2710706) q[2];
sx q[2];
rz(2.0588277) q[2];
rz(1.6735531) q[3];
sx q[3];
rz(-0.65410528) q[3];
sx q[3];
rz(0.42412529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
