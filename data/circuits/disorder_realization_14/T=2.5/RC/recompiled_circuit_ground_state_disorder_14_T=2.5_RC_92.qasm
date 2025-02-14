OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-3.0225585) q[0];
sx q[0];
rz(-1.1086388) q[0];
sx q[0];
rz(-2.570987) q[0];
rz(-1.5174436) q[1];
sx q[1];
rz(-2.1921506) q[1];
sx q[1];
rz(1.5631262) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1180029) q[0];
sx q[0];
rz(-0.45195107) q[0];
sx q[0];
rz(1.8382545) q[0];
rz(-pi) q[1];
rz(2.7907098) q[2];
sx q[2];
rz(-0.84215468) q[2];
sx q[2];
rz(-0.0082727783) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38754639) q[1];
sx q[1];
rz(-1.7472194) q[1];
sx q[1];
rz(0.37732084) q[1];
x q[2];
rz(-0.61042036) q[3];
sx q[3];
rz(-1.2812611) q[3];
sx q[3];
rz(0.60411835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7382875) q[2];
sx q[2];
rz(-0.83633509) q[2];
sx q[2];
rz(2.9264012) q[2];
rz(-3.134794) q[3];
sx q[3];
rz(-0.88383043) q[3];
sx q[3];
rz(-2.1782037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(2.6007518) q[0];
sx q[0];
rz(-0.90966666) q[0];
sx q[0];
rz(1.2700861) q[0];
rz(-2.0974244) q[1];
sx q[1];
rz(-2.7296598) q[1];
sx q[1];
rz(-2.5792436) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5836307) q[0];
sx q[0];
rz(-1.4578125) q[0];
sx q[0];
rz(-1.8336748) q[0];
rz(-pi) q[1];
x q[1];
rz(0.59074838) q[2];
sx q[2];
rz(-0.68522725) q[2];
sx q[2];
rz(2.6815989) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.54004242) q[1];
sx q[1];
rz(-1.6311495) q[1];
sx q[1];
rz(0.73436952) q[1];
rz(-pi) q[2];
rz(1.2581732) q[3];
sx q[3];
rz(-1.8653096) q[3];
sx q[3];
rz(-2.9972863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.0079785138) q[2];
sx q[2];
rz(-2.9972711) q[2];
sx q[2];
rz(-3.0001384) q[2];
rz(-2.134363) q[3];
sx q[3];
rz(-1.4814582) q[3];
sx q[3];
rz(3.0032515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36002845) q[0];
sx q[0];
rz(-1.6369632) q[0];
sx q[0];
rz(-0.052852782) q[0];
rz(1.3797034) q[1];
sx q[1];
rz(-2.3638937) q[1];
sx q[1];
rz(-2.8676829) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3223548) q[0];
sx q[0];
rz(-0.18095005) q[0];
sx q[0];
rz(1.7767679) q[0];
x q[1];
rz(-0.093294662) q[2];
sx q[2];
rz(-2.3226566) q[2];
sx q[2];
rz(-1.3064177) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9603091) q[1];
sx q[1];
rz(-0.77086721) q[1];
sx q[1];
rz(3.1032352) q[1];
x q[2];
rz(-1.307358) q[3];
sx q[3];
rz(-1.6826077) q[3];
sx q[3];
rz(-1.5105351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2770047) q[2];
sx q[2];
rz(-1.4952345) q[2];
sx q[2];
rz(2.541339) q[2];
rz(0.65435919) q[3];
sx q[3];
rz(-2.8221966) q[3];
sx q[3];
rz(-1.4734242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3510503) q[0];
sx q[0];
rz(-1.496614) q[0];
sx q[0];
rz(-0.45781621) q[0];
rz(-0.59902016) q[1];
sx q[1];
rz(-0.51282239) q[1];
sx q[1];
rz(-2.8694428) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2519855) q[0];
sx q[0];
rz(-2.4757644) q[0];
sx q[0];
rz(-0.4154887) q[0];
rz(0.53109269) q[2];
sx q[2];
rz(-1.102524) q[2];
sx q[2];
rz(2.7837769) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5065342) q[1];
sx q[1];
rz(-2.9218704) q[1];
sx q[1];
rz(-0.57172872) q[1];
rz(2.5013345) q[3];
sx q[3];
rz(-1.4319515) q[3];
sx q[3];
rz(-1.1306896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.51364556) q[2];
sx q[2];
rz(-2.7529035) q[2];
sx q[2];
rz(-1.909931) q[2];
rz(1.4218467) q[3];
sx q[3];
rz(-1.809241) q[3];
sx q[3];
rz(-0.96760145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.642639) q[0];
sx q[0];
rz(-0.068878219) q[0];
sx q[0];
rz(-2.4054085) q[0];
rz(0.63756293) q[1];
sx q[1];
rz(-1.2022737) q[1];
sx q[1];
rz(0.41086248) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7124568) q[0];
sx q[0];
rz(-0.68198689) q[0];
sx q[0];
rz(0.76759591) q[0];
rz(-0.21145282) q[2];
sx q[2];
rz(-2.1295068) q[2];
sx q[2];
rz(-1.3635456) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.84437925) q[1];
sx q[1];
rz(-2.2003268) q[1];
sx q[1];
rz(-1.5624863) q[1];
rz(-0.75848364) q[3];
sx q[3];
rz(-0.82421571) q[3];
sx q[3];
rz(3.0424524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.32834443) q[2];
sx q[2];
rz(-1.5466377) q[2];
sx q[2];
rz(0.52097875) q[2];
rz(3.0327435) q[3];
sx q[3];
rz(-2.1284926) q[3];
sx q[3];
rz(-2.525135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2206889) q[0];
sx q[0];
rz(-1.9158798) q[0];
sx q[0];
rz(-1.8826245) q[0];
rz(-2.068223) q[1];
sx q[1];
rz(-1.4732889) q[1];
sx q[1];
rz(0.64356709) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94144017) q[0];
sx q[0];
rz(-0.45108953) q[0];
sx q[0];
rz(-0.01157014) q[0];
x q[1];
rz(-2.6344755) q[2];
sx q[2];
rz(-2.7165453) q[2];
sx q[2];
rz(-3.0866462) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.7832598) q[1];
sx q[1];
rz(-1.4070916) q[1];
sx q[1];
rz(0.81061426) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57585133) q[3];
sx q[3];
rz(-1.1237877) q[3];
sx q[3];
rz(-0.67547638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.7966938) q[2];
sx q[2];
rz(-1.0641791) q[2];
sx q[2];
rz(1.3382443) q[2];
rz(-1.4517339) q[3];
sx q[3];
rz(-1.9386049) q[3];
sx q[3];
rz(0.10045997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.34375) q[0];
sx q[0];
rz(-1.2748953) q[0];
sx q[0];
rz(-2.1882958) q[0];
rz(-0.15996179) q[1];
sx q[1];
rz(-2.5123031) q[1];
sx q[1];
rz(0.86872855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9729444) q[0];
sx q[0];
rz(-2.4356844) q[0];
sx q[0];
rz(-1.9078518) q[0];
rz(-pi) q[1];
rz(0.73448278) q[2];
sx q[2];
rz(-2.1772263) q[2];
sx q[2];
rz(-0.36494985) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27426961) q[1];
sx q[1];
rz(-1.5918917) q[1];
sx q[1];
rz(-0.29550754) q[1];
rz(0.89959156) q[3];
sx q[3];
rz(-1.7760522) q[3];
sx q[3];
rz(-1.156607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4538883) q[2];
sx q[2];
rz(-2.0084461) q[2];
sx q[2];
rz(-2.647661) q[2];
rz(-1.4880873) q[3];
sx q[3];
rz(-2.2736214) q[3];
sx q[3];
rz(-0.021763703) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9138551) q[0];
sx q[0];
rz(-1.4233002) q[0];
sx q[0];
rz(-0.29105759) q[0];
rz(-0.13045467) q[1];
sx q[1];
rz(-2.7619402) q[1];
sx q[1];
rz(1.5694654) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0091189) q[0];
sx q[0];
rz(-1.557701) q[0];
sx q[0];
rz(1.520169) q[0];
rz(-2.9697598) q[2];
sx q[2];
rz(-0.90441695) q[2];
sx q[2];
rz(3.1058571) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3321905) q[1];
sx q[1];
rz(-0.6242663) q[1];
sx q[1];
rz(0.19849507) q[1];
rz(-pi) q[2];
rz(1.6984822) q[3];
sx q[3];
rz(-1.8771168) q[3];
sx q[3];
rz(-1.1088913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.36720413) q[2];
sx q[2];
rz(-0.17750326) q[2];
sx q[2];
rz(0.11167488) q[2];
rz(-0.81237927) q[3];
sx q[3];
rz(-1.2923391) q[3];
sx q[3];
rz(-1.5264548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-0.77686253) q[0];
sx q[0];
rz(-2.9193945) q[0];
sx q[0];
rz(-0.36866933) q[0];
rz(-0.61019507) q[1];
sx q[1];
rz(-1.5851603) q[1];
sx q[1];
rz(2.1048224) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42175779) q[0];
sx q[0];
rz(-1.0204691) q[0];
sx q[0];
rz(-1.3070413) q[0];
rz(-1.5794994) q[2];
sx q[2];
rz(-2.6285951) q[2];
sx q[2];
rz(2.7287116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1644124) q[1];
sx q[1];
rz(-2.3930139) q[1];
sx q[1];
rz(-1.776265) q[1];
x q[2];
rz(2.0572564) q[3];
sx q[3];
rz(-1.7170625) q[3];
sx q[3];
rz(1.676596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.6250352) q[2];
sx q[2];
rz(-0.80356193) q[2];
sx q[2];
rz(0.45300031) q[2];
rz(-2.1895444) q[3];
sx q[3];
rz(-2.0249849) q[3];
sx q[3];
rz(3.0208352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20671885) q[0];
sx q[0];
rz(-3.0937338) q[0];
sx q[0];
rz(-0.55024838) q[0];
rz(-2.0568559) q[1];
sx q[1];
rz(-2.5526498) q[1];
sx q[1];
rz(-1.3909856) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15875951) q[0];
sx q[0];
rz(-1.5313494) q[0];
sx q[0];
rz(2.0723913) q[0];
x q[1];
rz(-2.1065305) q[2];
sx q[2];
rz(-1.4260056) q[2];
sx q[2];
rz(-2.5404251) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8164313) q[1];
sx q[1];
rz(-1.0712938) q[1];
sx q[1];
rz(-2.4234467) q[1];
rz(-pi) q[2];
rz(2.3971294) q[3];
sx q[3];
rz(-3.0583707) q[3];
sx q[3];
rz(1.2017488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37798205) q[2];
sx q[2];
rz(-0.72145975) q[2];
sx q[2];
rz(1.9788431) q[2];
rz(-2.5871596) q[3];
sx q[3];
rz(-2.2576163) q[3];
sx q[3];
rz(-2.673431) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8155895) q[0];
sx q[0];
rz(-1.3023888) q[0];
sx q[0];
rz(1.6525255) q[0];
rz(-0.21492699) q[1];
sx q[1];
rz(-0.88354127) q[1];
sx q[1];
rz(3.0149928) q[1];
rz(0.28543546) q[2];
sx q[2];
rz(-0.68842059) q[2];
sx q[2];
rz(-0.60981087) q[2];
rz(-1.636408) q[3];
sx q[3];
rz(-0.8498522) q[3];
sx q[3];
rz(-0.1838799) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
