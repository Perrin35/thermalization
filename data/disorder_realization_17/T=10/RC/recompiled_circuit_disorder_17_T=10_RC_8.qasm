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
rz(1.1458122) q[0];
rz(-2.1687578) q[1];
sx q[1];
rz(-1.6701148) q[1];
sx q[1];
rz(2.8491128) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2746409) q[0];
sx q[0];
rz(-0.39817087) q[0];
sx q[0];
rz(1.7142332) q[0];
rz(-pi) q[1];
rz(3.0897065) q[2];
sx q[2];
rz(-2.3615026) q[2];
sx q[2];
rz(-0.29733301) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9521192) q[1];
sx q[1];
rz(-1.6166302) q[1];
sx q[1];
rz(0.97656753) q[1];
rz(-1.0969767) q[3];
sx q[3];
rz(-1.5690656) q[3];
sx q[3];
rz(2.2497821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.73074377) q[2];
sx q[2];
rz(-2.0480672) q[2];
sx q[2];
rz(2.6921819) q[2];
rz(-2.4959026) q[3];
sx q[3];
rz(-2.460545) q[3];
sx q[3];
rz(1.2908363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9280424) q[0];
sx q[0];
rz(-0.88590652) q[0];
sx q[0];
rz(0.74203062) q[0];
rz(-1.4713326) q[1];
sx q[1];
rz(-0.69683087) q[1];
sx q[1];
rz(1.3630294) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.610696) q[0];
sx q[0];
rz(-0.96141059) q[0];
sx q[0];
rz(-2.0531528) q[0];
x q[1];
rz(-2.734191) q[2];
sx q[2];
rz(-0.50440895) q[2];
sx q[2];
rz(0.87770578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0568309) q[1];
sx q[1];
rz(-2.8014604) q[1];
sx q[1];
rz(-3.0840193) q[1];
rz(-pi) q[2];
rz(-2.1915073) q[3];
sx q[3];
rz(-0.32334298) q[3];
sx q[3];
rz(-2.8107373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8403975) q[2];
sx q[2];
rz(-0.49673721) q[2];
sx q[2];
rz(1.2724686) q[2];
rz(2.8043591) q[3];
sx q[3];
rz(-1.971258) q[3];
sx q[3];
rz(0.01005323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0740046) q[0];
sx q[0];
rz(-2.8019866) q[0];
sx q[0];
rz(-0.2628251) q[0];
rz(0.35573959) q[1];
sx q[1];
rz(-1.4590615) q[1];
sx q[1];
rz(-1.9062818) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6567898) q[0];
sx q[0];
rz(-0.018253837) q[0];
sx q[0];
rz(-2.546893) q[0];
x q[1];
rz(-1.6791061) q[2];
sx q[2];
rz(-2.8998067) q[2];
sx q[2];
rz(2.0568648) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.77093) q[1];
sx q[1];
rz(-1.2491033) q[1];
sx q[1];
rz(-1.6817723) q[1];
rz(2.7945307) q[3];
sx q[3];
rz(-2.2357781) q[3];
sx q[3];
rz(-0.60225981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0170903) q[2];
sx q[2];
rz(-1.3829145) q[2];
sx q[2];
rz(2.9612605) q[2];
rz(-0.63594937) q[3];
sx q[3];
rz(-1.162581) q[3];
sx q[3];
rz(2.7699871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3746049) q[0];
sx q[0];
rz(-0.42414442) q[0];
sx q[0];
rz(-0.71907991) q[0];
rz(-2.6285697) q[1];
sx q[1];
rz(-1.9388371) q[1];
sx q[1];
rz(-1.0346574) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2186787) q[0];
sx q[0];
rz(-1.5905587) q[0];
sx q[0];
rz(0.058905525) q[0];
x q[1];
rz(-1.9070508) q[2];
sx q[2];
rz(-0.14306919) q[2];
sx q[2];
rz(2.2640995) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0532916) q[1];
sx q[1];
rz(-2.3846855) q[1];
sx q[1];
rz(-2.7234368) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33200522) q[3];
sx q[3];
rz(-0.49626353) q[3];
sx q[3];
rz(0.23976025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2085312) q[2];
sx q[2];
rz(-2.8716817) q[2];
sx q[2];
rz(-0.2362403) q[2];
rz(1.9832206) q[3];
sx q[3];
rz(-2.1393496) q[3];
sx q[3];
rz(0.85197824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1119969) q[0];
sx q[0];
rz(-2.1313666) q[0];
sx q[0];
rz(-0.90233666) q[0];
rz(2.2205655) q[1];
sx q[1];
rz(-2.5245456) q[1];
sx q[1];
rz(-0.62228084) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.92684) q[0];
sx q[0];
rz(-1.2878969) q[0];
sx q[0];
rz(-2.9515285) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0201966) q[2];
sx q[2];
rz(-0.65181323) q[2];
sx q[2];
rz(2.7898942) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2942549) q[1];
sx q[1];
rz(-1.6604742) q[1];
sx q[1];
rz(1.7404106) q[1];
rz(-2.4683687) q[3];
sx q[3];
rz(-2.554318) q[3];
sx q[3];
rz(2.4622963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.92419147) q[2];
sx q[2];
rz(-1.3240341) q[2];
sx q[2];
rz(3.1203111) q[2];
rz(-1.4422669) q[3];
sx q[3];
rz(-2.4271624) q[3];
sx q[3];
rz(-0.91059476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4845881) q[0];
sx q[0];
rz(-0.5963043) q[0];
sx q[0];
rz(3.0969627) q[0];
rz(-2.1394829) q[1];
sx q[1];
rz(-1.0005181) q[1];
sx q[1];
rz(-0.034428509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1761985) q[0];
sx q[0];
rz(-1.584443) q[0];
sx q[0];
rz(-1.7561046) q[0];
rz(1.351864) q[2];
sx q[2];
rz(-0.5133252) q[2];
sx q[2];
rz(-0.22253144) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-3.1019563) q[1];
sx q[1];
rz(-0.64412457) q[1];
sx q[1];
rz(2.1307751) q[1];
x q[2];
rz(-0.70528443) q[3];
sx q[3];
rz(-1.4449638) q[3];
sx q[3];
rz(-2.3182486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3520711) q[2];
sx q[2];
rz(-1.2723203) q[2];
sx q[2];
rz(2.1748523) q[2];
rz(-0.012185193) q[3];
sx q[3];
rz(-2.2120078) q[3];
sx q[3];
rz(-0.10908443) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0063342) q[0];
sx q[0];
rz(-2.874458) q[0];
sx q[0];
rz(0.99408856) q[0];
rz(3.0864691) q[1];
sx q[1];
rz(-1.8693285) q[1];
sx q[1];
rz(2.4023043) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2544884) q[0];
sx q[0];
rz(-1.4758037) q[0];
sx q[0];
rz(-1.7083005) q[0];
rz(-pi) q[1];
rz(0.57619608) q[2];
sx q[2];
rz(-1.155876) q[2];
sx q[2];
rz(-0.31838271) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1623605) q[1];
sx q[1];
rz(-1.3812961) q[1];
sx q[1];
rz(2.6081671) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39412734) q[3];
sx q[3];
rz(-2.4836802) q[3];
sx q[3];
rz(-3.1042299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.58549515) q[2];
sx q[2];
rz(-1.408351) q[2];
sx q[2];
rz(2.5778256) q[2];
rz(-0.24027763) q[3];
sx q[3];
rz(-2.6085745) q[3];
sx q[3];
rz(-1.7668004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69650841) q[0];
sx q[0];
rz(-2.962528) q[0];
sx q[0];
rz(0.58404303) q[0];
rz(-0.42075992) q[1];
sx q[1];
rz(-2.1003484) q[1];
sx q[1];
rz(2.8631794) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7137372) q[0];
sx q[0];
rz(-0.58361485) q[0];
sx q[0];
rz(1.7996656) q[0];
x q[1];
rz(0.39324795) q[2];
sx q[2];
rz(-1.3967782) q[2];
sx q[2];
rz(0.72088036) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.003493017) q[1];
sx q[1];
rz(-1.7894735) q[1];
sx q[1];
rz(-1.9298645) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1182546) q[3];
sx q[3];
rz(-1.3609144) q[3];
sx q[3];
rz(1.4048214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.86928308) q[2];
sx q[2];
rz(-0.81988207) q[2];
sx q[2];
rz(0.38796866) q[2];
rz(-1.7913473) q[3];
sx q[3];
rz(-0.46348059) q[3];
sx q[3];
rz(2.8665682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2824654) q[0];
sx q[0];
rz(-0.21331856) q[0];
sx q[0];
rz(-2.4097089) q[0];
rz(2.9751119) q[1];
sx q[1];
rz(-0.44547588) q[1];
sx q[1];
rz(-0.073721185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.026422231) q[0];
sx q[0];
rz(-1.1367043) q[0];
sx q[0];
rz(-3.0300481) q[0];
rz(-pi) q[1];
rz(0.080967112) q[2];
sx q[2];
rz(-1.4514187) q[2];
sx q[2];
rz(-0.93842426) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.32578711) q[1];
sx q[1];
rz(-0.98762883) q[1];
sx q[1];
rz(-1.5259685) q[1];
rz(1.681626) q[3];
sx q[3];
rz(-2.14114) q[3];
sx q[3];
rz(-3.0487206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5471197) q[2];
sx q[2];
rz(-0.23590817) q[2];
sx q[2];
rz(-2.7129042) q[2];
rz(-1.3171014) q[3];
sx q[3];
rz(-1.7842112) q[3];
sx q[3];
rz(-1.930442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73228943) q[0];
sx q[0];
rz(-1.4855054) q[0];
sx q[0];
rz(2.0987341) q[0];
rz(-1.6304784) q[1];
sx q[1];
rz(-1.379456) q[1];
sx q[1];
rz(-1.7932549) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74603426) q[0];
sx q[0];
rz(-1.4575882) q[0];
sx q[0];
rz(3.1157225) q[0];
rz(-pi) q[1];
rz(-0.83421631) q[2];
sx q[2];
rz(-1.5110821) q[2];
sx q[2];
rz(0.63876736) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9546982) q[1];
sx q[1];
rz(-2.8098626) q[1];
sx q[1];
rz(1.7897254) q[1];
rz(1.3799558) q[3];
sx q[3];
rz(-1.3340501) q[3];
sx q[3];
rz(-0.56425795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.80031359) q[2];
sx q[2];
rz(-0.51395243) q[2];
sx q[2];
rz(-0.094816118) q[2];
rz(1.173165) q[3];
sx q[3];
rz(-1.7404107) q[3];
sx q[3];
rz(-2.9828984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6288347) q[0];
sx q[0];
rz(-0.47369581) q[0];
sx q[0];
rz(1.0439903) q[0];
rz(1.6013153) q[1];
sx q[1];
rz(-1.5834783) q[1];
sx q[1];
rz(-1.6316354) q[1];
rz(1.9376471) q[2];
sx q[2];
rz(-0.22782142) q[2];
sx q[2];
rz(0.08624764) q[2];
rz(0.52900984) q[3];
sx q[3];
rz(-0.88376868) q[3];
sx q[3];
rz(0.39467011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
