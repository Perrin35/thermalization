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
rz(-1.834637) q[0];
sx q[0];
rz(-0.87943465) q[0];
sx q[0];
rz(-0.49361324) q[0];
rz(-1.2797132) q[1];
sx q[1];
rz(3.7682025) q[1];
sx q[1];
rz(11.087853) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077162077) q[0];
sx q[0];
rz(-1.0748942) q[0];
sx q[0];
rz(2.2999963) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8844312) q[2];
sx q[2];
rz(-2.2182121) q[2];
sx q[2];
rz(2.0055298) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2027758) q[1];
sx q[1];
rz(-2.4033053) q[1];
sx q[1];
rz(-2.8144005) q[1];
rz(-pi) q[2];
rz(-2.7731887) q[3];
sx q[3];
rz(-0.96782902) q[3];
sx q[3];
rz(1.6507698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1349858) q[2];
sx q[2];
rz(-2.0015494) q[2];
sx q[2];
rz(-1.1084278) q[2];
rz(-1.3125575) q[3];
sx q[3];
rz(-2.8969942) q[3];
sx q[3];
rz(0.31497064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38158622) q[0];
sx q[0];
rz(-0.8257603) q[0];
sx q[0];
rz(-0.015856892) q[0];
rz(-0.99041692) q[1];
sx q[1];
rz(-2.3742193) q[1];
sx q[1];
rz(-0.86318618) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1970438) q[0];
sx q[0];
rz(-1.4743351) q[0];
sx q[0];
rz(2.405557) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3349971) q[2];
sx q[2];
rz(-1.8827229) q[2];
sx q[2];
rz(0.27659135) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.51746115) q[1];
sx q[1];
rz(-1.7449813) q[1];
sx q[1];
rz(2.6427173) q[1];
rz(-pi) q[2];
rz(-0.14948577) q[3];
sx q[3];
rz(-2.2870697) q[3];
sx q[3];
rz(0.0090741875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4414759) q[2];
sx q[2];
rz(-0.2773383) q[2];
sx q[2];
rz(-2.6727943) q[2];
rz(-0.68764728) q[3];
sx q[3];
rz(-1.629849) q[3];
sx q[3];
rz(2.8414753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1402635) q[0];
sx q[0];
rz(-0.92324531) q[0];
sx q[0];
rz(2.4053251) q[0];
rz(0.34321347) q[1];
sx q[1];
rz(-2.2540269) q[1];
sx q[1];
rz(0.18377486) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90405078) q[0];
sx q[0];
rz(-1.5045368) q[0];
sx q[0];
rz(0.74426331) q[0];
rz(-pi) q[1];
rz(-0.60004514) q[2];
sx q[2];
rz(-0.79216865) q[2];
sx q[2];
rz(0.60952696) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0818834) q[1];
sx q[1];
rz(-2.3569111) q[1];
sx q[1];
rz(-0.42894026) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1322429) q[3];
sx q[3];
rz(-2.5261176) q[3];
sx q[3];
rz(1.6650247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3595714) q[2];
sx q[2];
rz(-0.51681334) q[2];
sx q[2];
rz(-2.6256631) q[2];
rz(-1.1082209) q[3];
sx q[3];
rz(-0.55086946) q[3];
sx q[3];
rz(-1.9903323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16294031) q[0];
sx q[0];
rz(-1.9300224) q[0];
sx q[0];
rz(-0.51280713) q[0];
rz(-0.45979744) q[1];
sx q[1];
rz(-1.5492946) q[1];
sx q[1];
rz(2.3777681) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1490246) q[0];
sx q[0];
rz(-1.4563174) q[0];
sx q[0];
rz(-0.13097288) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41641367) q[2];
sx q[2];
rz(-1.8685307) q[2];
sx q[2];
rz(2.1418051) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1000378) q[1];
sx q[1];
rz(-2.6943992) q[1];
sx q[1];
rz(2.2639422) q[1];
rz(2.7971036) q[3];
sx q[3];
rz(-1.7055344) q[3];
sx q[3];
rz(1.8014993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0189455) q[2];
sx q[2];
rz(-2.985869) q[2];
sx q[2];
rz(0.32749185) q[2];
rz(0.28199768) q[3];
sx q[3];
rz(-1.3982541) q[3];
sx q[3];
rz(-1.5544844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9649488) q[0];
sx q[0];
rz(-0.38589859) q[0];
sx q[0];
rz(-2.1642245) q[0];
rz(-2.7727959) q[1];
sx q[1];
rz(-2.2803523) q[1];
sx q[1];
rz(-1.7288953) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4810886) q[0];
sx q[0];
rz(-1.7605283) q[0];
sx q[0];
rz(-2.0697748) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0875889) q[2];
sx q[2];
rz(-1.5828504) q[2];
sx q[2];
rz(0.64753676) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8846795) q[1];
sx q[1];
rz(-1.2882782) q[1];
sx q[1];
rz(3.1090841) q[1];
x q[2];
rz(0.47980185) q[3];
sx q[3];
rz(-2.5069624) q[3];
sx q[3];
rz(3.0632085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29350027) q[2];
sx q[2];
rz(-2.3257181) q[2];
sx q[2];
rz(0.14981848) q[2];
rz(-2.7430429) q[3];
sx q[3];
rz(-1.6179251) q[3];
sx q[3];
rz(-0.15628763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3715816) q[0];
sx q[0];
rz(-3.0122029) q[0];
sx q[0];
rz(0.55321252) q[0];
rz(-0.54168701) q[1];
sx q[1];
rz(-2.0952416) q[1];
sx q[1];
rz(-0.4479301) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0174367) q[0];
sx q[0];
rz(-1.7154335) q[0];
sx q[0];
rz(0.38993024) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80336415) q[2];
sx q[2];
rz(-0.58310544) q[2];
sx q[2];
rz(1.8809812) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3975911) q[1];
sx q[1];
rz(-0.41644704) q[1];
sx q[1];
rz(1.2214425) q[1];
x q[2];
rz(0.16334857) q[3];
sx q[3];
rz(-2.8046097) q[3];
sx q[3];
rz(2.8548422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5003659) q[2];
sx q[2];
rz(-1.9772823) q[2];
sx q[2];
rz(2.3424303) q[2];
rz(2.7568119) q[3];
sx q[3];
rz(-2.81541) q[3];
sx q[3];
rz(-2.1414115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36050972) q[0];
sx q[0];
rz(-2.0979083) q[0];
sx q[0];
rz(2.5497896) q[0];
rz(-2.7599755) q[1];
sx q[1];
rz(-0.38002574) q[1];
sx q[1];
rz(0.15242481) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7922184) q[0];
sx q[0];
rz(-2.0830941) q[0];
sx q[0];
rz(-2.7854334) q[0];
rz(1.3954787) q[2];
sx q[2];
rz(-2.0421056) q[2];
sx q[2];
rz(-2.5702554) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.287474) q[1];
sx q[1];
rz(-1.2240613) q[1];
sx q[1];
rz(-1.5902014) q[1];
x q[2];
rz(-0.6432047) q[3];
sx q[3];
rz(-0.99785757) q[3];
sx q[3];
rz(-0.3518897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3064208) q[2];
sx q[2];
rz(-0.99205899) q[2];
sx q[2];
rz(1.6888118) q[2];
rz(-0.49550223) q[3];
sx q[3];
rz(-2.317954) q[3];
sx q[3];
rz(2.5775094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38967663) q[0];
sx q[0];
rz(-3.1041807) q[0];
sx q[0];
rz(-2.9801242) q[0];
rz(-0.031919315) q[1];
sx q[1];
rz(-2.4999764) q[1];
sx q[1];
rz(-1.2772824) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.204435) q[0];
sx q[0];
rz(-1.3124518) q[0];
sx q[0];
rz(0.72465557) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1487537) q[2];
sx q[2];
rz(-1.3598056) q[2];
sx q[2];
rz(-0.27560292) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9163852) q[1];
sx q[1];
rz(-0.3233094) q[1];
sx q[1];
rz(2.7955452) q[1];
rz(-pi) q[2];
rz(1.8648009) q[3];
sx q[3];
rz(-1.534003) q[3];
sx q[3];
rz(-0.5816254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45234597) q[2];
sx q[2];
rz(-0.9114868) q[2];
sx q[2];
rz(-2.7455043) q[2];
rz(-0.39997697) q[3];
sx q[3];
rz(-0.59498274) q[3];
sx q[3];
rz(2.2649435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9006627) q[0];
sx q[0];
rz(-0.018095896) q[0];
sx q[0];
rz(-0.15116365) q[0];
rz(2.1647272) q[1];
sx q[1];
rz(-2.7604389) q[1];
sx q[1];
rz(-0.74473286) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.789294) q[0];
sx q[0];
rz(-1.8478824) q[0];
sx q[0];
rz(-2.9267071) q[0];
rz(-pi) q[1];
rz(1.8877541) q[2];
sx q[2];
rz(-0.43192568) q[2];
sx q[2];
rz(-3.0639086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7691466) q[1];
sx q[1];
rz(-2.219104) q[1];
sx q[1];
rz(-1.5689956) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5517057) q[3];
sx q[3];
rz(-2.7253236) q[3];
sx q[3];
rz(1.2225162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.62290827) q[2];
sx q[2];
rz(-1.238995) q[2];
sx q[2];
rz(-0.94804478) q[2];
rz(1.0575804) q[3];
sx q[3];
rz(-1.6653929) q[3];
sx q[3];
rz(1.074056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.093216151) q[0];
sx q[0];
rz(-0.26234782) q[0];
sx q[0];
rz(-0.23396215) q[0];
rz(-1.1853064) q[1];
sx q[1];
rz(-0.92820853) q[1];
sx q[1];
rz(-1.8918461) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7980745) q[0];
sx q[0];
rz(-1.2815223) q[0];
sx q[0];
rz(-0.57855655) q[0];
x q[1];
rz(1.7777257) q[2];
sx q[2];
rz(-0.96038681) q[2];
sx q[2];
rz(2.0774942) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0873336) q[1];
sx q[1];
rz(-2.1076084) q[1];
sx q[1];
rz(2.3786503) q[1];
rz(-pi) q[2];
rz(1.3951282) q[3];
sx q[3];
rz(-0.89699706) q[3];
sx q[3];
rz(-1.2487981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.027111) q[2];
sx q[2];
rz(-1.1610843) q[2];
sx q[2];
rz(1.7005881) q[2];
rz(0.13127413) q[3];
sx q[3];
rz(-0.461853) q[3];
sx q[3];
rz(-0.24391267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.048653614) q[0];
sx q[0];
rz(-0.96554148) q[0];
sx q[0];
rz(1.1337793) q[0];
rz(-1.0406021) q[1];
sx q[1];
rz(-2.2941209) q[1];
sx q[1];
rz(2.6170731) q[1];
rz(1.9807627) q[2];
sx q[2];
rz(-0.66902918) q[2];
sx q[2];
rz(1.4794028) q[2];
rz(2.9374801) q[3];
sx q[3];
rz(-0.64328803) q[3];
sx q[3];
rz(-2.358123) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
