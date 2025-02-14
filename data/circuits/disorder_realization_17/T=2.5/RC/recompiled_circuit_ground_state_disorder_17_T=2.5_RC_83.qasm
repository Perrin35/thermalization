OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.390653) q[0];
sx q[0];
rz(-0.57354623) q[0];
sx q[0];
rz(3.1107922) q[0];
rz(0.38712328) q[1];
sx q[1];
rz(-2.9409599) q[1];
sx q[1];
rz(-0.90955847) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61197554) q[0];
sx q[0];
rz(-2.2155166) q[0];
sx q[0];
rz(-0.84512777) q[0];
rz(-pi) q[1];
x q[1];
rz(2.163226) q[2];
sx q[2];
rz(-0.9610976) q[2];
sx q[2];
rz(-2.040122) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2594257) q[1];
sx q[1];
rz(-1.7710238) q[1];
sx q[1];
rz(3.0012896) q[1];
x q[2];
rz(1.8424787) q[3];
sx q[3];
rz(-1.4702904) q[3];
sx q[3];
rz(-1.0043912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51547852) q[2];
sx q[2];
rz(-1.6823744) q[2];
sx q[2];
rz(0.98575753) q[2];
rz(1.23729) q[3];
sx q[3];
rz(-0.83043778) q[3];
sx q[3];
rz(-1.591709) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80376959) q[0];
sx q[0];
rz(-2.0812806) q[0];
sx q[0];
rz(-2.7983303) q[0];
rz(-0.23974165) q[1];
sx q[1];
rz(-0.38455757) q[1];
sx q[1];
rz(-3.1125617) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0615026) q[0];
sx q[0];
rz(-2.6250589) q[0];
sx q[0];
rz(0.37688984) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.042988792) q[2];
sx q[2];
rz(-2.413297) q[2];
sx q[2];
rz(-2.9790619) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.033054277) q[1];
sx q[1];
rz(-1.5863451) q[1];
sx q[1];
rz(-2.5485746) q[1];
rz(-0.47201158) q[3];
sx q[3];
rz(-2.4084512) q[3];
sx q[3];
rz(-2.5705119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3451738) q[2];
sx q[2];
rz(-0.56946218) q[2];
sx q[2];
rz(-0.84863895) q[2];
rz(0.38550115) q[3];
sx q[3];
rz(-1.4717088) q[3];
sx q[3];
rz(-2.8482385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8141396) q[0];
sx q[0];
rz(-1.3804945) q[0];
sx q[0];
rz(-3.0019548) q[0];
rz(-2.8756554) q[1];
sx q[1];
rz(-1.1775002) q[1];
sx q[1];
rz(-1.0136484) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8089897) q[0];
sx q[0];
rz(-2.4255891) q[0];
sx q[0];
rz(-2.1779275) q[0];
x q[1];
rz(-0.75516197) q[2];
sx q[2];
rz(-0.56729672) q[2];
sx q[2];
rz(0.082956359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2340548) q[1];
sx q[1];
rz(-0.50633303) q[1];
sx q[1];
rz(0.34409006) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4675713) q[3];
sx q[3];
rz(-1.6397392) q[3];
sx q[3];
rz(-0.063223039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2024112) q[2];
sx q[2];
rz(-2.1911669) q[2];
sx q[2];
rz(2.6964296) q[2];
rz(-2.4116481) q[3];
sx q[3];
rz(-1.4533307) q[3];
sx q[3];
rz(-1.9385519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3438123) q[0];
sx q[0];
rz(-0.36425632) q[0];
sx q[0];
rz(-1.9601747) q[0];
rz(-1.5161318) q[1];
sx q[1];
rz(-1.5189891) q[1];
sx q[1];
rz(2.6536062) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5833736) q[0];
sx q[0];
rz(-1.4130371) q[0];
sx q[0];
rz(-2.8789916) q[0];
x q[1];
rz(2.1806966) q[2];
sx q[2];
rz(-0.95003613) q[2];
sx q[2];
rz(1.4767978) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.5778325) q[1];
sx q[1];
rz(-1.8298733) q[1];
sx q[1];
rz(0.96430802) q[1];
rz(-pi) q[2];
rz(-1.2156297) q[3];
sx q[3];
rz(-1.9653335) q[3];
sx q[3];
rz(1.9882787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.92859149) q[2];
sx q[2];
rz(-1.4139621) q[2];
sx q[2];
rz(-0.50814) q[2];
rz(-2.7022434) q[3];
sx q[3];
rz(-2.5033247) q[3];
sx q[3];
rz(-2.5590069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4947434) q[0];
sx q[0];
rz(-1.9288394) q[0];
sx q[0];
rz(-2.5272722) q[0];
rz(0.99614227) q[1];
sx q[1];
rz(-0.83729815) q[1];
sx q[1];
rz(-3.0706629) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1114601) q[0];
sx q[0];
rz(-1.8778525) q[0];
sx q[0];
rz(2.2024406) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28648232) q[2];
sx q[2];
rz(-1.9989609) q[2];
sx q[2];
rz(-0.78543898) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48551895) q[1];
sx q[1];
rz(-0.89886256) q[1];
sx q[1];
rz(-0.80826064) q[1];
x q[2];
rz(0.67885375) q[3];
sx q[3];
rz(-2.1311789) q[3];
sx q[3];
rz(-2.7075775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7190242) q[2];
sx q[2];
rz(-1.7302128) q[2];
sx q[2];
rz(-3.0651429) q[2];
rz(3.0494173) q[3];
sx q[3];
rz(-2.01391) q[3];
sx q[3];
rz(2.5341212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344264) q[0];
sx q[0];
rz(-0.18944117) q[0];
sx q[0];
rz(-2.6763951) q[0];
rz(0.094296433) q[1];
sx q[1];
rz(-2.1262719) q[1];
sx q[1];
rz(-1.2634855) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91335653) q[0];
sx q[0];
rz(-2.7688518) q[0];
sx q[0];
rz(-1.2380225) q[0];
rz(-pi) q[1];
rz(0.0025778254) q[2];
sx q[2];
rz(-0.63226223) q[2];
sx q[2];
rz(-2.8272976) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.24119995) q[1];
sx q[1];
rz(-2.3896895) q[1];
sx q[1];
rz(-0.65691595) q[1];
rz(1.1912548) q[3];
sx q[3];
rz(-2.0058245) q[3];
sx q[3];
rz(-1.3635927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5491526) q[2];
sx q[2];
rz(-1.1391888) q[2];
sx q[2];
rz(2.450954) q[2];
rz(-1.6820827) q[3];
sx q[3];
rz(-3.1091318) q[3];
sx q[3];
rz(2.8800817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.9948147) q[0];
sx q[0];
rz(-2.1187466) q[0];
sx q[0];
rz(-0.47635517) q[0];
rz(-1.1789383) q[1];
sx q[1];
rz(-2.4504688) q[1];
sx q[1];
rz(-1.8144511) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.961273) q[0];
sx q[0];
rz(-1.1256582) q[0];
sx q[0];
rz(-2.0044786) q[0];
rz(-0.05169241) q[2];
sx q[2];
rz(-1.0982656) q[2];
sx q[2];
rz(-0.57628265) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.33846291) q[1];
sx q[1];
rz(-1.6410562) q[1];
sx q[1];
rz(-0.52609013) q[1];
rz(0.031257625) q[3];
sx q[3];
rz(-0.61504902) q[3];
sx q[3];
rz(0.94155967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2771161) q[2];
sx q[2];
rz(-1.4863374) q[2];
sx q[2];
rz(2.3790996) q[2];
rz(-1.3660733) q[3];
sx q[3];
rz(-1.450918) q[3];
sx q[3];
rz(-2.8097927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45258006) q[0];
sx q[0];
rz(-3.0735425) q[0];
sx q[0];
rz(-0.77823234) q[0];
rz(1.7315841) q[1];
sx q[1];
rz(-2.2739613) q[1];
sx q[1];
rz(2.5708139) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3616791) q[0];
sx q[0];
rz(-2.329149) q[0];
sx q[0];
rz(2.1606584) q[0];
rz(-1.2879175) q[2];
sx q[2];
rz(-1.6215417) q[2];
sx q[2];
rz(3.1078448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.135317) q[1];
sx q[1];
rz(-2.0915739) q[1];
sx q[1];
rz(-3.0515677) q[1];
rz(0.93511411) q[3];
sx q[3];
rz(-1.550727) q[3];
sx q[3];
rz(2.3754313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0290231) q[2];
sx q[2];
rz(-1.7971635) q[2];
sx q[2];
rz(-0.53978187) q[2];
rz(2.6953186) q[3];
sx q[3];
rz(-0.73791426) q[3];
sx q[3];
rz(0.81587434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096263252) q[0];
sx q[0];
rz(-0.49880609) q[0];
sx q[0];
rz(-2.2742284) q[0];
rz(-1.7034886) q[1];
sx q[1];
rz(-2.4080364) q[1];
sx q[1];
rz(-2.6677456) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4854914) q[0];
sx q[0];
rz(-0.37603077) q[0];
sx q[0];
rz(0.72304) q[0];
rz(-pi) q[1];
rz(-1.0269073) q[2];
sx q[2];
rz(-0.68008125) q[2];
sx q[2];
rz(0.56267738) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9349228) q[1];
sx q[1];
rz(-0.75074544) q[1];
sx q[1];
rz(-2.9331371) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2603277) q[3];
sx q[3];
rz(-1.9816958) q[3];
sx q[3];
rz(-1.9104065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4465686) q[2];
sx q[2];
rz(-0.5646255) q[2];
sx q[2];
rz(2.0107644) q[2];
rz(-1.1787777) q[3];
sx q[3];
rz(-1.7480353) q[3];
sx q[3];
rz(2.668837) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-0.21993615) q[0];
sx q[0];
rz(-0.49880767) q[0];
sx q[0];
rz(-0.98861277) q[0];
rz(-1.0722718) q[1];
sx q[1];
rz(-1.5360473) q[1];
sx q[1];
rz(-1.706121) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0611183) q[0];
sx q[0];
rz(-0.99493626) q[0];
sx q[0];
rz(3.1115576) q[0];
rz(0.033013952) q[2];
sx q[2];
rz(-2.9131914) q[2];
sx q[2];
rz(1.2778801) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.86462532) q[1];
sx q[1];
rz(-1.5706148) q[1];
sx q[1];
rz(1.5110043) q[1];
rz(-0.39809583) q[3];
sx q[3];
rz(-2.4090066) q[3];
sx q[3];
rz(1.7375515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84206218) q[2];
sx q[2];
rz(-0.64075035) q[2];
sx q[2];
rz(0.93471849) q[2];
rz(1.0457906) q[3];
sx q[3];
rz(-0.17242923) q[3];
sx q[3];
rz(2.893253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0228148) q[0];
sx q[0];
rz(-1.7775443) q[0];
sx q[0];
rz(-1.4015629) q[0];
rz(0.10764311) q[1];
sx q[1];
rz(-1.6247152) q[1];
sx q[1];
rz(-2.7693137) q[1];
rz(-1.2992181) q[2];
sx q[2];
rz(-2.3403946) q[2];
sx q[2];
rz(2.8190294) q[2];
rz(-2.6519898) q[3];
sx q[3];
rz(-1.0013178) q[3];
sx q[3];
rz(2.8538741) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
