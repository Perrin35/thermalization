OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(1.8811037) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(-1.7927875) q[1];
sx q[1];
rz(-0.92372149) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1200861) q[0];
sx q[0];
rz(-1.2354295) q[0];
sx q[0];
rz(2.04727) q[0];
rz(-1.2878296) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(-1.1793009) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.87286283) q[1];
sx q[1];
rz(-1.1464719) q[1];
sx q[1];
rz(-2.5904168) q[1];
rz(1.190891) q[3];
sx q[3];
rz(-1.1067179) q[3];
sx q[3];
rz(-2.3604148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(2.2062733) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(-0.67990047) q[1];
sx q[1];
rz(-0.49566832) q[1];
sx q[1];
rz(-1.4555567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50549492) q[0];
sx q[0];
rz(-0.60244766) q[0];
sx q[0];
rz(1.2732182) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37462072) q[2];
sx q[2];
rz(-1.6472367) q[2];
sx q[2];
rz(-2.3945216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93745366) q[1];
sx q[1];
rz(-2.0935645) q[1];
sx q[1];
rz(-0.015603113) q[1];
rz(-pi) q[2];
rz(1.5973813) q[3];
sx q[3];
rz(-1.1592602) q[3];
sx q[3];
rz(-0.4707903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(1.7896174) q[2];
rz(-0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(2.8296208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882554) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(1.172539) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3282933) q[0];
sx q[0];
rz(-1.0070224) q[0];
sx q[0];
rz(-1.203712) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.34611361) q[2];
sx q[2];
rz(-2.3556404) q[2];
sx q[2];
rz(-1.9217938) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8923924) q[1];
sx q[1];
rz(-1.8141659) q[1];
sx q[1];
rz(2.1057486) q[1];
rz(-pi) q[2];
rz(-3.0828589) q[3];
sx q[3];
rz(-1.8211094) q[3];
sx q[3];
rz(-0.55610031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0744434) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(2.2303936) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(-0.12810853) q[0];
rz(-0.076106636) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(0.52350837) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0353161) q[0];
sx q[0];
rz(-1.9395394) q[0];
sx q[0];
rz(0.62555255) q[0];
x q[1];
rz(-0.24387118) q[2];
sx q[2];
rz(-0.3393617) q[2];
sx q[2];
rz(1.7983758) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.6701811) q[1];
sx q[1];
rz(-1.3172611) q[1];
sx q[1];
rz(-0.86004911) q[1];
rz(-pi) q[2];
rz(-1.773049) q[3];
sx q[3];
rz(-0.60086717) q[3];
sx q[3];
rz(-0.7522538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(-2.5775487) q[2];
rz(-0.28856746) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600835) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(1.4915285) q[0];
rz(2.2619757) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-2.1496444) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0569699) q[0];
sx q[0];
rz(-1.3966494) q[0];
sx q[0];
rz(-1.6790381) q[0];
rz(-pi) q[1];
rz(1.8736585) q[2];
sx q[2];
rz(-1.5793243) q[2];
sx q[2];
rz(-1.0964583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.76584133) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(2.8326616) q[1];
rz(-pi) q[2];
rz(-0.81327849) q[3];
sx q[3];
rz(-1.4762029) q[3];
sx q[3];
rz(2.7300342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-0.36964881) q[2];
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
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.4145684) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(1.3487934) q[0];
rz(-2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(-1.7165002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3544918) q[0];
sx q[0];
rz(-2.0658501) q[0];
sx q[0];
rz(1.8255193) q[0];
rz(1.8259949) q[2];
sx q[2];
rz(-1.6591424) q[2];
sx q[2];
rz(-1.7905854) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.52275601) q[1];
sx q[1];
rz(-0.4546051) q[1];
sx q[1];
rz(1.7230117) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1720042) q[3];
sx q[3];
rz(-1.117327) q[3];
sx q[3];
rz(2.3525402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(-2.5777204) q[2];
rz(2.9610736) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(-0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5086223) q[0];
sx q[0];
rz(-2.9794725) q[0];
sx q[0];
rz(-0.41931835) q[0];
rz(-1.5527027) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(-2.3197876) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99188995) q[0];
sx q[0];
rz(-0.8677965) q[0];
sx q[0];
rz(2.4292612) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6452438) q[2];
sx q[2];
rz(-2.2349572) q[2];
sx q[2];
rz(1.6830483) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2591178) q[1];
sx q[1];
rz(-1.3438517) q[1];
sx q[1];
rz(-0.74101733) q[1];
rz(1.9415226) q[3];
sx q[3];
rz(-1.0154187) q[3];
sx q[3];
rz(0.45645082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3372779) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(3.0120567) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.725175) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(-0.30934632) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.8364505) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7628521) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(-0.61600323) q[0];
x q[1];
rz(-1.9403946) q[2];
sx q[2];
rz(-0.9005138) q[2];
sx q[2];
rz(-3.0073462) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.0933989) q[1];
sx q[1];
rz(-1.9462002) q[1];
sx q[1];
rz(2.4673389) q[1];
rz(1.5042138) q[3];
sx q[3];
rz(-0.3399907) q[3];
sx q[3];
rz(0.41223994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(-1.6513599) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(2.819678) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-0.70294356) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1212595) q[0];
sx q[0];
rz(-2.5998305) q[0];
sx q[0];
rz(-0.74777491) q[0];
x q[1];
rz(1.0614971) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(2.408037) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.96303899) q[1];
sx q[1];
rz(-0.56561618) q[1];
sx q[1];
rz(1.2600684) q[1];
x q[2];
rz(0.95548198) q[3];
sx q[3];
rz(-1.8468879) q[3];
sx q[3];
rz(2.626782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(1.8019603) q[2];
rz(-2.83589) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9777578) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(0.90351358) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(-0.46863619) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7077431) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(-1.7953403) q[0];
rz(-pi) q[1];
rz(-2.8154545) q[2];
sx q[2];
rz(-1.9378127) q[2];
sx q[2];
rz(1.067576) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.84506449) q[1];
sx q[1];
rz(-1.6390641) q[1];
sx q[1];
rz(-1.3045842) q[1];
rz(-pi) q[2];
rz(-2.7077984) q[3];
sx q[3];
rz(-1.1489831) q[3];
sx q[3];
rz(-2.1123561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.998385) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(-1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6464012) q[0];
sx q[0];
rz(-1.1947182) q[0];
sx q[0];
rz(2.4583046) q[0];
rz(0.62190965) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(-2.3399578) q[2];
sx q[2];
rz(-0.87052204) q[2];
sx q[2];
rz(-1.082765) q[2];
rz(0.078483742) q[3];
sx q[3];
rz(-2.2208636) q[3];
sx q[3];
rz(0.29490864) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
