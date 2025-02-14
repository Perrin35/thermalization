OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.2709687) q[0];
sx q[0];
rz(-0.55611098) q[0];
sx q[0];
rz(2.1882353) q[0];
rz(0.019729992) q[1];
sx q[1];
rz(-2.1058197) q[1];
sx q[1];
rz(0.9486202) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.079049) q[0];
sx q[0];
rz(-1.9808852) q[0];
sx q[0];
rz(-2.5586302) q[0];
rz(-3.1300053) q[2];
sx q[2];
rz(-0.23749781) q[2];
sx q[2];
rz(1.604014) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.030101209) q[1];
sx q[1];
rz(-1.1426569) q[1];
sx q[1];
rz(0.74049048) q[1];
rz(-2.7953732) q[3];
sx q[3];
rz(-1.0448714) q[3];
sx q[3];
rz(3.0478551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76401508) q[2];
sx q[2];
rz(-3.0674051) q[2];
sx q[2];
rz(-0.78262502) q[2];
rz(0.11893663) q[3];
sx q[3];
rz(-1.0340034) q[3];
sx q[3];
rz(2.9940166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19621944) q[0];
sx q[0];
rz(-2.0202899) q[0];
sx q[0];
rz(-2.8837606) q[0];
rz(-3.0505772) q[1];
sx q[1];
rz(-1.0992522) q[1];
sx q[1];
rz(-1.4965422) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.634406) q[0];
sx q[0];
rz(-1.0189462) q[0];
sx q[0];
rz(-3.0613012) q[0];
x q[1];
rz(-0.74380959) q[2];
sx q[2];
rz(-2.653476) q[2];
sx q[2];
rz(1.0783735) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61747256) q[1];
sx q[1];
rz(-1.5407498) q[1];
sx q[1];
rz(-1.9829795) q[1];
rz(-pi) q[2];
rz(-2.6378651) q[3];
sx q[3];
rz(-2.4085975) q[3];
sx q[3];
rz(-0.021857787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0195134) q[2];
sx q[2];
rz(-1.2747108) q[2];
sx q[2];
rz(-0.40461928) q[2];
rz(0.98958611) q[3];
sx q[3];
rz(-1.3000969) q[3];
sx q[3];
rz(-2.5185481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0290381) q[0];
sx q[0];
rz(-2.9974388) q[0];
sx q[0];
rz(-0.83830225) q[0];
rz(-0.84838947) q[1];
sx q[1];
rz(-1.219607) q[1];
sx q[1];
rz(0.99260509) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2783689) q[0];
sx q[0];
rz(-1.430817) q[0];
sx q[0];
rz(-0.72576688) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8072629) q[2];
sx q[2];
rz(-1.624776) q[2];
sx q[2];
rz(-2.8950429) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5287244) q[1];
sx q[1];
rz(-0.49002417) q[1];
sx q[1];
rz(-2.571066) q[1];
rz(-pi) q[2];
rz(-0.35723585) q[3];
sx q[3];
rz(-0.30870507) q[3];
sx q[3];
rz(1.5451408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.8695716) q[2];
sx q[2];
rz(-1.973899) q[2];
sx q[2];
rz(-2.0167548) q[2];
rz(0.62075067) q[3];
sx q[3];
rz(-0.97095942) q[3];
sx q[3];
rz(3.0264405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2426185) q[0];
sx q[0];
rz(-1.1599351) q[0];
sx q[0];
rz(-0.17380357) q[0];
rz(-2.4558892) q[1];
sx q[1];
rz(-1.655429) q[1];
sx q[1];
rz(2.3033843) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7156775) q[0];
sx q[0];
rz(-2.0260149) q[0];
sx q[0];
rz(1.0550189) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.051866) q[2];
sx q[2];
rz(-2.0806667) q[2];
sx q[2];
rz(3.1367658) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.3633109) q[1];
sx q[1];
rz(-1.4290591) q[1];
sx q[1];
rz(-2.4635386) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.066579) q[3];
sx q[3];
rz(-2.677644) q[3];
sx q[3];
rz(1.848965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.54245) q[2];
sx q[2];
rz(-1.6833545) q[2];
sx q[2];
rz(0.35169265) q[2];
rz(1.1951949) q[3];
sx q[3];
rz(-1.1870793) q[3];
sx q[3];
rz(-0.024141969) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.675932) q[0];
sx q[0];
rz(-2.6213578) q[0];
sx q[0];
rz(-2.7886673) q[0];
rz(-2.832761) q[1];
sx q[1];
rz(-2.1052723) q[1];
sx q[1];
rz(0.028506361) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7362979) q[0];
sx q[0];
rz(-0.60966821) q[0];
sx q[0];
rz(1.8788615) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1813004) q[2];
sx q[2];
rz(-1.2151698) q[2];
sx q[2];
rz(0.86330044) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9069017) q[1];
sx q[1];
rz(-2.3447835) q[1];
sx q[1];
rz(1.2654348) q[1];
x q[2];
rz(-1.1358374) q[3];
sx q[3];
rz(-1.7925781) q[3];
sx q[3];
rz(1.3835761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9286524) q[2];
sx q[2];
rz(-0.063455909) q[2];
sx q[2];
rz(2.1707936) q[2];
rz(2.094723) q[3];
sx q[3];
rz(-1.6887083) q[3];
sx q[3];
rz(0.48986062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8397119) q[0];
sx q[0];
rz(-1.3588926) q[0];
sx q[0];
rz(-0.090959892) q[0];
rz(-1.2184527) q[1];
sx q[1];
rz(-0.33088845) q[1];
sx q[1];
rz(-2.5023696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42933336) q[0];
sx q[0];
rz(-0.3540701) q[0];
sx q[0];
rz(2.7430915) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2940862) q[2];
sx q[2];
rz(-2.9581262) q[2];
sx q[2];
rz(-3.1364721) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5167349) q[1];
sx q[1];
rz(-1.0565041) q[1];
sx q[1];
rz(-0.74511294) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2442144) q[3];
sx q[3];
rz(-1.5726798) q[3];
sx q[3];
rz(2.7060946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.062139221) q[2];
sx q[2];
rz(-0.97525758) q[2];
sx q[2];
rz(-2.6677168) q[2];
rz(1.3300995) q[3];
sx q[3];
rz(-2.2606943) q[3];
sx q[3];
rz(-2.7661095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5652931) q[0];
sx q[0];
rz(-1.0813035) q[0];
sx q[0];
rz(2.9190049) q[0];
rz(-1.0635771) q[1];
sx q[1];
rz(-1.5779326) q[1];
sx q[1];
rz(1.1933914) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1656859) q[0];
sx q[0];
rz(-0.75729174) q[0];
sx q[0];
rz(2.971847) q[0];
rz(3.0203825) q[2];
sx q[2];
rz(-1.796399) q[2];
sx q[2];
rz(0.76162213) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8377922) q[1];
sx q[1];
rz(-1.1322316) q[1];
sx q[1];
rz(-1.720251) q[1];
rz(-pi) q[2];
rz(1.0289331) q[3];
sx q[3];
rz(-0.30908424) q[3];
sx q[3];
rz(-2.2268471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.1412763) q[2];
sx q[2];
rz(-0.87339425) q[2];
sx q[2];
rz(-2.5013962) q[2];
rz(-0.021942465) q[3];
sx q[3];
rz(-1.4098189) q[3];
sx q[3];
rz(-0.35023165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.491965) q[0];
sx q[0];
rz(-1.7438629) q[0];
sx q[0];
rz(-0.52038991) q[0];
rz(0.87977663) q[1];
sx q[1];
rz(-1.2326515) q[1];
sx q[1];
rz(-0.89281503) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3443106) q[0];
sx q[0];
rz(-1.6160674) q[0];
sx q[0];
rz(-3.0540953) q[0];
rz(-0.73142902) q[2];
sx q[2];
rz(-0.9584223) q[2];
sx q[2];
rz(2.2867672) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.5984447) q[1];
sx q[1];
rz(-2.3259427) q[1];
sx q[1];
rz(-1.1786908) q[1];
rz(-0.67075394) q[3];
sx q[3];
rz(-0.62707147) q[3];
sx q[3];
rz(-2.1334836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8799379) q[2];
sx q[2];
rz(-2.0377906) q[2];
sx q[2];
rz(-0.086816303) q[2];
rz(-2.1964729) q[3];
sx q[3];
rz(-0.34543959) q[3];
sx q[3];
rz(-2.0115578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18411186) q[0];
sx q[0];
rz(-0.090395398) q[0];
sx q[0];
rz(-0.064067319) q[0];
rz(-2.0059026) q[1];
sx q[1];
rz(-1.7013197) q[1];
sx q[1];
rz(0.51220977) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24146809) q[0];
sx q[0];
rz(-1.9808123) q[0];
sx q[0];
rz(1.2154237) q[0];
rz(0.86821235) q[2];
sx q[2];
rz(-1.3033452) q[2];
sx q[2];
rz(0.78154678) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.71609771) q[1];
sx q[1];
rz(-1.4852045) q[1];
sx q[1];
rz(-2.8817428) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5867611) q[3];
sx q[3];
rz(-1.1619688) q[3];
sx q[3];
rz(2.7328934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3247437) q[2];
sx q[2];
rz(-0.73255676) q[2];
sx q[2];
rz(2.0762439) q[2];
rz(-1.6522853) q[3];
sx q[3];
rz(-2.1842712) q[3];
sx q[3];
rz(-3.0491507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5866933) q[0];
sx q[0];
rz(-0.41189343) q[0];
sx q[0];
rz(2.8107585) q[0];
rz(-1.2384442) q[1];
sx q[1];
rz(-0.60791433) q[1];
sx q[1];
rz(-0.27473658) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3150644) q[0];
sx q[0];
rz(-1.2312268) q[0];
sx q[0];
rz(1.1134321) q[0];
rz(-pi) q[1];
rz(3.0255453) q[2];
sx q[2];
rz(-0.2751285) q[2];
sx q[2];
rz(-2.3809718) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.7995791) q[1];
sx q[1];
rz(-1.430961) q[1];
sx q[1];
rz(-1.0619771) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.072321691) q[3];
sx q[3];
rz(-2.157271) q[3];
sx q[3];
rz(-0.14432913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1642509) q[2];
sx q[2];
rz(-0.75426102) q[2];
sx q[2];
rz(-1.3573307) q[2];
rz(-0.50061289) q[3];
sx q[3];
rz(-1.5631792) q[3];
sx q[3];
rz(-1.4969426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5464583) q[0];
sx q[0];
rz(-1.5820137) q[0];
sx q[0];
rz(-0.59376846) q[0];
rz(-0.068269923) q[1];
sx q[1];
rz(-1.2208114) q[1];
sx q[1];
rz(0.023718871) q[1];
rz(-2.7958128) q[2];
sx q[2];
rz(-1.4152534) q[2];
sx q[2];
rz(-1.5243668) q[2];
rz(-0.81238644) q[3];
sx q[3];
rz(-2.8858622) q[3];
sx q[3];
rz(1.1271911) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
