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
rz(0.060184181) q[0];
sx q[0];
rz(4.2005778) q[0];
sx q[0];
rz(8.4146001) q[0];
rz(-0.8085568) q[1];
sx q[1];
rz(-0.29616907) q[1];
sx q[1];
rz(0.30997601) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55059075) q[0];
sx q[0];
rz(-2.7059116) q[0];
sx q[0];
rz(-2.8753619) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0344738) q[2];
sx q[2];
rz(-2.9351882) q[2];
sx q[2];
rz(0.53910461) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9604608) q[1];
sx q[1];
rz(-2.7560661) q[1];
sx q[1];
rz(2.1147644) q[1];
rz(1.3499851) q[3];
sx q[3];
rz(-2.3901) q[3];
sx q[3];
rz(1.3205075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4787204) q[2];
sx q[2];
rz(-2.8933849) q[2];
sx q[2];
rz(-0.09566801) q[2];
rz(0.11224789) q[3];
sx q[3];
rz(-0.87533689) q[3];
sx q[3];
rz(-1.7101425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9005168) q[0];
sx q[0];
rz(-2.2523585) q[0];
sx q[0];
rz(0.61479968) q[0];
rz(-0.88848937) q[1];
sx q[1];
rz(-1.588984) q[1];
sx q[1];
rz(-0.49957553) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0083778) q[0];
sx q[0];
rz(-1.7724123) q[0];
sx q[0];
rz(-1.9033543) q[0];
x q[1];
rz(0.26878727) q[2];
sx q[2];
rz(-0.54020451) q[2];
sx q[2];
rz(-2.7160783) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0608637) q[1];
sx q[1];
rz(-1.398096) q[1];
sx q[1];
rz(-0.91367803) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1533215) q[3];
sx q[3];
rz(-1.8559905) q[3];
sx q[3];
rz(2.5799283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2483612) q[2];
sx q[2];
rz(-1.2373368) q[2];
sx q[2];
rz(-3.1207747) q[2];
rz(-1.2402395) q[3];
sx q[3];
rz(-1.4695243) q[3];
sx q[3];
rz(1.9564691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6388539) q[0];
sx q[0];
rz(-2.8732712) q[0];
sx q[0];
rz(2.8079206) q[0];
rz(2.4036713) q[1];
sx q[1];
rz(-1.9155904) q[1];
sx q[1];
rz(-0.12942448) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20724587) q[0];
sx q[0];
rz(-0.6967623) q[0];
sx q[0];
rz(1.3214701) q[0];
rz(-2.1652075) q[2];
sx q[2];
rz(-0.98862851) q[2];
sx q[2];
rz(-0.79603031) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8379535) q[1];
sx q[1];
rz(-2.9898414) q[1];
sx q[1];
rz(0.93317731) q[1];
rz(-0.47862847) q[3];
sx q[3];
rz(-0.84175292) q[3];
sx q[3];
rz(-2.7117376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0176598) q[2];
sx q[2];
rz(-1.61597) q[2];
sx q[2];
rz(1.7714436) q[2];
rz(-2.395199) q[3];
sx q[3];
rz(-1.8236225) q[3];
sx q[3];
rz(-3.0588176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.133404) q[0];
sx q[0];
rz(-1.9318102) q[0];
sx q[0];
rz(-0.56513894) q[0];
rz(-0.42916974) q[1];
sx q[1];
rz(-2.4950835) q[1];
sx q[1];
rz(0.82829222) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1734742) q[0];
sx q[0];
rz(-0.83195639) q[0];
sx q[0];
rz(1.678874) q[0];
rz(-pi) q[1];
rz(-2.8548334) q[2];
sx q[2];
rz(-1.1017403) q[2];
sx q[2];
rz(0.048407528) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2678623) q[1];
sx q[1];
rz(-0.56247382) q[1];
sx q[1];
rz(0.51452325) q[1];
x q[2];
rz(2.8935562) q[3];
sx q[3];
rz(-2.0741077) q[3];
sx q[3];
rz(-0.65500427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68507489) q[2];
sx q[2];
rz(-1.061331) q[2];
sx q[2];
rz(1.7100517) q[2];
rz(-0.93959129) q[3];
sx q[3];
rz(-2.6288433) q[3];
sx q[3];
rz(0.00069869839) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13570304) q[0];
sx q[0];
rz(-2.8764184) q[0];
sx q[0];
rz(1.3294719) q[0];
rz(1.1520518) q[1];
sx q[1];
rz(-2.0538797) q[1];
sx q[1];
rz(0.00024814127) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5329929) q[0];
sx q[0];
rz(-1.8155351) q[0];
sx q[0];
rz(2.9518963) q[0];
x q[1];
rz(-2.1049961) q[2];
sx q[2];
rz(-1.4956258) q[2];
sx q[2];
rz(-0.49803621) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.21890404) q[1];
sx q[1];
rz(-2.5210025) q[1];
sx q[1];
rz(2.1922078) q[1];
rz(0.080623589) q[3];
sx q[3];
rz(-2.4084457) q[3];
sx q[3];
rz(2.6997363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0033215) q[2];
sx q[2];
rz(-1.2500117) q[2];
sx q[2];
rz(3.1357583) q[2];
rz(2.2516069) q[3];
sx q[3];
rz(-0.68166387) q[3];
sx q[3];
rz(-2.0148923) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24790813) q[0];
sx q[0];
rz(-1.4434781) q[0];
sx q[0];
rz(-1.4877315) q[0];
rz(-3.1112025) q[1];
sx q[1];
rz(-1.1266212) q[1];
sx q[1];
rz(-0.94246513) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1346325) q[0];
sx q[0];
rz(-2.83036) q[0];
sx q[0];
rz(0.65629058) q[0];
rz(1.8293657) q[2];
sx q[2];
rz(-1.8601918) q[2];
sx q[2];
rz(-1.0602151) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8854495) q[1];
sx q[1];
rz(-2.0383308) q[1];
sx q[1];
rz(-2.8199151) q[1];
rz(-pi) q[2];
rz(-2.0464155) q[3];
sx q[3];
rz(-1.3381357) q[3];
sx q[3];
rz(1.7847248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5421062) q[2];
sx q[2];
rz(-0.78080559) q[2];
sx q[2];
rz(-1.3055118) q[2];
rz(2.5944338) q[3];
sx q[3];
rz(-1.2055509) q[3];
sx q[3];
rz(-2.3356596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
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
rz(0.18043537) q[0];
sx q[0];
rz(-2.1077709) q[0];
sx q[0];
rz(-0.10051522) q[0];
rz(0.48249498) q[1];
sx q[1];
rz(-1.9827739) q[1];
sx q[1];
rz(-0.85711342) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.094497546) q[0];
sx q[0];
rz(-2.408354) q[0];
sx q[0];
rz(-0.49219699) q[0];
rz(-pi) q[1];
rz(2.4947462) q[2];
sx q[2];
rz(-1.610272) q[2];
sx q[2];
rz(2.0496862) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4883721) q[1];
sx q[1];
rz(-1.778942) q[1];
sx q[1];
rz(2.7339762) q[1];
rz(-1.1825652) q[3];
sx q[3];
rz(-2.534158) q[3];
sx q[3];
rz(1.5807815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7414005) q[2];
sx q[2];
rz(-0.97101784) q[2];
sx q[2];
rz(-0.3024438) q[2];
rz(2.1619469) q[3];
sx q[3];
rz(-2.1392348) q[3];
sx q[3];
rz(-1.2790595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.361146) q[0];
sx q[0];
rz(-0.13872153) q[0];
sx q[0];
rz(-3.1106023) q[0];
rz(0.57394761) q[1];
sx q[1];
rz(-1.6684883) q[1];
sx q[1];
rz(1.2773638) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63420682) q[0];
sx q[0];
rz(-1.9045826) q[0];
sx q[0];
rz(0.35315634) q[0];
x q[1];
rz(1.9979111) q[2];
sx q[2];
rz(-0.53887212) q[2];
sx q[2];
rz(0.32921916) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9121418) q[1];
sx q[1];
rz(-1.0849285) q[1];
sx q[1];
rz(-3.058601) q[1];
x q[2];
rz(1.7788309) q[3];
sx q[3];
rz(-0.29104189) q[3];
sx q[3];
rz(-2.7179949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.10451) q[2];
sx q[2];
rz(-3.0057378) q[2];
sx q[2];
rz(2.6541397) q[2];
rz(-0.69665748) q[3];
sx q[3];
rz(-0.928855) q[3];
sx q[3];
rz(-0.063974403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.071851991) q[0];
sx q[0];
rz(-0.47436473) q[0];
sx q[0];
rz(-2.790614) q[0];
rz(-2.9674496) q[1];
sx q[1];
rz(-1.6330279) q[1];
sx q[1];
rz(-1.4016271) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5183254) q[0];
sx q[0];
rz(-1.4484157) q[0];
sx q[0];
rz(-0.24952475) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0171595) q[2];
sx q[2];
rz(-2.2702262) q[2];
sx q[2];
rz(-2.1289189) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8216711) q[1];
sx q[1];
rz(-0.18583365) q[1];
sx q[1];
rz(1.7286848) q[1];
rz(-0.97934874) q[3];
sx q[3];
rz(-0.92100793) q[3];
sx q[3];
rz(-1.5929008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1104687) q[2];
sx q[2];
rz(-2.4565171) q[2];
sx q[2];
rz(-0.6366716) q[2];
rz(1.1104256) q[3];
sx q[3];
rz(-1.597155) q[3];
sx q[3];
rz(0.5184263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7745895) q[0];
sx q[0];
rz(-2.84802) q[0];
sx q[0];
rz(1.0585744) q[0];
rz(0.038657945) q[1];
sx q[1];
rz(-1.5267173) q[1];
sx q[1];
rz(-1.0640594) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2702613) q[0];
sx q[0];
rz(-0.0055905213) q[0];
sx q[0];
rz(-2.5020775) q[0];
x q[1];
rz(0.27711192) q[2];
sx q[2];
rz(-1.849035) q[2];
sx q[2];
rz(0.81467512) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.41760379) q[1];
sx q[1];
rz(-3.0688802) q[1];
sx q[1];
rz(-2.5515208) q[1];
x q[2];
rz(-2.6820716) q[3];
sx q[3];
rz(-1.6020613) q[3];
sx q[3];
rz(0.36324901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.71558636) q[2];
sx q[2];
rz(-0.49804372) q[2];
sx q[2];
rz(-3.0675724) q[2];
rz(0.38153875) q[3];
sx q[3];
rz(-1.3149202) q[3];
sx q[3];
rz(-2.7874302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.030180177) q[0];
sx q[0];
rz(-1.3127865) q[0];
sx q[0];
rz(2.4319613) q[0];
rz(-2.5297655) q[1];
sx q[1];
rz(-2.430293) q[1];
sx q[1];
rz(-1.5878955) q[1];
rz(1.189497) q[2];
sx q[2];
rz(-1.0972037) q[2];
sx q[2];
rz(0.9158132) q[2];
rz(-2.6388219) q[3];
sx q[3];
rz(-0.81450653) q[3];
sx q[3];
rz(-0.52342879) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
