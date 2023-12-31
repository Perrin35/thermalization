OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(1.5078478) q[0];
rz(0.047343407) q[1];
sx q[1];
rz(-2.3634057) q[1];
sx q[1];
rz(0.49931061) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9141657) q[0];
sx q[0];
rz(-1.6343071) q[0];
sx q[0];
rz(2.8660197) q[0];
rz(2.9688641) q[2];
sx q[2];
rz(-0.85731259) q[2];
sx q[2];
rz(1.8634862) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0807304) q[1];
sx q[1];
rz(-0.84134358) q[1];
sx q[1];
rz(-1.7727477) q[1];
rz(-1.6742168) q[3];
sx q[3];
rz(-1.5973583) q[3];
sx q[3];
rz(2.6288222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(-0.23400083) q[3];
sx q[3];
rz(-2.6205385) q[3];
sx q[3];
rz(-2.8570989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.6750591) q[0];
sx q[0];
rz(1.0043253) q[0];
rz(-1.6218119) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(1.0027592) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5402055) q[0];
sx q[0];
rz(-2.0048855) q[0];
sx q[0];
rz(0.65625221) q[0];
rz(-pi) q[1];
rz(-1.7582714) q[2];
sx q[2];
rz(-1.4008998) q[2];
sx q[2];
rz(1.0152917) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1042852) q[1];
sx q[1];
rz(-2.1530495) q[1];
sx q[1];
rz(-0.18584713) q[1];
rz(-pi) q[2];
rz(1.2172132) q[3];
sx q[3];
rz(-2.1401322) q[3];
sx q[3];
rz(-0.47488892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8721547) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(-0.41444591) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(0.088236563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(1.8931483) q[0];
sx q[0];
rz(-1.8284766) q[0];
sx q[0];
rz(-2.3625968) q[0];
rz(0.39930725) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-2.2580106) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1012264) q[0];
sx q[0];
rz(-1.7352312) q[0];
sx q[0];
rz(2.7357487) q[0];
rz(-pi) q[1];
rz(-2.8475464) q[2];
sx q[2];
rz(-1.9532734) q[2];
sx q[2];
rz(-1.8665714) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.59144943) q[1];
sx q[1];
rz(-1.5482229) q[1];
sx q[1];
rz(-1.8557465) q[1];
x q[2];
rz(1.4530573) q[3];
sx q[3];
rz(-0.91500926) q[3];
sx q[3];
rz(1.9829139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.3016004) q[2];
sx q[2];
rz(-1.8847382) q[2];
sx q[2];
rz(2.7123614) q[2];
rz(-2.1515576) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(2.278573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7097968) q[0];
sx q[0];
rz(-1.7683832) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(2.0344095) q[1];
sx q[1];
rz(-2.2829843) q[1];
sx q[1];
rz(2.591419) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.195939) q[0];
sx q[0];
rz(-1.494207) q[0];
sx q[0];
rz(-2.4038195) q[0];
rz(-2.8470464) q[2];
sx q[2];
rz(-1.6606332) q[2];
sx q[2];
rz(-1.3603269) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.78754567) q[1];
sx q[1];
rz(-0.58632942) q[1];
sx q[1];
rz(1.7802618) q[1];
rz(-pi) q[2];
rz(0.22927852) q[3];
sx q[3];
rz(-1.6212654) q[3];
sx q[3];
rz(-0.63427395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6198373) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(-2.8660529) q[2];
rz(0.11166212) q[3];
sx q[3];
rz(-1.9409981) q[3];
sx q[3];
rz(2.6707941) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6977285) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(2.2648947) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-2.2677939) q[1];
sx q[1];
rz(-0.91526389) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4782151) q[0];
sx q[0];
rz(-1.5236679) q[0];
sx q[0];
rz(-2.1992654) q[0];
rz(-pi) q[1];
rz(1.615633) q[2];
sx q[2];
rz(-2.0955288) q[2];
sx q[2];
rz(-2.5068138) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.96990578) q[1];
sx q[1];
rz(-1.415442) q[1];
sx q[1];
rz(3.1093662) q[1];
rz(2.4004585) q[3];
sx q[3];
rz(-0.94666615) q[3];
sx q[3];
rz(0.65628624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9203732) q[2];
sx q[2];
rz(-2.129014) q[2];
sx q[2];
rz(-0.4635703) q[2];
rz(2.5772337) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(2.213403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0267462) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(2.0625431) q[0];
rz(-0.59533978) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(0.39658305) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30124861) q[0];
sx q[0];
rz(-0.52919555) q[0];
sx q[0];
rz(-2.089414) q[0];
x q[1];
rz(2.5846892) q[2];
sx q[2];
rz(-0.54585005) q[2];
sx q[2];
rz(2.5715373) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0551758) q[1];
sx q[1];
rz(-0.91311087) q[1];
sx q[1];
rz(-3.0228826) q[1];
x q[2];
rz(-0.58432213) q[3];
sx q[3];
rz(-1.7752247) q[3];
sx q[3];
rz(-1.0014597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.98809272) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(1.5552103) q[2];
rz(-1.4554626) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(2.3144408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39847386) q[0];
sx q[0];
rz(-1.9819336) q[0];
sx q[0];
rz(-2.356785) q[0];
rz(1.8709042) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(2.0152337) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2130177) q[0];
sx q[0];
rz(-1.2733766) q[0];
sx q[0];
rz(-1.6638882) q[0];
rz(-pi) q[1];
rz(-1.8225841) q[2];
sx q[2];
rz(-2.6195824) q[2];
sx q[2];
rz(0.28282794) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3583402) q[1];
sx q[1];
rz(-1.1249152) q[1];
sx q[1];
rz(-1.3621484) q[1];
rz(1.8282312) q[3];
sx q[3];
rz(-0.35850393) q[3];
sx q[3];
rz(-3.0646216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.209098) q[2];
sx q[2];
rz(-2.035049) q[2];
sx q[2];
rz(-0.15360019) q[2];
rz(0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(1.37384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.7711733) q[0];
sx q[0];
rz(-2.6035247) q[0];
sx q[0];
rz(0.11288189) q[0];
rz(-1.0007292) q[1];
sx q[1];
rz(-2.395144) q[1];
sx q[1];
rz(3.1088366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5927133) q[0];
sx q[0];
rz(-1.2476876) q[0];
sx q[0];
rz(-0.26922853) q[0];
x q[1];
rz(-2.6997386) q[2];
sx q[2];
rz(-2.0896857) q[2];
sx q[2];
rz(-1.9372802) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36383648) q[1];
sx q[1];
rz(-1.6297852) q[1];
sx q[1];
rz(-1.9883518) q[1];
rz(-pi) q[2];
rz(1.350358) q[3];
sx q[3];
rz(-0.38808295) q[3];
sx q[3];
rz(-2.2409852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(0.62210554) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(-2.7511403) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0416097) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.6149678) q[0];
rz(2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(0.48535767) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1982614) q[0];
sx q[0];
rz(-1.7051464) q[0];
sx q[0];
rz(-1.6249379) q[0];
rz(-pi) q[1];
rz(-2.910026) q[2];
sx q[2];
rz(-2.5307557) q[2];
sx q[2];
rz(-1.4996741) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8977412) q[1];
sx q[1];
rz(-2.284639) q[1];
sx q[1];
rz(0.41390151) q[1];
rz(-pi) q[2];
rz(-0.58166196) q[3];
sx q[3];
rz(-1.186944) q[3];
sx q[3];
rz(1.4130842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(2.9821441) q[2];
rz(-1.4032646) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52206802) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(1.4341226) q[0];
rz(-1.8822949) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(-2.5433345) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3529417) q[0];
sx q[0];
rz(-1.3539679) q[0];
sx q[0];
rz(1.7822595) q[0];
x q[1];
rz(1.422545) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(1.460618) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4926589) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(1.2050864) q[1];
rz(-pi) q[2];
rz(2.9833097) q[3];
sx q[3];
rz(-1.7332819) q[3];
sx q[3];
rz(2.3815001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.8965315) q[2];
sx q[2];
rz(-2.4895978) q[2];
rz(2.5478798) q[3];
sx q[3];
rz(-1.1810602) q[3];
sx q[3];
rz(2.0098856) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.839529) q[0];
sx q[0];
rz(-0.65757127) q[0];
sx q[0];
rz(0.99304637) q[0];
rz(-1.9051753) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-1.4478222) q[2];
sx q[2];
rz(-0.95477827) q[2];
sx q[2];
rz(1.9670602) q[2];
rz(2.7491309) q[3];
sx q[3];
rz(-0.93467181) q[3];
sx q[3];
rz(-1.296464) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
