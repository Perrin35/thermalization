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
rz(0.33557284) q[0];
sx q[0];
rz(4.8290401) q[0];
sx q[0];
rz(7.3168559) q[0];
rz(1.4312862) q[1];
sx q[1];
rz(-0.55748504) q[1];
sx q[1];
rz(-1.0097591) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5752707) q[0];
sx q[0];
rz(-0.22804865) q[0];
sx q[0];
rz(0.98342498) q[0];
rz(-pi) q[1];
rz(0.0037438914) q[2];
sx q[2];
rz(-1.4883248) q[2];
sx q[2];
rz(-2.1090901) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.20045763) q[1];
sx q[1];
rz(-1.4054646) q[1];
sx q[1];
rz(1.453406) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5752154) q[3];
sx q[3];
rz(-1.7130245) q[3];
sx q[3];
rz(-2.9517236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7734163) q[2];
sx q[2];
rz(-2.2645617) q[2];
sx q[2];
rz(0.82484335) q[2];
rz(0.33251897) q[3];
sx q[3];
rz(-0.85086346) q[3];
sx q[3];
rz(-2.6578145) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9124209) q[0];
sx q[0];
rz(-0.084675463) q[0];
sx q[0];
rz(-2.605865) q[0];
rz(-0.76672211) q[1];
sx q[1];
rz(-1.7886536) q[1];
sx q[1];
rz(-1.3923233) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.509393) q[0];
sx q[0];
rz(-1.107111) q[0];
sx q[0];
rz(2.7576588) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5611982) q[2];
sx q[2];
rz(-1.9119153) q[2];
sx q[2];
rz(-2.7500602) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1050807) q[1];
sx q[1];
rz(-2.3859508) q[1];
sx q[1];
rz(-1.5393658) q[1];
rz(-pi) q[2];
rz(-1.4658607) q[3];
sx q[3];
rz(-1.4426878) q[3];
sx q[3];
rz(-0.48612938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9048189) q[2];
sx q[2];
rz(-0.83319131) q[2];
sx q[2];
rz(1.0986249) q[2];
rz(2.3084579) q[3];
sx q[3];
rz(-1.5435217) q[3];
sx q[3];
rz(-2.0875077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58040923) q[0];
sx q[0];
rz(-1.9920749) q[0];
sx q[0];
rz(0.68823254) q[0];
rz(1.8904103) q[1];
sx q[1];
rz(-0.65443188) q[1];
sx q[1];
rz(-1.2122663) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.040846856) q[0];
sx q[0];
rz(-2.0424574) q[0];
sx q[0];
rz(2.3673886) q[0];
rz(-pi) q[1];
rz(1.8463676) q[2];
sx q[2];
rz(-1.3418806) q[2];
sx q[2];
rz(2.6273397) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7292068) q[1];
sx q[1];
rz(-0.6589891) q[1];
sx q[1];
rz(-3.1324196) q[1];
x q[2];
rz(-3.1367407) q[3];
sx q[3];
rz(-1.650164) q[3];
sx q[3];
rz(-2.8341189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7324098) q[2];
sx q[2];
rz(-0.99486351) q[2];
sx q[2];
rz(2.691972) q[2];
rz(-0.85363394) q[3];
sx q[3];
rz(-2.4946404) q[3];
sx q[3];
rz(0.70931119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.7483826) q[0];
sx q[0];
rz(-2.125183) q[0];
sx q[0];
rz(-2.7384695) q[0];
rz(0.77278167) q[1];
sx q[1];
rz(-1.9085596) q[1];
sx q[1];
rz(0.82537878) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7220424) q[0];
sx q[0];
rz(-0.78762509) q[0];
sx q[0];
rz(1.1088637) q[0];
rz(-pi) q[1];
rz(-1.5947444) q[2];
sx q[2];
rz(-1.5976816) q[2];
sx q[2];
rz(-1.250911) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5063012) q[1];
sx q[1];
rz(-1.4951733) q[1];
sx q[1];
rz(0.98905321) q[1];
rz(-0.89572923) q[3];
sx q[3];
rz(-0.78215996) q[3];
sx q[3];
rz(-1.6767927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0394502) q[2];
sx q[2];
rz(-1.8723698) q[2];
sx q[2];
rz(-2.6329182) q[2];
rz(1.8968808) q[3];
sx q[3];
rz(-2.0679943) q[3];
sx q[3];
rz(1.8083474) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4077069) q[0];
sx q[0];
rz(-1.5579959) q[0];
sx q[0];
rz(-2.7852614) q[0];
rz(1.7286667) q[1];
sx q[1];
rz(-1.8316385) q[1];
sx q[1];
rz(1.8467356) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0644449) q[0];
sx q[0];
rz(-2.003621) q[0];
sx q[0];
rz(1.3291784) q[0];
rz(2.9239976) q[2];
sx q[2];
rz(-0.67309531) q[2];
sx q[2];
rz(2.4270647) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19902557) q[1];
sx q[1];
rz(-2.6771705) q[1];
sx q[1];
rz(-0.53408043) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6925907) q[3];
sx q[3];
rz(-1.0706361) q[3];
sx q[3];
rz(1.7412108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.20092189) q[2];
sx q[2];
rz(-1.385043) q[2];
sx q[2];
rz(-2.5015639) q[2];
rz(2.707543) q[3];
sx q[3];
rz(-0.95165747) q[3];
sx q[3];
rz(1.1205589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5956748) q[0];
sx q[0];
rz(-0.43788236) q[0];
sx q[0];
rz(2.3906999) q[0];
rz(2.3789876) q[1];
sx q[1];
rz(-1.4354939) q[1];
sx q[1];
rz(0.5336175) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5845244) q[0];
sx q[0];
rz(-1.8584588) q[0];
sx q[0];
rz(0.30867048) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3377588) q[2];
sx q[2];
rz(-1.2495923) q[2];
sx q[2];
rz(-2.8558232) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.92317947) q[1];
sx q[1];
rz(-2.235102) q[1];
sx q[1];
rz(-2.9914145) q[1];
x q[2];
rz(-0.814416) q[3];
sx q[3];
rz(-0.73668639) q[3];
sx q[3];
rz(-1.8131922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.077347191) q[2];
sx q[2];
rz(-0.40234819) q[2];
sx q[2];
rz(0.8858436) q[2];
rz(1.3028076) q[3];
sx q[3];
rz(-1.9582483) q[3];
sx q[3];
rz(-0.47811374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(1.9293905) q[0];
sx q[0];
rz(-0.92340702) q[0];
sx q[0];
rz(0.59342629) q[0];
rz(1.628283) q[1];
sx q[1];
rz(-1.1791041) q[1];
sx q[1];
rz(2.5679307) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57034111) q[0];
sx q[0];
rz(-1.9431337) q[0];
sx q[0];
rz(2.6254197) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5638177) q[2];
sx q[2];
rz(-1.7794357) q[2];
sx q[2];
rz(-1.8633757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.57854125) q[1];
sx q[1];
rz(-2.6669901) q[1];
sx q[1];
rz(2.5792349) q[1];
x q[2];
rz(0.078646831) q[3];
sx q[3];
rz(-1.6212731) q[3];
sx q[3];
rz(2.5223062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4297428) q[2];
sx q[2];
rz(-1.9387127) q[2];
sx q[2];
rz(1.708606) q[2];
rz(1.1687219) q[3];
sx q[3];
rz(-0.43787268) q[3];
sx q[3];
rz(-1.5168813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2699921) q[0];
sx q[0];
rz(-0.90684909) q[0];
sx q[0];
rz(-0.20371833) q[0];
rz(1.8261955) q[1];
sx q[1];
rz(-1.8351277) q[1];
sx q[1];
rz(1.3320097) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7527134) q[0];
sx q[0];
rz(-2.1915771) q[0];
sx q[0];
rz(0.78691532) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2526251) q[2];
sx q[2];
rz(-0.79304129) q[2];
sx q[2];
rz(-2.5630132) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.2329008) q[1];
sx q[1];
rz(-2.3756643) q[1];
sx q[1];
rz(0.58118622) q[1];
rz(-2.2409954) q[3];
sx q[3];
rz(-2.934747) q[3];
sx q[3];
rz(-1.2634009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.77067644) q[2];
sx q[2];
rz(-3.0597661) q[2];
sx q[2];
rz(-1.7601684) q[2];
rz(2.5904739) q[3];
sx q[3];
rz(-1.8128606) q[3];
sx q[3];
rz(-0.74530017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7119174) q[0];
sx q[0];
rz(-1.4497919) q[0];
sx q[0];
rz(1.2592738) q[0];
rz(-1.7970386) q[1];
sx q[1];
rz(-0.63251907) q[1];
sx q[1];
rz(1.9154027) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0622985) q[0];
sx q[0];
rz(-0.31656893) q[0];
sx q[0];
rz(-2.0189254) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0001423) q[2];
sx q[2];
rz(-1.9252021) q[2];
sx q[2];
rz(0.62858519) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.256168) q[1];
sx q[1];
rz(-0.82535911) q[1];
sx q[1];
rz(-1.481545) q[1];
x q[2];
rz(-2.8060447) q[3];
sx q[3];
rz(-2.0651544) q[3];
sx q[3];
rz(-1.7691607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5598477) q[2];
sx q[2];
rz(-0.86220828) q[2];
sx q[2];
rz(1.5398514) q[2];
rz(1.3567443) q[3];
sx q[3];
rz(-0.53250766) q[3];
sx q[3];
rz(-1.2342359) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27720472) q[0];
sx q[0];
rz(-1.5437523) q[0];
sx q[0];
rz(-3.0371015) q[0];
rz(-1.3970207) q[1];
sx q[1];
rz(-1.3518159) q[1];
sx q[1];
rz(1.6207961) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45986555) q[0];
sx q[0];
rz(-2.2978373) q[0];
sx q[0];
rz(-1.8521502) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7303329) q[2];
sx q[2];
rz(-1.9070093) q[2];
sx q[2];
rz(-1.6415063) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4222153) q[1];
sx q[1];
rz(-1.6746192) q[1];
sx q[1];
rz(2.8524998) q[1];
x q[2];
rz(0.057985882) q[3];
sx q[3];
rz(-2.1515905) q[3];
sx q[3];
rz(-0.58422663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4349159) q[2];
sx q[2];
rz(-0.84182635) q[2];
sx q[2];
rz(-1.8527156) q[2];
rz(1.0051109) q[3];
sx q[3];
rz(-1.0594599) q[3];
sx q[3];
rz(3.03481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2224251) q[0];
sx q[0];
rz(-1.5821624) q[0];
sx q[0];
rz(-0.17056175) q[0];
rz(-0.87986058) q[1];
sx q[1];
rz(-0.51645551) q[1];
sx q[1];
rz(-2.5194306) q[1];
rz(-2.6420319) q[2];
sx q[2];
rz(-1.5141664) q[2];
sx q[2];
rz(0.92922633) q[2];
rz(-2.5414657) q[3];
sx q[3];
rz(-1.5059581) q[3];
sx q[3];
rz(-3.0653421) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
