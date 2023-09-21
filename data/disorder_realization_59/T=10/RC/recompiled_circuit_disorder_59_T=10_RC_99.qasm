OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3988848) q[0];
sx q[0];
rz(-2.3595915) q[0];
sx q[0];
rz(-1.8703823) q[0];
rz(3.4186163) q[1];
sx q[1];
rz(3.613598) q[1];
sx q[1];
rz(9.4233905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1479552) q[0];
sx q[0];
rz(-1.7339098) q[0];
sx q[0];
rz(-1.1975343) q[0];
x q[1];
rz(-0.087287993) q[2];
sx q[2];
rz(-2.6929571) q[2];
sx q[2];
rz(-2.0729614) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.063349799) q[1];
sx q[1];
rz(-1.6410876) q[1];
sx q[1];
rz(2.6438144) q[1];
rz(-pi) q[2];
rz(1.9418342) q[3];
sx q[3];
rz(-1.3317809) q[3];
sx q[3];
rz(-1.0306851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15443054) q[2];
sx q[2];
rz(-0.61750948) q[2];
sx q[2];
rz(-0.74938613) q[2];
rz(2.1253712) q[3];
sx q[3];
rz(-1.1775492) q[3];
sx q[3];
rz(-2.7367676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7063023) q[0];
sx q[0];
rz(-0.82536936) q[0];
sx q[0];
rz(-2.170927) q[0];
rz(2.1043815) q[1];
sx q[1];
rz(-1.7036006) q[1];
sx q[1];
rz(0.81545365) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52112752) q[0];
sx q[0];
rz(-1.5107811) q[0];
sx q[0];
rz(-0.64882664) q[0];
rz(0.72950659) q[2];
sx q[2];
rz(-1.365005) q[2];
sx q[2];
rz(-1.4480928) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.38824575) q[1];
sx q[1];
rz(-1.8716295) q[1];
sx q[1];
rz(1.6080329) q[1];
x q[2];
rz(-0.29173298) q[3];
sx q[3];
rz(-2.4749304) q[3];
sx q[3];
rz(-3.047903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6796391) q[2];
sx q[2];
rz(-1.5755499) q[2];
sx q[2];
rz(-2.5088076) q[2];
rz(-1.9880382) q[3];
sx q[3];
rz(-0.76806918) q[3];
sx q[3];
rz(2.8320584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(0.84045029) q[0];
sx q[0];
rz(-1.9294894) q[0];
sx q[0];
rz(0.87483037) q[0];
rz(-1.3300928) q[1];
sx q[1];
rz(-1.4346088) q[1];
sx q[1];
rz(2.1420746) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8130428) q[0];
sx q[0];
rz(-1.8694436) q[0];
sx q[0];
rz(2.9850328) q[0];
rz(-0.91018422) q[2];
sx q[2];
rz(-2.0816457) q[2];
sx q[2];
rz(2.3597033) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7175908) q[1];
sx q[1];
rz(-2.058299) q[1];
sx q[1];
rz(1.8600149) q[1];
rz(-pi) q[2];
rz(-1.0619034) q[3];
sx q[3];
rz(-1.5012) q[3];
sx q[3];
rz(0.059046179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.53753608) q[2];
sx q[2];
rz(-0.92210046) q[2];
sx q[2];
rz(-2.5615454) q[2];
rz(0.81702685) q[3];
sx q[3];
rz(-1.7592808) q[3];
sx q[3];
rz(1.9918611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76628768) q[0];
sx q[0];
rz(-1.5505318) q[0];
sx q[0];
rz(-0.91039175) q[0];
rz(-2.6903649) q[1];
sx q[1];
rz(-1.5463566) q[1];
sx q[1];
rz(-2.8667563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66514689) q[0];
sx q[0];
rz(-3.0229212) q[0];
sx q[0];
rz(-0.84400405) q[0];
rz(1.574013) q[2];
sx q[2];
rz(-2.6811757) q[2];
sx q[2];
rz(0.38052961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6064925) q[1];
sx q[1];
rz(-0.98787687) q[1];
sx q[1];
rz(0.54775723) q[1];
x q[2];
rz(1.4196017) q[3];
sx q[3];
rz(-0.70221838) q[3];
sx q[3];
rz(-0.37213009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.794902) q[2];
sx q[2];
rz(-1.1179504) q[2];
sx q[2];
rz(1.6332731) q[2];
rz(1.1446965) q[3];
sx q[3];
rz(-2.4016524) q[3];
sx q[3];
rz(0.16170734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4836924) q[0];
sx q[0];
rz(-1.9265441) q[0];
sx q[0];
rz(2.143798) q[0];
rz(-0.18355852) q[1];
sx q[1];
rz(-1.4869556) q[1];
sx q[1];
rz(-1.516974) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.957513) q[0];
sx q[0];
rz(-2.3947869) q[0];
sx q[0];
rz(1.0190796) q[0];
rz(2.5324608) q[2];
sx q[2];
rz(-0.62437781) q[2];
sx q[2];
rz(2.6513211) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.5114054) q[1];
sx q[1];
rz(-1.1668219) q[1];
sx q[1];
rz(1.1955111) q[1];
rz(-2.3715641) q[3];
sx q[3];
rz(-0.75776811) q[3];
sx q[3];
rz(-3.1185574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3395485) q[2];
sx q[2];
rz(-0.99207726) q[2];
sx q[2];
rz(-0.3240164) q[2];
rz(-1.8185395) q[3];
sx q[3];
rz(-2.3855305) q[3];
sx q[3];
rz(-1.5312622) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3762387) q[0];
sx q[0];
rz(-2.1026251) q[0];
sx q[0];
rz(-1.2639686) q[0];
rz(0.91066796) q[1];
sx q[1];
rz(-1.940454) q[1];
sx q[1];
rz(2.8009159) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9522889) q[0];
sx q[0];
rz(-1.5049107) q[0];
sx q[0];
rz(0.031866372) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8315115) q[2];
sx q[2];
rz(-0.78044621) q[2];
sx q[2];
rz(2.9030637) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.50960474) q[1];
sx q[1];
rz(-2.0443516) q[1];
sx q[1];
rz(-1.8297086) q[1];
x q[2];
rz(2.2716899) q[3];
sx q[3];
rz(-0.60855908) q[3];
sx q[3];
rz(0.22939798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1910151) q[2];
sx q[2];
rz(-0.58379972) q[2];
sx q[2];
rz(2.3699956) q[2];
rz(2.5937882) q[3];
sx q[3];
rz(-2.1717725) q[3];
sx q[3];
rz(-0.56345338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9724378) q[0];
sx q[0];
rz(-1.6945524) q[0];
sx q[0];
rz(-2.7959438) q[0];
rz(-0.06282839) q[1];
sx q[1];
rz(-2.6627916) q[1];
sx q[1];
rz(-2.6766434) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3195254) q[0];
sx q[0];
rz(-1.7626581) q[0];
sx q[0];
rz(1.096154) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2642235) q[2];
sx q[2];
rz(-0.69958985) q[2];
sx q[2];
rz(0.96044651) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3685776) q[1];
sx q[1];
rz(-1.7891208) q[1];
sx q[1];
rz(1.6792084) q[1];
x q[2];
rz(0.025860272) q[3];
sx q[3];
rz(-1.6169294) q[3];
sx q[3];
rz(0.24054724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0039625) q[2];
sx q[2];
rz(-1.5045065) q[2];
sx q[2];
rz(-0.31759343) q[2];
rz(2.5701304) q[3];
sx q[3];
rz(-2.1025889) q[3];
sx q[3];
rz(2.8542744) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6222318) q[0];
sx q[0];
rz(-1.2943635) q[0];
sx q[0];
rz(2.8572594) q[0];
rz(0.55150664) q[1];
sx q[1];
rz(-0.14177828) q[1];
sx q[1];
rz(0.078358738) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4975472) q[0];
sx q[0];
rz(-0.80701485) q[0];
sx q[0];
rz(-0.09597309) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9044754) q[2];
sx q[2];
rz(-1.8097005) q[2];
sx q[2];
rz(0.6616) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6137177) q[1];
sx q[1];
rz(-2.4637239) q[1];
sx q[1];
rz(-1.2916958) q[1];
rz(-pi) q[2];
rz(-3.1073242) q[3];
sx q[3];
rz(-1.4141603) q[3];
sx q[3];
rz(0.53965118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4006965) q[2];
sx q[2];
rz(-0.90831465) q[2];
sx q[2];
rz(0.25137869) q[2];
rz(-2.5583983) q[3];
sx q[3];
rz(-2.0299032) q[3];
sx q[3];
rz(-1.73197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79779977) q[0];
sx q[0];
rz(-0.082158953) q[0];
sx q[0];
rz(-3.0902241) q[0];
rz(-0.92357606) q[1];
sx q[1];
rz(-2.4802465) q[1];
sx q[1];
rz(2.267568) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6459991) q[0];
sx q[0];
rz(-2.3388303) q[0];
sx q[0];
rz(-1.5828703) q[0];
x q[1];
rz(-1.9782412) q[2];
sx q[2];
rz(-1.904084) q[2];
sx q[2];
rz(-2.8826706) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8585426) q[1];
sx q[1];
rz(-2.0701323) q[1];
sx q[1];
rz(0.081375558) q[1];
x q[2];
rz(2.5328818) q[3];
sx q[3];
rz(-1.2864283) q[3];
sx q[3];
rz(1.9539208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.41436568) q[2];
sx q[2];
rz(-0.7545158) q[2];
sx q[2];
rz(2.5218463) q[2];
rz(1.184458) q[3];
sx q[3];
rz(-1.8871566) q[3];
sx q[3];
rz(1.363389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0062362) q[0];
sx q[0];
rz(-1.0422491) q[0];
sx q[0];
rz(-2.4172879) q[0];
rz(-2.9528217) q[1];
sx q[1];
rz(-0.17938463) q[1];
sx q[1];
rz(1.1788517) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41335426) q[0];
sx q[0];
rz(-1.2984707) q[0];
sx q[0];
rz(-0.21270919) q[0];
x q[1];
rz(-1.4775425) q[2];
sx q[2];
rz(-1.4443195) q[2];
sx q[2];
rz(1.7699514) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9465543) q[1];
sx q[1];
rz(-0.95609162) q[1];
sx q[1];
rz(-0.22919319) q[1];
rz(-pi) q[2];
rz(0.3801109) q[3];
sx q[3];
rz(-0.5553402) q[3];
sx q[3];
rz(-1.7172608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.05802352) q[2];
sx q[2];
rz(-2.100779) q[2];
sx q[2];
rz(0.89938346) q[2];
rz(-0.87456885) q[3];
sx q[3];
rz(-2.7159297) q[3];
sx q[3];
rz(1.4609059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(-2.5205004) q[0];
sx q[0];
rz(-2.7324471) q[0];
sx q[0];
rz(0.24656217) q[0];
rz(-0.75795603) q[1];
sx q[1];
rz(-1.4823722) q[1];
sx q[1];
rz(-1.4588251) q[1];
rz(-1.1307217) q[2];
sx q[2];
rz(-1.6025087) q[2];
sx q[2];
rz(-0.58379731) q[2];
rz(-0.48537985) q[3];
sx q[3];
rz(-2.831922) q[3];
sx q[3];
rz(0.61556863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
