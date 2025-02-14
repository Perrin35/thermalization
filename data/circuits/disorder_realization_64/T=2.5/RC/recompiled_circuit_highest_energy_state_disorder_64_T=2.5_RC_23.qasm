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
rz(3.1177899) q[0];
sx q[0];
rz(-1.1060214) q[0];
sx q[0];
rz(-0.7769146) q[0];
rz(1.39224) q[1];
sx q[1];
rz(-1.3148146) q[1];
sx q[1];
rz(-0.97631747) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1084749) q[0];
sx q[0];
rz(-2.7078848) q[0];
sx q[0];
rz(1.2483622) q[0];
x q[1];
rz(0.56371477) q[2];
sx q[2];
rz(-1.6275121) q[2];
sx q[2];
rz(-2.4617755) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3562505) q[1];
sx q[1];
rz(-1.9948729) q[1];
sx q[1];
rz(0.29845684) q[1];
x q[2];
rz(1.5209274) q[3];
sx q[3];
rz(-2.5679776) q[3];
sx q[3];
rz(0.22465868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6355847) q[2];
sx q[2];
rz(-3.0878461) q[2];
sx q[2];
rz(-2.7347943) q[2];
rz(-0.16945101) q[3];
sx q[3];
rz(-2.6120766) q[3];
sx q[3];
rz(1.0725526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88103831) q[0];
sx q[0];
rz(-0.23242234) q[0];
sx q[0];
rz(-3.1378003) q[0];
rz(-0.077839851) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(-0.30581623) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3681889) q[0];
sx q[0];
rz(-2.0819252) q[0];
sx q[0];
rz(-2.9298733) q[0];
rz(2.9980837) q[2];
sx q[2];
rz(-2.3845551) q[2];
sx q[2];
rz(0.17998634) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.21286035) q[1];
sx q[1];
rz(-1.1747735) q[1];
sx q[1];
rz(-0.81800445) q[1];
x q[2];
rz(-2.3313794) q[3];
sx q[3];
rz(-0.088220291) q[3];
sx q[3];
rz(-0.42217964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1514312) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(2.4988417) q[2];
rz(3.0691872) q[3];
sx q[3];
rz(-1.0591155) q[3];
sx q[3];
rz(1.7806627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.088242315) q[0];
sx q[0];
rz(-2.998816) q[0];
sx q[0];
rz(2.9852168) q[0];
rz(-0.0414255) q[1];
sx q[1];
rz(-0.62774575) q[1];
sx q[1];
rz(-1.5511537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5807626) q[0];
sx q[0];
rz(-1.7770045) q[0];
sx q[0];
rz(2.5872562) q[0];
rz(-pi) q[1];
rz(1.5029491) q[2];
sx q[2];
rz(-2.2257651) q[2];
sx q[2];
rz(0.52815765) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52122766) q[1];
sx q[1];
rz(-1.4583734) q[1];
sx q[1];
rz(-2.6112154) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87167344) q[3];
sx q[3];
rz(-2.1101885) q[3];
sx q[3];
rz(-1.9339069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7330043) q[2];
sx q[2];
rz(-2.6027347) q[2];
sx q[2];
rz(-0.020922529) q[2];
rz(2.9544592) q[3];
sx q[3];
rz(-2.9360866) q[3];
sx q[3];
rz(3.0278897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1784096) q[0];
sx q[0];
rz(-2.8091176) q[0];
sx q[0];
rz(-2.7030429) q[0];
rz(1.6167538) q[1];
sx q[1];
rz(-2.8068145) q[1];
sx q[1];
rz(-0.24510342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.995718) q[0];
sx q[0];
rz(-0.85332131) q[0];
sx q[0];
rz(-3.0307253) q[0];
rz(-pi) q[1];
rz(-1.7282053) q[2];
sx q[2];
rz(-1.9345967) q[2];
sx q[2];
rz(2.8828893) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0673128) q[1];
sx q[1];
rz(-1.7099705) q[1];
sx q[1];
rz(-2.1140631) q[1];
rz(-pi) q[2];
rz(-0.33632261) q[3];
sx q[3];
rz(-2.3592279) q[3];
sx q[3];
rz(-0.7366283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(2.8098246) q[2];
rz(2.6541384) q[3];
sx q[3];
rz(-2.09477) q[3];
sx q[3];
rz(-2.148518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496534) q[0];
sx q[0];
rz(-1.4537469) q[0];
sx q[0];
rz(-2.3680903) q[0];
rz(-1.1812814) q[1];
sx q[1];
rz(-2.9995194) q[1];
sx q[1];
rz(-1.389651) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4163602) q[0];
sx q[0];
rz(-1.2073887) q[0];
sx q[0];
rz(-2.711722) q[0];
rz(-pi) q[1];
rz(-2.2821065) q[2];
sx q[2];
rz(-2.0469249) q[2];
sx q[2];
rz(-1.7993594) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.70848318) q[1];
sx q[1];
rz(-0.481284) q[1];
sx q[1];
rz(-2.9631056) q[1];
rz(-2.4944011) q[3];
sx q[3];
rz(-0.29509896) q[3];
sx q[3];
rz(1.6012675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2191849) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(0.47214559) q[2];
rz(-1.2989429) q[3];
sx q[3];
rz(-1.8483714) q[3];
sx q[3];
rz(2.3310272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.9206813) q[0];
sx q[0];
rz(-2.8784316) q[0];
sx q[0];
rz(0.26350185) q[0];
rz(2.0384516) q[1];
sx q[1];
rz(-1.3145072) q[1];
sx q[1];
rz(-0.37364328) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220396) q[0];
sx q[0];
rz(-3.0470938) q[0];
sx q[0];
rz(-2.2631133) q[0];
rz(-0.12219001) q[2];
sx q[2];
rz(-0.8475248) q[2];
sx q[2];
rz(-3.0360589) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3195575) q[1];
sx q[1];
rz(-0.6614092) q[1];
sx q[1];
rz(-1.0378077) q[1];
rz(-pi) q[2];
rz(-2.5534036) q[3];
sx q[3];
rz(-1.08687) q[3];
sx q[3];
rz(-1.5671135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.11334795) q[2];
sx q[2];
rz(-0.17013203) q[2];
sx q[2];
rz(-2.587758) q[2];
rz(1.3977741) q[3];
sx q[3];
rz(-2.5388986) q[3];
sx q[3];
rz(-2.8288614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5572307) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(-1.9118017) q[0];
rz(2.9025485) q[1];
sx q[1];
rz(-1.5115279) q[1];
sx q[1];
rz(-0.30034932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8883483) q[0];
sx q[0];
rz(-1.7333366) q[0];
sx q[0];
rz(-1.6881315) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4392468) q[2];
sx q[2];
rz(-1.9611729) q[2];
sx q[2];
rz(-0.74884383) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.21619851) q[1];
sx q[1];
rz(-1.5860737) q[1];
sx q[1];
rz(3.141204) q[1];
rz(-pi) q[2];
rz(-0.024552931) q[3];
sx q[3];
rz(-0.85806393) q[3];
sx q[3];
rz(-2.7803382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.131669) q[2];
sx q[2];
rz(-1.4574304) q[2];
sx q[2];
rz(0.24492502) q[2];
rz(0.51982546) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(0.68827099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35172611) q[0];
sx q[0];
rz(-0.34857294) q[0];
sx q[0];
rz(-1.9785731) q[0];
rz(3.0746958) q[1];
sx q[1];
rz(-1.4935378) q[1];
sx q[1];
rz(-1.012872) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71105768) q[0];
sx q[0];
rz(-1.0037046) q[0];
sx q[0];
rz(2.5398769) q[0];
rz(-2.9628721) q[2];
sx q[2];
rz(-2.131049) q[2];
sx q[2];
rz(2.5541039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0433181) q[1];
sx q[1];
rz(-1.9008844) q[1];
sx q[1];
rz(-1.2374452) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4550405) q[3];
sx q[3];
rz(-1.2408537) q[3];
sx q[3];
rz(-0.8984962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0248727) q[2];
sx q[2];
rz(-0.97815424) q[2];
sx q[2];
rz(2.8187974) q[2];
rz(0.60574496) q[3];
sx q[3];
rz(-2.3492458) q[3];
sx q[3];
rz(2.7927223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054319687) q[0];
sx q[0];
rz(-2.9948586) q[0];
sx q[0];
rz(-3.1291381) q[0];
rz(-0.74673486) q[1];
sx q[1];
rz(-0.92339271) q[1];
sx q[1];
rz(2.8616203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50282798) q[0];
sx q[0];
rz(-1.6831213) q[0];
sx q[0];
rz(0.051044271) q[0];
x q[1];
rz(1.1367646) q[2];
sx q[2];
rz(-1.869259) q[2];
sx q[2];
rz(-0.091574319) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1713531) q[1];
sx q[1];
rz(-0.84988028) q[1];
sx q[1];
rz(1.9848787) q[1];
rz(-pi) q[2];
rz(0.98484184) q[3];
sx q[3];
rz(-2.7142314) q[3];
sx q[3];
rz(-2.7387184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3296457) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(-2.9296056) q[2];
rz(0.82344615) q[3];
sx q[3];
rz(-1.4444838) q[3];
sx q[3];
rz(-0.25920355) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14281808) q[0];
sx q[0];
rz(-0.057567216) q[0];
sx q[0];
rz(2.4488191) q[0];
rz(-2.5686) q[1];
sx q[1];
rz(-1.3447821) q[1];
sx q[1];
rz(2.7105892) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8739024) q[0];
sx q[0];
rz(-0.55492102) q[0];
sx q[0];
rz(-3.0303427) q[0];
x q[1];
rz(-2.0153322) q[2];
sx q[2];
rz(-2.1241786) q[2];
sx q[2];
rz(-2.0972507) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1173874) q[1];
sx q[1];
rz(-0.75953249) q[1];
sx q[1];
rz(2.7717436) q[1];
rz(-0.30461664) q[3];
sx q[3];
rz(-1.0416789) q[3];
sx q[3];
rz(1.9068789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4219605) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-0.56023041) q[2];
rz(2.6719921) q[3];
sx q[3];
rz(-0.40265366) q[3];
sx q[3];
rz(0.71389055) q[3];
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
rz(-pi) q[3];
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
rz(-1.8568759) q[0];
sx q[0];
rz(-1.4061883) q[0];
sx q[0];
rz(-1.0549369) q[0];
rz(0.84125413) q[1];
sx q[1];
rz(-1.1019191) q[1];
sx q[1];
rz(-0.046774653) q[1];
rz(2.2711783) q[2];
sx q[2];
rz(-2.6725548) q[2];
sx q[2];
rz(-0.15147333) q[2];
rz(3.1145949) q[3];
sx q[3];
rz(-0.91201966) q[3];
sx q[3];
rz(-2.1569679) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
