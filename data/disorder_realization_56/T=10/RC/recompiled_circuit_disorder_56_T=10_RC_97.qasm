OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874415) q[0];
sx q[0];
rz(3.677877) q[0];
sx q[0];
rz(10.372547) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(-1.0277494) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0004814) q[0];
sx q[0];
rz(-1.9465465) q[0];
sx q[0];
rz(1.5953338) q[0];
rz(-pi) q[1];
rz(2.920354) q[2];
sx q[2];
rz(-2.049963) q[2];
sx q[2];
rz(2.0146807) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7386908) q[1];
sx q[1];
rz(-1.9539023) q[1];
sx q[1];
rz(-0.72523592) q[1];
rz(-pi) q[2];
rz(-2.0184228) q[3];
sx q[3];
rz(-2.1493704) q[3];
sx q[3];
rz(2.6920126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.4346069) q[2];
sx q[2];
rz(-2.091308) q[2];
rz(1.1132647) q[3];
sx q[3];
rz(-0.89171019) q[3];
sx q[3];
rz(-3.0734857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0691836) q[0];
sx q[0];
rz(-1.8962815) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(0.61966664) q[1];
sx q[1];
rz(-1.0071808) q[1];
sx q[1];
rz(-2.0334977) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8544234) q[0];
sx q[0];
rz(-2.4872635) q[0];
sx q[0];
rz(-1.9186583) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1531395) q[2];
sx q[2];
rz(-2.4465912) q[2];
sx q[2];
rz(2.7671791) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.88156466) q[1];
sx q[1];
rz(-1.2177055) q[1];
sx q[1];
rz(2.7295223) q[1];
rz(-1.3012582) q[3];
sx q[3];
rz(-0.88629913) q[3];
sx q[3];
rz(-1.0052296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(-2.1662946) q[2];
rz(-2.2235928) q[3];
sx q[3];
rz(-1.2851597) q[3];
sx q[3];
rz(-2.8454034) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96238962) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(-0.54291022) q[0];
rz(0.88223282) q[1];
sx q[1];
rz(-1.135332) q[1];
sx q[1];
rz(-2.1767445) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1436413) q[0];
sx q[0];
rz(-0.43295857) q[0];
sx q[0];
rz(2.4233682) q[0];
rz(-pi) q[1];
rz(-2.8446571) q[2];
sx q[2];
rz(-1.529971) q[2];
sx q[2];
rz(0.21619913) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1988586) q[1];
sx q[1];
rz(-3.004289) q[1];
sx q[1];
rz(-1.0820461) q[1];
rz(1.996702) q[3];
sx q[3];
rz(-2.689038) q[3];
sx q[3];
rz(-1.2061314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0470011) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(2.0084521) q[2];
rz(2.9099693) q[3];
sx q[3];
rz(-1.8685721) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27947458) q[0];
sx q[0];
rz(-3.1311488) q[0];
sx q[0];
rz(-1.3765155) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.8771749) q[1];
sx q[1];
rz(2.8994697) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0771368) q[0];
sx q[0];
rz(-2.4479439) q[0];
sx q[0];
rz(-0.25607381) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5838727) q[2];
sx q[2];
rz(-1.1411975) q[2];
sx q[2];
rz(-1.0193046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8230069) q[1];
sx q[1];
rz(-1.7404403) q[1];
sx q[1];
rz(0.45348788) q[1];
x q[2];
rz(1.4297156) q[3];
sx q[3];
rz(-1.9851079) q[3];
sx q[3];
rz(2.8431818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.018192856) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(3.0467395) q[2];
rz(1.3421966) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(2.7024787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78385329) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(2.7807996) q[0];
rz(1.7533253) q[1];
sx q[1];
rz(-1.3307945) q[1];
sx q[1];
rz(-2.0070019) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73571262) q[0];
sx q[0];
rz(-0.94731936) q[0];
sx q[0];
rz(-2.1924125) q[0];
x q[1];
rz(1.1410494) q[2];
sx q[2];
rz(-1.2417214) q[2];
sx q[2];
rz(0.82440257) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1281631) q[1];
sx q[1];
rz(-2.014782) q[1];
sx q[1];
rz(-1.7432937) q[1];
rz(-pi) q[2];
rz(1.6547104) q[3];
sx q[3];
rz(-2.0062371) q[3];
sx q[3];
rz(0.30121379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0052884) q[2];
sx q[2];
rz(-1.4207999) q[2];
sx q[2];
rz(-2.5689382) q[2];
rz(-2.2128361) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(-1.9740392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3271493) q[0];
sx q[0];
rz(-1.0961908) q[0];
sx q[0];
rz(-1.8922528) q[0];
rz(-1.2231187) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(1.9893601) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0056861176) q[0];
sx q[0];
rz(-1.583477) q[0];
sx q[0];
rz(-1.6760875) q[0];
rz(-0.50243369) q[2];
sx q[2];
rz(-1.9366493) q[2];
sx q[2];
rz(2.1937214) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8730948) q[1];
sx q[1];
rz(-2.9151306) q[1];
sx q[1];
rz(-0.48392673) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1762987) q[3];
sx q[3];
rz(-0.98635841) q[3];
sx q[3];
rz(2.8402929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46889177) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(2.1506298) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-0.9655374) q[3];
sx q[3];
rz(-0.34753862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-2.8836442) q[0];
sx q[0];
rz(-0.22709665) q[0];
sx q[0];
rz(-3.0793072) q[0];
rz(2.9557872) q[1];
sx q[1];
rz(-1.4567016) q[1];
sx q[1];
rz(-0.3947765) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5236854) q[0];
sx q[0];
rz(-0.46645188) q[0];
sx q[0];
rz(-2.594069) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4418418) q[2];
sx q[2];
rz(-0.47669461) q[2];
sx q[2];
rz(1.2303908) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9638622) q[1];
sx q[1];
rz(-1.9964255) q[1];
sx q[1];
rz(-1.8427986) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9881334) q[3];
sx q[3];
rz(-1.6073265) q[3];
sx q[3];
rz(0.032144459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4930967) q[2];
sx q[2];
rz(-2.811537) q[2];
sx q[2];
rz(0.27302343) q[2];
rz(-1.8388883) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(0.31204143) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39772314) q[0];
sx q[0];
rz(-1.1345154) q[0];
sx q[0];
rz(1.2851108) q[0];
rz(-1.6400736) q[1];
sx q[1];
rz(-1.39095) q[1];
sx q[1];
rz(-1.8008908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9467981) q[0];
sx q[0];
rz(-1.9614944) q[0];
sx q[0];
rz(2.7892053) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69581823) q[2];
sx q[2];
rz(-1.5317481) q[2];
sx q[2];
rz(2.2829636) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50780523) q[1];
sx q[1];
rz(-1.9314737) q[1];
sx q[1];
rz(-0.37640576) q[1];
x q[2];
rz(-2.826346) q[3];
sx q[3];
rz(-1.1479706) q[3];
sx q[3];
rz(1.8002312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(-1.0236053) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-0.39026323) q[3];
sx q[3];
rz(-0.24188724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-0.99656314) q[0];
sx q[0];
rz(-1.6212844) q[0];
rz(2.8114491) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(-2.3044589) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5569816) q[0];
sx q[0];
rz(-0.71174445) q[0];
sx q[0];
rz(-0.31194375) q[0];
rz(-pi) q[1];
x q[1];
rz(0.43948549) q[2];
sx q[2];
rz(-2.5186933) q[2];
sx q[2];
rz(-2.5650131) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7652119) q[1];
sx q[1];
rz(-1.4188671) q[1];
sx q[1];
rz(-3.0357749) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0765452) q[3];
sx q[3];
rz(-2.5878083) q[3];
sx q[3];
rz(-1.6365285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5403486) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(0.17962757) q[2];
rz(-2.1458697) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(-1.3109591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3778465) q[0];
sx q[0];
rz(-0.34559956) q[0];
sx q[0];
rz(-2.0843704) q[0];
rz(-3.0341042) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(0.9799788) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7769023) q[0];
sx q[0];
rz(-1.6216462) q[0];
sx q[0];
rz(1.6965673) q[0];
x q[1];
rz(-0.51771848) q[2];
sx q[2];
rz(-2.0132408) q[2];
sx q[2];
rz(2.8874318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.53303888) q[1];
sx q[1];
rz(-1.8095008) q[1];
sx q[1];
rz(-2.0134316) q[1];
rz(-pi) q[2];
rz(-3.0409056) q[3];
sx q[3];
rz(-1.2486613) q[3];
sx q[3];
rz(1.611657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8367299) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(1.0277964) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.8324865) q[3];
sx q[3];
rz(-2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2733611) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(2.1297395) q[2];
sx q[2];
rz(-1.4197299) q[2];
sx q[2];
rz(1.1001669) q[2];
rz(1.1888614) q[3];
sx q[3];
rz(-1.100913) q[3];
sx q[3];
rz(2.037896) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
