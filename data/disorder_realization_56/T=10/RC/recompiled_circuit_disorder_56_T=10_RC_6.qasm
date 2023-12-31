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
rz(-2.6053083) q[0];
sx q[0];
rz(0.94776881) q[0];
rz(1.8127958) q[1];
sx q[1];
rz(-1.2674018) q[1];
sx q[1];
rz(2.1138432) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7209137) q[0];
sx q[0];
rz(-1.5936216) q[0];
sx q[0];
rz(2.7657397) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.920354) q[2];
sx q[2];
rz(-1.0916296) q[2];
sx q[2];
rz(2.0146807) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4029018) q[1];
sx q[1];
rz(-1.9539023) q[1];
sx q[1];
rz(-2.4163567) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0184228) q[3];
sx q[3];
rz(-0.99222224) q[3];
sx q[3];
rz(-2.6920126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0119005) q[2];
sx q[2];
rz(-1.7069858) q[2];
sx q[2];
rz(-1.0502846) q[2];
rz(2.0283279) q[3];
sx q[3];
rz(-2.2498825) q[3];
sx q[3];
rz(0.068107001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.072409078) q[0];
sx q[0];
rz(-1.2453112) q[0];
sx q[0];
rz(-0.29775277) q[0];
rz(-0.61966664) q[1];
sx q[1];
rz(-2.1344118) q[1];
sx q[1];
rz(-2.0334977) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0034804) q[0];
sx q[0];
rz(-1.3618042) q[0];
sx q[0];
rz(-0.94603993) q[0];
rz(-pi) q[1];
x q[1];
rz(0.98845311) q[2];
sx q[2];
rz(-2.4465912) q[2];
sx q[2];
rz(-0.37441355) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.88156466) q[1];
sx q[1];
rz(-1.2177055) q[1];
sx q[1];
rz(2.7295223) q[1];
x q[2];
rz(2.439019) q[3];
sx q[3];
rz(-1.3630023) q[3];
sx q[3];
rz(0.39263615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6790598) q[2];
sx q[2];
rz(-0.97683895) q[2];
sx q[2];
rz(2.1662946) q[2];
rz(0.9179999) q[3];
sx q[3];
rz(-1.8564329) q[3];
sx q[3];
rz(2.8454034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.179203) q[0];
sx q[0];
rz(-0.85150349) q[0];
sx q[0];
rz(-2.5986824) q[0];
rz(2.2593598) q[1];
sx q[1];
rz(-2.0062607) q[1];
sx q[1];
rz(0.96484819) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.097792625) q[0];
sx q[0];
rz(-1.8505197) q[0];
sx q[0];
rz(-0.33491896) q[0];
x q[1];
rz(-1.5281048) q[2];
sx q[2];
rz(-1.8674769) q[2];
sx q[2];
rz(-1.3670849) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.94273401) q[1];
sx q[1];
rz(-0.13730362) q[1];
sx q[1];
rz(2.0595466) q[1];
rz(-2.9433555) q[3];
sx q[3];
rz(-1.9803515) q[3];
sx q[3];
rz(-2.4026681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.09459153) q[2];
sx q[2];
rz(-2.5307405) q[2];
sx q[2];
rz(-2.0084521) q[2];
rz(0.23162332) q[3];
sx q[3];
rz(-1.2730205) q[3];
sx q[3];
rz(2.384322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27947458) q[0];
sx q[0];
rz(-0.010443895) q[0];
sx q[0];
rz(1.3765155) q[0];
rz(-0.51849413) q[1];
sx q[1];
rz(-1.2644178) q[1];
sx q[1];
rz(0.24212295) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69232363) q[0];
sx q[0];
rz(-1.7334492) q[0];
sx q[0];
rz(-0.67740324) q[0];
rz(0.42963117) q[2];
sx q[2];
rz(-1.5826844) q[2];
sx q[2];
rz(-0.55693835) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3185858) q[1];
sx q[1];
rz(-1.7404403) q[1];
sx q[1];
rz(-0.45348788) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8321213) q[3];
sx q[3];
rz(-0.43635338) q[3];
sx q[3];
rz(-2.5040124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1233998) q[2];
sx q[2];
rz(-2.172956) q[2];
sx q[2];
rz(-0.094853178) q[2];
rz(1.799396) q[3];
sx q[3];
rz(-1.3972524) q[3];
sx q[3];
rz(0.43911394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3577394) q[0];
sx q[0];
rz(-1.8308324) q[0];
sx q[0];
rz(-0.36079303) q[0];
rz(1.7533253) q[1];
sx q[1];
rz(-1.8107982) q[1];
sx q[1];
rz(-1.1345908) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2312233) q[0];
sx q[0];
rz(-1.0783505) q[0];
sx q[0];
rz(2.4173196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35933944) q[2];
sx q[2];
rz(-1.1655072) q[2];
sx q[2];
rz(-0.59935024) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1281631) q[1];
sx q[1];
rz(-2.014782) q[1];
sx q[1];
rz(1.7432937) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.43679045) q[3];
sx q[3];
rz(-1.646864) q[3];
sx q[3];
rz(1.9074744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0052884) q[2];
sx q[2];
rz(-1.7207928) q[2];
sx q[2];
rz(-2.5689382) q[2];
rz(-2.2128361) q[3];
sx q[3];
rz(-0.52162617) q[3];
sx q[3];
rz(1.1675534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3271493) q[0];
sx q[0];
rz(-2.0454018) q[0];
sx q[0];
rz(-1.2493398) q[0];
rz(1.918474) q[1];
sx q[1];
rz(-1.5253116) q[1];
sx q[1];
rz(-1.1522326) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.457068) q[0];
sx q[0];
rz(-0.10604924) q[0];
sx q[0];
rz(-1.690879) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50243369) q[2];
sx q[2];
rz(-1.2049434) q[2];
sx q[2];
rz(0.94787129) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.17105477) q[1];
sx q[1];
rz(-1.4661403) q[1];
sx q[1];
rz(-2.9403951) q[1];
rz(-pi) q[2];
x q[2];
rz(0.67752083) q[3];
sx q[3];
rz(-2.065425) q[3];
sx q[3];
rz(0.90466162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6727009) q[2];
sx q[2];
rz(-1.2142618) q[2];
sx q[2];
rz(0.99096283) q[2];
rz(2.4937566) q[3];
sx q[3];
rz(-2.1760553) q[3];
sx q[3];
rz(-2.794054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(2.7468162) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6900401) q[0];
sx q[0];
rz(-1.3344904) q[0];
sx q[0];
rz(0.4060181) q[0];
rz(1.2497181) q[2];
sx q[2];
rz(-1.9294538) q[2];
sx q[2];
rz(1.9888339) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3697606) q[1];
sx q[1];
rz(-2.6409915) q[1];
sx q[1];
rz(2.6066149) q[1];
rz(-pi) q[2];
rz(-1.6077605) q[3];
sx q[3];
rz(-1.4174403) q[3];
sx q[3];
rz(-1.6085898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.64849598) q[2];
sx q[2];
rz(-0.33005565) q[2];
sx q[2];
rz(-0.27302343) q[2];
rz(-1.3027044) q[3];
sx q[3];
rz(-1.3132934) q[3];
sx q[3];
rz(-0.31204143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7438695) q[0];
sx q[0];
rz(-2.0070772) q[0];
sx q[0];
rz(1.8564818) q[0];
rz(1.5015191) q[1];
sx q[1];
rz(-1.7506426) q[1];
sx q[1];
rz(1.8008908) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96269755) q[0];
sx q[0];
rz(-0.51998752) q[0];
sx q[0];
rz(-2.2682701) q[0];
x q[1];
rz(0.69581823) q[2];
sx q[2];
rz(-1.5317481) q[2];
sx q[2];
rz(-0.85862904) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.50780523) q[1];
sx q[1];
rz(-1.2101189) q[1];
sx q[1];
rz(2.7651869) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.826346) q[3];
sx q[3];
rz(-1.1479706) q[3];
sx q[3];
rz(-1.3413615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0354707) q[2];
sx q[2];
rz(-2.3351228) q[2];
sx q[2];
rz(2.1179874) q[2];
rz(-2.9566531) q[3];
sx q[3];
rz(-2.7513294) q[3];
sx q[3];
rz(-2.8997054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0446562) q[0];
sx q[0];
rz(-2.1450295) q[0];
sx q[0];
rz(-1.5203083) q[0];
rz(-0.3301436) q[1];
sx q[1];
rz(-1.2076999) q[1];
sx q[1];
rz(0.83713371) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22568208) q[0];
sx q[0];
rz(-1.7726232) q[0];
sx q[0];
rz(-0.68737824) q[0];
x q[1];
rz(-0.57640055) q[2];
sx q[2];
rz(-1.8216368) q[2];
sx q[2];
rz(0.62945156) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7652119) q[1];
sx q[1];
rz(-1.4188671) q[1];
sx q[1];
rz(-0.10581776) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0749531) q[3];
sx q[3];
rz(-1.3131724) q[3];
sx q[3];
rz(2.6356634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.60124406) q[2];
sx q[2];
rz(-1.0849755) q[2];
sx q[2];
rz(2.9619651) q[2];
rz(0.99572292) q[3];
sx q[3];
rz(-1.2487753) q[3];
sx q[3];
rz(1.8306336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-2.3778465) q[0];
sx q[0];
rz(-2.7959931) q[0];
sx q[0];
rz(-2.0843704) q[0];
rz(-3.0341042) q[1];
sx q[1];
rz(-1.2534393) q[1];
sx q[1];
rz(-2.1616139) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5532593) q[0];
sx q[0];
rz(-0.13561121) q[0];
sx q[0];
rz(1.1853663) q[0];
x q[1];
rz(1.0716295) q[2];
sx q[2];
rz(-1.1071148) q[2];
sx q[2];
rz(-2.0641363) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.500463) q[1];
sx q[1];
rz(-2.6424721) q[1];
sx q[1];
rz(-2.0874546) q[1];
x q[2];
rz(-1.2471334) q[3];
sx q[3];
rz(-1.666288) q[3];
sx q[3];
rz(-3.132706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8367299) q[2];
sx q[2];
rz(-1.9367846) q[2];
sx q[2];
rz(-1.0277964) q[2];
rz(-1.7547539) q[3];
sx q[3];
rz(-1.3091062) q[3];
sx q[3];
rz(2.8579779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86823157) q[0];
sx q[0];
rz(-2.1049451) q[0];
sx q[0];
rz(2.0275397) q[0];
rz(-0.83203075) q[1];
sx q[1];
rz(-2.6770626) q[1];
sx q[1];
rz(-2.4774036) q[1];
rz(0.17765799) q[2];
sx q[2];
rz(-2.1226317) q[2];
sx q[2];
rz(2.5771099) q[2];
rz(1.9527312) q[3];
sx q[3];
rz(-2.0406796) q[3];
sx q[3];
rz(-1.1036967) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
