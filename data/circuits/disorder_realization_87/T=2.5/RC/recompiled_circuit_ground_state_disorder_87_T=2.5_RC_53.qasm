OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7423695) q[0];
sx q[0];
rz(2.7288781) q[0];
sx q[0];
rz(5.4469845) q[0];
rz(-2.3669481) q[1];
sx q[1];
rz(-3.0795842) q[1];
sx q[1];
rz(1.0083415) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3140393) q[0];
sx q[0];
rz(-1.9441105) q[0];
sx q[0];
rz(-3.0568152) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3353592) q[2];
sx q[2];
rz(-1.2004235) q[2];
sx q[2];
rz(0.78200227) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.28532449) q[1];
sx q[1];
rz(-0.50438499) q[1];
sx q[1];
rz(-2.8125202) q[1];
rz(-pi) q[2];
rz(-1.6761682) q[3];
sx q[3];
rz(-2.947315) q[3];
sx q[3];
rz(0.36997488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8334373) q[2];
sx q[2];
rz(-0.94352949) q[2];
sx q[2];
rz(0.65307871) q[2];
rz(3.0270882) q[3];
sx q[3];
rz(-1.3936035) q[3];
sx q[3];
rz(-2.0792927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(0.36650518) q[0];
sx q[0];
rz(-1.0307642) q[0];
sx q[0];
rz(-1.3436226) q[0];
rz(-1.9812745) q[1];
sx q[1];
rz(-0.41028816) q[1];
sx q[1];
rz(0.62058273) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7337429) q[0];
sx q[0];
rz(-1.3792999) q[0];
sx q[0];
rz(-1.922185) q[0];
rz(-pi) q[1];
rz(0.36850117) q[2];
sx q[2];
rz(-0.93612367) q[2];
sx q[2];
rz(2.2919185) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0479792) q[1];
sx q[1];
rz(-1.2097712) q[1];
sx q[1];
rz(-2.4314636) q[1];
x q[2];
rz(-0.95497953) q[3];
sx q[3];
rz(-1.8035676) q[3];
sx q[3];
rz(1.7580838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4429984) q[2];
sx q[2];
rz(-1.6220762) q[2];
sx q[2];
rz(-0.2300187) q[2];
rz(1.5551785) q[3];
sx q[3];
rz(-1.14862) q[3];
sx q[3];
rz(-0.27568278) q[3];
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
rz(2.8860633) q[0];
sx q[0];
rz(-2.7524188) q[0];
sx q[0];
rz(-1.080876) q[0];
rz(-0.64322645) q[1];
sx q[1];
rz(-0.79935646) q[1];
sx q[1];
rz(0.08531514) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9179419) q[0];
sx q[0];
rz(-2.1408014) q[0];
sx q[0];
rz(2.355769) q[0];
rz(-pi) q[1];
rz(-1.7833461) q[2];
sx q[2];
rz(-2.1324799) q[2];
sx q[2];
rz(-1.4457653) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.18135611) q[1];
sx q[1];
rz(-1.4048368) q[1];
sx q[1];
rz(-1.7039429) q[1];
rz(-pi) q[2];
rz(3.078104) q[3];
sx q[3];
rz(-2.5194019) q[3];
sx q[3];
rz(2.6486462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9366511) q[2];
sx q[2];
rz(-3.1320269) q[2];
sx q[2];
rz(2.9452475) q[2];
rz(2.8593072) q[3];
sx q[3];
rz(-1.117027) q[3];
sx q[3];
rz(2.096368) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1333756) q[0];
sx q[0];
rz(-0.5468002) q[0];
sx q[0];
rz(-0.38544449) q[0];
rz(1.4133981) q[1];
sx q[1];
rz(-1.2047267) q[1];
sx q[1];
rz(-2.8916496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1544558) q[0];
sx q[0];
rz(-2.3780795) q[0];
sx q[0];
rz(-2.187243) q[0];
rz(-0.59993292) q[2];
sx q[2];
rz(-0.10641185) q[2];
sx q[2];
rz(2.0787042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39482597) q[1];
sx q[1];
rz(-1.4020846) q[1];
sx q[1];
rz(-0.49612237) q[1];
rz(-pi) q[2];
rz(-2.945845) q[3];
sx q[3];
rz(-2.3382814) q[3];
sx q[3];
rz(2.3565528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.46127737) q[2];
sx q[2];
rz(-1.2948493) q[2];
sx q[2];
rz(-2.9052022) q[2];
rz(-1.2605028) q[3];
sx q[3];
rz(-0.25006306) q[3];
sx q[3];
rz(-2.4642956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33609718) q[0];
sx q[0];
rz(-1.6039055) q[0];
sx q[0];
rz(0.48511037) q[0];
rz(-1.9612034) q[1];
sx q[1];
rz(-1.5943269) q[1];
sx q[1];
rz(2.3050883) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1331951) q[0];
sx q[0];
rz(-0.48497619) q[0];
sx q[0];
rz(0.72686813) q[0];
x q[1];
rz(2.2274687) q[2];
sx q[2];
rz(-0.49490041) q[2];
sx q[2];
rz(-0.47593853) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4662469) q[1];
sx q[1];
rz(-2.2398178) q[1];
sx q[1];
rz(-3.0803842) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6338634) q[3];
sx q[3];
rz(-1.9708037) q[3];
sx q[3];
rz(2.0162752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2430719) q[2];
sx q[2];
rz(-0.73183766) q[2];
sx q[2];
rz(2.3720429) q[2];
rz(-2.2122993) q[3];
sx q[3];
rz(-1.1309036) q[3];
sx q[3];
rz(0.23269674) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9418075) q[0];
sx q[0];
rz(-0.48518825) q[0];
sx q[0];
rz(-2.3727681) q[0];
rz(-1.7898412) q[1];
sx q[1];
rz(-0.47305802) q[1];
sx q[1];
rz(0.88968712) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2425491) q[0];
sx q[0];
rz(-0.43936037) q[0];
sx q[0];
rz(-2.0813294) q[0];
x q[1];
rz(-0.0013498505) q[2];
sx q[2];
rz(-1.5264411) q[2];
sx q[2];
rz(-2.7173017) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5701712) q[1];
sx q[1];
rz(-1.716905) q[1];
sx q[1];
rz(0.04695462) q[1];
rz(-0.23402782) q[3];
sx q[3];
rz(-0.46189538) q[3];
sx q[3];
rz(1.6017101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49030226) q[2];
sx q[2];
rz(-1.4375261) q[2];
sx q[2];
rz(1.5865405) q[2];
rz(-1.2706612) q[3];
sx q[3];
rz(-2.7226166) q[3];
sx q[3];
rz(-0.13203013) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9689869) q[0];
sx q[0];
rz(-1.1382599) q[0];
sx q[0];
rz(-1.9975115) q[0];
rz(0.83089685) q[1];
sx q[1];
rz(-1.7832719) q[1];
sx q[1];
rz(0.78713083) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68987304) q[0];
sx q[0];
rz(-1.4362808) q[0];
sx q[0];
rz(-3.0089507) q[0];
rz(-pi) q[1];
rz(-2.2550341) q[2];
sx q[2];
rz(-1.1349003) q[2];
sx q[2];
rz(2.2615711) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5980144) q[1];
sx q[1];
rz(-2.0435963) q[1];
sx q[1];
rz(2.6101357) q[1];
x q[2];
rz(-0.55946405) q[3];
sx q[3];
rz(-1.9064404) q[3];
sx q[3];
rz(-2.1223202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0234915) q[2];
sx q[2];
rz(-0.96994895) q[2];
sx q[2];
rz(-3.1084295) q[2];
rz(-1.5432594) q[3];
sx q[3];
rz(-2.4256746) q[3];
sx q[3];
rz(1.6395114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24932662) q[0];
sx q[0];
rz(-1.5557657) q[0];
sx q[0];
rz(1.3931042) q[0];
rz(2.6847367) q[1];
sx q[1];
rz(-1.1445878) q[1];
sx q[1];
rz(1.1005864) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7129214) q[0];
sx q[0];
rz(-1.572519) q[0];
sx q[0];
rz(3.0655906) q[0];
rz(-pi) q[1];
x q[1];
rz(0.7347091) q[2];
sx q[2];
rz(-0.8438973) q[2];
sx q[2];
rz(-1.5555842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7010569) q[1];
sx q[1];
rz(-0.66796366) q[1];
sx q[1];
rz(-1.6255767) q[1];
x q[2];
rz(-0.32842095) q[3];
sx q[3];
rz(-0.41116086) q[3];
sx q[3];
rz(2.823165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0988203) q[2];
sx q[2];
rz(-0.49751147) q[2];
sx q[2];
rz(1.5607321) q[2];
rz(2.3780195) q[3];
sx q[3];
rz(-1.6288501) q[3];
sx q[3];
rz(1.0360576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7324657) q[0];
sx q[0];
rz(-0.75508535) q[0];
sx q[0];
rz(0.29148802) q[0];
rz(1.2358707) q[1];
sx q[1];
rz(-1.9449077) q[1];
sx q[1];
rz(-3.018697) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9366347) q[0];
sx q[0];
rz(-2.7911148) q[0];
sx q[0];
rz(2.9913783) q[0];
x q[1];
rz(0.9083304) q[2];
sx q[2];
rz(-2.789302) q[2];
sx q[2];
rz(-0.34991821) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0258515) q[1];
sx q[1];
rz(-2.9065955) q[1];
sx q[1];
rz(2.0629289) q[1];
rz(1.6121665) q[3];
sx q[3];
rz(-0.14474104) q[3];
sx q[3];
rz(-0.54171692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5128532) q[2];
sx q[2];
rz(-1.3924007) q[2];
sx q[2];
rz(2.0780308) q[2];
rz(-0.17744803) q[3];
sx q[3];
rz(-2.1833503) q[3];
sx q[3];
rz(-2.1232846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.43750957) q[0];
sx q[0];
rz(-2.6884485) q[0];
sx q[0];
rz(1.4310687) q[0];
rz(-0.60072947) q[1];
sx q[1];
rz(-2.6917916) q[1];
sx q[1];
rz(-2.5240555) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8059313) q[0];
sx q[0];
rz(-1.7863723) q[0];
sx q[0];
rz(-0.12713253) q[0];
x q[1];
rz(2.3138758) q[2];
sx q[2];
rz(-1.9736145) q[2];
sx q[2];
rz(-0.043443505) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.85166288) q[1];
sx q[1];
rz(-1.549822) q[1];
sx q[1];
rz(2.1037654) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57281877) q[3];
sx q[3];
rz(-1.8316934) q[3];
sx q[3];
rz(0.12810055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2403229) q[2];
sx q[2];
rz(-2.1412886) q[2];
sx q[2];
rz(-0.97736248) q[2];
rz(-2.0591002) q[3];
sx q[3];
rz(-2.2668224) q[3];
sx q[3];
rz(0.055518363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33901535) q[0];
sx q[0];
rz(-0.76887283) q[0];
sx q[0];
rz(-0.73706891) q[0];
rz(-3.0958685) q[1];
sx q[1];
rz(-0.8174236) q[1];
sx q[1];
rz(1.5541706) q[1];
rz(0.81188079) q[2];
sx q[2];
rz(-2.4627081) q[2];
sx q[2];
rz(-2.2344786) q[2];
rz(0.65683881) q[3];
sx q[3];
rz(-1.7056864) q[3];
sx q[3];
rz(1.0581072) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
