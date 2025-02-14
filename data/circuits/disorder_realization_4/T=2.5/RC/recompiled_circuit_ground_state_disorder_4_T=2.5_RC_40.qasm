OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.5623915) q[0];
sx q[0];
rz(-2.3810823) q[0];
sx q[0];
rz(-1.0150681) q[0];
rz(-1.1633582) q[1];
sx q[1];
rz(-2.7203163) q[1];
sx q[1];
rz(-1.2383229) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9559798) q[0];
sx q[0];
rz(-0.85991128) q[0];
sx q[0];
rz(0.86998633) q[0];
rz(0.36739393) q[2];
sx q[2];
rz(-1.9695373) q[2];
sx q[2];
rz(-0.90786394) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4290966) q[1];
sx q[1];
rz(-1.1084021) q[1];
sx q[1];
rz(-2.7346947) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7389033) q[3];
sx q[3];
rz(-2.1283177) q[3];
sx q[3];
rz(1.3489189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.51332659) q[2];
sx q[2];
rz(-1.790975) q[2];
sx q[2];
rz(-2.4386621) q[2];
rz(0.85567307) q[3];
sx q[3];
rz(-0.84508768) q[3];
sx q[3];
rz(-1.2969016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.6913476) q[0];
sx q[0];
rz(-2.0424728) q[0];
sx q[0];
rz(2.3968089) q[0];
rz(-1.567747) q[1];
sx q[1];
rz(-2.4796922) q[1];
sx q[1];
rz(-0.94211284) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73977913) q[0];
sx q[0];
rz(-1.3688068) q[0];
sx q[0];
rz(0.76633472) q[0];
rz(0.81824192) q[2];
sx q[2];
rz(-1.628792) q[2];
sx q[2];
rz(2.4596283) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.86650634) q[1];
sx q[1];
rz(-0.77436354) q[1];
sx q[1];
rz(1.7656209) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5553988) q[3];
sx q[3];
rz(-1.6209416) q[3];
sx q[3];
rz(-1.3231089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2524903) q[2];
sx q[2];
rz(-0.95688755) q[2];
sx q[2];
rz(-2.0460184) q[2];
rz(-1.4787632) q[3];
sx q[3];
rz(-0.79083276) q[3];
sx q[3];
rz(-1.3862632) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0304994) q[0];
sx q[0];
rz(-1.4249304) q[0];
sx q[0];
rz(1.7362562) q[0];
rz(1.3966712) q[1];
sx q[1];
rz(-1.7878572) q[1];
sx q[1];
rz(0.90000802) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42381313) q[0];
sx q[0];
rz(-0.85067816) q[0];
sx q[0];
rz(-0.49085842) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58061231) q[2];
sx q[2];
rz(-1.0913278) q[2];
sx q[2];
rz(1.978385) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.69086466) q[1];
sx q[1];
rz(-2.303586) q[1];
sx q[1];
rz(-1.3370675) q[1];
rz(-pi) q[2];
rz(2.4803512) q[3];
sx q[3];
rz(-2.0783278) q[3];
sx q[3];
rz(0.37933168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46093837) q[2];
sx q[2];
rz(-0.21549455) q[2];
sx q[2];
rz(-2.1920965) q[2];
rz(-2.5203868) q[3];
sx q[3];
rz(-2.216279) q[3];
sx q[3];
rz(1.9882103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42156521) q[0];
sx q[0];
rz(-1.8864487) q[0];
sx q[0];
rz(3.0928639) q[0];
rz(-0.85982927) q[1];
sx q[1];
rz(-2.8099334) q[1];
sx q[1];
rz(-2.8299832) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6831419) q[0];
sx q[0];
rz(-1.299017) q[0];
sx q[0];
rz(-2.0608451) q[0];
rz(1.0428869) q[2];
sx q[2];
rz(-1.2706869) q[2];
sx q[2];
rz(-1.8703415) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5587363) q[1];
sx q[1];
rz(-1.8355825) q[1];
sx q[1];
rz(-0.89465054) q[1];
rz(2.0120088) q[3];
sx q[3];
rz(-1.8407243) q[3];
sx q[3];
rz(2.2306311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.17778808) q[2];
sx q[2];
rz(-2.9298941) q[2];
sx q[2];
rz(1.6507899) q[2];
rz(-2.7866411) q[3];
sx q[3];
rz(-1.4070516) q[3];
sx q[3];
rz(1.2995592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0519003) q[0];
sx q[0];
rz(-2.7841452) q[0];
sx q[0];
rz(-1.2667013) q[0];
rz(3.025324) q[1];
sx q[1];
rz(-2.1513042) q[1];
sx q[1];
rz(-1.6273392) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8775692) q[0];
sx q[0];
rz(-2.2248587) q[0];
sx q[0];
rz(-2.7436849) q[0];
x q[1];
rz(-1.9562938) q[2];
sx q[2];
rz(-0.23699871) q[2];
sx q[2];
rz(2.1182107) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0872495) q[1];
sx q[1];
rz(-1.6218446) q[1];
sx q[1];
rz(-2.1150949) q[1];
x q[2];
rz(-1.5290401) q[3];
sx q[3];
rz(-2.4126518) q[3];
sx q[3];
rz(1.6870013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2221471) q[2];
sx q[2];
rz(-2.6667892) q[2];
sx q[2];
rz(1.9661281) q[2];
rz(0.90421024) q[3];
sx q[3];
rz(-1.892482) q[3];
sx q[3];
rz(2.1632532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-1.0670052) q[0];
sx q[0];
rz(-2.8103204) q[0];
sx q[0];
rz(-2.7401155) q[0];
rz(1.1270771) q[1];
sx q[1];
rz(-1.3104442) q[1];
sx q[1];
rz(2.3289767) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82040374) q[0];
sx q[0];
rz(-2.2049453) q[0];
sx q[0];
rz(-3.1364596) q[0];
rz(-pi) q[1];
rz(-0.27490669) q[2];
sx q[2];
rz(-1.3454352) q[2];
sx q[2];
rz(-1.7079086) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5548812) q[1];
sx q[1];
rz(-2.5062392) q[1];
sx q[1];
rz(-0.34907267) q[1];
x q[2];
rz(0.62536247) q[3];
sx q[3];
rz(-0.84983045) q[3];
sx q[3];
rz(1.687885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.20069417) q[2];
sx q[2];
rz(-0.55018598) q[2];
sx q[2];
rz(-1.5852488) q[2];
rz(2.9511792) q[3];
sx q[3];
rz(-2.3868581) q[3];
sx q[3];
rz(1.93369) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82454005) q[0];
sx q[0];
rz(-2.7523478) q[0];
sx q[0];
rz(-3.0188766) q[0];
rz(1.1931194) q[1];
sx q[1];
rz(-0.90843186) q[1];
sx q[1];
rz(-2.3741123) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.240399) q[0];
sx q[0];
rz(-1.6220105) q[0];
sx q[0];
rz(-2.2177494) q[0];
rz(-1.6558455) q[2];
sx q[2];
rz(-1.951393) q[2];
sx q[2];
rz(-2.2589661) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2042522) q[1];
sx q[1];
rz(-2.2240727) q[1];
sx q[1];
rz(-0.055257992) q[1];
rz(0.010482739) q[3];
sx q[3];
rz(-0.40626486) q[3];
sx q[3];
rz(0.83028136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7439338) q[2];
sx q[2];
rz(-0.021641061) q[2];
sx q[2];
rz(-1.5519315) q[2];
rz(1.4895561) q[3];
sx q[3];
rz(-1.7347615) q[3];
sx q[3];
rz(1.1154729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3234696) q[0];
sx q[0];
rz(-0.36190811) q[0];
sx q[0];
rz(-0.37539151) q[0];
rz(0.88343945) q[1];
sx q[1];
rz(-1.6049623) q[1];
sx q[1];
rz(0.72867957) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7405072) q[0];
sx q[0];
rz(-2.2621763) q[0];
sx q[0];
rz(-2.2702433) q[0];
rz(-2.905718) q[2];
sx q[2];
rz(-1.7348924) q[2];
sx q[2];
rz(0.99170384) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9364215) q[1];
sx q[1];
rz(-1.9269619) q[1];
sx q[1];
rz(2.2486399) q[1];
x q[2];
rz(-0.75499423) q[3];
sx q[3];
rz(-2.0557311) q[3];
sx q[3];
rz(2.274853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4484619) q[2];
sx q[2];
rz(-1.5784266) q[2];
sx q[2];
rz(-0.7473839) q[2];
rz(-1.735431) q[3];
sx q[3];
rz(-1.2179255) q[3];
sx q[3];
rz(1.5860484) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2324227) q[0];
sx q[0];
rz(-2.2012043) q[0];
sx q[0];
rz(2.4160093) q[0];
rz(-1.8046509) q[1];
sx q[1];
rz(-1.2533816) q[1];
sx q[1];
rz(-0.57428378) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5919898) q[0];
sx q[0];
rz(-1.4583602) q[0];
sx q[0];
rz(-1.4103343) q[0];
x q[1];
rz(-2.1663675) q[2];
sx q[2];
rz(-0.65845602) q[2];
sx q[2];
rz(2.0582046) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7567819) q[1];
sx q[1];
rz(-0.73721209) q[1];
sx q[1];
rz(-1.4016777) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58135191) q[3];
sx q[3];
rz(-0.87331088) q[3];
sx q[3];
rz(0.70840981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9743222) q[2];
sx q[2];
rz(-1.3784626) q[2];
sx q[2];
rz(2.4794225) q[2];
rz(-0.76464379) q[3];
sx q[3];
rz(-1.5352826) q[3];
sx q[3];
rz(1.3242599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26226703) q[0];
sx q[0];
rz(-0.58795324) q[0];
sx q[0];
rz(1.3978488) q[0];
rz(2.0715879) q[1];
sx q[1];
rz(-1.1281697) q[1];
sx q[1];
rz(-1.3319344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4797213) q[0];
sx q[0];
rz(-0.1537696) q[0];
sx q[0];
rz(-1.3698306) q[0];
x q[1];
rz(1.0584303) q[2];
sx q[2];
rz(-2.67423) q[2];
sx q[2];
rz(-0.99603727) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3146914) q[1];
sx q[1];
rz(-1.3413635) q[1];
sx q[1];
rz(-1.8073842) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2666107) q[3];
sx q[3];
rz(-2.0379645) q[3];
sx q[3];
rz(2.9206729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3146882) q[2];
sx q[2];
rz(-2.0411699) q[2];
sx q[2];
rz(2.9617214) q[2];
rz(1.7449069) q[3];
sx q[3];
rz(-1.7161918) q[3];
sx q[3];
rz(0.21656187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3708645) q[0];
sx q[0];
rz(-0.48038078) q[0];
sx q[0];
rz(0.90850716) q[0];
rz(2.6188359) q[1];
sx q[1];
rz(-1.1154543) q[1];
sx q[1];
rz(1.9784068) q[1];
rz(2.592587) q[2];
sx q[2];
rz(-2.4291951) q[2];
sx q[2];
rz(0.090428314) q[2];
rz(1.9369851) q[3];
sx q[3];
rz(-2.0156751) q[3];
sx q[3];
rz(-1.0734476) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
