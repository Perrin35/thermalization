OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.81340462) q[0];
sx q[0];
rz(-0.60941154) q[0];
sx q[0];
rz(-0.038490064) q[0];
rz(2.6961532) q[1];
sx q[1];
rz(-0.69094509) q[1];
sx q[1];
rz(3.0145187) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7434397) q[0];
sx q[0];
rz(-2.2619704) q[0];
sx q[0];
rz(-1.0831867) q[0];
x q[1];
rz(-1.5928361) q[2];
sx q[2];
rz(-0.0052050455) q[2];
sx q[2];
rz(1.5798868) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5161612) q[1];
sx q[1];
rz(-1.8998084) q[1];
sx q[1];
rz(2.8686868) q[1];
x q[2];
rz(-2.1640763) q[3];
sx q[3];
rz(-1.9818056) q[3];
sx q[3];
rz(-2.0899977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.27158633) q[2];
sx q[2];
rz(-2.8903457) q[2];
sx q[2];
rz(2.0269488) q[2];
rz(-2.8031269) q[3];
sx q[3];
rz(-1.7287799) q[3];
sx q[3];
rz(-2.4131405) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.912643) q[0];
sx q[0];
rz(-0.27353188) q[0];
sx q[0];
rz(1.9563414) q[0];
rz(-1.395902) q[1];
sx q[1];
rz(-1.4811367) q[1];
sx q[1];
rz(-3.0286068) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4596371) q[0];
sx q[0];
rz(-1.5779061) q[0];
sx q[0];
rz(3.0856641) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.127203) q[2];
sx q[2];
rz(-1.796486) q[2];
sx q[2];
rz(2.0114102) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.13571067) q[1];
sx q[1];
rz(-0.23488472) q[1];
sx q[1];
rz(-1.2967111) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2295938) q[3];
sx q[3];
rz(-0.56345255) q[3];
sx q[3];
rz(0.4215301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.2286171) q[2];
sx q[2];
rz(-1.3591707) q[2];
sx q[2];
rz(-1.011298) q[2];
rz(-3.0806115) q[3];
sx q[3];
rz(-1.1119548) q[3];
sx q[3];
rz(0.27420592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36534742) q[0];
sx q[0];
rz(-1.341935) q[0];
sx q[0];
rz(1.9976529) q[0];
rz(1.2155608) q[1];
sx q[1];
rz(-2.2823915) q[1];
sx q[1];
rz(-2.5124195) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0191285) q[0];
sx q[0];
rz(-2.1681684) q[0];
sx q[0];
rz(0.040014223) q[0];
rz(2.0780656) q[2];
sx q[2];
rz(-1.4794297) q[2];
sx q[2];
rz(0.12742119) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6110826) q[1];
sx q[1];
rz(-1.2092918) q[1];
sx q[1];
rz(-1.481249) q[1];
rz(2.8820417) q[3];
sx q[3];
rz(-1.6208315) q[3];
sx q[3];
rz(-0.52122859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4048142) q[2];
sx q[2];
rz(-1.7958769) q[2];
sx q[2];
rz(0.2573615) q[2];
rz(-1.8836053) q[3];
sx q[3];
rz(-1.4895997) q[3];
sx q[3];
rz(-0.45132288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314972) q[0];
sx q[0];
rz(-0.74084145) q[0];
sx q[0];
rz(1.6175057) q[0];
rz(1.357088) q[1];
sx q[1];
rz(-2.4391104) q[1];
sx q[1];
rz(1.9125028) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1462018) q[0];
sx q[0];
rz(-2.6553391) q[0];
sx q[0];
rz(2.3381837) q[0];
rz(0.88987197) q[2];
sx q[2];
rz(-2.2069815) q[2];
sx q[2];
rz(3.0539587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9214529) q[1];
sx q[1];
rz(-1.8760257) q[1];
sx q[1];
rz(1.5120927) q[1];
rz(-pi) q[2];
rz(-0.92023682) q[3];
sx q[3];
rz(-0.71862513) q[3];
sx q[3];
rz(-2.7176734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7288397) q[2];
sx q[2];
rz(-1.3430026) q[2];
sx q[2];
rz(-0.084224852) q[2];
rz(-1.6526875) q[3];
sx q[3];
rz(-0.050857734) q[3];
sx q[3];
rz(0.97676718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65234891) q[0];
sx q[0];
rz(-0.69559613) q[0];
sx q[0];
rz(-3.1074281) q[0];
rz(-1.661352) q[1];
sx q[1];
rz(-2.7378597) q[1];
sx q[1];
rz(1.2108948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2171777) q[0];
sx q[0];
rz(-2.6519288) q[0];
sx q[0];
rz(-2.2368466) q[0];
rz(-pi) q[1];
rz(1.239085) q[2];
sx q[2];
rz(-1.0704652) q[2];
sx q[2];
rz(2.1608888) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8029866) q[1];
sx q[1];
rz(-2.698773) q[1];
sx q[1];
rz(-1.7426331) q[1];
rz(-pi) q[2];
rz(-0.072739425) q[3];
sx q[3];
rz(-1.3499198) q[3];
sx q[3];
rz(-1.8386724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5898798) q[2];
sx q[2];
rz(-2.5248542) q[2];
sx q[2];
rz(1.7390772) q[2];
rz(-1.9139404) q[3];
sx q[3];
rz(-0.66870767) q[3];
sx q[3];
rz(-0.68784586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8703406) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(-1.2649076) q[0];
rz(1.2875617) q[1];
sx q[1];
rz(-2.2076905) q[1];
sx q[1];
rz(-2.5087779) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098645602) q[0];
sx q[0];
rz(-1.8967231) q[0];
sx q[0];
rz(-0.21975369) q[0];
x q[1];
rz(0.44356029) q[2];
sx q[2];
rz(-2.38667) q[2];
sx q[2];
rz(0.62204966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9196284) q[1];
sx q[1];
rz(-1.1884513) q[1];
sx q[1];
rz(-1.4098806) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1232093) q[3];
sx q[3];
rz(-1.2515178) q[3];
sx q[3];
rz(2.3553987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5878933) q[2];
sx q[2];
rz(-2.9554695) q[2];
sx q[2];
rz(2.5659989) q[2];
rz(-2.0268188) q[3];
sx q[3];
rz(-1.9284748) q[3];
sx q[3];
rz(2.7093844) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1681528) q[0];
sx q[0];
rz(-1.7519209) q[0];
sx q[0];
rz(2.6218276) q[0];
rz(1.686056) q[1];
sx q[1];
rz(-2.3958903) q[1];
sx q[1];
rz(-2.0533452) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35991372) q[0];
sx q[0];
rz(-2.3895411) q[0];
sx q[0];
rz(-0.52668026) q[0];
x q[1];
rz(1.1920405) q[2];
sx q[2];
rz(-0.57504762) q[2];
sx q[2];
rz(-1.4894007) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.89061531) q[1];
sx q[1];
rz(-0.979662) q[1];
sx q[1];
rz(2.0059523) q[1];
x q[2];
rz(1.3081) q[3];
sx q[3];
rz(-2.1593931) q[3];
sx q[3];
rz(-1.3610507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.18152848) q[2];
sx q[2];
rz(-1.1390511) q[2];
sx q[2];
rz(0.10614928) q[2];
rz(-1.6541121) q[3];
sx q[3];
rz(-2.5663576) q[3];
sx q[3];
rz(0.15644786) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625668) q[0];
sx q[0];
rz(-1.8667969) q[0];
sx q[0];
rz(-1.4803192) q[0];
rz(1.5777499) q[1];
sx q[1];
rz(-2.5778975) q[1];
sx q[1];
rz(-3.1196583) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2527538) q[0];
sx q[0];
rz(-1.6710179) q[0];
sx q[0];
rz(-2.1747111) q[0];
rz(-1.1925132) q[2];
sx q[2];
rz(-1.6334051) q[2];
sx q[2];
rz(0.68990842) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1749461) q[1];
sx q[1];
rz(-1.7632329) q[1];
sx q[1];
rz(0.92144664) q[1];
rz(-pi) q[2];
rz(1.3335732) q[3];
sx q[3];
rz(-2.7586474) q[3];
sx q[3];
rz(2.4701729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1070626) q[2];
sx q[2];
rz(-1.7835534) q[2];
sx q[2];
rz(-3.0800842) q[2];
rz(-2.6574078) q[3];
sx q[3];
rz(-1.8248841) q[3];
sx q[3];
rz(0.22783247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43559647) q[0];
sx q[0];
rz(-1.0974925) q[0];
sx q[0];
rz(-2.66658) q[0];
rz(-1.2844405) q[1];
sx q[1];
rz(-0.78452763) q[1];
sx q[1];
rz(0.49819836) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6817956) q[0];
sx q[0];
rz(-1.7234992) q[0];
sx q[0];
rz(2.5891749) q[0];
rz(-pi) q[1];
rz(1.5412386) q[2];
sx q[2];
rz(-0.65612462) q[2];
sx q[2];
rz(-1.6610638) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0908689) q[1];
sx q[1];
rz(-1.5962287) q[1];
sx q[1];
rz(2.9821616) q[1];
rz(-pi) q[2];
rz(2.5352802) q[3];
sx q[3];
rz(-0.91266914) q[3];
sx q[3];
rz(0.50901124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66607928) q[2];
sx q[2];
rz(-0.84045118) q[2];
sx q[2];
rz(1.5290414) q[2];
rz(0.1651925) q[3];
sx q[3];
rz(-1.2673204) q[3];
sx q[3];
rz(-0.40273777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2861479) q[0];
sx q[0];
rz(-2.5832472) q[0];
sx q[0];
rz(-0.98854351) q[0];
rz(-0.63829154) q[1];
sx q[1];
rz(-2.0604362) q[1];
sx q[1];
rz(-0.036103006) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7571018) q[0];
sx q[0];
rz(-0.57390139) q[0];
sx q[0];
rz(-1.4623653) q[0];
rz(-pi) q[1];
rz(-2.2746207) q[2];
sx q[2];
rz(-0.72153202) q[2];
sx q[2];
rz(1.7448278) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3152781) q[1];
sx q[1];
rz(-2.3857496) q[1];
sx q[1];
rz(-0.83725137) q[1];
x q[2];
rz(2.3309541) q[3];
sx q[3];
rz(-2.8346857) q[3];
sx q[3];
rz(-1.7561654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.710076) q[2];
sx q[2];
rz(-0.71147951) q[2];
sx q[2];
rz(-2.9284488) q[2];
rz(1.3171116) q[3];
sx q[3];
rz(-0.20753838) q[3];
sx q[3];
rz(-2.3402787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31981907) q[0];
sx q[0];
rz(-1.4052916) q[0];
sx q[0];
rz(-0.92394335) q[0];
rz(-0.79628235) q[1];
sx q[1];
rz(-0.84838198) q[1];
sx q[1];
rz(2.7774245) q[1];
rz(0.71103617) q[2];
sx q[2];
rz(-2.0954676) q[2];
sx q[2];
rz(-2.3116805) q[2];
rz(-2.2880461) q[3];
sx q[3];
rz(-0.86409909) q[3];
sx q[3];
rz(1.7649337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
