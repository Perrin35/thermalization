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
rz(-0.44543946) q[1];
sx q[1];
rz(-2.4506476) q[1];
sx q[1];
rz(-3.0145187) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7434397) q[0];
sx q[0];
rz(-0.87962224) q[0];
sx q[0];
rz(2.0584059) q[0];
rz(-pi) q[1];
rz(-1.5760001) q[2];
sx q[2];
rz(-1.570911) q[2];
sx q[2];
rz(3.1104627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5161612) q[1];
sx q[1];
rz(-1.8998084) q[1];
sx q[1];
rz(0.27290588) q[1];
rz(-pi) q[2];
rz(-0.48396516) q[3];
sx q[3];
rz(-2.1089156) q[3];
sx q[3];
rz(-2.3592268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.27158633) q[2];
sx q[2];
rz(-0.25124696) q[2];
sx q[2];
rz(-1.1146438) q[2];
rz(2.8031269) q[3];
sx q[3];
rz(-1.7287799) q[3];
sx q[3];
rz(-0.72845212) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.912643) q[0];
sx q[0];
rz(-0.27353188) q[0];
sx q[0];
rz(1.9563414) q[0];
rz(-1.7456906) q[1];
sx q[1];
rz(-1.4811367) q[1];
sx q[1];
rz(3.0286068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6819555) q[0];
sx q[0];
rz(-1.5779061) q[0];
sx q[0];
rz(0.05592859) q[0];
rz(1.0795679) q[2];
sx q[2];
rz(-2.6473111) q[2];
sx q[2];
rz(3.1411067) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.005882) q[1];
sx q[1];
rz(-0.23488472) q[1];
sx q[1];
rz(1.2967111) q[1];
rz(-pi) q[2];
rz(2.7725622) q[3];
sx q[3];
rz(-1.1347767) q[3];
sx q[3];
rz(-0.31991239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.9129755) q[2];
sx q[2];
rz(-1.7824219) q[2];
sx q[2];
rz(1.011298) q[2];
rz(-3.0806115) q[3];
sx q[3];
rz(-2.0296378) q[3];
sx q[3];
rz(2.8673867) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36534742) q[0];
sx q[0];
rz(-1.341935) q[0];
sx q[0];
rz(1.9976529) q[0];
rz(1.9260319) q[1];
sx q[1];
rz(-0.85920119) q[1];
sx q[1];
rz(-2.5124195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1224642) q[0];
sx q[0];
rz(-0.97342426) q[0];
sx q[0];
rz(-0.040014223) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0371527) q[2];
sx q[2];
rz(-2.0757489) q[2];
sx q[2];
rz(-1.4940408) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8596955) q[1];
sx q[1];
rz(-0.37195656) q[1];
sx q[1];
rz(0.23223784) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9488936) q[3];
sx q[3];
rz(-0.26422285) q[3];
sx q[3];
rz(1.2357354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4048142) q[2];
sx q[2];
rz(-1.3457158) q[2];
sx q[2];
rz(0.2573615) q[2];
rz(-1.8836053) q[3];
sx q[3];
rz(-1.6519929) q[3];
sx q[3];
rz(-2.6902698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.314972) q[0];
sx q[0];
rz(-2.4007512) q[0];
sx q[0];
rz(1.524087) q[0];
rz(1.7845047) q[1];
sx q[1];
rz(-2.4391104) q[1];
sx q[1];
rz(-1.9125028) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1462018) q[0];
sx q[0];
rz(-0.48625356) q[0];
sx q[0];
rz(0.803409) q[0];
x q[1];
rz(-0.76009373) q[2];
sx q[2];
rz(-2.1017535) q[2];
sx q[2];
rz(2.107055) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.22013979) q[1];
sx q[1];
rz(-1.8760257) q[1];
sx q[1];
rz(1.6294999) q[1];
rz(-pi) q[2];
rz(2.1788254) q[3];
sx q[3];
rz(-1.1606779) q[3];
sx q[3];
rz(-1.6671163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.412753) q[2];
sx q[2];
rz(-1.3430026) q[2];
sx q[2];
rz(0.084224852) q[2];
rz(-1.6526875) q[3];
sx q[3];
rz(-3.0907349) q[3];
sx q[3];
rz(-0.97676718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65234891) q[0];
sx q[0];
rz(-2.4459965) q[0];
sx q[0];
rz(0.034164567) q[0];
rz(-1.661352) q[1];
sx q[1];
rz(-0.40373293) q[1];
sx q[1];
rz(-1.2108948) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1015625) q[0];
sx q[0];
rz(-1.8656601) q[0];
sx q[0];
rz(-1.1739789) q[0];
rz(-1.239085) q[2];
sx q[2];
rz(-2.0711275) q[2];
sx q[2];
rz(-0.98070383) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3386061) q[1];
sx q[1];
rz(-0.44281964) q[1];
sx q[1];
rz(-1.7426331) q[1];
rz(-1.3493531) q[3];
sx q[3];
rz(-1.6417656) q[3];
sx q[3];
rz(0.25191316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5898798) q[2];
sx q[2];
rz(-0.61673841) q[2];
sx q[2];
rz(1.7390772) q[2];
rz(-1.2276522) q[3];
sx q[3];
rz(-0.66870767) q[3];
sx q[3];
rz(0.68784586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8703406) q[0];
sx q[0];
rz(-1.7043865) q[0];
sx q[0];
rz(1.8766851) q[0];
rz(1.2875617) q[1];
sx q[1];
rz(-2.2076905) q[1];
sx q[1];
rz(0.63281473) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4007551) q[0];
sx q[0];
rz(-1.3627865) q[0];
sx q[0];
rz(-1.9041787) q[0];
x q[1];
rz(-0.70437141) q[2];
sx q[2];
rz(-1.8692817) q[2];
sx q[2];
rz(0.61559144) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2883451) q[1];
sx q[1];
rz(-1.4215905) q[1];
sx q[1];
rz(0.38686727) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7711677) q[3];
sx q[3];
rz(-2.0923449) q[3];
sx q[3];
rz(-0.5934779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5878933) q[2];
sx q[2];
rz(-0.18612315) q[2];
sx q[2];
rz(0.57559377) q[2];
rz(2.0268188) q[3];
sx q[3];
rz(-1.9284748) q[3];
sx q[3];
rz(0.43220821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9734398) q[0];
sx q[0];
rz(-1.7519209) q[0];
sx q[0];
rz(2.6218276) q[0];
rz(-1.686056) q[1];
sx q[1];
rz(-2.3958903) q[1];
sx q[1];
rz(2.0533452) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0323419) q[0];
sx q[0];
rz(-2.2025488) q[0];
sx q[0];
rz(-2.0103309) q[0];
rz(-pi) q[1];
rz(-2.9063792) q[2];
sx q[2];
rz(-2.1005582) q[2];
sx q[2];
rz(1.0465682) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.89061531) q[1];
sx q[1];
rz(-2.1619307) q[1];
sx q[1];
rz(-1.1356403) q[1];
rz(-pi) q[2];
rz(-1.3081) q[3];
sx q[3];
rz(-0.98219959) q[3];
sx q[3];
rz(-1.3610507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9600642) q[2];
sx q[2];
rz(-2.0025415) q[2];
sx q[2];
rz(3.0354434) q[2];
rz(1.4874805) q[3];
sx q[3];
rz(-2.5663576) q[3];
sx q[3];
rz(0.15644786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079025896) q[0];
sx q[0];
rz(-1.2747958) q[0];
sx q[0];
rz(1.6612735) q[0];
rz(1.5777499) q[1];
sx q[1];
rz(-0.5636951) q[1];
sx q[1];
rz(3.1196583) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53793478) q[0];
sx q[0];
rz(-0.61115038) q[0];
sx q[0];
rz(-1.7460599) q[0];
rz(0.067358067) q[2];
sx q[2];
rz(-1.948301) q[2];
sx q[2];
rz(0.8560271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74833588) q[1];
sx q[1];
rz(-0.93538688) q[1];
sx q[1];
rz(2.901668) q[1];
x q[2];
rz(-0.094385967) q[3];
sx q[3];
rz(-1.1991074) q[3];
sx q[3];
rz(-0.92640141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.03453001) q[2];
sx q[2];
rz(-1.7835534) q[2];
sx q[2];
rz(-3.0800842) q[2];
rz(-2.6574078) q[3];
sx q[3];
rz(-1.8248841) q[3];
sx q[3];
rz(-2.9137602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7059962) q[0];
sx q[0];
rz(-1.0974925) q[0];
sx q[0];
rz(-2.66658) q[0];
rz(1.2844405) q[1];
sx q[1];
rz(-2.357065) q[1];
sx q[1];
rz(0.49819836) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9371004) q[0];
sx q[0];
rz(-2.1160567) q[0];
sx q[0];
rz(-1.3919361) q[0];
rz(2.2267098) q[2];
sx q[2];
rz(-1.5527662) q[2];
sx q[2];
rz(-3.0747482) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47583844) q[1];
sx q[1];
rz(-1.7301754) q[1];
sx q[1];
rz(1.5450374) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60631246) q[3];
sx q[3];
rz(-2.2289235) q[3];
sx q[3];
rz(-2.6325814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66607928) q[2];
sx q[2];
rz(-0.84045118) q[2];
sx q[2];
rz(-1.6125512) q[2];
rz(0.1651925) q[3];
sx q[3];
rz(-1.2673204) q[3];
sx q[3];
rz(2.7388549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85544473) q[0];
sx q[0];
rz(-0.5583455) q[0];
sx q[0];
rz(0.98854351) q[0];
rz(-0.63829154) q[1];
sx q[1];
rz(-2.0604362) q[1];
sx q[1];
rz(-0.036103006) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7571018) q[0];
sx q[0];
rz(-2.5676913) q[0];
sx q[0];
rz(1.6792273) q[0];
rz(-pi) q[1];
rz(2.1616012) q[2];
sx q[2];
rz(-2.0124751) q[2];
sx q[2];
rz(0.74143386) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3152781) q[1];
sx q[1];
rz(-2.3857496) q[1];
sx q[1];
rz(2.3043413) q[1];
x q[2];
rz(-0.21499459) q[3];
sx q[3];
rz(-1.3500596) q[3];
sx q[3];
rz(-2.1695015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.710076) q[2];
sx q[2];
rz(-2.4301131) q[2];
sx q[2];
rz(-2.9284488) q[2];
rz(-1.3171116) q[3];
sx q[3];
rz(-2.9340543) q[3];
sx q[3];
rz(-2.3402787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
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
rz(2.4160855) q[2];
sx q[2];
rz(-0.85559597) q[2];
sx q[2];
rz(-1.2676452) q[2];
rz(2.4854971) q[3];
sx q[3];
rz(-2.1811206) q[3];
sx q[3];
rz(2.695786) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
