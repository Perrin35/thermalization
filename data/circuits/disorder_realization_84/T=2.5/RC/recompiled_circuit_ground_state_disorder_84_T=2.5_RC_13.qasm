OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.328188) q[0];
sx q[0];
rz(-2.5321811) q[0];
sx q[0];
rz(-3.1031026) q[0];
rz(2.6961532) q[1];
sx q[1];
rz(-0.69094509) q[1];
sx q[1];
rz(-0.12707392) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0920757) q[0];
sx q[0];
rz(-2.319515) q[0];
sx q[0];
rz(0.51527937) q[0];
x q[1];
rz(1.5655925) q[2];
sx q[2];
rz(-1.570911) q[2];
sx q[2];
rz(3.1104627) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.03555402) q[1];
sx q[1];
rz(-1.8287204) q[1];
sx q[1];
rz(1.2300806) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1640763) q[3];
sx q[3];
rz(-1.9818056) q[3];
sx q[3];
rz(-1.051595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.27158633) q[2];
sx q[2];
rz(-0.25124696) q[2];
sx q[2];
rz(1.1146438) q[2];
rz(-0.33846578) q[3];
sx q[3];
rz(-1.4128127) q[3];
sx q[3];
rz(-2.4131405) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2289497) q[0];
sx q[0];
rz(-2.8680608) q[0];
sx q[0];
rz(1.9563414) q[0];
rz(-1.395902) q[1];
sx q[1];
rz(-1.4811367) q[1];
sx q[1];
rz(-3.0286068) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015152935) q[0];
sx q[0];
rz(-3.0852144) q[0];
sx q[0];
rz(3.0150816) q[0];
x q[1];
rz(-1.0795679) q[2];
sx q[2];
rz(-2.6473111) q[2];
sx q[2];
rz(-3.1411067) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14571602) q[1];
sx q[1];
rz(-1.7967581) q[1];
sx q[1];
rz(0.06468062) q[1];
rz(-pi) q[2];
rz(2.2295938) q[3];
sx q[3];
rz(-0.56345255) q[3];
sx q[3];
rz(-2.7200626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9129755) q[2];
sx q[2];
rz(-1.7824219) q[2];
sx q[2];
rz(-1.011298) q[2];
rz(0.060981123) q[3];
sx q[3];
rz(-2.0296378) q[3];
sx q[3];
rz(-0.27420592) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7762452) q[0];
sx q[0];
rz(-1.341935) q[0];
sx q[0];
rz(1.9976529) q[0];
rz(-1.9260319) q[1];
sx q[1];
rz(-2.2823915) q[1];
sx q[1];
rz(0.62917319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1935223) q[0];
sx q[0];
rz(-2.5430449) q[0];
sx q[0];
rz(-1.6295327) q[0];
rz(-pi) q[1];
rz(0.10444) q[2];
sx q[2];
rz(-1.0658437) q[2];
sx q[2];
rz(-1.6475519) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.28189714) q[1];
sx q[1];
rz(-2.7696361) q[1];
sx q[1];
rz(2.9093548) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6225624) q[3];
sx q[3];
rz(-1.8300149) q[3];
sx q[3];
rz(1.0362877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4048142) q[2];
sx q[2];
rz(-1.3457158) q[2];
sx q[2];
rz(-0.2573615) q[2];
rz(1.2579873) q[3];
sx q[3];
rz(-1.6519929) q[3];
sx q[3];
rz(-2.6902698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8266206) q[0];
sx q[0];
rz(-2.4007512) q[0];
sx q[0];
rz(-1.524087) q[0];
rz(1.357088) q[1];
sx q[1];
rz(-0.70248228) q[1];
sx q[1];
rz(-1.9125028) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1462018) q[0];
sx q[0];
rz(-0.48625356) q[0];
sx q[0];
rz(2.3381837) q[0];
rz(-pi) q[1];
rz(-0.76009373) q[2];
sx q[2];
rz(-1.0398391) q[2];
sx q[2];
rz(1.0345376) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1145852) q[1];
sx q[1];
rz(-0.31064992) q[1];
sx q[1];
rz(0.18410054) q[1];
rz(-pi) q[2];
rz(-0.92023682) q[3];
sx q[3];
rz(-2.4229675) q[3];
sx q[3];
rz(2.7176734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.412753) q[2];
sx q[2];
rz(-1.3430026) q[2];
sx q[2];
rz(-3.0573678) q[2];
rz(1.6526875) q[3];
sx q[3];
rz(-3.0907349) q[3];
sx q[3];
rz(-2.1648255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4892437) q[0];
sx q[0];
rz(-2.4459965) q[0];
sx q[0];
rz(-3.1074281) q[0];
rz(-1.661352) q[1];
sx q[1];
rz(-0.40373293) q[1];
sx q[1];
rz(1.9306978) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6519484) q[0];
sx q[0];
rz(-1.9495954) q[0];
sx q[0];
rz(2.8234711) q[0];
x q[1];
rz(-2.6043774) q[2];
sx q[2];
rz(-0.59241931) q[2];
sx q[2];
rz(2.7835961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9927499) q[1];
sx q[1];
rz(-2.0066432) q[1];
sx q[1];
rz(0.080912605) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7922395) q[3];
sx q[3];
rz(-1.4998271) q[3];
sx q[3];
rz(0.25191316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5898798) q[2];
sx q[2];
rz(-2.5248542) q[2];
sx q[2];
rz(1.7390772) q[2];
rz(1.9139404) q[3];
sx q[3];
rz(-0.66870767) q[3];
sx q[3];
rz(0.68784586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(0.27125204) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(1.2649076) q[0];
rz(-1.854031) q[1];
sx q[1];
rz(-2.2076905) q[1];
sx q[1];
rz(-2.5087779) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0429471) q[0];
sx q[0];
rz(-1.2448695) q[0];
sx q[0];
rz(-0.21975369) q[0];
rz(2.6980324) q[2];
sx q[2];
rz(-0.75492263) q[2];
sx q[2];
rz(-2.519543) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.63232374) q[1];
sx q[1];
rz(-0.41328584) q[1];
sx q[1];
rz(0.37915165) q[1];
rz(-1.0086083) q[3];
sx q[3];
rz(-2.5119532) q[3];
sx q[3];
rz(-1.255577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5536993) q[2];
sx q[2];
rz(-0.18612315) q[2];
sx q[2];
rz(0.57559377) q[2];
rz(-1.1147739) q[3];
sx q[3];
rz(-1.2131178) q[3];
sx q[3];
rz(-0.43220821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9734398) q[0];
sx q[0];
rz(-1.7519209) q[0];
sx q[0];
rz(2.6218276) q[0];
rz(1.4555367) q[1];
sx q[1];
rz(-0.74570233) q[1];
sx q[1];
rz(1.0882475) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3322743) q[0];
sx q[0];
rz(-1.921321) q[0];
sx q[0];
rz(2.4615898) q[0];
rz(-pi) q[1];
rz(1.0287766) q[2];
sx q[2];
rz(-1.7732829) q[2];
sx q[2];
rz(2.7378792) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.42667761) q[1];
sx q[1];
rz(-1.2132056) q[1];
sx q[1];
rz(2.5044051) q[1];
rz(-pi) q[2];
rz(0.37101118) q[3];
sx q[3];
rz(-2.5034222) q[3];
sx q[3];
rz(1.3294892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9600642) q[2];
sx q[2];
rz(-1.1390511) q[2];
sx q[2];
rz(-3.0354434) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0625668) q[0];
sx q[0];
rz(-1.8667969) q[0];
sx q[0];
rz(-1.6612735) q[0];
rz(-1.5638428) q[1];
sx q[1];
rz(-0.5636951) q[1];
sx q[1];
rz(-0.021934358) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88883884) q[0];
sx q[0];
rz(-1.4705747) q[0];
sx q[0];
rz(-0.96688156) q[0];
x q[1];
rz(-1.1925132) q[2];
sx q[2];
rz(-1.6334051) q[2];
sx q[2];
rz(-2.4516842) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1749461) q[1];
sx q[1];
rz(-1.3783598) q[1];
sx q[1];
rz(2.220146) q[1];
rz(-1.1975953) q[3];
sx q[3];
rz(-1.4828728) q[3];
sx q[3];
rz(0.67876354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.03453001) q[2];
sx q[2];
rz(-1.7835534) q[2];
sx q[2];
rz(0.061508451) q[2];
rz(-2.6574078) q[3];
sx q[3];
rz(-1.8248841) q[3];
sx q[3];
rz(0.22783247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43559647) q[0];
sx q[0];
rz(-1.0974925) q[0];
sx q[0];
rz(2.66658) q[0];
rz(-1.2844405) q[1];
sx q[1];
rz(-0.78452763) q[1];
sx q[1];
rz(0.49819836) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2044922) q[0];
sx q[0];
rz(-1.0255359) q[0];
sx q[0];
rz(1.7496566) q[0];
x q[1];
rz(3.118843) q[2];
sx q[2];
rz(-2.2265847) q[2];
sx q[2];
rz(1.5178258) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5048154) q[1];
sx q[1];
rz(-2.9801629) q[1];
sx q[1];
rz(-2.982711) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.60631246) q[3];
sx q[3];
rz(-2.2289235) q[3];
sx q[3];
rz(2.6325814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4755134) q[2];
sx q[2];
rz(-0.84045118) q[2];
sx q[2];
rz(-1.6125512) q[2];
rz(2.9764002) q[3];
sx q[3];
rz(-1.2673204) q[3];
sx q[3];
rz(-2.7388549) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2861479) q[0];
sx q[0];
rz(-0.5583455) q[0];
sx q[0];
rz(-0.98854351) q[0];
rz(0.63829154) q[1];
sx q[1];
rz(-2.0604362) q[1];
sx q[1];
rz(0.036103006) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7571018) q[0];
sx q[0];
rz(-0.57390139) q[0];
sx q[0];
rz(-1.6792273) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1616012) q[2];
sx q[2];
rz(-2.0124751) q[2];
sx q[2];
rz(-2.4001588) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8263146) q[1];
sx q[1];
rz(-0.75584303) q[1];
sx q[1];
rz(-2.3043413) q[1];
rz(-2.3309541) q[3];
sx q[3];
rz(-0.30690696) q[3];
sx q[3];
rz(-1.7561654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4315167) q[2];
sx q[2];
rz(-2.4301131) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8217736) q[0];
sx q[0];
rz(-1.736301) q[0];
sx q[0];
rz(2.2176493) q[0];
rz(-2.3453103) q[1];
sx q[1];
rz(-2.2932107) q[1];
sx q[1];
rz(-0.36416818) q[1];
rz(2.4160855) q[2];
sx q[2];
rz(-0.85559597) q[2];
sx q[2];
rz(-1.2676452) q[2];
rz(-0.85354652) q[3];
sx q[3];
rz(-2.2774936) q[3];
sx q[3];
rz(-1.3766589) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
