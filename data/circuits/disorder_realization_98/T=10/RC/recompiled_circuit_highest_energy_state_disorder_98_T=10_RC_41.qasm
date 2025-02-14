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
rz(0.75818169) q[0];
sx q[0];
rz(-2.7739006) q[0];
sx q[0];
rz(-2.0453069) q[0];
rz(0.81049377) q[1];
sx q[1];
rz(-0.23263045) q[1];
sx q[1];
rz(-2.0356324) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61845022) q[0];
sx q[0];
rz(-1.9621657) q[0];
sx q[0];
rz(2.6871293) q[0];
rz(-pi) q[1];
rz(0.26768406) q[2];
sx q[2];
rz(-1.4776728) q[2];
sx q[2];
rz(2.6665319) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.80820647) q[1];
sx q[1];
rz(-1.574834) q[1];
sx q[1];
rz(1.5334849) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8572122) q[3];
sx q[3];
rz(-1.577498) q[3];
sx q[3];
rz(2.0363765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7652863) q[2];
sx q[2];
rz(-0.47069612) q[2];
sx q[2];
rz(0.37962309) q[2];
rz(1.5073353) q[3];
sx q[3];
rz(-2.0416656) q[3];
sx q[3];
rz(-0.96989337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64604243) q[0];
sx q[0];
rz(-2.8046785) q[0];
sx q[0];
rz(2.2406793) q[0];
rz(0.13149978) q[1];
sx q[1];
rz(-1.200518) q[1];
sx q[1];
rz(0.44201717) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64009066) q[0];
sx q[0];
rz(-1.8953865) q[0];
sx q[0];
rz(-1.6571655) q[0];
rz(3.1369741) q[2];
sx q[2];
rz(-1.5717662) q[2];
sx q[2];
rz(-2.0624954) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.17168853) q[1];
sx q[1];
rz(-1.101788) q[1];
sx q[1];
rz(0.63290872) q[1];
rz(-1.6824746) q[3];
sx q[3];
rz(-1.968921) q[3];
sx q[3];
rz(0.83440895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3758731) q[2];
sx q[2];
rz(-1.4292382) q[2];
sx q[2];
rz(0.6089375) q[2];
rz(0.40306148) q[3];
sx q[3];
rz(-1.9654704) q[3];
sx q[3];
rz(-2.6398931) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.32560638) q[0];
sx q[0];
rz(-0.41929647) q[0];
sx q[0];
rz(-2.0035279) q[0];
rz(-1.552938) q[1];
sx q[1];
rz(-0.76858968) q[1];
sx q[1];
rz(0.80702153) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.606756) q[0];
sx q[0];
rz(-1.5436158) q[0];
sx q[0];
rz(0.099221275) q[0];
rz(2.8093178) q[2];
sx q[2];
rz(-2.6746779) q[2];
sx q[2];
rz(-0.98540598) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0453607) q[1];
sx q[1];
rz(-1.7590471) q[1];
sx q[1];
rz(-1.4493311) q[1];
rz(1.9495973) q[3];
sx q[3];
rz(-1.3037494) q[3];
sx q[3];
rz(-2.7821845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1300065) q[2];
sx q[2];
rz(-2.5277972) q[2];
sx q[2];
rz(-1.2137132) q[2];
rz(0.064149292) q[3];
sx q[3];
rz(-1.6545273) q[3];
sx q[3];
rz(0.00042644342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59378687) q[0];
sx q[0];
rz(-1.2603899) q[0];
sx q[0];
rz(-2.2547145) q[0];
rz(-0.15874323) q[1];
sx q[1];
rz(-0.49247772) q[1];
sx q[1];
rz(1.278272) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87065164) q[0];
sx q[0];
rz(-0.88490153) q[0];
sx q[0];
rz(0.57620462) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4492369) q[2];
sx q[2];
rz(-2.468363) q[2];
sx q[2];
rz(3.0557951) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6311196) q[1];
sx q[1];
rz(-2.2620631) q[1];
sx q[1];
rz(0.078029836) q[1];
x q[2];
rz(0.15110357) q[3];
sx q[3];
rz(-1.0211111) q[3];
sx q[3];
rz(-2.0696039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.13901237) q[2];
sx q[2];
rz(-0.85550344) q[2];
sx q[2];
rz(0.95376897) q[2];
rz(1.9468797) q[3];
sx q[3];
rz(-0.84269968) q[3];
sx q[3];
rz(2.6660624) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139528) q[0];
sx q[0];
rz(-0.71641818) q[0];
sx q[0];
rz(-2.5415976) q[0];
rz(0.68296877) q[1];
sx q[1];
rz(-2.3536847) q[1];
sx q[1];
rz(2.8111828) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70008695) q[0];
sx q[0];
rz(-1.7528025) q[0];
sx q[0];
rz(-2.299593) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1938964) q[2];
sx q[2];
rz(-2.8665339) q[2];
sx q[2];
rz(2.3291921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8176387) q[1];
sx q[1];
rz(-1.7637225) q[1];
sx q[1];
rz(2.0617538) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0903301) q[3];
sx q[3];
rz(-2.3805729) q[3];
sx q[3];
rz(0.42945592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43634513) q[2];
sx q[2];
rz(-1.0785495) q[2];
sx q[2];
rz(-0.62502965) q[2];
rz(0.69508067) q[3];
sx q[3];
rz(-1.86097) q[3];
sx q[3];
rz(-3.0888016) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92796749) q[0];
sx q[0];
rz(-1.1518421) q[0];
sx q[0];
rz(1.3731765) q[0];
rz(-1.3881418) q[1];
sx q[1];
rz(-2.2699247) q[1];
sx q[1];
rz(1.53481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4563303) q[0];
sx q[0];
rz(-1.7182516) q[0];
sx q[0];
rz(-0.035650226) q[0];
rz(2.0892145) q[2];
sx q[2];
rz(-0.3568584) q[2];
sx q[2];
rz(-1.269358) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.36027137) q[1];
sx q[1];
rz(-2.2151673) q[1];
sx q[1];
rz(-1.0553947) q[1];
rz(-pi) q[2];
rz(-1.2535048) q[3];
sx q[3];
rz(-2.7219238) q[3];
sx q[3];
rz(-2.9515226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9390949) q[2];
sx q[2];
rz(-1.1700583) q[2];
sx q[2];
rz(1.6004174) q[2];
rz(2.1206858) q[3];
sx q[3];
rz(-1.7424135) q[3];
sx q[3];
rz(-0.30103621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.130126) q[0];
sx q[0];
rz(-0.13700329) q[0];
sx q[0];
rz(-0.89114183) q[0];
rz(-2.4961684) q[1];
sx q[1];
rz(-1.3963457) q[1];
sx q[1];
rz(-0.83271629) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1306097) q[0];
sx q[0];
rz(-2.0314706) q[0];
sx q[0];
rz(-1.4185262) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4786167) q[2];
sx q[2];
rz(-2.5155009) q[2];
sx q[2];
rz(-2.7220059) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.558185) q[1];
sx q[1];
rz(-1.3761569) q[1];
sx q[1];
rz(0.50048142) q[1];
rz(1.9644968) q[3];
sx q[3];
rz(-1.0648921) q[3];
sx q[3];
rz(-0.33123744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5523395) q[2];
sx q[2];
rz(-2.4746042) q[2];
sx q[2];
rz(1.5935295) q[2];
rz(-0.42406905) q[3];
sx q[3];
rz(-1.4196906) q[3];
sx q[3];
rz(-2.8467395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9645204) q[0];
sx q[0];
rz(-1.1966713) q[0];
sx q[0];
rz(-0.82343423) q[0];
rz(1.9442762) q[1];
sx q[1];
rz(-1.4875393) q[1];
sx q[1];
rz(-1.4607325) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7410008) q[0];
sx q[0];
rz(-1.680458) q[0];
sx q[0];
rz(-1.2506668) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9731198) q[2];
sx q[2];
rz(-2.1572067) q[2];
sx q[2];
rz(2.7975067) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7825111) q[1];
sx q[1];
rz(-2.4788279) q[1];
sx q[1];
rz(-3.0726391) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36007215) q[3];
sx q[3];
rz(-1.8723893) q[3];
sx q[3];
rz(-0.87226471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0921649) q[2];
sx q[2];
rz(-0.66114134) q[2];
sx q[2];
rz(-3.0180422) q[2];
rz(-2.3030247) q[3];
sx q[3];
rz(-1.1127915) q[3];
sx q[3];
rz(-2.6113094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11414828) q[0];
sx q[0];
rz(-2.1338978) q[0];
sx q[0];
rz(-1.7281519) q[0];
rz(0.89883262) q[1];
sx q[1];
rz(-0.90679589) q[1];
sx q[1];
rz(-2.5487505) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7031738) q[0];
sx q[0];
rz(-2.1577564) q[0];
sx q[0];
rz(-0.86721768) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4845092) q[2];
sx q[2];
rz(-1.3893428) q[2];
sx q[2];
rz(1.726651) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1072126) q[1];
sx q[1];
rz(-1.980956) q[1];
sx q[1];
rz(-0.35663794) q[1];
x q[2];
rz(-0.96372202) q[3];
sx q[3];
rz(-2.1603909) q[3];
sx q[3];
rz(2.251925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8574519) q[2];
sx q[2];
rz(-2.393674) q[2];
sx q[2];
rz(0.22135529) q[2];
rz(1.0434693) q[3];
sx q[3];
rz(-0.9404434) q[3];
sx q[3];
rz(-0.38844696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8429883) q[0];
sx q[0];
rz(-0.28268155) q[0];
sx q[0];
rz(-2.2931732) q[0];
rz(-0.09659718) q[1];
sx q[1];
rz(-1.074147) q[1];
sx q[1];
rz(0.75327795) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4254508) q[0];
sx q[0];
rz(-1.4708859) q[0];
sx q[0];
rz(1.7990756) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6069372) q[2];
sx q[2];
rz(-0.74191626) q[2];
sx q[2];
rz(-2.9336172) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.2836783) q[1];
sx q[1];
rz(-1.1925624) q[1];
sx q[1];
rz(-2.0973732) q[1];
x q[2];
rz(-2.4270699) q[3];
sx q[3];
rz(-1.9252637) q[3];
sx q[3];
rz(-0.14795854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0640556) q[2];
sx q[2];
rz(-0.83751837) q[2];
sx q[2];
rz(2.688664) q[2];
rz(0.60106599) q[3];
sx q[3];
rz(-1.6744924) q[3];
sx q[3];
rz(-2.1452904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30408981) q[0];
sx q[0];
rz(-1.6927728) q[0];
sx q[0];
rz(0.46020831) q[0];
rz(-2.9651463) q[1];
sx q[1];
rz(-0.23527589) q[1];
sx q[1];
rz(-1.489524) q[1];
rz(1.5714191) q[2];
sx q[2];
rz(-1.2221619) q[2];
sx q[2];
rz(-2.6978343) q[2];
rz(2.4963958) q[3];
sx q[3];
rz(-0.83958902) q[3];
sx q[3];
rz(-2.9506172) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
