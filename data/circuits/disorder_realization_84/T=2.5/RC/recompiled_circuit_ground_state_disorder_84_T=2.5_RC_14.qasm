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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.049517) q[0];
sx q[0];
rz(-2.319515) q[0];
sx q[0];
rz(0.51527937) q[0];
rz(1.5487566) q[2];
sx q[2];
rz(-0.0052050455) q[2];
sx q[2];
rz(1.5798868) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.03555402) q[1];
sx q[1];
rz(-1.3128723) q[1];
sx q[1];
rz(-1.9115121) q[1];
rz(-pi) q[2];
x q[2];
rz(0.97751632) q[3];
sx q[3];
rz(-1.1597871) q[3];
sx q[3];
rz(2.0899977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8700063) q[2];
sx q[2];
rz(-0.25124696) q[2];
sx q[2];
rz(-1.1146438) q[2];
rz(-0.33846578) q[3];
sx q[3];
rz(-1.4128127) q[3];
sx q[3];
rz(-2.4131405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2289497) q[0];
sx q[0];
rz(-2.8680608) q[0];
sx q[0];
rz(1.9563414) q[0];
rz(1.7456906) q[1];
sx q[1];
rz(-1.4811367) q[1];
sx q[1];
rz(0.11298583) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1264397) q[0];
sx q[0];
rz(-0.056378214) q[0];
sx q[0];
rz(-3.0150816) q[0];
rz(-pi) q[1];
rz(2.0143896) q[2];
sx q[2];
rz(-1.3451066) q[2];
sx q[2];
rz(-2.0114102) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.13571067) q[1];
sx q[1];
rz(-0.23488472) q[1];
sx q[1];
rz(-1.2967111) q[1];
rz(1.107502) q[3];
sx q[3];
rz(-1.2377081) q[3];
sx q[3];
rz(-1.4127914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.2286171) q[2];
sx q[2];
rz(-1.3591707) q[2];
sx q[2];
rz(1.011298) q[2];
rz(-0.060981123) q[3];
sx q[3];
rz(-2.0296378) q[3];
sx q[3];
rz(-2.8673867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7762452) q[0];
sx q[0];
rz(-1.7996576) q[0];
sx q[0];
rz(1.1439398) q[0];
rz(-1.2155608) q[1];
sx q[1];
rz(-0.85920119) q[1];
sx q[1];
rz(0.62917319) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0191285) q[0];
sx q[0];
rz(-2.1681684) q[0];
sx q[0];
rz(-0.040014223) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0371527) q[2];
sx q[2];
rz(-2.0757489) q[2];
sx q[2];
rz(-1.4940408) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.28189714) q[1];
sx q[1];
rz(-2.7696361) q[1];
sx q[1];
rz(-0.23223784) q[1];
x q[2];
rz(2.8820417) q[3];
sx q[3];
rz(-1.5207612) q[3];
sx q[3];
rz(0.52122859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.7367785) q[2];
sx q[2];
rz(-1.7958769) q[2];
sx q[2];
rz(0.2573615) q[2];
rz(-1.2579873) q[3];
sx q[3];
rz(-1.6519929) q[3];
sx q[3];
rz(2.6902698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8266206) q[0];
sx q[0];
rz(-2.4007512) q[0];
sx q[0];
rz(-1.524087) q[0];
rz(1.7845047) q[1];
sx q[1];
rz(-2.4391104) q[1];
sx q[1];
rz(-1.9125028) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8600677) q[0];
sx q[0];
rz(-1.2403786) q[0];
sx q[0];
rz(1.9343253) q[0];
rz(-2.3814989) q[2];
sx q[2];
rz(-1.0398391) q[2];
sx q[2];
rz(2.107055) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.3683161) q[1];
sx q[1];
rz(-1.6267836) q[1];
sx q[1];
rz(-2.8358688) q[1];
x q[2];
rz(-2.2213558) q[3];
sx q[3];
rz(-0.71862513) q[3];
sx q[3];
rz(-0.42391923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.412753) q[2];
sx q[2];
rz(-1.7985901) q[2];
sx q[2];
rz(3.0573678) q[2];
rz(1.6526875) q[3];
sx q[3];
rz(-0.050857734) q[3];
sx q[3];
rz(-0.97676718) q[3];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65234891) q[0];
sx q[0];
rz(-2.4459965) q[0];
sx q[0];
rz(3.1074281) q[0];
rz(1.4802406) q[1];
sx q[1];
rz(-2.7378597) q[1];
sx q[1];
rz(1.2108948) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1015625) q[0];
sx q[0];
rz(-1.8656601) q[0];
sx q[0];
rz(1.9676137) q[0];
rz(-pi) q[1];
x q[1];
rz(1.239085) q[2];
sx q[2];
rz(-1.0704652) q[2];
sx q[2];
rz(2.1608888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.38773266) q[1];
sx q[1];
rz(-1.6441303) q[1];
sx q[1];
rz(2.0078986) q[1];
rz(-pi) q[2];
rz(1.7922395) q[3];
sx q[3];
rz(-1.4998271) q[3];
sx q[3];
rz(2.8896795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5898798) q[2];
sx q[2];
rz(-2.5248542) q[2];
sx q[2];
rz(-1.4025154) q[2];
rz(-1.9139404) q[3];
sx q[3];
rz(-2.472885) q[3];
sx q[3];
rz(-2.4537468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27125204) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(-1.8766851) q[0];
rz(-1.854031) q[1];
sx q[1];
rz(-2.2076905) q[1];
sx q[1];
rz(-2.5087779) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7408376) q[0];
sx q[0];
rz(-1.3627865) q[0];
sx q[0];
rz(1.9041787) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1870474) q[2];
sx q[2];
rz(-2.2381659) q[2];
sx q[2];
rz(-1.9415034) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9196284) q[1];
sx q[1];
rz(-1.9531413) q[1];
sx q[1];
rz(1.4098806) q[1];
rz(-pi) q[2];
rz(-2.7711677) q[3];
sx q[3];
rz(-2.0923449) q[3];
sx q[3];
rz(0.5934779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.5536993) q[2];
sx q[2];
rz(-0.18612315) q[2];
sx q[2];
rz(-0.57559377) q[2];
rz(1.1147739) q[3];
sx q[3];
rz(-1.2131178) q[3];
sx q[3];
rz(0.43220821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9734398) q[0];
sx q[0];
rz(-1.3896717) q[0];
sx q[0];
rz(-2.6218276) q[0];
rz(1.686056) q[1];
sx q[1];
rz(-2.3958903) q[1];
sx q[1];
rz(-2.0533452) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1092508) q[0];
sx q[0];
rz(-0.93904385) q[0];
sx q[0];
rz(2.0103309) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9495522) q[2];
sx q[2];
rz(-0.57504762) q[2];
sx q[2];
rz(-1.4894007) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.714915) q[1];
sx q[1];
rz(-1.928387) q[1];
sx q[1];
rz(0.63718759) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5367706) q[3];
sx q[3];
rz(-1.3530952) q[3];
sx q[3];
rz(3.0800501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.18152848) q[2];
sx q[2];
rz(-2.0025415) q[2];
sx q[2];
rz(-3.0354434) q[2];
rz(-1.4874805) q[3];
sx q[3];
rz(-2.5663576) q[3];
sx q[3];
rz(-0.15644786) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079025896) q[0];
sx q[0];
rz(-1.8667969) q[0];
sx q[0];
rz(-1.6612735) q[0];
rz(-1.5638428) q[1];
sx q[1];
rz(-2.5778975) q[1];
sx q[1];
rz(0.021934358) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2527538) q[0];
sx q[0];
rz(-1.4705747) q[0];
sx q[0];
rz(0.96688156) q[0];
rz(-pi) q[1];
rz(1.9490795) q[2];
sx q[2];
rz(-1.6334051) q[2];
sx q[2];
rz(-2.4516842) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3932568) q[1];
sx q[1];
rz(-2.2062058) q[1];
sx q[1];
rz(-2.901668) q[1];
rz(-pi) q[2];
rz(0.094385967) q[3];
sx q[3];
rz(-1.1991074) q[3];
sx q[3];
rz(-2.2151912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.03453001) q[2];
sx q[2];
rz(-1.3580393) q[2];
sx q[2];
rz(0.061508451) q[2];
rz(-0.48418489) q[3];
sx q[3];
rz(-1.8248841) q[3];
sx q[3];
rz(2.9137602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7059962) q[0];
sx q[0];
rz(-2.0441002) q[0];
sx q[0];
rz(0.47501269) q[0];
rz(-1.8571521) q[1];
sx q[1];
rz(-2.357065) q[1];
sx q[1];
rz(0.49819836) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2044922) q[0];
sx q[0];
rz(-1.0255359) q[0];
sx q[0];
rz(1.3919361) q[0];
rz(-2.2267098) q[2];
sx q[2];
rz(-1.5888264) q[2];
sx q[2];
rz(-3.0747482) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0908689) q[1];
sx q[1];
rz(-1.5453639) q[1];
sx q[1];
rz(2.9821616) q[1];
rz(-pi) q[2];
rz(2.3257019) q[3];
sx q[3];
rz(-1.1031086) q[3];
sx q[3];
rz(-0.66064686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66607928) q[2];
sx q[2];
rz(-0.84045118) q[2];
sx q[2];
rz(-1.5290414) q[2];
rz(-2.9764002) q[3];
sx q[3];
rz(-1.8742722) q[3];
sx q[3];
rz(0.40273777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2861479) q[0];
sx q[0];
rz(-0.5583455) q[0];
sx q[0];
rz(-0.98854351) q[0];
rz(2.5033011) q[1];
sx q[1];
rz(-1.0811564) q[1];
sx q[1];
rz(-3.1054896) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38449088) q[0];
sx q[0];
rz(-2.5676913) q[0];
sx q[0];
rz(1.6792273) q[0];
x q[1];
rz(0.51757054) q[2];
sx q[2];
rz(-2.0985275) q[2];
sx q[2];
rz(0.55014683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.42369595) q[1];
sx q[1];
rz(-1.0362018) q[1];
sx q[1];
rz(2.5786492) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3309541) q[3];
sx q[3];
rz(-0.30690696) q[3];
sx q[3];
rz(1.7561654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4315167) q[2];
sx q[2];
rz(-2.4301131) q[2];
sx q[2];
rz(-0.21314387) q[2];
rz(1.824481) q[3];
sx q[3];
rz(-0.20753838) q[3];
sx q[3];
rz(2.3402787) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8217736) q[0];
sx q[0];
rz(-1.4052916) q[0];
sx q[0];
rz(-0.92394335) q[0];
rz(0.79628235) q[1];
sx q[1];
rz(-2.2932107) q[1];
sx q[1];
rz(-0.36416818) q[1];
rz(2.2231215) q[2];
sx q[2];
rz(-0.9705636) q[2];
sx q[2];
rz(2.8080429) q[2];
rz(2.2939479) q[3];
sx q[3];
rz(-2.0942735) q[3];
sx q[3];
rz(0.70944722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
