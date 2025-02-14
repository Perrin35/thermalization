OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.0319808) q[0];
sx q[0];
rz(4.1469753) q[0];
sx q[0];
rz(8.4561705) q[0];
rz(0.80945102) q[1];
sx q[1];
rz(-1.5341772) q[1];
sx q[1];
rz(3.1117575) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063970621) q[0];
sx q[0];
rz(-2.1093624) q[0];
sx q[0];
rz(-1.2993815) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0709463) q[2];
sx q[2];
rz(-0.81977568) q[2];
sx q[2];
rz(1.9762612) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3181495) q[1];
sx q[1];
rz(-1.9572721) q[1];
sx q[1];
rz(2.3837368) q[1];
rz(-pi) q[2];
rz(-1.6601059) q[3];
sx q[3];
rz(-1.7062643) q[3];
sx q[3];
rz(-0.40218807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5030824) q[2];
sx q[2];
rz(-2.5677887) q[2];
sx q[2];
rz(-2.8199675) q[2];
rz(0.62421787) q[3];
sx q[3];
rz(-1.3892684) q[3];
sx q[3];
rz(1.6577087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.6931848) q[0];
sx q[0];
rz(-0.25633651) q[0];
sx q[0];
rz(-2.2053027) q[0];
rz(-1.2721277) q[1];
sx q[1];
rz(-1.5446168) q[1];
sx q[1];
rz(1.7535271) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17860912) q[0];
sx q[0];
rz(-0.22679193) q[0];
sx q[0];
rz(1.3139477) q[0];
rz(-pi) q[1];
x q[1];
rz(0.77735591) q[2];
sx q[2];
rz(-1.1749975) q[2];
sx q[2];
rz(-2.0344025) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.85390988) q[1];
sx q[1];
rz(-1.0398932) q[1];
sx q[1];
rz(-2.4413051) q[1];
rz(2.2514492) q[3];
sx q[3];
rz(-0.33484866) q[3];
sx q[3];
rz(-2.072842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0442514) q[2];
sx q[2];
rz(-2.1772431) q[2];
sx q[2];
rz(2.468289) q[2];
rz(1.8246957) q[3];
sx q[3];
rz(-1.8484867) q[3];
sx q[3];
rz(2.3015658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8116542) q[0];
sx q[0];
rz(-1.543777) q[0];
sx q[0];
rz(3.0918308) q[0];
rz(0.20637575) q[1];
sx q[1];
rz(-1.3866837) q[1];
sx q[1];
rz(-1.6516364) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0127209) q[0];
sx q[0];
rz(-1.6908592) q[0];
sx q[0];
rz(0.61693667) q[0];
rz(-pi) q[1];
rz(-2.5875638) q[2];
sx q[2];
rz(-1.1893335) q[2];
sx q[2];
rz(-1.9395525) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94585431) q[1];
sx q[1];
rz(-0.39195132) q[1];
sx q[1];
rz(2.3183104) q[1];
rz(-0.95155119) q[3];
sx q[3];
rz(-0.54446044) q[3];
sx q[3];
rz(1.8306554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9516307) q[2];
sx q[2];
rz(-2.5962574) q[2];
sx q[2];
rz(-3.087431) q[2];
rz(1.1545898) q[3];
sx q[3];
rz(-1.2879939) q[3];
sx q[3];
rz(-2.0258928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40256777) q[0];
sx q[0];
rz(-1.4250647) q[0];
sx q[0];
rz(-0.85790747) q[0];
rz(0.67445451) q[1];
sx q[1];
rz(-2.1916788) q[1];
sx q[1];
rz(-1.8106921) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75174369) q[0];
sx q[0];
rz(-1.3350335) q[0];
sx q[0];
rz(-1.430473) q[0];
rz(-pi) q[1];
rz(1.2670008) q[2];
sx q[2];
rz(-1.4080321) q[2];
sx q[2];
rz(3.0050584) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1787599) q[1];
sx q[1];
rz(-0.58738999) q[1];
sx q[1];
rz(-1.2580832) q[1];
rz(-pi) q[2];
rz(-3.1382124) q[3];
sx q[3];
rz(-1.9900609) q[3];
sx q[3];
rz(1.5730891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7669507) q[2];
sx q[2];
rz(-0.16051126) q[2];
sx q[2];
rz(-3.0373419) q[2];
rz(0.57779622) q[3];
sx q[3];
rz(-1.0775074) q[3];
sx q[3];
rz(-0.92208636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1330971) q[0];
sx q[0];
rz(-1.47559) q[0];
sx q[0];
rz(-0.9444899) q[0];
rz(1.3635483) q[1];
sx q[1];
rz(-1.3282158) q[1];
sx q[1];
rz(1.4469226) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3730225) q[0];
sx q[0];
rz(-1.3813586) q[0];
sx q[0];
rz(0.85590881) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7124921) q[2];
sx q[2];
rz(-2.8887914) q[2];
sx q[2];
rz(-1.6293874) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.31476281) q[1];
sx q[1];
rz(-2.3896425) q[1];
sx q[1];
rz(0.58677267) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.99926104) q[3];
sx q[3];
rz(-1.0715228) q[3];
sx q[3];
rz(1.5443791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1030582) q[2];
sx q[2];
rz(-0.47163042) q[2];
sx q[2];
rz(-0.72059694) q[2];
rz(-3.0125812) q[3];
sx q[3];
rz(-1.3611662) q[3];
sx q[3];
rz(-1.2989929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1593889) q[0];
sx q[0];
rz(-0.99697462) q[0];
sx q[0];
rz(2.5250315) q[0];
rz(0.95245862) q[1];
sx q[1];
rz(-1.1943694) q[1];
sx q[1];
rz(-1.8985101) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12981249) q[0];
sx q[0];
rz(-1.8300548) q[0];
sx q[0];
rz(-3.0002322) q[0];
rz(-pi) q[1];
rz(-0.41105777) q[2];
sx q[2];
rz(-2.2349544) q[2];
sx q[2];
rz(1.039584) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9437517) q[1];
sx q[1];
rz(-1.5103522) q[1];
sx q[1];
rz(-0.72058679) q[1];
rz(-pi) q[2];
rz(-2.1922621) q[3];
sx q[3];
rz(-1.5104745) q[3];
sx q[3];
rz(0.5746791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.2359961) q[2];
sx q[2];
rz(-0.2874898) q[2];
sx q[2];
rz(3.0671885) q[2];
rz(1.0807886) q[3];
sx q[3];
rz(-1.6483043) q[3];
sx q[3];
rz(1.0409482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3022795) q[0];
sx q[0];
rz(-1.9957207) q[0];
sx q[0];
rz(-3.0795414) q[0];
rz(-1.7265559) q[1];
sx q[1];
rz(-2.3686385) q[1];
sx q[1];
rz(-3.0516023) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71848559) q[0];
sx q[0];
rz(-2.2540917) q[0];
sx q[0];
rz(-2.8124269) q[0];
x q[1];
rz(-1.3700468) q[2];
sx q[2];
rz(-0.69789574) q[2];
sx q[2];
rz(-3.0332886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9518826) q[1];
sx q[1];
rz(-1.9196577) q[1];
sx q[1];
rz(-0.93367082) q[1];
rz(-0.84109938) q[3];
sx q[3];
rz(-2.3537618) q[3];
sx q[3];
rz(2.6177518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.1354847) q[2];
sx q[2];
rz(-1.8718448) q[2];
sx q[2];
rz(-0.6832068) q[2];
rz(-1.8190544) q[3];
sx q[3];
rz(-0.95219487) q[3];
sx q[3];
rz(-1.3873842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4676062) q[0];
sx q[0];
rz(-1.4961996) q[0];
sx q[0];
rz(2.6737387) q[0];
rz(2.4521008) q[1];
sx q[1];
rz(-2.4263224) q[1];
sx q[1];
rz(-0.47324866) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.015206451) q[0];
sx q[0];
rz(-1.6412982) q[0];
sx q[0];
rz(-1.7470933) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32903657) q[2];
sx q[2];
rz(-2.6758782) q[2];
sx q[2];
rz(0.28458111) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0508214) q[1];
sx q[1];
rz(-1.5863401) q[1];
sx q[1];
rz(-1.5883716) q[1];
rz(-pi) q[2];
rz(2.4062633) q[3];
sx q[3];
rz(-0.37041723) q[3];
sx q[3];
rz(0.35740023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0143934) q[2];
sx q[2];
rz(-1.8202123) q[2];
sx q[2];
rz(2.7180706) q[2];
rz(-1.8262156) q[3];
sx q[3];
rz(-0.49153057) q[3];
sx q[3];
rz(-1.8486842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.2143283) q[0];
sx q[0];
rz(-0.94896999) q[0];
sx q[0];
rz(1.8073136) q[0];
rz(1.8233874) q[1];
sx q[1];
rz(-2.2315836) q[1];
sx q[1];
rz(-3.1277025) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2023808) q[0];
sx q[0];
rz(-1.9343543) q[0];
sx q[0];
rz(-0.027287539) q[0];
x q[1];
rz(-2.7330796) q[2];
sx q[2];
rz(-2.7613104) q[2];
sx q[2];
rz(-2.399596) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1760343) q[1];
sx q[1];
rz(-1.3234884) q[1];
sx q[1];
rz(-0.97977389) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4327496) q[3];
sx q[3];
rz(-2.4073753) q[3];
sx q[3];
rz(-0.54121298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5615329) q[2];
sx q[2];
rz(-1.1416898) q[2];
sx q[2];
rz(-1.053099) q[2];
rz(-0.92755353) q[3];
sx q[3];
rz(-1.901123) q[3];
sx q[3];
rz(-2.659306) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3574361) q[0];
sx q[0];
rz(-2.744839) q[0];
sx q[0];
rz(2.5411153) q[0];
rz(2.8720169) q[1];
sx q[1];
rz(-2.4495864) q[1];
sx q[1];
rz(1.3788266) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10256448) q[0];
sx q[0];
rz(-1.2192654) q[0];
sx q[0];
rz(2.8536056) q[0];
rz(-pi) q[1];
rz(-1.8767886) q[2];
sx q[2];
rz(-1.9255133) q[2];
sx q[2];
rz(1.6911094) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.32374183) q[1];
sx q[1];
rz(-1.3881589) q[1];
sx q[1];
rz(-0.43437965) q[1];
rz(-pi) q[2];
rz(2.5559101) q[3];
sx q[3];
rz(-2.437311) q[3];
sx q[3];
rz(-1.7862939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1158585) q[2];
sx q[2];
rz(-0.33418843) q[2];
sx q[2];
rz(-2.1791229) q[2];
rz(2.1879503) q[3];
sx q[3];
rz(-1.0735268) q[3];
sx q[3];
rz(1.7474705) q[3];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7643323) q[0];
sx q[0];
rz(-0.88484103) q[0];
sx q[0];
rz(-1.9405889) q[0];
rz(1.3537021) q[1];
sx q[1];
rz(-1.9722912) q[1];
sx q[1];
rz(-0.73696662) q[1];
rz(-0.66966343) q[2];
sx q[2];
rz(-2.1047932) q[2];
sx q[2];
rz(1.3303458) q[2];
rz(-1.4312454) q[3];
sx q[3];
rz(-1.0679676) q[3];
sx q[3];
rz(-2.3936489) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
