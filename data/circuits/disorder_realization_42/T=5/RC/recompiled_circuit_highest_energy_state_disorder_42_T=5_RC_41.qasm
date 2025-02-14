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
rz(0.49405721) q[0];
sx q[0];
rz(2.9394975) q[0];
sx q[0];
rz(9.1261368) q[0];
rz(2.5477297) q[1];
sx q[1];
rz(-1.9057823) q[1];
sx q[1];
rz(-1.2716582) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49079417) q[0];
sx q[0];
rz(-1.4795924) q[0];
sx q[0];
rz(-1.5861938) q[0];
rz(-pi) q[1];
rz(3.0230672) q[2];
sx q[2];
rz(-1.6942319) q[2];
sx q[2];
rz(-0.92230421) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.16954452) q[1];
sx q[1];
rz(-1.3005592) q[1];
sx q[1];
rz(-2.1317905) q[1];
rz(-pi) q[2];
rz(0.94489495) q[3];
sx q[3];
rz(-2.6941819) q[3];
sx q[3];
rz(-1.5571896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9876081) q[2];
sx q[2];
rz(-0.5927023) q[2];
sx q[2];
rz(0.1565557) q[2];
rz(-0.17078677) q[3];
sx q[3];
rz(-0.6905061) q[3];
sx q[3];
rz(-1.4120215) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60468948) q[0];
sx q[0];
rz(-1.5645744) q[0];
sx q[0];
rz(2.1942196) q[0];
rz(-2.0715711) q[1];
sx q[1];
rz(-1.0173631) q[1];
sx q[1];
rz(-1.0191466) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4812093) q[0];
sx q[0];
rz(-1.0374746) q[0];
sx q[0];
rz(2.9571499) q[0];
rz(-pi) q[1];
x q[1];
rz(1.972946) q[2];
sx q[2];
rz(-1.1199335) q[2];
sx q[2];
rz(-0.27489812) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.85657287) q[1];
sx q[1];
rz(-1.1390705) q[1];
sx q[1];
rz(1.747471) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2484577) q[3];
sx q[3];
rz(-0.97468978) q[3];
sx q[3];
rz(-2.346938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3715839) q[2];
sx q[2];
rz(-1.2086478) q[2];
sx q[2];
rz(-0.51119512) q[2];
rz(2.8080071) q[3];
sx q[3];
rz(-0.89190069) q[3];
sx q[3];
rz(-2.3489478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0148575) q[0];
sx q[0];
rz(-0.52913409) q[0];
sx q[0];
rz(2.1281309) q[0];
rz(0.86499372) q[1];
sx q[1];
rz(-1.2428872) q[1];
sx q[1];
rz(1.5369044) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8415926) q[0];
sx q[0];
rz(-0.16491297) q[0];
sx q[0];
rz(2.4191036) q[0];
rz(-pi) q[1];
rz(-2.2584707) q[2];
sx q[2];
rz(-0.97849023) q[2];
sx q[2];
rz(-1.2186714) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9158244) q[1];
sx q[1];
rz(-1.9729923) q[1];
sx q[1];
rz(0.56176995) q[1];
x q[2];
rz(-2.7551094) q[3];
sx q[3];
rz(-1.0253275) q[3];
sx q[3];
rz(2.2487933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4889739) q[2];
sx q[2];
rz(-1.7964615) q[2];
sx q[2];
rz(-2.6441003) q[2];
rz(-0.8688212) q[3];
sx q[3];
rz(-2.7610064) q[3];
sx q[3];
rz(-0.07180056) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82197613) q[0];
sx q[0];
rz(-3.138534) q[0];
sx q[0];
rz(-0.36913607) q[0];
rz(-0.71040756) q[1];
sx q[1];
rz(-2.0901168) q[1];
sx q[1];
rz(-1.8416454) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6746691) q[0];
sx q[0];
rz(-2.3116768) q[0];
sx q[0];
rz(-3.039917) q[0];
rz(1.7252015) q[2];
sx q[2];
rz(-1.1325628) q[2];
sx q[2];
rz(1.1296425) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0365606) q[1];
sx q[1];
rz(-1.2249882) q[1];
sx q[1];
rz(1.2086091) q[1];
rz(-3.0336712) q[3];
sx q[3];
rz(-1.125959) q[3];
sx q[3];
rz(0.51714424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9909782) q[2];
sx q[2];
rz(-1.6852448) q[2];
sx q[2];
rz(0.72188226) q[2];
rz(2.7041096) q[3];
sx q[3];
rz(-0.36967725) q[3];
sx q[3];
rz(0.30445254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40780145) q[0];
sx q[0];
rz(-2.3423539) q[0];
sx q[0];
rz(-2.9275628) q[0];
rz(-0.80263823) q[1];
sx q[1];
rz(-1.1624348) q[1];
sx q[1];
rz(-0.42073694) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3823809) q[0];
sx q[0];
rz(-1.7457643) q[0];
sx q[0];
rz(0.079240327) q[0];
rz(-pi) q[1];
rz(-0.56511648) q[2];
sx q[2];
rz(-1.775913) q[2];
sx q[2];
rz(-0.8520593) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.72374684) q[1];
sx q[1];
rz(-1.6366948) q[1];
sx q[1];
rz(1.4562277) q[1];
x q[2];
rz(0.84263148) q[3];
sx q[3];
rz(-1.6602605) q[3];
sx q[3];
rz(0.35736978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.78571) q[2];
sx q[2];
rz(-2.7648401) q[2];
sx q[2];
rz(2.9316087) q[2];
rz(-2.0436132) q[3];
sx q[3];
rz(-2.4185601) q[3];
sx q[3];
rz(-0.22600225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(1.8534773) q[0];
sx q[0];
rz(-2.2324039) q[0];
sx q[0];
rz(-3.1262681) q[0];
rz(0.76459908) q[1];
sx q[1];
rz(-1.8761643) q[1];
sx q[1];
rz(0.24078029) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2368589) q[0];
sx q[0];
rz(-0.52182996) q[0];
sx q[0];
rz(-0.93469341) q[0];
rz(-pi) q[1];
rz(2.639995) q[2];
sx q[2];
rz(-1.7946417) q[2];
sx q[2];
rz(-2.1783592) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2041596) q[1];
sx q[1];
rz(-1.6361409) q[1];
sx q[1];
rz(-1.3974299) q[1];
rz(-pi) q[2];
rz(2.094029) q[3];
sx q[3];
rz(-1.2480471) q[3];
sx q[3];
rz(-2.9595295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8814055) q[2];
sx q[2];
rz(-2.0892961) q[2];
sx q[2];
rz(2.6336199) q[2];
rz(-0.14902614) q[3];
sx q[3];
rz(-0.38145724) q[3];
sx q[3];
rz(-1.6833359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7070865) q[0];
sx q[0];
rz(-2.2266946) q[0];
sx q[0];
rz(-1.1900505) q[0];
rz(3.1022364) q[1];
sx q[1];
rz(-1.0330361) q[1];
sx q[1];
rz(-2.8796223) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7527237) q[0];
sx q[0];
rz(-1.6181728) q[0];
sx q[0];
rz(-0.33426826) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2798669) q[2];
sx q[2];
rz(-0.48469886) q[2];
sx q[2];
rz(-2.16207) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6050755) q[1];
sx q[1];
rz(-1.9585155) q[1];
sx q[1];
rz(1.8907029) q[1];
rz(-pi) q[2];
x q[2];
rz(0.39851578) q[3];
sx q[3];
rz(-0.57999963) q[3];
sx q[3];
rz(2.2452584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1956341) q[2];
sx q[2];
rz(-2.8527263) q[2];
sx q[2];
rz(1.757901) q[2];
rz(-1.1727928) q[3];
sx q[3];
rz(-1.5044418) q[3];
sx q[3];
rz(0.25243944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70343542) q[0];
sx q[0];
rz(-0.96980888) q[0];
sx q[0];
rz(-2.5954212) q[0];
rz(-0.81900412) q[1];
sx q[1];
rz(-1.9592229) q[1];
sx q[1];
rz(2.6738653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87238246) q[0];
sx q[0];
rz(-2.8765686) q[0];
sx q[0];
rz(1.3779968) q[0];
rz(1.1557577) q[2];
sx q[2];
rz(-1.9318244) q[2];
sx q[2];
rz(2.4280809) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.88087) q[1];
sx q[1];
rz(-2.4308009) q[1];
sx q[1];
rz(2.3265762) q[1];
x q[2];
rz(1.3276576) q[3];
sx q[3];
rz(-1.4068651) q[3];
sx q[3];
rz(-1.0019913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.24432448) q[2];
sx q[2];
rz(-1.6017598) q[2];
sx q[2];
rz(0.77678219) q[2];
rz(1.6952093) q[3];
sx q[3];
rz(-1.8518238) q[3];
sx q[3];
rz(-2.9608534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1553797) q[0];
sx q[0];
rz(-1.7387094) q[0];
sx q[0];
rz(-3.0631995) q[0];
rz(2.8502803) q[1];
sx q[1];
rz(-0.71543175) q[1];
sx q[1];
rz(1.327347) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3122897) q[0];
sx q[0];
rz(-0.96594496) q[0];
sx q[0];
rz(0.74440794) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7545596) q[2];
sx q[2];
rz(-3.0245028) q[2];
sx q[2];
rz(2.4486604) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1256333) q[1];
sx q[1];
rz(-1.2984097) q[1];
sx q[1];
rz(0.65954882) q[1];
x q[2];
rz(2.6895608) q[3];
sx q[3];
rz(-0.31637438) q[3];
sx q[3];
rz(-1.8282426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3751117) q[2];
sx q[2];
rz(-1.8152619) q[2];
sx q[2];
rz(-2.7074936) q[2];
rz(-0.33657524) q[3];
sx q[3];
rz(-1.429129) q[3];
sx q[3];
rz(-0.3955287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3340988) q[0];
sx q[0];
rz(-2.8039126) q[0];
sx q[0];
rz(0.39921528) q[0];
rz(2.4576064) q[1];
sx q[1];
rz(-1.8268879) q[1];
sx q[1];
rz(-1.4113873) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.08049649) q[0];
sx q[0];
rz(-2.7640124) q[0];
sx q[0];
rz(0.49712531) q[0];
rz(-pi) q[1];
rz(-0.26587395) q[2];
sx q[2];
rz(-0.69593799) q[2];
sx q[2];
rz(1.5741866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.70727247) q[1];
sx q[1];
rz(-2.1967683) q[1];
sx q[1];
rz(-1.1202588) q[1];
rz(-pi) q[2];
rz(2.3606922) q[3];
sx q[3];
rz(-1.077543) q[3];
sx q[3];
rz(0.79003143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4513678) q[2];
sx q[2];
rz(-1.6121612) q[2];
sx q[2];
rz(-0.31768793) q[2];
rz(1.0090656) q[3];
sx q[3];
rz(-1.3721481) q[3];
sx q[3];
rz(-0.38869977) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7379363) q[0];
sx q[0];
rz(-1.6725578) q[0];
sx q[0];
rz(-2.9148711) q[0];
rz(0.33568385) q[1];
sx q[1];
rz(-0.93135584) q[1];
sx q[1];
rz(2.9109536) q[1];
rz(3.1138218) q[2];
sx q[2];
rz(-2.4774629) q[2];
sx q[2];
rz(-3.1240875) q[2];
rz(0.54343358) q[3];
sx q[3];
rz(-0.68250485) q[3];
sx q[3];
rz(-1.8547081) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
