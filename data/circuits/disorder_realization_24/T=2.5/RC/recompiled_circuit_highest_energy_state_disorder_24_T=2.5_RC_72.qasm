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
rz(2.6296122) q[0];
sx q[0];
rz(2.5706302) q[0];
sx q[0];
rz(9.6124967) q[0];
rz(-2.3277148) q[1];
sx q[1];
rz(-1.4717646) q[1];
sx q[1];
rz(-1.2893113) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9780802) q[0];
sx q[0];
rz(-2.8239282) q[0];
sx q[0];
rz(-0.72084041) q[0];
rz(1.5411096) q[2];
sx q[2];
rz(-2.1412537) q[2];
sx q[2];
rz(0.96889673) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.493452) q[1];
sx q[1];
rz(-1.8187742) q[1];
sx q[1];
rz(0.47391717) q[1];
rz(-0.23203316) q[3];
sx q[3];
rz(-1.9153144) q[3];
sx q[3];
rz(-1.0330878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2671555) q[2];
sx q[2];
rz(-1.1031373) q[2];
sx q[2];
rz(-0.77679408) q[2];
rz(1.3022425) q[3];
sx q[3];
rz(-1.4812508) q[3];
sx q[3];
rz(1.2031901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5341107) q[0];
sx q[0];
rz(-2.6597839) q[0];
sx q[0];
rz(2.1433461) q[0];
rz(-0.36969319) q[1];
sx q[1];
rz(-1.0414711) q[1];
sx q[1];
rz(1.1179771) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3415754) q[0];
sx q[0];
rz(-1.9722) q[0];
sx q[0];
rz(2.163351) q[0];
rz(-pi) q[1];
rz(-0.68685617) q[2];
sx q[2];
rz(-1.1026898) q[2];
sx q[2];
rz(-1.3324225) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.79827944) q[1];
sx q[1];
rz(-1.4237397) q[1];
sx q[1];
rz(1.1215854) q[1];
rz(-pi) q[2];
rz(-2.8727628) q[3];
sx q[3];
rz(-2.1850366) q[3];
sx q[3];
rz(-1.5419095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2758241) q[2];
sx q[2];
rz(-1.3752702) q[2];
sx q[2];
rz(0.38530525) q[2];
rz(-0.038711874) q[3];
sx q[3];
rz(-2.7888515) q[3];
sx q[3];
rz(2.501131) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63967079) q[0];
sx q[0];
rz(-0.47698912) q[0];
sx q[0];
rz(-1.3013526) q[0];
rz(-3.0838857) q[1];
sx q[1];
rz(-2.6951908) q[1];
sx q[1];
rz(-1.3267964) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8463752) q[0];
sx q[0];
rz(-2.1070988) q[0];
sx q[0];
rz(0.55732507) q[0];
rz(-0.71452629) q[2];
sx q[2];
rz(-2.9305612) q[2];
sx q[2];
rz(-2.2843366) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.51242764) q[1];
sx q[1];
rz(-1.8574134) q[1];
sx q[1];
rz(2.5932339) q[1];
rz(2.2928724) q[3];
sx q[3];
rz(-1.3009138) q[3];
sx q[3];
rz(0.11972846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3441299) q[2];
sx q[2];
rz(-0.8351438) q[2];
sx q[2];
rz(0.75378913) q[2];
rz(1.7720743) q[3];
sx q[3];
rz(-1.6794208) q[3];
sx q[3];
rz(0.65521017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264483) q[0];
sx q[0];
rz(-1.0855874) q[0];
sx q[0];
rz(1.7468859) q[0];
rz(-1.1766379) q[1];
sx q[1];
rz(-1.5086915) q[1];
sx q[1];
rz(-2.7746157) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94818753) q[0];
sx q[0];
rz(-0.087554878) q[0];
sx q[0];
rz(0.8669391) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.835145) q[2];
sx q[2];
rz(-1.2472787) q[2];
sx q[2];
rz(2.8217962) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6527892) q[1];
sx q[1];
rz(-1.3769883) q[1];
sx q[1];
rz(0.36313063) q[1];
x q[2];
rz(-2.7744297) q[3];
sx q[3];
rz(-1.490574) q[3];
sx q[3];
rz(1.0368766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7299812) q[2];
sx q[2];
rz(-1.9133762) q[2];
sx q[2];
rz(1.8709987) q[2];
rz(0.96108428) q[3];
sx q[3];
rz(-1.4189439) q[3];
sx q[3];
rz(0.050475033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5807895) q[0];
sx q[0];
rz(-2.2152948) q[0];
sx q[0];
rz(1.8286937) q[0];
rz(-1.1543697) q[1];
sx q[1];
rz(-1.1416953) q[1];
sx q[1];
rz(-0.55783522) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9168222) q[0];
sx q[0];
rz(-1.0616515) q[0];
sx q[0];
rz(2.1780685) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16061546) q[2];
sx q[2];
rz(-1.6470419) q[2];
sx q[2];
rz(-0.74197021) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8164929) q[1];
sx q[1];
rz(-1.8728647) q[1];
sx q[1];
rz(2.3343071) q[1];
rz(0.29739012) q[3];
sx q[3];
rz(-2.6618086) q[3];
sx q[3];
rz(2.9088717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4981726) q[2];
sx q[2];
rz(-1.8489445) q[2];
sx q[2];
rz(-0.84929973) q[2];
rz(-2.9863206) q[3];
sx q[3];
rz(-2.1657491) q[3];
sx q[3];
rz(2.7707905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3575386) q[0];
sx q[0];
rz(-0.58670601) q[0];
sx q[0];
rz(-0.46911711) q[0];
rz(0.47438374) q[1];
sx q[1];
rz(-2.3553039) q[1];
sx q[1];
rz(-0.75538409) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3123388) q[0];
sx q[0];
rz(-2.8702998) q[0];
sx q[0];
rz(-1.5872699) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9017436) q[2];
sx q[2];
rz(-0.96815434) q[2];
sx q[2];
rz(1.2801666) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1052221) q[1];
sx q[1];
rz(-1.3579988) q[1];
sx q[1];
rz(3.0119621) q[1];
x q[2];
rz(0.83306082) q[3];
sx q[3];
rz(-1.2445881) q[3];
sx q[3];
rz(1.8859832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0581806) q[2];
sx q[2];
rz(-1.3231134) q[2];
sx q[2];
rz(-0.18784909) q[2];
rz(-1.2457054) q[3];
sx q[3];
rz(-1.3789504) q[3];
sx q[3];
rz(-2.0511621) q[3];
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
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.93194) q[0];
sx q[0];
rz(-1.2670452) q[0];
sx q[0];
rz(-0.39091045) q[0];
rz(0.93005013) q[1];
sx q[1];
rz(-1.7101733) q[1];
sx q[1];
rz(-1.3040868) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5374889) q[0];
sx q[0];
rz(-1.2199645) q[0];
sx q[0];
rz(3.119191) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6565336) q[2];
sx q[2];
rz(-1.2303599) q[2];
sx q[2];
rz(-0.96735937) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3871349) q[1];
sx q[1];
rz(-1.4031193) q[1];
sx q[1];
rz(2.0281591) q[1];
rz(2.845302) q[3];
sx q[3];
rz(-1.5888831) q[3];
sx q[3];
rz(0.91218218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8730674) q[2];
sx q[2];
rz(-0.26669845) q[2];
sx q[2];
rz(-0.11338691) q[2];
rz(3.125627) q[3];
sx q[3];
rz(-1.4957875) q[3];
sx q[3];
rz(-0.16460831) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3575344) q[0];
sx q[0];
rz(-1.4736195) q[0];
sx q[0];
rz(3.0175324) q[0];
rz(0.1217753) q[1];
sx q[1];
rz(-0.75141326) q[1];
sx q[1];
rz(1.6990936) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0589813) q[0];
sx q[0];
rz(-1.1878723) q[0];
sx q[0];
rz(0.10848606) q[0];
rz(-2.997918) q[2];
sx q[2];
rz(-2.8468003) q[2];
sx q[2];
rz(0.93570342) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.6119163) q[1];
sx q[1];
rz(-1.9135336) q[1];
sx q[1];
rz(0.14118282) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8201039) q[3];
sx q[3];
rz(-2.2852906) q[3];
sx q[3];
rz(2.6282603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.96559912) q[2];
sx q[2];
rz(-1.3600574) q[2];
sx q[2];
rz(-0.74481258) q[2];
rz(-0.92721573) q[3];
sx q[3];
rz(-1.5085446) q[3];
sx q[3];
rz(-2.0492699) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48831478) q[0];
sx q[0];
rz(-1.4396242) q[0];
sx q[0];
rz(2.9162245) q[0];
rz(-2.0291746) q[1];
sx q[1];
rz(-0.77443361) q[1];
sx q[1];
rz(0.89881277) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8169237) q[0];
sx q[0];
rz(-2.3153439) q[0];
sx q[0];
rz(0.81314317) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1583011) q[2];
sx q[2];
rz(-2.0480683) q[2];
sx q[2];
rz(2.0606747) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4087993) q[1];
sx q[1];
rz(-2.2299754) q[1];
sx q[1];
rz(-0.2562457) q[1];
rz(-1.2779253) q[3];
sx q[3];
rz(-1.5808269) q[3];
sx q[3];
rz(-2.9603279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.044346873) q[2];
sx q[2];
rz(-1.711742) q[2];
sx q[2];
rz(-0.35935768) q[2];
rz(0.25092956) q[3];
sx q[3];
rz(-2.0693306) q[3];
sx q[3];
rz(-2.3738764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.3971685) q[0];
sx q[0];
rz(-0.1305307) q[0];
sx q[0];
rz(0.8771483) q[0];
rz(1.5769222) q[1];
sx q[1];
rz(-1.6404057) q[1];
sx q[1];
rz(-2.3559779) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4920542) q[0];
sx q[0];
rz(-1.7888594) q[0];
sx q[0];
rz(-1.3278924) q[0];
rz(-0.65746515) q[2];
sx q[2];
rz(-0.6415002) q[2];
sx q[2];
rz(2.7335707) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4686376) q[1];
sx q[1];
rz(-1.6412927) q[1];
sx q[1];
rz(-2.7967909) q[1];
rz(-pi) q[2];
rz(2.8749171) q[3];
sx q[3];
rz(-0.52191041) q[3];
sx q[3];
rz(-0.44346186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2239573) q[2];
sx q[2];
rz(-0.79006299) q[2];
sx q[2];
rz(-2.4033974) q[2];
rz(-2.2186642) q[3];
sx q[3];
rz(-1.3083369) q[3];
sx q[3];
rz(0.2963399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34687635) q[0];
sx q[0];
rz(-1.3608169) q[0];
sx q[0];
rz(0.97186744) q[0];
rz(2.234266) q[1];
sx q[1];
rz(-2.0569888) q[1];
sx q[1];
rz(2.8203698) q[1];
rz(2.1512866) q[2];
sx q[2];
rz(-1.2377501) q[2];
sx q[2];
rz(-0.35828423) q[2];
rz(-0.54013822) q[3];
sx q[3];
rz(-2.028699) q[3];
sx q[3];
rz(1.5069458) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
