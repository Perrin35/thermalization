OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5443213) q[0];
sx q[0];
rz(3.5452329) q[0];
sx q[0];
rz(9.7950254) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(-1.4189781) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0231853) q[0];
sx q[0];
rz(-1.7405131) q[0];
sx q[0];
rz(-0.24766185) q[0];
x q[1];
rz(-0.45067694) q[2];
sx q[2];
rz(-2.5833231) q[2];
sx q[2];
rz(-2.6435341) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8542006) q[1];
sx q[1];
rz(-2.5890447) q[1];
sx q[1];
rz(1.6709177) q[1];
rz(1.8889514) q[3];
sx q[3];
rz(-2.3080024) q[3];
sx q[3];
rz(-0.51154256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7449164) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(-0.35787004) q[2];
rz(0.19168028) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382004) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(-0.068280846) q[0];
rz(-2.1349019) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(2.4904747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90322813) q[0];
sx q[0];
rz(-2.4698386) q[0];
sx q[0];
rz(-2.8147459) q[0];
x q[1];
rz(2.6158193) q[2];
sx q[2];
rz(-1.5915046) q[2];
sx q[2];
rz(2.6244342) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5755641) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(0.50599392) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49115527) q[3];
sx q[3];
rz(-1.3782116) q[3];
sx q[3];
rz(2.0302949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(2.2382656) q[2];
rz(-2.3870758) q[3];
sx q[3];
rz(-1.6468331) q[3];
sx q[3];
rz(-0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2925401) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-2.1123871) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2632336) q[0];
sx q[0];
rz(-1.2326816) q[0];
sx q[0];
rz(2.6792206) q[0];
rz(0.35691397) q[2];
sx q[2];
rz(-1.1733574) q[2];
sx q[2];
rz(-0.39427653) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5268847) q[1];
sx q[1];
rz(-1.3612559) q[1];
sx q[1];
rz(0.31063147) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1251827) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(-0.5806877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2091973) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(-2.5629937) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349174) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(2.7329965) q[0];
rz(1.714255) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(0.88159195) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5170167) q[0];
sx q[0];
rz(-1.3991742) q[0];
sx q[0];
rz(-2.4044328) q[0];
rz(-pi) q[1];
rz(-0.84076968) q[2];
sx q[2];
rz(-2.2731014) q[2];
sx q[2];
rz(2.4384987) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2243005) q[1];
sx q[1];
rz(-0.91887337) q[1];
sx q[1];
rz(0.93445458) q[1];
rz(-pi) q[2];
rz(0.14933417) q[3];
sx q[3];
rz(-2.7200326) q[3];
sx q[3];
rz(1.4967407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0435698) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(-2.0289452) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(-1.5295193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1495789) q[0];
sx q[0];
rz(-1.9861789) q[0];
sx q[0];
rz(-1.3647112) q[0];
rz(0.31750202) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(-3.024335) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011615959) q[0];
sx q[0];
rz(-0.59966171) q[0];
sx q[0];
rz(2.9090803) q[0];
rz(-2.8357382) q[2];
sx q[2];
rz(-0.62539414) q[2];
sx q[2];
rz(-1.2139699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.14355625) q[1];
sx q[1];
rz(-1.7245502) q[1];
sx q[1];
rz(-2.9424332) q[1];
x q[2];
rz(-3.0646867) q[3];
sx q[3];
rz(-0.98499005) q[3];
sx q[3];
rz(-2.2211071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7964898) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(1.7618746) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538486) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(2.4353819) q[0];
rz(-2.0300991) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(0.10791735) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4554493) q[0];
sx q[0];
rz(-2.3685072) q[0];
sx q[0];
rz(-0.42868741) q[0];
rz(-pi) q[1];
rz(-0.27643398) q[2];
sx q[2];
rz(-0.041639608) q[2];
sx q[2];
rz(0.68539884) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2006827) q[1];
sx q[1];
rz(-0.74146491) q[1];
sx q[1];
rz(0.42592589) q[1];
x q[2];
rz(0.88166724) q[3];
sx q[3];
rz(-1.9534743) q[3];
sx q[3];
rz(-1.6430287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8191021) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-1.1090013) q[2];
rz(-1.2935982) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(-3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7727707) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(0.014904508) q[0];
rz(0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(0.79963911) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945471) q[0];
sx q[0];
rz(-2.4106815) q[0];
sx q[0];
rz(-0.38696179) q[0];
x q[1];
rz(2.4104426) q[2];
sx q[2];
rz(-0.41229782) q[2];
sx q[2];
rz(1.9129802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9866242) q[1];
sx q[1];
rz(-1.5151549) q[1];
sx q[1];
rz(0.14543047) q[1];
rz(-pi) q[2];
rz(1.3560489) q[3];
sx q[3];
rz(-0.575892) q[3];
sx q[3];
rz(1.4240595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5853167) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(0.9220534) q[2];
rz(-1.8317892) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17469445) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(2.877537) q[0];
rz(-1.7407725) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(-0.31731269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50002977) q[0];
sx q[0];
rz(-0.39425685) q[0];
sx q[0];
rz(-1.3803598) q[0];
rz(1.9192341) q[2];
sx q[2];
rz(-1.4197883) q[2];
sx q[2];
rz(-2.3210971) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4006153) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(-3.0269347) q[1];
rz(-pi) q[2];
rz(2.92225) q[3];
sx q[3];
rz(-2.4845124) q[3];
sx q[3];
rz(-0.59857644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.6395817) q[2];
rz(0.25029415) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(-1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.2450927) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(1.4244351) q[0];
rz(-2.6622488) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(0.11553484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6087225) q[0];
sx q[0];
rz(-1.0057797) q[0];
sx q[0];
rz(0.67824058) q[0];
rz(-pi) q[1];
x q[1];
rz(0.3165452) q[2];
sx q[2];
rz(-2.4015421) q[2];
sx q[2];
rz(-0.22063247) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.301831) q[1];
sx q[1];
rz(-1.2823732) q[1];
sx q[1];
rz(1.1327869) q[1];
rz(-pi) q[2];
rz(-0.1435189) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(1.6539751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-2.0319735) q[2];
sx q[2];
rz(-1.4257365) q[2];
rz(1.4005631) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(2.8588296) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(2.8724331) q[0];
rz(-1.0956988) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(-1.3795308) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9513959) q[0];
sx q[0];
rz(-1.405389) q[0];
sx q[0];
rz(0.10683807) q[0];
rz(-pi) q[1];
rz(-1.1534259) q[2];
sx q[2];
rz(-1.4977314) q[2];
sx q[2];
rz(3.0550115) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3456612) q[1];
sx q[1];
rz(-2.965651) q[1];
sx q[1];
rz(-1.3962586) q[1];
rz(-pi) q[2];
rz(-1.8809324) q[3];
sx q[3];
rz(-1.2536067) q[3];
sx q[3];
rz(2.1807501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.049008869) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(1.5157549) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.2172858) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(-0.45463195) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(2.0438647) q[2];
sx q[2];
rz(-1.7251144) q[2];
sx q[2];
rz(-0.27380064) q[2];
rz(-0.061108246) q[3];
sx q[3];
rz(-1.8295049) q[3];
sx q[3];
rz(-2.6873333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
