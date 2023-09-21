OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5972714) q[0];
sx q[0];
rz(-0.40364021) q[0];
sx q[0];
rz(2.7713452) q[0];
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1365294) q[0];
sx q[0];
rz(-2.8423474) q[0];
sx q[0];
rz(2.5314999) q[0];
rz(0.45067694) q[2];
sx q[2];
rz(-0.55826954) q[2];
sx q[2];
rz(-2.6435341) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8542006) q[1];
sx q[1];
rz(-0.55254793) q[1];
sx q[1];
rz(1.470675) q[1];
rz(0.3317823) q[3];
sx q[3];
rz(-0.79091573) q[3];
sx q[3];
rz(3.0856109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7449164) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(2.7837226) q[2];
rz(0.19168028) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0033922694) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(3.0733118) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(-2.4904747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8296482) q[0];
sx q[0];
rz(-0.94046578) q[0];
sx q[0];
rz(-1.8207361) q[0];
rz(-0.5257734) q[2];
sx q[2];
rz(-1.5500881) q[2];
sx q[2];
rz(-2.6244342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.880289) q[1];
sx q[1];
rz(-2.0118879) q[1];
sx q[1];
rz(2.1217568) q[1];
rz(-pi) q[2];
x q[2];
rz(0.49115527) q[3];
sx q[3];
rz(-1.7633811) q[3];
sx q[3];
rz(1.1112978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(0.90332705) q[2];
rz(-0.7545169) q[3];
sx q[3];
rz(-1.6468331) q[3];
sx q[3];
rz(0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
x q[1];
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
rz(1.0292056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27965901) q[0];
sx q[0];
rz(-2.5760981) q[0];
sx q[0];
rz(0.66753597) q[0];
rz(-pi) q[1];
rz(0.35691397) q[2];
sx q[2];
rz(-1.1733574) q[2];
sx q[2];
rz(-0.39427653) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.88941105) q[1];
sx q[1];
rz(-1.2671789) q[1];
sx q[1];
rz(1.3510515) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0164099) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(2.560905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9323953) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(-1.8445245) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(-0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2066752) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-0.40859616) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(-0.88159195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5170167) q[0];
sx q[0];
rz(-1.7424184) q[0];
sx q[0];
rz(0.73715985) q[0];
rz(-0.84882952) q[2];
sx q[2];
rz(-1.0366882) q[2];
sx q[2];
rz(1.7497077) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9095163) q[1];
sx q[1];
rz(-2.0627923) q[1];
sx q[1];
rz(2.3823649) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6374171) q[3];
sx q[3];
rz(-1.9873709) q[3];
sx q[3];
rz(1.6601603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0980229) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(-1.1126474) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.6120733) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(-1.7768815) q[0];
rz(2.8240906) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(3.024335) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011615959) q[0];
sx q[0];
rz(-2.5419309) q[0];
sx q[0];
rz(0.2325124) q[0];
rz(-1.7848894) q[2];
sx q[2];
rz(-2.1630686) q[2];
sx q[2];
rz(-1.5562197) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9980364) q[1];
sx q[1];
rz(-1.4170425) q[1];
sx q[1];
rz(2.9424332) q[1];
x q[2];
rz(-0.98362555) q[3];
sx q[3];
rz(-1.5067325) q[3];
sx q[3];
rz(0.6928882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7964898) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(1.7618746) q[2];
rz(1.1896677) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18774408) q[0];
sx q[0];
rz(-1.6402316) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(1.1114936) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(0.10791735) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7100922) q[0];
sx q[0];
rz(-1.2762696) q[0];
sx q[0];
rz(-0.72580238) q[0];
rz(-pi) q[1];
x q[1];
rz(0.040060476) q[2];
sx q[2];
rz(-1.5594348) q[2];
sx q[2];
rz(0.60919112) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.38938746) q[1];
sx q[1];
rz(-2.2332193) q[1];
sx q[1];
rz(-1.2090769) q[1];
rz(-pi) q[2];
rz(-2.2599254) q[3];
sx q[3];
rz(-1.1881184) q[3];
sx q[3];
rz(1.6430287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-2.0325913) q[2];
rz(-1.8479944) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-2.368822) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(-3.1266881) q[0];
rz(-2.7203454) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(-0.79963911) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22916238) q[0];
sx q[0];
rz(-1.3161356) q[0];
sx q[0];
rz(2.4486662) q[0];
rz(-pi) q[1];
rz(-1.2866576) q[2];
sx q[2];
rz(-1.8737027) q[2];
sx q[2];
rz(-2.6877833) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0529424) q[1];
sx q[1];
rz(-0.155641) q[1];
sx q[1];
rz(2.774653) q[1];
x q[2];
rz(-1.0054672) q[3];
sx q[3];
rz(-1.4544832) q[3];
sx q[3];
rz(-0.034193589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(0.9220534) q[2];
rz(-1.3098035) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9668982) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(-2.877537) q[0];
rz(-1.4008201) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(2.82428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4357294) q[0];
sx q[0];
rz(-1.9575435) q[0];
sx q[0];
rz(3.0630037) q[0];
rz(-pi) q[1];
rz(1.9192341) q[2];
sx q[2];
rz(-1.7218044) q[2];
sx q[2];
rz(2.3210971) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.322284) q[1];
sx q[1];
rz(-1.6849736) q[1];
sx q[1];
rz(-1.6627922) q[1];
rz(-0.64538892) q[3];
sx q[3];
rz(-1.4374975) q[3];
sx q[3];
rz(-0.79750878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.6395817) q[2];
rz(-0.25029415) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(1.7171575) q[0];
rz(0.4793438) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-0.11553484) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53287017) q[0];
sx q[0];
rz(-1.0057797) q[0];
sx q[0];
rz(0.67824058) q[0];
rz(-pi) q[1];
rz(-1.2938415) q[2];
sx q[2];
rz(-0.87522725) q[2];
sx q[2];
rz(-2.9447174) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.1852473) q[1];
sx q[1];
rz(-2.6223409) q[1];
sx q[1];
rz(-0.96038702) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4453836) q[3];
sx q[3];
rz(-1.4783825) q[3];
sx q[3];
rz(0.19314167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66403786) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.4257365) q[2];
rz(-1.7410295) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7857159) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(0.26915959) q[0];
rz(-1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(-1.7620618) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3740765) q[0];
sx q[0];
rz(-0.19664581) q[0];
sx q[0];
rz(-2.1392512) q[0];
x q[1];
rz(-0.079897957) q[2];
sx q[2];
rz(-1.986984) q[2];
sx q[2];
rz(-1.4518567) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.60296842) q[1];
sx q[1];
rz(-1.6011964) q[1];
sx q[1];
rz(1.3974742) q[1];
x q[2];
rz(-1.8809324) q[3];
sx q[3];
rz(-1.2536067) q[3];
sx q[3];
rz(-0.96084259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.049008869) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(-1.6258378) q[2];
rz(1.9745291) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(-2.1102171) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2172858) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(-2.6869607) q[1];
sx q[1];
rz(-2.0352719) q[1];
sx q[1];
rz(-0.24771053) q[1];
rz(0.17299962) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
rz(-1.3439988) q[3];
sx q[3];
rz(-2.8759225) q[3];
sx q[3];
rz(0.21951036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];