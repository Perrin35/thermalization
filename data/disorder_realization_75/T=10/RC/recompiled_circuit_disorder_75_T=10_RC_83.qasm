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
rz(1.8887853) q[1];
sx q[1];
rz(10.843756) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1184074) q[0];
sx q[0];
rz(-1.7405131) q[0];
sx q[0];
rz(-2.8939308) q[0];
rz(-pi) q[1];
x q[1];
rz(0.51214829) q[2];
sx q[2];
rz(-1.8036267) q[2];
sx q[2];
rz(1.4621967) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7367243) q[1];
sx q[1];
rz(-2.1202592) q[1];
sx q[1];
rz(3.0800372) q[1];
x q[2];
rz(-1.2526413) q[3];
sx q[3];
rz(-0.83359026) q[3];
sx q[3];
rz(0.51154256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(-0.35787004) q[2];
rz(2.9499124) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(2.5884957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382004) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(3.0733118) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(0.65111792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2383645) q[0];
sx q[0];
rz(-0.67175409) q[0];
sx q[0];
rz(0.3268468) q[0];
rz(-1.5947371) q[2];
sx q[2];
rz(-1.0451473) q[2];
sx q[2];
rz(-2.0999694) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2613036) q[1];
sx q[1];
rz(-2.0118879) q[1];
sx q[1];
rz(-2.1217568) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7495456) q[3];
sx q[3];
rz(-0.52467504) q[3];
sx q[3];
rz(-0.80313659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(0.90332705) q[2];
rz(-0.7545169) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(0.49355155) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-2.1123871) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4714067) q[0];
sx q[0];
rz(-1.1364511) q[0];
sx q[0];
rz(-1.1964519) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87666338) q[2];
sx q[2];
rz(-2.6138517) q[2];
sx q[2];
rz(1.1610247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5268847) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(2.8309612) q[1];
x q[2];
rz(-1.0263222) q[3];
sx q[3];
rz(-1.4635651) q[3];
sx q[3];
rz(-0.925392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2091973) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(-1.8445245) q[2];
rz(-0.57859892) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(-0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(1.9349174) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(0.40859616) q[0];
rz(1.714255) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(-2.2600007) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13215412) q[0];
sx q[0];
rz(-0.75320019) q[0];
sx q[0];
rz(-0.25235812) q[0];
rz(-pi) q[1];
rz(2.2927631) q[2];
sx q[2];
rz(-1.0366882) q[2];
sx q[2];
rz(1.7497077) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2320764) q[1];
sx q[1];
rz(-1.0788003) q[1];
sx q[1];
rz(-0.75922774) q[1];
rz(-pi) q[2];
rz(-2.9922585) q[3];
sx q[3];
rz(-2.7200326) q[3];
sx q[3];
rz(-1.644852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0435698) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(-1.1126474) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(0.31750202) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(0.11725765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8506354) q[0];
sx q[0];
rz(-0.98941776) q[0];
sx q[0];
rz(1.7270389) q[0];
rz(-pi) q[1];
rz(-2.8357382) q[2];
sx q[2];
rz(-0.62539414) q[2];
sx q[2];
rz(1.9276227) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77814233) q[1];
sx q[1];
rz(-0.25097978) q[1];
sx q[1];
rz(-0.66448934) q[1];
rz(2.1579671) q[3];
sx q[3];
rz(-1.5067325) q[3];
sx q[3];
rz(-2.4487045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7964898) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.3797181) q[2];
rz(1.1896677) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538486) q[0];
sx q[0];
rz(-1.6402316) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(2.0300991) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(-0.10791735) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2544125) q[0];
sx q[0];
rz(-0.88258703) q[0];
sx q[0];
rz(1.1855017) q[0];
x q[1];
rz(-1.5594257) q[2];
sx q[2];
rz(-1.6108542) q[2];
sx q[2];
rz(-0.96206059) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7522052) q[1];
sx q[1];
rz(-2.2332193) q[1];
sx q[1];
rz(1.2090769) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48072731) q[3];
sx q[3];
rz(-2.2016964) q[3];
sx q[3];
rz(-2.7709099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(2.0325913) q[2];
rz(-1.8479944) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(-0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7727707) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(3.1266881) q[0];
rz(2.7203454) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(-2.3419535) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124303) q[0];
sx q[0];
rz(-1.3161356) q[0];
sx q[0];
rz(2.4486662) q[0];
rz(-pi) q[1];
rz(1.854935) q[2];
sx q[2];
rz(-1.8737027) q[2];
sx q[2];
rz(0.45380935) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0886503) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(-0.36693962) q[1];
rz(3.0040967) q[3];
sx q[3];
rz(-1.0097479) q[3];
sx q[3];
rz(1.4631127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5853167) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(2.2195393) q[2];
rz(-1.3098035) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(2.183765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17469445) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(-0.2640557) q[0];
rz(1.7407725) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(2.82428) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6415629) q[0];
sx q[0];
rz(-0.39425685) q[0];
sx q[0];
rz(1.7612329) q[0];
rz(-pi) q[1];
rz(-0.1605026) q[2];
sx q[2];
rz(-1.9151033) q[2];
sx q[2];
rz(0.80489327) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.74097733) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(3.0269347) q[1];
rz(-pi) q[2];
rz(-2.92225) q[3];
sx q[3];
rz(-2.4845124) q[3];
sx q[3];
rz(-2.5430162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7060966) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(1.502011) q[2];
rz(2.8912985) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(-1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89649993) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(1.4244351) q[0];
rz(-0.4793438) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(-0.11553484) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45142052) q[0];
sx q[0];
rz(-0.85309404) q[0];
sx q[0];
rz(-0.79057981) q[0];
x q[1];
rz(-1.2938415) q[2];
sx q[2];
rz(-2.2663654) q[2];
sx q[2];
rz(2.9447174) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2781196) q[1];
sx q[1];
rz(-1.9895456) q[1];
sx q[1];
rz(-0.31660415) q[1];
rz(-pi) q[2];
rz(1.6910016) q[3];
sx q[3];
rz(-2.2634441) q[3];
sx q[3];
rz(1.3006749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-2.0319735) q[2];
sx q[2];
rz(-1.4257365) q[2];
rz(-1.4005631) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(-0.28276309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.7857159) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(0.26915959) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(1.7620618) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7675161) q[0];
sx q[0];
rz(-0.19664581) q[0];
sx q[0];
rz(2.1392512) q[0];
rz(-1.3921521) q[2];
sx q[2];
rz(-0.42334712) q[2];
sx q[2];
rz(1.64738) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.79593149) q[1];
sx q[1];
rz(-0.17594166) q[1];
sx q[1];
rz(-1.7453341) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8096301) q[3];
sx q[3];
rz(-1.276607) q[3];
sx q[3];
rz(2.4320137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0925838) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(-1.6258378) q[2];
rz(1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-2.1102171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(0.17299962) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
rz(3.0804844) q[3];
sx q[3];
rz(-1.8295049) q[3];
sx q[3];
rz(-2.6873333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];