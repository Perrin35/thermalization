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
rz(-0.37024745) q[0];
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(1.7226146) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1365294) q[0];
sx q[0];
rz(-2.8423474) q[0];
sx q[0];
rz(-2.5314999) q[0];
rz(1.3051891) q[2];
sx q[2];
rz(-2.0678389) q[2];
sx q[2];
rz(0.020393919) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4048684) q[1];
sx q[1];
rz(-2.1202592) q[1];
sx q[1];
rz(3.0800372) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2526413) q[3];
sx q[3];
rz(-0.83359026) q[3];
sx q[3];
rz(-2.6300501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7449164) q[2];
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
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382004) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-3.0733118) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(-0.65111792) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40819528) q[0];
sx q[0];
rz(-1.3696284) q[0];
sx q[0];
rz(-2.4961619) q[0];
rz(-pi) q[1];
rz(-1.5947371) q[2];
sx q[2];
rz(-1.0451473) q[2];
sx q[2];
rz(1.0416232) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.843833) q[1];
sx q[1];
rz(-2.4503772) q[1];
sx q[1];
rz(2.3046231) q[1];
rz(-1.3531559) q[3];
sx q[3];
rz(-2.0520891) q[3];
sx q[3];
rz(-2.7841115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.45289257) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(0.90332705) q[2];
rz(0.7545169) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-0.49355155) q[0];
rz(-1.6429398) q[1];
sx q[1];
rz(-0.40619266) q[1];
sx q[1];
rz(2.1123871) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4714067) q[0];
sx q[0];
rz(-2.0051415) q[0];
sx q[0];
rz(-1.9451408) q[0];
x q[1];
rz(1.9919954) q[2];
sx q[2];
rz(-1.8987978) q[2];
sx q[2];
rz(-2.1084107) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.61470795) q[1];
sx q[1];
rz(-1.3612559) q[1];
sx q[1];
rz(0.31063147) q[1];
rz(-2.1152705) q[3];
sx q[3];
rz(-1.6780276) q[3];
sx q[3];
rz(-0.925392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(1.8445245) q[2];
rz(-0.57859892) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(-2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2066752) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(2.7329965) q[0];
rz(1.714255) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(-2.2600007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9339869) q[0];
sx q[0];
rz(-0.84689492) q[0];
sx q[0];
rz(-1.8007604) q[0];
rz(-pi) q[1];
rz(0.84882952) q[2];
sx q[2];
rz(-2.1049044) q[2];
sx q[2];
rz(-1.3918849) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.2320764) q[1];
sx q[1];
rz(-2.0627923) q[1];
sx q[1];
rz(0.75922774) q[1];
rz(-pi) q[2];
rz(2.7241957) q[3];
sx q[3];
rz(-1.5098803) q[3];
sx q[3];
rz(0.062373769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0980229) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(2.0289452) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9920138) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(-0.31750202) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(3.024335) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8506354) q[0];
sx q[0];
rz(-2.1521749) q[0];
sx q[0];
rz(1.7270389) q[0];
rz(-pi) q[1];
rz(-1.3567032) q[2];
sx q[2];
rz(-2.1630686) q[2];
sx q[2];
rz(1.5562197) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.3634503) q[1];
sx q[1];
rz(-0.25097978) q[1];
sx q[1];
rz(-2.4771033) q[1];
rz(-pi) q[2];
rz(1.4555143) q[3];
sx q[3];
rz(-0.59024631) q[3];
sx q[3];
rz(0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(1.3797181) q[2];
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
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(2.4353819) q[0];
rz(2.0300991) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(-0.10791735) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7100922) q[0];
sx q[0];
rz(-1.8653231) q[0];
sx q[0];
rz(0.72580238) q[0];
rz(-1.582167) q[2];
sx q[2];
rz(-1.6108542) q[2];
sx q[2];
rz(-2.1795321) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1888215) q[1];
sx q[1];
rz(-1.8535887) q[1];
sx q[1];
rz(-2.4464843) q[1];
x q[2];
rz(0.48072731) q[3];
sx q[3];
rz(-2.2016964) q[3];
sx q[3];
rz(2.7709099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8191021) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(-2.0325913) q[2];
rz(1.8479944) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7727707) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(3.1266881) q[0];
rz(-2.7203454) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(0.79963911) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0470456) q[0];
sx q[0];
rz(-2.4106815) q[0];
sx q[0];
rz(0.38696179) q[0];
x q[1];
rz(-0.31475474) q[2];
sx q[2];
rz(-1.2999279) q[2];
sx q[2];
rz(1.0300919) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0529424) q[1];
sx q[1];
rz(-0.155641) q[1];
sx q[1];
rz(0.36693962) q[1];
x q[2];
rz(1.3560489) q[3];
sx q[3];
rz(-2.5657006) q[3];
sx q[3];
rz(1.7175331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(2.2195393) q[2];
rz(1.8317892) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(-0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17469445) q[0];
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
rz(0.50002977) q[0];
sx q[0];
rz(-2.7473358) q[0];
sx q[0];
rz(1.7612329) q[0];
rz(1.1515456) q[2];
sx q[2];
rz(-0.37852415) q[2];
sx q[2];
rz(0.35767698) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.81930868) q[1];
sx q[1];
rz(-1.6849736) q[1];
sx q[1];
rz(-1.6627922) q[1];
rz(-pi) q[2];
rz(-2.92225) q[3];
sx q[3];
rz(-0.6570802) q[3];
sx q[3];
rz(-0.59857644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43549609) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(0.25029415) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450927) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(-1.4244351) q[0];
rz(0.4793438) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(0.11553484) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45142052) q[0];
sx q[0];
rz(-0.85309404) q[0];
sx q[0];
rz(2.3510128) q[0];
rz(-pi) q[1];
rz(2.4268482) q[2];
sx q[2];
rz(-1.359316) q[2];
sx q[2];
rz(-1.5541058) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2781196) q[1];
sx q[1];
rz(-1.1520471) q[1];
sx q[1];
rz(2.8249885) q[1];
rz(2.9980738) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(-1.4876175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-2.0319735) q[2];
sx q[2];
rz(1.7158562) q[2];
rz(-1.7410295) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(-0.28276309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7857159) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(-2.8724331) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(1.3795308) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7675161) q[0];
sx q[0];
rz(-2.9449468) q[0];
sx q[0];
rz(2.1392512) q[0];
rz(-pi) q[1];
rz(-3.0616947) q[2];
sx q[2];
rz(-1.1546087) q[2];
sx q[2];
rz(-1.4518567) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5386242) q[1];
sx q[1];
rz(-1.5403962) q[1];
sx q[1];
rz(-1.7441185) q[1];
rz(2.3926211) q[3];
sx q[3];
rz(-0.43991551) q[3];
sx q[3];
rz(0.16187748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(1.6258378) q[2];
rz(1.1670636) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(-1.0313755) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9243069) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(0.45463195) q[1];
sx q[1];
rz(-2.0352719) q[1];
sx q[1];
rz(-0.24771053) q[1];
rz(-2.968593) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
rz(-3.0804844) q[3];
sx q[3];
rz(-1.3120878) q[3];
sx q[3];
rz(0.45425934) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
