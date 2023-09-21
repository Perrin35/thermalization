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
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6366633) q[0];
sx q[0];
rz(-1.3267656) q[0];
sx q[0];
rz(1.3958449) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3051891) q[2];
sx q[2];
rz(-1.0737537) q[2];
sx q[2];
rz(0.020393919) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4048684) q[1];
sx q[1];
rz(-1.0213335) q[1];
sx q[1];
rz(0.0615555) q[1];
rz(-pi) q[2];
rz(2.3787093) q[3];
sx q[3];
rz(-1.3370822) q[3];
sx q[3];
rz(1.8644615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(-2.7837226) q[2];
rz(2.9499124) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(-0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1382004) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(3.0733118) q[0];
rz(-2.1349019) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(0.65111792) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90322813) q[0];
sx q[0];
rz(-2.4698386) q[0];
sx q[0];
rz(-0.3268468) q[0];
rz(-2.6158193) q[2];
sx q[2];
rz(-1.5500881) q[2];
sx q[2];
rz(2.6244342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.843833) q[1];
sx q[1];
rz(-0.69121541) q[1];
sx q[1];
rz(-0.83696951) q[1];
x q[2];
rz(-1.7884368) q[3];
sx q[3];
rz(-1.0895035) q[3];
sx q[3];
rz(0.35748112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(-0.90332705) q[2];
rz(0.7545169) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(2.6480411) q[0];
rz(-1.6429398) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-2.1123871) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.670186) q[0];
sx q[0];
rz(-1.1364511) q[0];
sx q[0];
rz(1.1964519) q[0];
rz(-2.2649293) q[2];
sx q[2];
rz(-0.52774094) q[2];
sx q[2];
rz(-1.1610247) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.88941105) q[1];
sx q[1];
rz(-1.2671789) q[1];
sx q[1];
rz(-1.7905411) q[1];
x q[2];
rz(-0.1251827) q[3];
sx q[3];
rz(-2.1117961) q[3];
sx q[3];
rz(-0.5806877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9323953) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(1.2970682) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2066752) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(0.40859616) q[0];
rz(-1.4273377) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-0.88159195) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0094385) q[0];
sx q[0];
rz(-0.75320019) q[0];
sx q[0];
rz(0.25235812) q[0];
rz(-pi) q[1];
rz(-2.4741715) q[2];
sx q[2];
rz(-2.1758658) q[2];
sx q[2];
rz(2.8990926) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8005912) q[1];
sx q[1];
rz(-0.87716555) q[1];
sx q[1];
rz(0.6615521) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9922585) q[3];
sx q[3];
rz(-0.42156005) q[3];
sx q[3];
rz(-1.4967407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0435698) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(-2.0289452) q[2];
rz(-0.55551314) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.5295193) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9920138) q[0];
sx q[0];
rz(-1.9861789) q[0];
sx q[0];
rz(-1.3647112) q[0];
rz(0.31750202) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(3.024335) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7754585) q[0];
sx q[0];
rz(-1.701208) q[0];
sx q[0];
rz(2.5545757) q[0];
rz(0.30585441) q[2];
sx q[2];
rz(-0.62539414) q[2];
sx q[2];
rz(1.9276227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14355625) q[1];
sx q[1];
rz(-1.4170425) q[1];
sx q[1];
rz(-2.9424332) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4555143) q[3];
sx q[3];
rz(-2.5513463) q[3];
sx q[3];
rz(2.3595927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7964898) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(-1.3797181) q[2];
rz(-1.951925) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9538486) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(-0.10791735) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2544125) q[0];
sx q[0];
rz(-2.2590056) q[0];
sx q[0];
rz(1.956091) q[0];
rz(-1.5594257) q[2];
sx q[2];
rz(-1.6108542) q[2];
sx q[2];
rz(-0.96206059) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.94090998) q[1];
sx q[1];
rz(-0.74146491) q[1];
sx q[1];
rz(0.42592589) q[1];
rz(2.6608653) q[3];
sx q[3];
rz(-2.2016964) q[3];
sx q[3];
rz(0.37068278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8191021) q[2];
sx q[2];
rz(-0.97444797) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368822) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(3.1266881) q[0];
rz(-0.42124721) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(0.79963911) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5478202) q[0];
sx q[0];
rz(-0.9043588) q[0];
sx q[0];
rz(-1.244546) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4104426) q[2];
sx q[2];
rz(-2.7292948) q[2];
sx q[2];
rz(1.9129802) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7176197) q[1];
sx q[1];
rz(-1.4255925) q[1];
sx q[1];
rz(-1.6270301) q[1];
x q[2];
rz(-1.3560489) q[3];
sx q[3];
rz(-0.575892) q[3];
sx q[3];
rz(1.7175331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(-0.9220534) q[2];
rz(-1.3098035) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(2.183765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17469445) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(2.877537) q[0];
rz(1.4008201) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(2.82428) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70586328) q[0];
sx q[0];
rz(-1.9575435) q[0];
sx q[0];
rz(-3.0630037) q[0];
rz(2.9810901) q[2];
sx q[2];
rz(-1.2264894) q[2];
sx q[2];
rz(-0.80489327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.74097733) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(0.11465794) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64538892) q[3];
sx q[3];
rz(-1.7040952) q[3];
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
rz(-1.7264265) q[3];
sx q[3];
rz(1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89649993) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(1.4244351) q[0];
rz(0.4793438) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(-3.0260578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4452267) q[0];
sx q[0];
rz(-1.0123024) q[0];
sx q[0];
rz(0.88748705) q[0];
rz(1.8477511) q[2];
sx q[2];
rz(-2.2663654) q[2];
sx q[2];
rz(-0.19687523) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.301831) q[1];
sx q[1];
rz(-1.8592195) q[1];
sx q[1];
rz(2.0088058) q[1];
rz(-pi) q[2];
rz(1.6910016) q[3];
sx q[3];
rz(-0.87814858) q[3];
sx q[3];
rz(1.8409178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.66403786) q[2];
sx q[2];
rz(-2.0319735) q[2];
sx q[2];
rz(1.4257365) q[2];
rz(1.4005631) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(-0.28276309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7857159) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(2.8724331) q[0];
rz(-1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(-1.7620618) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1901967) q[0];
sx q[0];
rz(-1.7362036) q[0];
sx q[0];
rz(0.10683807) q[0];
rz(-pi) q[1];
rz(-3.0616947) q[2];
sx q[2];
rz(-1.986984) q[2];
sx q[2];
rz(-1.6897359) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1684432) q[1];
sx q[1];
rz(-1.7440376) q[1];
sx q[1];
rz(-3.1107305) q[1];
x q[2];
rz(1.2606603) q[3];
sx q[3];
rz(-1.2536067) q[3];
sx q[3];
rz(2.1807501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(1.5157549) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9243069) q[0];
sx q[0];
rz(-1.8368245) q[0];
sx q[0];
rz(0.40689847) q[0];
rz(2.6869607) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(-1.8998014) q[2];
sx q[2];
rz(-0.49578373) q[2];
sx q[2];
rz(1.0052581) q[2];
rz(-1.8299673) q[3];
sx q[3];
rz(-1.5117241) q[3];
sx q[3];
rz(-1.1321887) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
