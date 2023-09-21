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
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1184074) q[0];
sx q[0];
rz(-1.7405131) q[0];
sx q[0];
rz(0.24766185) q[0];
x q[1];
rz(0.45067694) q[2];
sx q[2];
rz(-2.5833231) q[2];
sx q[2];
rz(-0.4980586) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7367243) q[1];
sx q[1];
rz(-2.1202592) q[1];
sx q[1];
rz(3.0800372) q[1];
rz(-1.8889514) q[3];
sx q[3];
rz(-2.3080024) q[3];
sx q[3];
rz(0.51154256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7449164) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(0.35787004) q[2];
rz(2.9499124) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0033922694) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(0.068280846) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(-0.65111792) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7333974) q[0];
sx q[0];
rz(-1.3696284) q[0];
sx q[0];
rz(-2.4961619) q[0];
x q[1];
rz(-0.5257734) q[2];
sx q[2];
rz(-1.5500881) q[2];
sx q[2];
rz(0.5171585) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29775961) q[1];
sx q[1];
rz(-0.69121541) q[1];
sx q[1];
rz(-2.3046231) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7495456) q[3];
sx q[3];
rz(-0.52467504) q[3];
sx q[3];
rz(-2.3384561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6887001) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(0.90332705) q[2];
rz(-2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(-2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(-1.6429398) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(1.0292056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4714067) q[0];
sx q[0];
rz(-2.0051415) q[0];
sx q[0];
rz(1.1964519) q[0];
rz(-0.35691397) q[2];
sx q[2];
rz(-1.9682353) q[2];
sx q[2];
rz(-0.39427653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2521816) q[1];
sx q[1];
rz(-1.8744138) q[1];
sx q[1];
rz(1.7905411) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7756996) q[3];
sx q[3];
rz(-0.5538867) q[3];
sx q[3];
rz(-2.3212471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2091973) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.8445245) q[2];
rz(-0.57859892) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9349174) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(2.7329965) q[0];
rz(1.714255) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-0.88159195) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0094385) q[0];
sx q[0];
rz(-0.75320019) q[0];
sx q[0];
rz(-0.25235812) q[0];
rz(-pi) q[1];
rz(-0.84076968) q[2];
sx q[2];
rz(-2.2731014) q[2];
sx q[2];
rz(2.4384987) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.2320764) q[1];
sx q[1];
rz(-1.0788003) q[1];
sx q[1];
rz(2.3823649) q[1];
x q[2];
rz(1.5041755) q[3];
sx q[3];
rz(-1.9873709) q[3];
sx q[3];
rz(-1.4814324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0435698) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(1.1126474) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-2.1798539) q[1];
sx q[1];
rz(-3.024335) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011615959) q[0];
sx q[0];
rz(-0.59966171) q[0];
sx q[0];
rz(-0.2325124) q[0];
rz(-2.5385777) q[2];
sx q[2];
rz(-1.7479959) q[2];
sx q[2];
rz(-0.10620968) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.14355625) q[1];
sx q[1];
rz(-1.7245502) q[1];
sx q[1];
rz(2.9424332) q[1];
x q[2];
rz(0.98362555) q[3];
sx q[3];
rz(-1.5067325) q[3];
sx q[3];
rz(-0.6928882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(-1.7618746) q[2];
rz(1.951925) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(-3.0670847) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(2.4353819) q[0];
rz(1.1114936) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(0.10791735) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4554493) q[0];
sx q[0];
rz(-2.3685072) q[0];
sx q[0];
rz(-0.42868741) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5594257) q[2];
sx q[2];
rz(-1.5307384) q[2];
sx q[2];
rz(-0.96206059) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-pi) q[2];
x q[2];
rz(2.6608653) q[3];
sx q[3];
rz(-0.93989621) q[3];
sx q[3];
rz(-0.37068278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8191021) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(1.1090013) q[2];
rz(-1.2935982) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(0.7727707) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(-3.1266881) q[0];
rz(-0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(2.3419535) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22916238) q[0];
sx q[0];
rz(-1.3161356) q[0];
sx q[0];
rz(-0.69292648) q[0];
rz(-0.73115007) q[2];
sx q[2];
rz(-2.7292948) q[2];
sx q[2];
rz(-1.9129802) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7176197) q[1];
sx q[1];
rz(-1.7160001) q[1];
sx q[1];
rz(-1.6270301) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1361254) q[3];
sx q[3];
rz(-1.6871095) q[3];
sx q[3];
rz(-3.1073991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(2.2195393) q[2];
rz(1.3098035) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9668982) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(-0.2640557) q[0];
rz(1.4008201) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(0.31731269) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.7839157) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4006153) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(-3.0269347) q[1];
rz(-1.7371014) q[3];
sx q[3];
rz(-0.93207031) q[3];
sx q[3];
rz(-0.8730264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(-1.6395817) q[2];
rz(-2.8912985) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(-1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(1.4244351) q[0];
rz(2.6622488) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(0.11553484) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4452267) q[0];
sx q[0];
rz(-1.0123024) q[0];
sx q[0];
rz(-2.2541056) q[0];
rz(1.8477511) q[2];
sx q[2];
rz(-0.87522725) q[2];
sx q[2];
rz(0.19687523) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9563453) q[1];
sx q[1];
rz(-2.6223409) q[1];
sx q[1];
rz(0.96038702) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9980738) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(-1.4876175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(1.7158562) q[2];
rz(-1.7410295) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(-2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(0.26915959) q[0];
rz(-1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(-1.7620618) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7675161) q[0];
sx q[0];
rz(-0.19664581) q[0];
sx q[0];
rz(1.0023414) q[0];
x q[1];
rz(-1.3921521) q[2];
sx q[2];
rz(-2.7182455) q[2];
sx q[2];
rz(-1.64738) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79593149) q[1];
sx q[1];
rz(-2.965651) q[1];
sx q[1];
rz(-1.7453341) q[1];
rz(-2.8096301) q[3];
sx q[3];
rz(-1.8649857) q[3];
sx q[3];
rz(2.4320137) q[3];
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
rz(1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-2.1102171) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.17299962) q[2];
sx q[2];
rz(-2.0377918) q[2];
sx q[2];
rz(1.3755058) q[2];
rz(-1.7975939) q[3];
sx q[3];
rz(-0.26567017) q[3];
sx q[3];
rz(-2.9220823) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];