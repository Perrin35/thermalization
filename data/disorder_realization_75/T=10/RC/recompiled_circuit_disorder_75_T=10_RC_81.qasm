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
rz(0.20180841) q[1];
sx q[1];
rz(1.8887853) q[1];
sx q[1];
rz(10.843756) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1184074) q[0];
sx q[0];
rz(-1.4010795) q[0];
sx q[0];
rz(-2.8939308) q[0];
x q[1];
rz(-1.3051891) q[2];
sx q[2];
rz(-1.0737537) q[2];
sx q[2];
rz(0.020393919) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.19810361) q[1];
sx q[1];
rz(-1.6232821) q[1];
sx q[1];
rz(-1.0204888) q[1];
rz(-pi) q[2];
rz(1.8889514) q[3];
sx q[3];
rz(-0.83359026) q[3];
sx q[3];
rz(0.51154256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(2.7837226) q[2];
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
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0033922694) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(0.068280846) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(2.4904747) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3119445) q[0];
sx q[0];
rz(-0.94046578) q[0];
sx q[0];
rz(-1.8207361) q[0];
rz(-pi) q[1];
rz(-3.100349) q[2];
sx q[2];
rz(-2.6154499) q[2];
sx q[2];
rz(-1.0893084) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5755641) q[1];
sx q[1];
rz(-1.0776507) q[1];
sx q[1];
rz(2.6355987) q[1];
x q[2];
rz(0.49115527) q[3];
sx q[3];
rz(-1.3782116) q[3];
sx q[3];
rz(-1.1112978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.45289257) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(0.90332705) q[2];
rz(0.7545169) q[3];
sx q[3];
rz(-1.6468331) q[3];
sx q[3];
rz(2.8675458) q[3];
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
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84905255) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(0.49355155) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-2.1123871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4714067) q[0];
sx q[0];
rz(-2.0051415) q[0];
sx q[0];
rz(-1.9451408) q[0];
x q[1];
rz(2.2649293) q[2];
sx q[2];
rz(-0.52774094) q[2];
sx q[2];
rz(-1.9805679) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5268847) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(-0.31063147) q[1];
x q[2];
rz(-1.0263222) q[3];
sx q[3];
rz(-1.6780276) q[3];
sx q[3];
rz(0.925392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2091973) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(1.2970682) q[2];
rz(-2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(-2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-1.2066752) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(-2.7329965) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-2.2600007) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13215412) q[0];
sx q[0];
rz(-2.3883925) q[0];
sx q[0];
rz(0.25235812) q[0];
x q[1];
rz(-2.300823) q[2];
sx q[2];
rz(-2.2731014) q[2];
sx q[2];
rz(0.70309397) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.2320764) q[1];
sx q[1];
rz(-1.0788003) q[1];
sx q[1];
rz(0.75922774) q[1];
rz(-pi) q[2];
rz(-0.14933417) q[3];
sx q[3];
rz(-0.42156005) q[3];
sx q[3];
rz(-1.644852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0980229) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(2.0289452) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(-1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9920138) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(-2.8240906) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(-3.024335) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2909572) q[0];
sx q[0];
rz(-2.1521749) q[0];
sx q[0];
rz(1.4145538) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8357382) q[2];
sx q[2];
rz(-2.5161985) q[2];
sx q[2];
rz(-1.9276227) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.77814233) q[1];
sx q[1];
rz(-2.8906129) q[1];
sx q[1];
rz(-2.4771033) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0646867) q[3];
sx q[3];
rz(-0.98499005) q[3];
sx q[3];
rz(-0.92048551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(-1.3797181) q[2];
rz(1.1896677) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(-0.074507944) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18774408) q[0];
sx q[0];
rz(-1.6402316) q[0];
sx q[0];
rz(2.4353819) q[0];
rz(1.1114936) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(0.10791735) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6861434) q[0];
sx q[0];
rz(-2.3685072) q[0];
sx q[0];
rz(2.7129052) q[0];
x q[1];
rz(1.5594257) q[2];
sx q[2];
rz(-1.6108542) q[2];
sx q[2];
rz(0.96206059) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2006827) q[1];
sx q[1];
rz(-2.4001277) q[1];
sx q[1];
rz(-2.7156668) q[1];
rz(-pi) q[2];
rz(0.48072731) q[3];
sx q[3];
rz(-0.93989621) q[3];
sx q[3];
rz(-2.7709099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8191021) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(1.1090013) q[2];
rz(1.2935982) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.368822) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(-0.014904508) q[0];
rz(-0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(2.3419535) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5478202) q[0];
sx q[0];
rz(-2.2372338) q[0];
sx q[0];
rz(1.8970467) q[0];
rz(-pi) q[1];
rz(-2.4104426) q[2];
sx q[2];
rz(-0.41229782) q[2];
sx q[2];
rz(-1.9129802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0529424) q[1];
sx q[1];
rz(-0.155641) q[1];
sx q[1];
rz(2.774653) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0054672) q[3];
sx q[3];
rz(-1.6871095) q[3];
sx q[3];
rz(-3.1073991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(-0.9220534) q[2];
rz(1.3098035) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17469445) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(0.2640557) q[0];
rz(-1.4008201) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(0.31731269) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50002977) q[0];
sx q[0];
rz(-2.7473358) q[0];
sx q[0];
rz(-1.7612329) q[0];
rz(-pi) q[1];
rz(-2.9810901) q[2];
sx q[2];
rz(-1.9151033) q[2];
sx q[2];
rz(-0.80489327) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74097733) q[1];
sx q[1];
rz(-1.6621915) q[1];
sx q[1];
rz(-3.0269347) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7371014) q[3];
sx q[3];
rz(-2.2095223) q[3];
sx q[3];
rz(2.2685662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.502011) q[2];
rz(0.25029415) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89649993) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(1.4244351) q[0];
rz(-0.4793438) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-3.0260578) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45142052) q[0];
sx q[0];
rz(-0.85309404) q[0];
sx q[0];
rz(-2.3510128) q[0];
x q[1];
rz(1.8477511) q[2];
sx q[2];
rz(-2.2663654) q[2];
sx q[2];
rz(-0.19687523) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.83976165) q[1];
sx q[1];
rz(-1.8592195) q[1];
sx q[1];
rz(1.1327869) q[1];
rz(-pi) q[2];
rz(-1.6910016) q[3];
sx q[3];
rz(-0.87814858) q[3];
sx q[3];
rz(-1.8409178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.7158562) q[2];
rz(1.7410295) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(-0.28276309) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(-2.8724331) q[0];
rz(-1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(-1.7620618) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36294286) q[0];
sx q[0];
rz(-1.6761707) q[0];
sx q[0];
rz(1.7371348) q[0];
rz(-pi) q[1];
rz(0.079897957) q[2];
sx q[2];
rz(-1.1546087) q[2];
sx q[2];
rz(-1.4518567) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1684432) q[1];
sx q[1];
rz(-1.397555) q[1];
sx q[1];
rz(-3.1107305) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8096301) q[3];
sx q[3];
rz(-1.276607) q[3];
sx q[3];
rz(-0.709579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.049008869) q[2];
sx q[2];
rz(-2.9276586) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(1.2417912) q[2];
sx q[2];
rz(-0.49578373) q[2];
sx q[2];
rz(1.0052581) q[2];
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