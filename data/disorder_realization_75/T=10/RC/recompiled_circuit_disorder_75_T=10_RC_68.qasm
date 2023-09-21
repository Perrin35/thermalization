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
rz(-2.7379524) q[0];
sx q[0];
rz(0.37024745) q[0];
rz(-2.9397842) q[1];
sx q[1];
rz(-1.8887853) q[1];
sx q[1];
rz(-1.7226146) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50492935) q[0];
sx q[0];
rz(-1.3267656) q[0];
sx q[0];
rz(-1.3958449) q[0];
rz(-pi) q[1];
rz(-1.3051891) q[2];
sx q[2];
rz(-2.0678389) q[2];
sx q[2];
rz(3.1211987) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.19810361) q[1];
sx q[1];
rz(-1.5183105) q[1];
sx q[1];
rz(-2.1211038) q[1];
rz(-pi) q[2];
rz(2.3787093) q[3];
sx q[3];
rz(-1.8045104) q[3];
sx q[3];
rz(1.2771311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7449164) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(0.35787004) q[2];
rz(-2.9499124) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-3.1382004) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-0.068280846) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(0.65111792) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7333974) q[0];
sx q[0];
rz(-1.7719643) q[0];
sx q[0];
rz(-0.64543076) q[0];
x q[1];
rz(0.041243677) q[2];
sx q[2];
rz(-0.52614279) q[2];
sx q[2];
rz(-2.0522842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5755641) q[1];
sx q[1];
rz(-2.063942) q[1];
sx q[1];
rz(2.6355987) q[1];
rz(-pi) q[2];
rz(2.6504374) q[3];
sx q[3];
rz(-1.7633811) q[3];
sx q[3];
rz(2.0302949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45289257) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(-2.2382656) q[2];
rz(-2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84905255) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-0.49355155) q[0];
rz(-1.4986528) q[1];
sx q[1];
rz(-0.40619266) q[1];
sx q[1];
rz(-2.1123871) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4714067) q[0];
sx q[0];
rz(-1.1364511) q[0];
sx q[0];
rz(-1.9451408) q[0];
rz(-2.2649293) q[2];
sx q[2];
rz(-0.52774094) q[2];
sx q[2];
rz(-1.1610247) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61470795) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(-0.31063147) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3658931) q[3];
sx q[3];
rz(-0.5538867) q[3];
sx q[3];
rz(2.3212471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(-2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349174) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-2.7329965) q[0];
rz(1.4273377) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(2.2600007) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20760575) q[0];
sx q[0];
rz(-0.84689492) q[0];
sx q[0];
rz(1.3408322) q[0];
rz(2.4741715) q[2];
sx q[2];
rz(-0.9657269) q[2];
sx q[2];
rz(2.8990926) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.2320764) q[1];
sx q[1];
rz(-2.0627923) q[1];
sx q[1];
rz(-0.75922774) q[1];
x q[2];
rz(2.7241957) q[3];
sx q[3];
rz(-1.6317123) q[3];
sx q[3];
rz(3.0792189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0435698) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(-1.1126474) q[2];
rz(-0.55551314) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(-1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1495789) q[0];
sx q[0];
rz(-1.9861789) q[0];
sx q[0];
rz(1.3647112) q[0];
rz(0.31750202) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(3.024335) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2909572) q[0];
sx q[0];
rz(-0.98941776) q[0];
sx q[0];
rz(1.7270389) q[0];
rz(1.7848894) q[2];
sx q[2];
rz(-0.97852409) q[2];
sx q[2];
rz(-1.5562197) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9980364) q[1];
sx q[1];
rz(-1.7245502) q[1];
sx q[1];
rz(0.19915944) q[1];
rz(-1.4555143) q[3];
sx q[3];
rz(-0.59024631) q[3];
sx q[3];
rz(-0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[1];
rz(pi/2) q[1];
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
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88718016) q[0];
sx q[0];
rz(-0.88258703) q[0];
sx q[0];
rz(-1.956091) q[0];
rz(-0.27643398) q[2];
sx q[2];
rz(-3.099953) q[2];
sx q[2];
rz(-0.68539884) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.38938746) q[1];
sx q[1];
rz(-2.2332193) q[1];
sx q[1];
rz(1.9325158) q[1];
x q[2];
rz(0.48072731) q[3];
sx q[3];
rz(-2.2016964) q[3];
sx q[3];
rz(2.7709099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.32249054) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(-1.1090013) q[2];
rz(1.8479944) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
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
rz(-0.7727707) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(-0.014904508) q[0];
rz(0.42124721) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(2.3419535) q[1];
rz(-pi) q[2];
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
rz(2.8268379) q[2];
sx q[2];
rz(-1.2999279) q[2];
sx q[2];
rz(1.0300919) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0529424) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(-0.36693962) q[1];
rz(-2.1361254) q[3];
sx q[3];
rz(-1.6871095) q[3];
sx q[3];
rz(3.1073991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17469445) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(2.877537) q[0];
rz(1.7407725) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(2.82428) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70586328) q[0];
sx q[0];
rz(-1.1840491) q[0];
sx q[0];
rz(3.0630037) q[0];
rz(1.2223586) q[2];
sx q[2];
rz(-1.4197883) q[2];
sx q[2];
rz(2.3210971) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.322284) q[1];
sx q[1];
rz(-1.6849736) q[1];
sx q[1];
rz(-1.6627922) q[1];
x q[2];
rz(1.4044912) q[3];
sx q[3];
rz(-2.2095223) q[3];
sx q[3];
rz(-2.2685662) q[3];
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
rz(-0.25029415) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(1.1423473) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89649993) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.7171575) q[0];
rz(-0.4793438) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(3.0260578) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6901721) q[0];
sx q[0];
rz(-0.85309404) q[0];
sx q[0];
rz(2.3510128) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3165452) q[2];
sx q[2];
rz(-0.74005055) q[2];
sx q[2];
rz(-0.22063247) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9563453) q[1];
sx q[1];
rz(-0.51925175) q[1];
sx q[1];
rz(2.1812056) q[1];
rz(-pi) q[2];
rz(2.9980738) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(1.6539751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(1.4257365) q[2];
rz(-1.7410295) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
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
rz(1.3795308) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7786498) q[0];
sx q[0];
rz(-1.4654219) q[0];
sx q[0];
rz(-1.7371348) q[0];
rz(0.079897957) q[2];
sx q[2];
rz(-1.1546087) q[2];
sx q[2];
rz(-1.4518567) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.60296842) q[1];
sx q[1];
rz(-1.6011964) q[1];
sx q[1];
rz(1.7441185) q[1];
rz(-pi) q[2];
rz(0.33196253) q[3];
sx q[3];
rz(-1.276607) q[3];
sx q[3];
rz(-2.4320137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.049008869) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(1.6258378) q[2];
rz(1.9745291) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(-2.1102171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2172858) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(-0.45463195) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(2.968593) q[2];
sx q[2];
rz(-2.0377918) q[2];
sx q[2];
rz(1.3755058) q[2];
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
