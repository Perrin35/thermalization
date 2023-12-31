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
rz(1.7226146) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6366633) q[0];
sx q[0];
rz(-1.814827) q[0];
sx q[0];
rz(-1.3958449) q[0];
rz(-pi) q[1];
rz(2.6909157) q[2];
sx q[2];
rz(-2.5833231) q[2];
sx q[2];
rz(-2.6435341) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4048684) q[1];
sx q[1];
rz(-1.0213335) q[1];
sx q[1];
rz(-0.0615555) q[1];
x q[2];
rz(-2.8098104) q[3];
sx q[3];
rz(-0.79091573) q[3];
sx q[3];
rz(-0.055981759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3966763) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(-2.7837226) q[2];
rz(0.19168028) q[3];
sx q[3];
rz(-2.7087757) q[3];
sx q[3];
rz(-0.55309692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-0.0033922694) q[0];
sx q[0];
rz(-1.1509742) q[0];
sx q[0];
rz(3.0733118) q[0];
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
rz(1.3119445) q[0];
sx q[0];
rz(-2.2011269) q[0];
sx q[0];
rz(-1.8207361) q[0];
x q[1];
rz(2.6158193) q[2];
sx q[2];
rz(-1.5500881) q[2];
sx q[2];
rz(-2.6244342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.56602851) q[1];
sx q[1];
rz(-1.0776507) q[1];
sx q[1];
rz(0.50599392) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7884368) q[3];
sx q[3];
rz(-1.0895035) q[3];
sx q[3];
rz(-0.35748112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(-2.2382656) q[2];
rz(2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(-0.27404684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.40619266) q[1];
sx q[1];
rz(2.1123871) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8619336) q[0];
sx q[0];
rz(-2.5760981) q[0];
sx q[0];
rz(2.4740567) q[0];
x q[1];
rz(2.7846787) q[2];
sx q[2];
rz(-1.9682353) q[2];
sx q[2];
rz(2.7473161) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.610565) q[1];
sx q[1];
rz(-0.37279168) q[1];
sx q[1];
rz(0.60786604) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3658931) q[3];
sx q[3];
rz(-2.587706) q[3];
sx q[3];
rz(-2.3212471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2091973) q[2];
sx q[2];
rz(-2.5019427) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(-0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9349174) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(2.7329965) q[0];
rz(-1.4273377) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(0.88159195) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13215412) q[0];
sx q[0];
rz(-2.3883925) q[0];
sx q[0];
rz(2.8892345) q[0];
x q[1];
rz(-0.84882952) q[2];
sx q[2];
rz(-2.1049044) q[2];
sx q[2];
rz(1.3918849) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2320764) q[1];
sx q[1];
rz(-1.0788003) q[1];
sx q[1];
rz(-0.75922774) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6374171) q[3];
sx q[3];
rz(-1.1542218) q[3];
sx q[3];
rz(1.6601603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0980229) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(-1.1126474) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(-1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(2.8240906) q[1];
sx q[1];
rz(-0.96173871) q[1];
sx q[1];
rz(-0.11725765) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1299767) q[0];
sx q[0];
rz(-0.59966171) q[0];
sx q[0];
rz(0.2325124) q[0];
rz(-0.60301493) q[2];
sx q[2];
rz(-1.3935967) q[2];
sx q[2];
rz(3.035383) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.458141) q[1];
sx q[1];
rz(-1.7675753) q[1];
sx q[1];
rz(1.4139921) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6860784) q[3];
sx q[3];
rz(-0.59024631) q[3];
sx q[3];
rz(0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(1.3797181) q[2];
rz(1.951925) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(-0.074507944) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18774408) q[0];
sx q[0];
rz(-1.501361) q[0];
sx q[0];
rz(2.4353819) q[0];
rz(-1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(-3.0336753) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(1.582167) q[2];
sx q[2];
rz(-1.5307384) q[2];
sx q[2];
rz(-2.1795321) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2006827) q[1];
sx q[1];
rz(-2.4001277) q[1];
sx q[1];
rz(0.42592589) q[1];
rz(-pi) q[2];
rz(1.0064441) q[3];
sx q[3];
rz(-0.77277771) q[3];
sx q[3];
rz(2.7882862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8191021) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(-2.0325913) q[2];
rz(1.2935982) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368822) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(0.014904508) q[0];
rz(2.7203454) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(-0.79963911) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9124303) q[0];
sx q[0];
rz(-1.3161356) q[0];
sx q[0];
rz(-0.69292648) q[0];
rz(-0.73115007) q[2];
sx q[2];
rz(-2.7292948) q[2];
sx q[2];
rz(-1.9129802) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7176197) q[1];
sx q[1];
rz(-1.7160001) q[1];
sx q[1];
rz(1.5145626) q[1];
rz(-pi) q[2];
rz(3.0040967) q[3];
sx q[3];
rz(-2.1318448) q[3];
sx q[3];
rz(1.6784799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.556276) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(-0.9220534) q[2];
rz(1.3098035) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(2.183765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9668982) q[0];
sx q[0];
rz(-2.0860724) q[0];
sx q[0];
rz(0.2640557) q[0];
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
rz(2.4357294) q[0];
sx q[0];
rz(-1.1840491) q[0];
sx q[0];
rz(-0.078588967) q[0];
rz(-pi) q[1];
rz(1.1515456) q[2];
sx q[2];
rz(-2.7630685) q[2];
sx q[2];
rz(-0.35767698) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4006153) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(3.0269347) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(2.7060966) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.6395817) q[2];
rz(0.25029415) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.2450927) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(-1.7171575) q[0];
rz(-0.4793438) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-3.0260578) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4452267) q[0];
sx q[0];
rz(-1.0123024) q[0];
sx q[0];
rz(-2.2541056) q[0];
rz(-pi) q[1];
rz(-2.8250474) q[2];
sx q[2];
rz(-2.4015421) q[2];
sx q[2];
rz(-0.22063247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.1852473) q[1];
sx q[1];
rz(-0.51925175) q[1];
sx q[1];
rz(2.1812056) q[1];
rz(-pi) q[2];
rz(0.1435189) q[3];
sx q[3];
rz(-0.70129881) q[3];
sx q[3];
rz(1.4876175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.4257365) q[2];
rz(1.7410295) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(-0.26915959) q[0];
rz(-1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(-1.7620618) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9513959) q[0];
sx q[0];
rz(-1.405389) q[0];
sx q[0];
rz(3.0347546) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3921521) q[2];
sx q[2];
rz(-0.42334712) q[2];
sx q[2];
rz(1.64738) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79593149) q[1];
sx q[1];
rz(-0.17594166) q[1];
sx q[1];
rz(-1.7453341) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8809324) q[3];
sx q[3];
rz(-1.8879859) q[3];
sx q[3];
rz(-0.96084259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(1.5157549) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(2.1102171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9243069) q[0];
sx q[0];
rz(-1.3047682) q[0];
sx q[0];
rz(-2.7346942) q[0];
rz(-2.6869607) q[1];
sx q[1];
rz(-2.0352719) q[1];
sx q[1];
rz(-0.24771053) q[1];
rz(1.2417912) q[2];
sx q[2];
rz(-0.49578373) q[2];
sx q[2];
rz(1.0052581) q[2];
rz(1.8299673) q[3];
sx q[3];
rz(-1.6298686) q[3];
sx q[3];
rz(2.009404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
