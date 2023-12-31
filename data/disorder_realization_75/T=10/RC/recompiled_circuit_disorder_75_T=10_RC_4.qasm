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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1184074) q[0];
sx q[0];
rz(-1.4010795) q[0];
sx q[0];
rz(0.24766185) q[0];
x q[1];
rz(0.45067694) q[2];
sx q[2];
rz(-2.5833231) q[2];
sx q[2];
rz(-0.4980586) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.19810361) q[1];
sx q[1];
rz(-1.6232821) q[1];
sx q[1];
rz(-1.0204888) q[1];
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
rz(-pi) q[1];
rz(1.3966763) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(-2.7837226) q[2];
rz(0.19168028) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(-2.5884957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0033922694) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-3.0733118) q[0];
rz(-1.0066907) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(0.65111792) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2383645) q[0];
sx q[0];
rz(-2.4698386) q[0];
sx q[0];
rz(-2.8147459) q[0];
x q[1];
rz(1.5468555) q[2];
sx q[2];
rz(-2.0964453) q[2];
sx q[2];
rz(-1.0416232) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.843833) q[1];
sx q[1];
rz(-2.4503772) q[1];
sx q[1];
rz(0.83696951) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6504374) q[3];
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
rz(-2.2382656) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.2925401) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-2.6480411) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-0.40619266) q[1];
sx q[1];
rz(-1.0292056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27965901) q[0];
sx q[0];
rz(-2.5760981) q[0];
sx q[0];
rz(-0.66753597) q[0];
rz(-pi) q[1];
rz(2.7846787) q[2];
sx q[2];
rz(-1.1733574) q[2];
sx q[2];
rz(0.39427653) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5310276) q[1];
sx q[1];
rz(-0.37279168) q[1];
sx q[1];
rz(0.60786604) q[1];
x q[2];
rz(3.0164099) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(-2.560905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(1.8445245) q[2];
rz(2.5629937) q[3];
sx q[3];
rz(-1.2207299) q[3];
sx q[3];
rz(-2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2066752) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(-2.7329965) q[0];
rz(-1.714255) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(0.88159195) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0094385) q[0];
sx q[0];
rz(-0.75320019) q[0];
sx q[0];
rz(2.8892345) q[0];
x q[1];
rz(0.66742113) q[2];
sx q[2];
rz(-2.1758658) q[2];
sx q[2];
rz(-0.24250008) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2243005) q[1];
sx q[1];
rz(-2.2227193) q[1];
sx q[1];
rz(-0.93445458) q[1];
rz(-pi) q[2];
rz(-1.6374171) q[3];
sx q[3];
rz(-1.1542218) q[3];
sx q[3];
rz(1.4814324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0435698) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(-2.0289452) q[2];
rz(-0.55551314) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(-1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1495789) q[0];
sx q[0];
rz(-1.9861789) q[0];
sx q[0];
rz(-1.3647112) q[0];
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
rz(1.7848894) q[2];
sx q[2];
rz(-2.1630686) q[2];
sx q[2];
rz(1.5562197) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6834516) q[1];
sx q[1];
rz(-1.7675753) q[1];
sx q[1];
rz(1.7276006) q[1];
x q[2];
rz(-0.076905964) q[3];
sx q[3];
rz(-2.1566026) q[3];
sx q[3];
rz(-2.2211071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7964898) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.3797181) q[2];
rz(1.951925) q[3];
sx q[3];
rz(-0.15933557) q[3];
sx q[3];
rz(3.0670847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88718016) q[0];
sx q[0];
rz(-0.88258703) q[0];
sx q[0];
rz(-1.1855017) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8651587) q[2];
sx q[2];
rz(-3.099953) q[2];
sx q[2];
rz(-0.68539884) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7522052) q[1];
sx q[1];
rz(-0.90837332) q[1];
sx q[1];
rz(1.2090769) q[1];
x q[2];
rz(-1.0064441) q[3];
sx q[3];
rz(-0.77277771) q[3];
sx q[3];
rz(-2.7882862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.32249054) q[2];
sx q[2];
rz(-0.97444797) q[2];
sx q[2];
rz(2.0325913) q[2];
rz(-1.2935982) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368822) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(-3.1266881) q[0];
rz(-0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(2.3419535) q[1];
rz(pi/2) q[2];
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
rz(2.4104426) q[2];
sx q[2];
rz(-2.7292948) q[2];
sx q[2];
rz(-1.9129802) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0529424) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(-2.774653) q[1];
rz(-1.0054672) q[3];
sx q[3];
rz(-1.6871095) q[3];
sx q[3];
rz(-3.1073991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5853167) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(-0.9220534) q[2];
rz(-1.8317892) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17469445) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(-0.2640557) q[0];
rz(1.4008201) q[1];
sx q[1];
rz(-1.6903279) q[1];
sx q[1];
rz(-2.82428) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89462751) q[0];
sx q[0];
rz(-1.4980226) q[0];
sx q[0];
rz(-1.9586246) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9900471) q[2];
sx q[2];
rz(-0.37852415) q[2];
sx q[2];
rz(-0.35767698) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.322284) q[1];
sx q[1];
rz(-1.6849736) q[1];
sx q[1];
rz(-1.4788005) q[1];
rz(-pi) q[2];
rz(1.7371014) q[3];
sx q[3];
rz(-0.93207031) q[3];
sx q[3];
rz(0.8730264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7060966) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(-0.25029415) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(-1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450927) q[0];
sx q[0];
rz(-2.8332906) q[0];
sx q[0];
rz(1.4244351) q[0];
rz(2.6622488) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-3.0260578) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45142052) q[0];
sx q[0];
rz(-0.85309404) q[0];
sx q[0];
rz(0.79057981) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8250474) q[2];
sx q[2];
rz(-2.4015421) q[2];
sx q[2];
rz(0.22063247) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.83976165) q[1];
sx q[1];
rz(-1.2823732) q[1];
sx q[1];
rz(2.0088058) q[1];
x q[2];
rz(-2.9980738) q[3];
sx q[3];
rz(-2.4402938) q[3];
sx q[3];
rz(1.6539751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.7158562) q[2];
rz(-1.7410295) q[3];
sx q[3];
rz(-1.0481342) q[3];
sx q[3];
rz(2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.35587674) q[0];
sx q[0];
rz(-0.73369217) q[0];
sx q[0];
rz(-2.8724331) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-0.91047374) q[1];
sx q[1];
rz(1.7620618) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3740765) q[0];
sx q[0];
rz(-2.9449468) q[0];
sx q[0];
rz(1.0023414) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3921521) q[2];
sx q[2];
rz(-2.7182455) q[2];
sx q[2];
rz(1.4942126) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5386242) q[1];
sx q[1];
rz(-1.6011964) q[1];
sx q[1];
rz(1.7441185) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33196253) q[3];
sx q[3];
rz(-1.276607) q[3];
sx q[3];
rz(0.709579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0925838) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(-1.6258378) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(1.0313755) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
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
rz(2.6869607) q[1];
sx q[1];
rz(-1.1063207) q[1];
sx q[1];
rz(2.8938821) q[1];
rz(-2.968593) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
rz(-1.3116253) q[3];
sx q[3];
rz(-1.6298686) q[3];
sx q[3];
rz(2.009404) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
