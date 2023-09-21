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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1184074) q[0];
sx q[0];
rz(-1.4010795) q[0];
sx q[0];
rz(-2.8939308) q[0];
rz(0.51214829) q[2];
sx q[2];
rz(-1.3379659) q[2];
sx q[2];
rz(1.6793959) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7367243) q[1];
sx q[1];
rz(-2.1202592) q[1];
sx q[1];
rz(-3.0800372) q[1];
x q[2];
rz(0.76288335) q[3];
sx q[3];
rz(-1.3370822) q[3];
sx q[3];
rz(1.2771311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7449164) q[2];
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
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0033922694) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-3.0733118) q[0];
rz(2.1349019) q[1];
sx q[1];
rz(-3.0196562) q[1];
sx q[1];
rz(2.4904747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40819528) q[0];
sx q[0];
rz(-1.3696284) q[0];
sx q[0];
rz(0.64543076) q[0];
rz(3.100349) q[2];
sx q[2];
rz(-0.52614279) q[2];
sx q[2];
rz(-1.0893084) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56602851) q[1];
sx q[1];
rz(-1.0776507) q[1];
sx q[1];
rz(-0.50599392) q[1];
rz(-0.49115527) q[3];
sx q[3];
rz(-1.3782116) q[3];
sx q[3];
rz(1.1112978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.45289257) q[2];
sx q[2];
rz(-2.9850027) q[2];
sx q[2];
rz(-2.2382656) q[2];
rz(2.3870758) q[3];
sx q[3];
rz(-1.4947596) q[3];
sx q[3];
rz(2.8675458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2925401) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(2.6480411) q[0];
rz(-1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(-1.0292056) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8619336) q[0];
sx q[0];
rz(-0.56549457) q[0];
sx q[0];
rz(-0.66753597) q[0];
rz(-pi) q[1];
x q[1];
rz(0.35691397) q[2];
sx q[2];
rz(-1.9682353) q[2];
sx q[2];
rz(-2.7473161) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61470795) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(0.31063147) q[1];
rz(-pi) q[2];
rz(1.3658931) q[3];
sx q[3];
rz(-2.587706) q[3];
sx q[3];
rz(-2.3212471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2066752) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-0.40859616) q[0];
rz(-1.714255) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(-0.88159195) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20760575) q[0];
sx q[0];
rz(-2.2946977) q[0];
sx q[0];
rz(1.3408322) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.84882952) q[2];
sx q[2];
rz(-2.1049044) q[2];
sx q[2];
rz(-1.7497077) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3410014) q[1];
sx q[1];
rz(-0.87716555) q[1];
sx q[1];
rz(-2.4800406) q[1];
rz(-pi) q[2];
rz(-0.14933417) q[3];
sx q[3];
rz(-0.42156005) q[3];
sx q[3];
rz(1.4967407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0980229) q[2];
sx q[2];
rz(-1.6514401) q[2];
sx q[2];
rz(1.1126474) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(-1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9920138) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(-1.7768815) q[0];
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
rz(0.011615959) q[0];
sx q[0];
rz(-2.5419309) q[0];
sx q[0];
rz(0.2325124) q[0];
rz(-pi) q[1];
rz(-2.5385777) q[2];
sx q[2];
rz(-1.7479959) q[2];
sx q[2];
rz(3.035383) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3634503) q[1];
sx q[1];
rz(-0.25097978) q[1];
sx q[1];
rz(-0.66448934) q[1];
x q[2];
rz(1.6860784) q[3];
sx q[3];
rz(-0.59024631) q[3];
sx q[3];
rz(-0.78199996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34510288) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(-1.1114936) q[1];
sx q[1];
rz(-2.5128384) q[1];
sx q[1];
rz(3.0336753) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2544125) q[0];
sx q[0];
rz(-0.88258703) q[0];
sx q[0];
rz(-1.1855017) q[0];
rz(-pi) q[1];
rz(-0.27643398) q[2];
sx q[2];
rz(-0.041639608) q[2];
sx q[2];
rz(-2.4561938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.95277112) q[1];
sx q[1];
rz(-1.2880039) q[1];
sx q[1];
rz(2.4464843) q[1];
x q[2];
rz(-2.1351486) q[3];
sx q[3];
rz(-0.77277771) q[3];
sx q[3];
rz(2.7882862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8191021) q[2];
sx q[2];
rz(-2.1671447) q[2];
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
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7727707) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(-3.1266881) q[0];
rz(0.42124721) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(0.79963911) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945471) q[0];
sx q[0];
rz(-2.4106815) q[0];
sx q[0];
rz(2.7546309) q[0];
rz(-pi) q[1];
rz(-1.854935) q[2];
sx q[2];
rz(-1.8737027) q[2];
sx q[2];
rz(2.6877833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0886503) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(-2.774653) q[1];
rz(3.0040967) q[3];
sx q[3];
rz(-2.1318448) q[3];
sx q[3];
rz(-1.4631127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(2.2195393) q[2];
rz(-1.8317892) q[3];
sx q[3];
rz(-1.9366879) q[3];
sx q[3];
rz(2.183765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-0.31731269) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89462751) q[0];
sx q[0];
rz(-1.4980226) q[0];
sx q[0];
rz(-1.9586246) q[0];
x q[1];
rz(-0.1605026) q[2];
sx q[2];
rz(-1.9151033) q[2];
sx q[2];
rz(-2.3366994) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4999769) q[1];
sx q[1];
rz(-0.14650211) q[1];
sx q[1];
rz(0.67540692) q[1];
rz(-1.7371014) q[3];
sx q[3];
rz(-0.93207031) q[3];
sx q[3];
rz(2.2685662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7060966) q[2];
sx q[2];
rz(-0.93402445) q[2];
sx q[2];
rz(-1.502011) q[2];
rz(0.25029415) q[3];
sx q[3];
rz(-1.7264265) q[3];
sx q[3];
rz(1.9992453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.7171575) q[0];
rz(0.4793438) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-0.11553484) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.696366) q[0];
sx q[0];
rz(-2.1292902) q[0];
sx q[0];
rz(0.88748705) q[0];
x q[1];
rz(1.8477511) q[2];
sx q[2];
rz(-2.2663654) q[2];
sx q[2];
rz(2.9447174) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.83976165) q[1];
sx q[1];
rz(-1.2823732) q[1];
sx q[1];
rz(1.1327869) q[1];
rz(-pi) q[2];
rz(1.450591) q[3];
sx q[3];
rz(-2.2634441) q[3];
sx q[3];
rz(1.8409178) q[3];
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
rz(-2.0934584) q[3];
sx q[3];
rz(0.28276309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
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
rz(-2.2311189) q[1];
sx q[1];
rz(-1.3795308) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3740765) q[0];
sx q[0];
rz(-0.19664581) q[0];
sx q[0];
rz(1.0023414) q[0];
x q[1];
rz(1.3921521) q[2];
sx q[2];
rz(-2.7182455) q[2];
sx q[2];
rz(-1.4942126) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1684432) q[1];
sx q[1];
rz(-1.7440376) q[1];
sx q[1];
rz(0.030862191) q[1];
x q[2];
rz(-0.33196253) q[3];
sx q[3];
rz(-1.8649857) q[3];
sx q[3];
rz(-2.4320137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.0925838) q[2];
sx q[2];
rz(-2.9276586) q[2];
sx q[2];
rz(-1.5157549) q[2];
rz(-1.1670636) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(-1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-2.0438647) q[2];
sx q[2];
rz(-1.4164783) q[2];
sx q[2];
rz(2.867792) q[2];
rz(1.7975939) q[3];
sx q[3];
rz(-2.8759225) q[3];
sx q[3];
rz(0.21951036) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];