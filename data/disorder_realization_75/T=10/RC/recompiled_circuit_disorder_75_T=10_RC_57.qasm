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
rz(0.20180841) q[1];
sx q[1];
rz(-1.2528074) q[1];
sx q[1];
rz(1.7226146) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0231853) q[0];
sx q[0];
rz(-1.4010795) q[0];
sx q[0];
rz(-2.8939308) q[0];
x q[1];
rz(1.8364036) q[2];
sx q[2];
rz(-2.0678389) q[2];
sx q[2];
rz(3.1211987) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7367243) q[1];
sx q[1];
rz(-1.0213335) q[1];
sx q[1];
rz(0.0615555) q[1];
rz(1.8889514) q[3];
sx q[3];
rz(-2.3080024) q[3];
sx q[3];
rz(2.6300501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3966763) q[2];
sx q[2];
rz(-1.4531206) q[2];
sx q[2];
rz(-0.35787004) q[2];
rz(-2.9499124) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(-2.5884957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1382004) q[0];
sx q[0];
rz(-1.9906185) q[0];
sx q[0];
rz(-0.068280846) q[0];
rz(2.1349019) q[1];
sx q[1];
rz(-0.12193646) q[1];
sx q[1];
rz(-2.4904747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8296482) q[0];
sx q[0];
rz(-2.2011269) q[0];
sx q[0];
rz(-1.3208566) q[0];
rz(0.041243677) q[2];
sx q[2];
rz(-0.52614279) q[2];
sx q[2];
rz(1.0893084) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.56602851) q[1];
sx q[1];
rz(-1.0776507) q[1];
sx q[1];
rz(2.6355987) q[1];
rz(-pi) q[2];
rz(-0.39204709) q[3];
sx q[3];
rz(-2.6169176) q[3];
sx q[3];
rz(-2.3384561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2925401) q[0];
sx q[0];
rz(-2.4981869) q[0];
sx q[0];
rz(-0.49355155) q[0];
rz(1.4986528) q[1];
sx q[1];
rz(-2.7354) q[1];
sx q[1];
rz(1.0292056) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27965901) q[0];
sx q[0];
rz(-2.5760981) q[0];
sx q[0];
rz(2.4740567) q[0];
rz(-pi) q[1];
rz(-0.87666338) q[2];
sx q[2];
rz(-0.52774094) q[2];
sx q[2];
rz(-1.9805679) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.610565) q[1];
sx q[1];
rz(-0.37279168) q[1];
sx q[1];
rz(2.5337266) q[1];
rz(-pi) q[2];
rz(-2.1152705) q[3];
sx q[3];
rz(-1.4635651) q[3];
sx q[3];
rz(-2.2162007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2091973) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.2970682) q[2];
rz(0.57859892) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(-2.4659618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2066752) q[0];
sx q[0];
rz(-2.4949555) q[0];
sx q[0];
rz(-2.7329965) q[0];
rz(1.714255) q[1];
sx q[1];
rz(-1.4988377) q[1];
sx q[1];
rz(0.88159195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9339869) q[0];
sx q[0];
rz(-2.2946977) q[0];
sx q[0];
rz(1.3408322) q[0];
rz(-pi) q[1];
rz(0.66742113) q[2];
sx q[2];
rz(-2.1758658) q[2];
sx q[2];
rz(-0.24250008) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.91729212) q[1];
sx q[1];
rz(-2.2227193) q[1];
sx q[1];
rz(-2.2071381) q[1];
rz(-pi) q[2];
rz(2.7241957) q[3];
sx q[3];
rz(-1.5098803) q[3];
sx q[3];
rz(-3.0792189) q[3];
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
rz(1.1126474) q[2];
rz(-2.5860795) q[3];
sx q[3];
rz(-1.3534618) q[3];
sx q[3];
rz(-1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1495789) q[0];
sx q[0];
rz(-1.1554138) q[0];
sx q[0];
rz(1.3647112) q[0];
rz(-0.31750202) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(0.11725765) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1299767) q[0];
sx q[0];
rz(-0.59966171) q[0];
sx q[0];
rz(-2.9090803) q[0];
rz(-0.60301493) q[2];
sx q[2];
rz(-1.7479959) q[2];
sx q[2];
rz(-3.035383) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.77814233) q[1];
sx q[1];
rz(-0.25097978) q[1];
sx q[1];
rz(2.4771033) q[1];
x q[2];
rz(3.0646867) q[3];
sx q[3];
rz(-2.1566026) q[3];
sx q[3];
rz(-2.2211071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.34510288) q[2];
sx q[2];
rz(-1.3619276) q[2];
sx q[2];
rz(1.7618746) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(-3.0670847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538486) q[0];
sx q[0];
rz(-1.6402316) q[0];
sx q[0];
rz(-2.4353819) q[0];
rz(-1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(0.10791735) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2544125) q[0];
sx q[0];
rz(-0.88258703) q[0];
sx q[0];
rz(1.956091) q[0];
x q[1];
rz(-0.27643398) q[2];
sx q[2];
rz(-0.041639608) q[2];
sx q[2];
rz(-2.4561938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7522052) q[1];
sx q[1];
rz(-0.90837332) q[1];
sx q[1];
rz(1.9325158) q[1];
rz(-2.2599254) q[3];
sx q[3];
rz(-1.1881184) q[3];
sx q[3];
rz(1.6430287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.32249054) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(2.0325913) q[2];
rz(1.2935982) q[3];
sx q[3];
rz(-1.785991) q[3];
sx q[3];
rz(-0.10722815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.368822) q[0];
sx q[0];
rz(-1.6596376) q[0];
sx q[0];
rz(-0.014904508) q[0];
rz(0.42124721) q[1];
sx q[1];
rz(-1.0528456) q[1];
sx q[1];
rz(2.3419535) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0470456) q[0];
sx q[0];
rz(-2.4106815) q[0];
sx q[0];
rz(2.7546309) q[0];
rz(-pi) q[1];
rz(2.8268379) q[2];
sx q[2];
rz(-1.8416648) q[2];
sx q[2];
rz(2.1115007) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0886503) q[1];
sx q[1];
rz(-2.9859516) q[1];
sx q[1];
rz(-2.774653) q[1];
rz(-3.0040967) q[3];
sx q[3];
rz(-1.0097479) q[3];
sx q[3];
rz(-1.4631127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5853167) q[2];
sx q[2];
rz(-2.5053146) q[2];
sx q[2];
rz(-2.2195393) q[2];
rz(-1.8317892) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(-2.183765) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9668982) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(2.877537) q[0];
rz(-1.7407725) q[1];
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
rz(-0.89462751) q[0];
sx q[0];
rz(-1.6435701) q[0];
sx q[0];
rz(-1.1829681) q[0];
rz(-2.9810901) q[2];
sx q[2];
rz(-1.9151033) q[2];
sx q[2];
rz(-0.80489327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4006153) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(-0.11465794) q[1];
rz(2.4962037) q[3];
sx q[3];
rz(-1.4374975) q[3];
sx q[3];
rz(-0.79750878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.43549609) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.502011) q[2];
rz(-0.25029415) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(-1.1423473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.4244351) q[0];
rz(-0.4793438) q[1];
sx q[1];
rz(-1.4762676) q[1];
sx q[1];
rz(-0.11553484) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4452267) q[0];
sx q[0];
rz(-2.1292902) q[0];
sx q[0];
rz(2.2541056) q[0];
rz(1.2938415) q[2];
sx q[2];
rz(-2.2663654) q[2];
sx q[2];
rz(-2.9447174) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.301831) q[1];
sx q[1];
rz(-1.2823732) q[1];
sx q[1];
rz(-2.0088058) q[1];
x q[2];
rz(-1.6910016) q[3];
sx q[3];
rz(-0.87814858) q[3];
sx q[3];
rz(-1.8409178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.66403786) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(1.4257365) q[2];
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
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35587674) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(-2.8724331) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(1.3795308) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7786498) q[0];
sx q[0];
rz(-1.6761707) q[0];
sx q[0];
rz(-1.4044579) q[0];
rz(-pi) q[1];
rz(0.079897957) q[2];
sx q[2];
rz(-1.1546087) q[2];
sx q[2];
rz(-1.4518567) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.79593149) q[1];
sx q[1];
rz(-2.965651) q[1];
sx q[1];
rz(-1.3962586) q[1];
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
rz(-0.2139341) q[2];
sx q[2];
rz(-1.5157549) q[2];
rz(1.9745291) q[3];
sx q[3];
rz(-1.5853106) q[3];
sx q[3];
rz(2.1102171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-1.2417912) q[2];
sx q[2];
rz(-2.6458089) q[2];
sx q[2];
rz(-2.1363346) q[2];
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
