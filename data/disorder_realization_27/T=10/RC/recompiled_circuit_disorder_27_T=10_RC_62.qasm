OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.0032089631) q[0];
sx q[0];
rz(-0.15455833) q[0];
sx q[0];
rz(0.69252339) q[0];
rz(5.0737557) q[1];
sx q[1];
rz(4.3901246) q[1];
sx q[1];
rz(7.6683383) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6801075) q[0];
sx q[0];
rz(-2.2930817) q[0];
sx q[0];
rz(1.1397584) q[0];
rz(-pi) q[1];
rz(-1.295747) q[2];
sx q[2];
rz(-1.0359456) q[2];
sx q[2];
rz(1.5616852) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.79917158) q[1];
sx q[1];
rz(-2.85596) q[1];
sx q[1];
rz(2.530453) q[1];
x q[2];
rz(2.5296506) q[3];
sx q[3];
rz(-0.7512593) q[3];
sx q[3];
rz(-2.2236168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2549071) q[2];
sx q[2];
rz(-2.343785) q[2];
sx q[2];
rz(0.20516667) q[2];
rz(-2.3702879) q[3];
sx q[3];
rz(-2.3588534) q[3];
sx q[3];
rz(-2.0390959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40760621) q[0];
sx q[0];
rz(-0.74626958) q[0];
sx q[0];
rz(2.6876887) q[0];
rz(-2.1167963) q[1];
sx q[1];
rz(-2.7339934) q[1];
sx q[1];
rz(1.227238) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42392143) q[0];
sx q[0];
rz(-2.9917891) q[0];
sx q[0];
rz(-1.040209) q[0];
rz(-pi) q[1];
rz(1.111755) q[2];
sx q[2];
rz(-1.0978205) q[2];
sx q[2];
rz(-2.3943261) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2482359) q[1];
sx q[1];
rz(-1.3125988) q[1];
sx q[1];
rz(-0.30171079) q[1];
rz(-pi) q[2];
x q[2];
rz(0.20083986) q[3];
sx q[3];
rz(-1.4174995) q[3];
sx q[3];
rz(0.64985819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.1002905) q[2];
sx q[2];
rz(-1.9561448) q[2];
sx q[2];
rz(-2.5741637) q[2];
rz(-0.36519095) q[3];
sx q[3];
rz(-1.7285715) q[3];
sx q[3];
rz(2.173483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.48297468) q[0];
sx q[0];
rz(-0.56476074) q[0];
sx q[0];
rz(-2.2429402) q[0];
rz(0.99575106) q[1];
sx q[1];
rz(-1.5834705) q[1];
sx q[1];
rz(2.8083037) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62045287) q[0];
sx q[0];
rz(-0.95690173) q[0];
sx q[0];
rz(2.1472473) q[0];
x q[1];
rz(1.3938445) q[2];
sx q[2];
rz(-1.4743917) q[2];
sx q[2];
rz(-2.5071438) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8771613) q[1];
sx q[1];
rz(-0.96211551) q[1];
sx q[1];
rz(-0.18552893) q[1];
rz(-2.1786147) q[3];
sx q[3];
rz(-2.4993948) q[3];
sx q[3];
rz(-1.9574788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68625346) q[2];
sx q[2];
rz(-1.7909966) q[2];
sx q[2];
rz(-1.9906445) q[2];
rz(0.84093705) q[3];
sx q[3];
rz(-1.1154113) q[3];
sx q[3];
rz(1.9870728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-0.44160098) q[0];
sx q[0];
rz(-1.5804407) q[0];
sx q[0];
rz(-0.68471318) q[0];
rz(1.0355863) q[1];
sx q[1];
rz(-0.5077478) q[1];
sx q[1];
rz(-1.9365786) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8367774) q[0];
sx q[0];
rz(-2.4495821) q[0];
sx q[0];
rz(-1.6230323) q[0];
x q[1];
rz(2.2150546) q[2];
sx q[2];
rz(-1.4069923) q[2];
sx q[2];
rz(-1.0545727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9884831) q[1];
sx q[1];
rz(-1.8433237) q[1];
sx q[1];
rz(-2.2793798) q[1];
rz(-pi) q[2];
rz(-0.86592309) q[3];
sx q[3];
rz(-1.5161627) q[3];
sx q[3];
rz(-1.6161402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.24923199) q[2];
sx q[2];
rz(-1.426733) q[2];
sx q[2];
rz(-0.37115804) q[2];
rz(-1.7403729) q[3];
sx q[3];
rz(-2.4818082) q[3];
sx q[3];
rz(1.1192809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054759653) q[0];
sx q[0];
rz(-2.355447) q[0];
sx q[0];
rz(-0.13312419) q[0];
rz(0.99331028) q[1];
sx q[1];
rz(-1.3860093) q[1];
sx q[1];
rz(-0.55508074) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.426794) q[0];
sx q[0];
rz(-1.5580651) q[0];
sx q[0];
rz(-1.8861594) q[0];
rz(3.073408) q[2];
sx q[2];
rz(-1.1567332) q[2];
sx q[2];
rz(-0.33472862) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.811378) q[1];
sx q[1];
rz(-2.0159617) q[1];
sx q[1];
rz(0.021275714) q[1];
x q[2];
rz(-2.5693232) q[3];
sx q[3];
rz(-3.045445) q[3];
sx q[3];
rz(0.26086807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.30620265) q[2];
sx q[2];
rz(-2.1357048) q[2];
sx q[2];
rz(0.13892697) q[2];
rz(-0.94240087) q[3];
sx q[3];
rz(-1.5709632) q[3];
sx q[3];
rz(-0.55148235) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54365629) q[0];
sx q[0];
rz(-2.5718226) q[0];
sx q[0];
rz(0.57975769) q[0];
rz(0.12750164) q[1];
sx q[1];
rz(-1.189905) q[1];
sx q[1];
rz(1.5396083) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.126639) q[0];
sx q[0];
rz(-0.74288988) q[0];
sx q[0];
rz(-1.8735621) q[0];
x q[1];
rz(-1.6164262) q[2];
sx q[2];
rz(-1.4403733) q[2];
sx q[2];
rz(1.4554731) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.08338883) q[1];
sx q[1];
rz(-1.997588) q[1];
sx q[1];
rz(-0.31174387) q[1];
rz(-pi) q[2];
rz(1.2195915) q[3];
sx q[3];
rz(-0.71577365) q[3];
sx q[3];
rz(0.90482611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5876028) q[2];
sx q[2];
rz(-2.8911399) q[2];
sx q[2];
rz(0.26947752) q[2];
rz(-2.907471) q[3];
sx q[3];
rz(-0.52474371) q[3];
sx q[3];
rz(-0.060119303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44678974) q[0];
sx q[0];
rz(-2.5233874) q[0];
sx q[0];
rz(-0.68429464) q[0];
rz(-0.11958312) q[1];
sx q[1];
rz(-1.2922829) q[1];
sx q[1];
rz(2.6228242) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51661922) q[0];
sx q[0];
rz(-0.74665753) q[0];
sx q[0];
rz(1.3730511) q[0];
rz(-pi) q[1];
rz(1.7905551) q[2];
sx q[2];
rz(-0.81232729) q[2];
sx q[2];
rz(-2.1213558) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59473945) q[1];
sx q[1];
rz(-2.2097128) q[1];
sx q[1];
rz(-0.98313318) q[1];
rz(3.1154237) q[3];
sx q[3];
rz(-1.9713638) q[3];
sx q[3];
rz(-2.8699584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8873022) q[2];
sx q[2];
rz(-1.7742949) q[2];
sx q[2];
rz(-0.56345144) q[2];
rz(0.051579483) q[3];
sx q[3];
rz(-2.0023465) q[3];
sx q[3];
rz(-0.95190597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96034399) q[0];
sx q[0];
rz(-2.612817) q[0];
sx q[0];
rz(-1.7425591) q[0];
rz(2.3545806) q[1];
sx q[1];
rz(-1.1287289) q[1];
sx q[1];
rz(2.3972437) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58957419) q[0];
sx q[0];
rz(-0.15757832) q[0];
sx q[0];
rz(-2.7872859) q[0];
rz(1.0561981) q[2];
sx q[2];
rz(-1.4375763) q[2];
sx q[2];
rz(-1.0366057) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.0155322) q[1];
sx q[1];
rz(-2.0704401) q[1];
sx q[1];
rz(1.3152895) q[1];
x q[2];
rz(2.6870071) q[3];
sx q[3];
rz(-1.5496407) q[3];
sx q[3];
rz(-1.3190312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7156334) q[2];
sx q[2];
rz(-1.7461494) q[2];
sx q[2];
rz(-2.5320833) q[2];
rz(-2.4842747) q[3];
sx q[3];
rz(-2.4980563) q[3];
sx q[3];
rz(0.26143423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9534) q[0];
sx q[0];
rz(-0.094390079) q[0];
sx q[0];
rz(1.5040065) q[0];
rz(-1.9001182) q[1];
sx q[1];
rz(-1.991661) q[1];
sx q[1];
rz(-0.77493587) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0477714) q[0];
sx q[0];
rz(-1.1822961) q[0];
sx q[0];
rz(-0.85431487) q[0];
rz(-pi) q[1];
rz(-2.0140892) q[2];
sx q[2];
rz(-1.7121676) q[2];
sx q[2];
rz(3.1183426) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.7271125) q[1];
sx q[1];
rz(-0.4948805) q[1];
sx q[1];
rz(2.7369569) q[1];
rz(-pi) q[2];
rz(-1.3516515) q[3];
sx q[3];
rz(-1.431576) q[3];
sx q[3];
rz(-0.62300357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.187414) q[2];
sx q[2];
rz(-0.22970197) q[2];
sx q[2];
rz(2.9837218) q[2];
rz(-1.9291417) q[3];
sx q[3];
rz(-1.0586497) q[3];
sx q[3];
rz(-1.8635748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0697486) q[0];
sx q[0];
rz(-2.1691515) q[0];
sx q[0];
rz(0.21433314) q[0];
rz(0.65746039) q[1];
sx q[1];
rz(-2.9174556) q[1];
sx q[1];
rz(1.0459895) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8911274) q[0];
sx q[0];
rz(-1.9292826) q[0];
sx q[0];
rz(-1.4573218) q[0];
rz(-pi) q[1];
rz(-1.5418391) q[2];
sx q[2];
rz(-0.98000295) q[2];
sx q[2];
rz(2.7383885) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6194832) q[1];
sx q[1];
rz(-1.932718) q[1];
sx q[1];
rz(2.1543909) q[1];
rz(-pi) q[2];
rz(-3.0797327) q[3];
sx q[3];
rz(-0.40611551) q[3];
sx q[3];
rz(1.5554242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41632286) q[2];
sx q[2];
rz(-3.0815093) q[2];
sx q[2];
rz(2.0521169) q[2];
rz(1.5754835) q[3];
sx q[3];
rz(-1.8982866) q[3];
sx q[3];
rz(0.65264788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8626704) q[0];
sx q[0];
rz(-0.60386064) q[0];
sx q[0];
rz(0.84558564) q[0];
rz(1.5325585) q[1];
sx q[1];
rz(-1.6747723) q[1];
sx q[1];
rz(2.0369045) q[1];
rz(-0.40907787) q[2];
sx q[2];
rz(-1.8553875) q[2];
sx q[2];
rz(2.6864048) q[2];
rz(-0.55003008) q[3];
sx q[3];
rz(-1.6246272) q[3];
sx q[3];
rz(0.18679242) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];