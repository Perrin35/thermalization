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
rz(1.8887853) q[1];
sx q[1];
rz(10.843756) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1365294) q[0];
sx q[0];
rz(-2.8423474) q[0];
sx q[0];
rz(-2.5314999) q[0];
rz(-pi) q[1];
rz(-2.6909157) q[2];
sx q[2];
rz(-2.5833231) q[2];
sx q[2];
rz(2.6435341) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8542006) q[1];
sx q[1];
rz(-0.55254793) q[1];
sx q[1];
rz(-1.470675) q[1];
rz(1.2526413) q[3];
sx q[3];
rz(-2.3080024) q[3];
sx q[3];
rz(-2.6300501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7449164) q[2];
sx q[2];
rz(-1.6884721) q[2];
sx q[2];
rz(-2.7837226) q[2];
rz(2.9499124) q[3];
sx q[3];
rz(-0.43281698) q[3];
sx q[3];
rz(2.5884957) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-0.12193646) q[1];
sx q[1];
rz(0.65111792) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40819528) q[0];
sx q[0];
rz(-1.7719643) q[0];
sx q[0];
rz(-2.4961619) q[0];
rz(1.5947371) q[2];
sx q[2];
rz(-1.0451473) q[2];
sx q[2];
rz(2.0999694) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.56602851) q[1];
sx q[1];
rz(-1.0776507) q[1];
sx q[1];
rz(0.50599392) q[1];
x q[2];
rz(-0.39204709) q[3];
sx q[3];
rz(-0.52467504) q[3];
sx q[3];
rz(-0.80313659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6887001) q[2];
sx q[2];
rz(-0.15658997) q[2];
sx q[2];
rz(-2.2382656) q[2];
rz(2.3870758) q[3];
sx q[3];
rz(-1.6468331) q[3];
sx q[3];
rz(-2.8675458) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2925401) q[0];
sx q[0];
rz(-0.64340574) q[0];
sx q[0];
rz(0.49355155) q[0];
rz(-1.6429398) q[1];
sx q[1];
rz(-0.40619266) q[1];
sx q[1];
rz(-1.0292056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2632336) q[0];
sx q[0];
rz(-1.2326816) q[0];
sx q[0];
rz(-2.6792206) q[0];
rz(-pi) q[1];
rz(-0.35691397) q[2];
sx q[2];
rz(-1.9682353) q[2];
sx q[2];
rz(-0.39427653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.61470795) q[1];
sx q[1];
rz(-1.7803368) q[1];
sx q[1];
rz(-2.8309612) q[1];
rz(-pi) q[2];
rz(0.1251827) q[3];
sx q[3];
rz(-1.0297965) q[3];
sx q[3];
rz(2.560905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9323953) q[2];
sx q[2];
rz(-0.63964996) q[2];
sx q[2];
rz(-1.8445245) q[2];
rz(-2.5629937) q[3];
sx q[3];
rz(-1.9208627) q[3];
sx q[3];
rz(0.67563081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(-1.2066752) q[0];
sx q[0];
rz(-0.6466372) q[0];
sx q[0];
rz(-0.40859616) q[0];
rz(1.714255) q[1];
sx q[1];
rz(-1.6427549) q[1];
sx q[1];
rz(-0.88159195) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6245759) q[0];
sx q[0];
rz(-1.3991742) q[0];
sx q[0];
rz(-0.73715985) q[0];
x q[1];
rz(-2.4741715) q[2];
sx q[2];
rz(-2.1758658) q[2];
sx q[2];
rz(-0.24250008) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3410014) q[1];
sx q[1];
rz(-2.2644271) q[1];
sx q[1];
rz(0.6615521) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41739695) q[3];
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
rz(1.0980229) q[2];
sx q[2];
rz(-1.4901525) q[2];
sx q[2];
rz(-1.1126474) q[2];
rz(0.55551314) q[3];
sx q[3];
rz(-1.7881309) q[3];
sx q[3];
rz(1.6120733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1495789) q[0];
sx q[0];
rz(-1.9861789) q[0];
sx q[0];
rz(1.7768815) q[0];
rz(-0.31750202) q[1];
sx q[1];
rz(-2.1798539) q[1];
sx q[1];
rz(0.11725765) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1299767) q[0];
sx q[0];
rz(-0.59966171) q[0];
sx q[0];
rz(2.9090803) q[0];
rz(-pi) q[1];
rz(2.8357382) q[2];
sx q[2];
rz(-2.5161985) q[2];
sx q[2];
rz(-1.2139699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.458141) q[1];
sx q[1];
rz(-1.7675753) q[1];
sx q[1];
rz(-1.7276006) q[1];
rz(3.0646867) q[3];
sx q[3];
rz(-2.1566026) q[3];
sx q[3];
rz(-2.2211071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.34510288) q[2];
sx q[2];
rz(-1.779665) q[2];
sx q[2];
rz(1.7618746) q[2];
rz(-1.1896677) q[3];
sx q[3];
rz(-2.9822571) q[3];
sx q[3];
rz(0.074507944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9538486) q[0];
sx q[0];
rz(-1.6402316) q[0];
sx q[0];
rz(0.70621079) q[0];
rz(1.1114936) q[1];
sx q[1];
rz(-0.62875426) q[1];
sx q[1];
rz(-0.10791735) q[1];
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
rz(-pi) q[1];
x q[1];
rz(0.040060476) q[2];
sx q[2];
rz(-1.5594348) q[2];
sx q[2];
rz(-2.5324015) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7522052) q[1];
sx q[1];
rz(-2.2332193) q[1];
sx q[1];
rz(1.2090769) q[1];
rz(-pi) q[2];
x q[2];
rz(0.88166724) q[3];
sx q[3];
rz(-1.9534743) q[3];
sx q[3];
rz(1.4985639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8191021) q[2];
sx q[2];
rz(-2.1671447) q[2];
sx q[2];
rz(-2.0325913) q[2];
rz(1.8479944) q[3];
sx q[3];
rz(-1.3556017) q[3];
sx q[3];
rz(3.0343645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7727707) q[0];
sx q[0];
rz(-1.4819551) q[0];
sx q[0];
rz(0.014904508) q[0];
rz(-2.7203454) q[1];
sx q[1];
rz(-2.0887471) q[1];
sx q[1];
rz(0.79963911) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5478202) q[0];
sx q[0];
rz(-0.9043588) q[0];
sx q[0];
rz(-1.8970467) q[0];
x q[1];
rz(-0.31475474) q[2];
sx q[2];
rz(-1.2999279) q[2];
sx q[2];
rz(-2.1115007) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9866242) q[1];
sx q[1];
rz(-1.6264377) q[1];
sx q[1];
rz(2.9961622) q[1];
rz(-pi) q[2];
rz(-2.1361254) q[3];
sx q[3];
rz(-1.4544832) q[3];
sx q[3];
rz(-3.1073991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.556276) q[2];
sx q[2];
rz(-0.63627807) q[2];
sx q[2];
rz(0.9220534) q[2];
rz(1.3098035) q[3];
sx q[3];
rz(-1.2049048) q[3];
sx q[3];
rz(0.95782763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17469445) q[0];
sx q[0];
rz(-1.0555203) q[0];
sx q[0];
rz(2.877537) q[0];
rz(-1.4008201) q[1];
sx q[1];
rz(-1.4512647) q[1];
sx q[1];
rz(-2.82428) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6415629) q[0];
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
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.74097733) q[1];
sx q[1];
rz(-1.4794011) q[1];
sx q[1];
rz(-3.0269347) q[1];
rz(0.64538892) q[3];
sx q[3];
rz(-1.4374975) q[3];
sx q[3];
rz(-2.3440839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43549609) q[2];
sx q[2];
rz(-2.2075682) q[2];
sx q[2];
rz(1.502011) q[2];
rz(-2.8912985) q[3];
sx q[3];
rz(-1.4151662) q[3];
sx q[3];
rz(1.1423473) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2450927) q[0];
sx q[0];
rz(-0.30830202) q[0];
sx q[0];
rz(-1.7171575) q[0];
rz(2.6622488) q[1];
sx q[1];
rz(-1.665325) q[1];
sx q[1];
rz(-3.0260578) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6901721) q[0];
sx q[0];
rz(-2.2884986) q[0];
sx q[0];
rz(2.3510128) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.3165452) q[2];
sx q[2];
rz(-0.74005055) q[2];
sx q[2];
rz(-0.22063247) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.1852473) q[1];
sx q[1];
rz(-0.51925175) q[1];
sx q[1];
rz(-2.1812056) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4453836) q[3];
sx q[3];
rz(-1.4783825) q[3];
sx q[3];
rz(-2.948451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4775548) q[2];
sx q[2];
rz(-1.1096191) q[2];
sx q[2];
rz(-1.4257365) q[2];
rz(1.4005631) q[3];
sx q[3];
rz(-2.0934584) q[3];
sx q[3];
rz(-2.8588296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7857159) q[0];
sx q[0];
rz(-2.4079005) q[0];
sx q[0];
rz(0.26915959) q[0];
rz(1.0956988) q[1];
sx q[1];
rz(-2.2311189) q[1];
sx q[1];
rz(-1.7620618) q[1];
rz(-pi/2) q[2];
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
x q[1];
rz(1.7494406) q[2];
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
rz(-pi/2) q[0];
rz(2.3456612) q[1];
sx q[1];
rz(-2.965651) q[1];
sx q[1];
rz(1.3962586) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8809324) q[3];
sx q[3];
rz(-1.2536067) q[3];
sx q[3];
rz(-0.96084259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0925838) q[2];
sx q[2];
rz(-0.2139341) q[2];
sx q[2];
rz(-1.5157549) q[2];
rz(1.9745291) q[3];
sx q[3];
rz(-1.556282) q[3];
sx q[3];
rz(1.0313755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-2.6869607) q[1];
sx q[1];
rz(-2.0352719) q[1];
sx q[1];
rz(-0.24771053) q[1];
rz(0.17299962) q[2];
sx q[2];
rz(-1.1038008) q[2];
sx q[2];
rz(-1.7660869) q[2];
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
