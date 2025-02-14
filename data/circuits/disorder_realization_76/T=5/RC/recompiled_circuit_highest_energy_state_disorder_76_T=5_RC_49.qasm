OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4401378) q[0];
sx q[0];
rz(-0.48818809) q[0];
sx q[0];
rz(-2.7421283) q[0];
rz(5.2108555) q[1];
sx q[1];
rz(4.0401754) q[1];
sx q[1];
rz(5.9730692) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7061569) q[0];
sx q[0];
rz(-2.6118267) q[0];
sx q[0];
rz(2.4166773) q[0];
rz(-pi) q[1];
rz(-2.1793876) q[2];
sx q[2];
rz(-2.1748389) q[2];
sx q[2];
rz(-2.5317328) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.063608544) q[1];
sx q[1];
rz(-2.0493174) q[1];
sx q[1];
rz(1.6850745) q[1];
rz(-pi) q[2];
rz(0.10914321) q[3];
sx q[3];
rz(-1.1880842) q[3];
sx q[3];
rz(-0.34527895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.255379) q[2];
sx q[2];
rz(-1.1905406) q[2];
sx q[2];
rz(-1.0211241) q[2];
rz(1.6518263) q[3];
sx q[3];
rz(-1.3083873) q[3];
sx q[3];
rz(-2.3511353) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0806231) q[0];
sx q[0];
rz(-2.1014093) q[0];
sx q[0];
rz(2.8676497) q[0];
rz(-0.25320369) q[1];
sx q[1];
rz(-2.4174523) q[1];
sx q[1];
rz(0.23987548) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1254933) q[0];
sx q[0];
rz(-2.3224127) q[0];
sx q[0];
rz(-2.3045834) q[0];
rz(-0.48440964) q[2];
sx q[2];
rz(-0.61867117) q[2];
sx q[2];
rz(1.8687488) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1170711) q[1];
sx q[1];
rz(-1.2327984) q[1];
sx q[1];
rz(1.3282969) q[1];
rz(-pi) q[2];
rz(-1.3680458) q[3];
sx q[3];
rz(-0.79907954) q[3];
sx q[3];
rz(-2.5917724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6720752) q[2];
sx q[2];
rz(-1.8263475) q[2];
sx q[2];
rz(1.7496656) q[2];
rz(2.4972656) q[3];
sx q[3];
rz(-1.143012) q[3];
sx q[3];
rz(-2.2279975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7740087) q[0];
sx q[0];
rz(-0.82355654) q[0];
sx q[0];
rz(-2.9665663) q[0];
rz(-1.7063829) q[1];
sx q[1];
rz(-1.274546) q[1];
sx q[1];
rz(-2.3451436) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1131409) q[0];
sx q[0];
rz(-2.7243834) q[0];
sx q[0];
rz(0.88930486) q[0];
rz(0.59993125) q[2];
sx q[2];
rz(-1.33432) q[2];
sx q[2];
rz(-1.8452871) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0092344) q[1];
sx q[1];
rz(-1.9446484) q[1];
sx q[1];
rz(-2.641885) q[1];
rz(1.4723331) q[3];
sx q[3];
rz(-1.0659127) q[3];
sx q[3];
rz(-1.7652924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1548057) q[2];
sx q[2];
rz(-0.044175819) q[2];
sx q[2];
rz(3.0955443) q[2];
rz(1.1151399) q[3];
sx q[3];
rz(-2.1348848) q[3];
sx q[3];
rz(-1.0435638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1107156) q[0];
sx q[0];
rz(-2.5512888) q[0];
sx q[0];
rz(-0.16511551) q[0];
rz(1.1394399) q[1];
sx q[1];
rz(-2.20729) q[1];
sx q[1];
rz(-1.6868235) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.063140537) q[0];
sx q[0];
rz(-1.8179896) q[0];
sx q[0];
rz(-1.4726223) q[0];
rz(-pi) q[1];
rz(-2.7023849) q[2];
sx q[2];
rz(-1.3813218) q[2];
sx q[2];
rz(-0.89292347) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36697179) q[1];
sx q[1];
rz(-0.72155385) q[1];
sx q[1];
rz(-1.6375035) q[1];
x q[2];
rz(-0.99253486) q[3];
sx q[3];
rz(-0.49697486) q[3];
sx q[3];
rz(-0.24028905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4399061) q[2];
sx q[2];
rz(-0.88752282) q[2];
sx q[2];
rz(-2.9627724) q[2];
rz(1.5170826) q[3];
sx q[3];
rz(-2.8330018) q[3];
sx q[3];
rz(1.6592244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6521859) q[0];
sx q[0];
rz(-1.711373) q[0];
sx q[0];
rz(1.7328316) q[0];
rz(0.57213712) q[1];
sx q[1];
rz(-1.1502384) q[1];
sx q[1];
rz(2.9170654) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.029789797) q[0];
sx q[0];
rz(-1.1953314) q[0];
sx q[0];
rz(2.445604) q[0];
rz(1.9064181) q[2];
sx q[2];
rz(-0.76424176) q[2];
sx q[2];
rz(-2.3847876) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.71345879) q[1];
sx q[1];
rz(-1.9947707) q[1];
sx q[1];
rz(-2.3522373) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5314674) q[3];
sx q[3];
rz(-0.49558276) q[3];
sx q[3];
rz(-0.15935824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.81395927) q[2];
sx q[2];
rz(-1.2322793) q[2];
sx q[2];
rz(0.28780469) q[2];
rz(-0.52590251) q[3];
sx q[3];
rz(-1.9764427) q[3];
sx q[3];
rz(2.5743918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7228058) q[0];
sx q[0];
rz(-2.4961508) q[0];
sx q[0];
rz(0.01471113) q[0];
rz(-2.9171433) q[1];
sx q[1];
rz(-1.6818455) q[1];
sx q[1];
rz(0.86123484) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5435641) q[0];
sx q[0];
rz(-0.042994067) q[0];
sx q[0];
rz(-0.046648101) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9056896) q[2];
sx q[2];
rz(-1.5318724) q[2];
sx q[2];
rz(0.30920497) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3726681) q[1];
sx q[1];
rz(-0.75809352) q[1];
sx q[1];
rz(-0.43129403) q[1];
rz(0.16448824) q[3];
sx q[3];
rz(-0.93014088) q[3];
sx q[3];
rz(-1.3817087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.33846778) q[2];
sx q[2];
rz(-2.4305692) q[2];
sx q[2];
rz(3.0306446) q[2];
rz(1.126193) q[3];
sx q[3];
rz(-1.6323099) q[3];
sx q[3];
rz(0.66669983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.8768537) q[0];
sx q[0];
rz(-2.5781317) q[0];
sx q[0];
rz(1.2891084) q[0];
rz(1.8076757) q[1];
sx q[1];
rz(-1.8218808) q[1];
sx q[1];
rz(1.2510373) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97020075) q[0];
sx q[0];
rz(-1.3193865) q[0];
sx q[0];
rz(-0.89604495) q[0];
rz(-pi) q[1];
rz(-1.8952605) q[2];
sx q[2];
rz(-1.0203385) q[2];
sx q[2];
rz(0.40263995) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.93537718) q[1];
sx q[1];
rz(-0.70997974) q[1];
sx q[1];
rz(-1.1919271) q[1];
rz(-pi) q[2];
rz(0.23836191) q[3];
sx q[3];
rz(-1.2796448) q[3];
sx q[3];
rz(-1.0358159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.053607792) q[2];
sx q[2];
rz(-1.3259652) q[2];
sx q[2];
rz(-1.717022) q[2];
rz(2.9234431) q[3];
sx q[3];
rz(-1.4616707) q[3];
sx q[3];
rz(0.48733369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3394534) q[0];
sx q[0];
rz(-0.55291432) q[0];
sx q[0];
rz(2.6735027) q[0];
rz(-0.83174902) q[1];
sx q[1];
rz(-2.2305198) q[1];
sx q[1];
rz(2.1591878) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8196682) q[0];
sx q[0];
rz(-1.5005932) q[0];
sx q[0];
rz(1.5062259) q[0];
x q[1];
rz(1.8562819) q[2];
sx q[2];
rz(-2.4997556) q[2];
sx q[2];
rz(1.8571719) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1156128) q[1];
sx q[1];
rz(-1.8048688) q[1];
sx q[1];
rz(1.267141) q[1];
x q[2];
rz(0.053388552) q[3];
sx q[3];
rz(-1.777919) q[3];
sx q[3];
rz(0.46536756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.19782797) q[2];
sx q[2];
rz(-0.48913893) q[2];
sx q[2];
rz(0.93734199) q[2];
rz(-1.7249974) q[3];
sx q[3];
rz(-1.2920486) q[3];
sx q[3];
rz(0.14911252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26405239) q[0];
sx q[0];
rz(-2.6314681) q[0];
sx q[0];
rz(-0.75793761) q[0];
rz(2.0856048) q[1];
sx q[1];
rz(-0.8989369) q[1];
sx q[1];
rz(-0.24076375) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7997076) q[0];
sx q[0];
rz(-2.1043964) q[0];
sx q[0];
rz(-1.207229) q[0];
x q[1];
rz(0.82334907) q[2];
sx q[2];
rz(-1.2687131) q[2];
sx q[2];
rz(0.27034098) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.71320885) q[1];
sx q[1];
rz(-1.2782103) q[1];
sx q[1];
rz(-0.87134113) q[1];
rz(-pi) q[2];
rz(1.2077755) q[3];
sx q[3];
rz(-1.2524093) q[3];
sx q[3];
rz(-2.8045079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9869438) q[2];
sx q[2];
rz(-1.9754396) q[2];
sx q[2];
rz(1.1692952) q[2];
rz(-2.5985006) q[3];
sx q[3];
rz(-1.8648059) q[3];
sx q[3];
rz(2.5663466) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52694046) q[0];
sx q[0];
rz(-1.4450547) q[0];
sx q[0];
rz(2.8618983) q[0];
rz(-1.1138629) q[1];
sx q[1];
rz(-1.056347) q[1];
sx q[1];
rz(2.9934771) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.516) q[0];
sx q[0];
rz(-0.91004386) q[0];
sx q[0];
rz(1.8126677) q[0];
rz(0.19377562) q[2];
sx q[2];
rz(-2.4420718) q[2];
sx q[2];
rz(-1.7278838) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.95947) q[1];
sx q[1];
rz(-2.39678) q[1];
sx q[1];
rz(-2.4306137) q[1];
rz(2.8612265) q[3];
sx q[3];
rz(-2.4383645) q[3];
sx q[3];
rz(0.91152465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1088045) q[2];
sx q[2];
rz(-0.84454909) q[2];
sx q[2];
rz(0.27842251) q[2];
rz(-1.6935211) q[3];
sx q[3];
rz(-1.672594) q[3];
sx q[3];
rz(1.9765123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096810452) q[0];
sx q[0];
rz(-0.7881931) q[0];
sx q[0];
rz(-1.0544554) q[0];
rz(1.0983559) q[1];
sx q[1];
rz(-1.337468) q[1];
sx q[1];
rz(-2.1877847) q[1];
rz(1.4020709) q[2];
sx q[2];
rz(-1.0015971) q[2];
sx q[2];
rz(0.88605273) q[2];
rz(-1.9270509) q[3];
sx q[3];
rz(-0.14394017) q[3];
sx q[3];
rz(-2.3059358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
