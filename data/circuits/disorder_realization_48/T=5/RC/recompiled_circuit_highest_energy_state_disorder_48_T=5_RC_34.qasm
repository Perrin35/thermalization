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
rz(0.91520619) q[0];
sx q[0];
rz(2.6829166) q[0];
sx q[0];
rz(12.248326) q[0];
rz(7.6492352) q[1];
sx q[1];
rz(0.69861689) q[1];
sx q[1];
rz(3.0787992) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1284769) q[0];
sx q[0];
rz(-2.4139691) q[0];
sx q[0];
rz(0.81650556) q[0];
rz(2.096296) q[2];
sx q[2];
rz(-0.7514779) q[2];
sx q[2];
rz(-2.1906302) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5972728) q[1];
sx q[1];
rz(-1.8493127) q[1];
sx q[1];
rz(2.0803948) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1476521) q[3];
sx q[3];
rz(-1.2341043) q[3];
sx q[3];
rz(-1.5722881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.723145) q[2];
sx q[2];
rz(-1.8258839) q[2];
sx q[2];
rz(-1.0203699) q[2];
rz(2.4127035) q[3];
sx q[3];
rz(-2.708669) q[3];
sx q[3];
rz(1.6093904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.51434022) q[0];
sx q[0];
rz(-1.9753375) q[0];
sx q[0];
rz(-0.038272055) q[0];
rz(-3.017784) q[1];
sx q[1];
rz(-2.6688711) q[1];
sx q[1];
rz(1.5708057) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7109315) q[0];
sx q[0];
rz(-1.5607906) q[0];
sx q[0];
rz(1.5902993) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10199396) q[2];
sx q[2];
rz(-2.3677459) q[2];
sx q[2];
rz(2.4255861) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1473737) q[1];
sx q[1];
rz(-3.1404853) q[1];
sx q[1];
rz(1.963041) q[1];
rz(-pi) q[2];
rz(2.6255076) q[3];
sx q[3];
rz(-1.4654963) q[3];
sx q[3];
rz(-0.95296958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.52246419) q[2];
sx q[2];
rz(-1.8121441) q[2];
sx q[2];
rz(2.5627356) q[2];
rz(2.0959496) q[3];
sx q[3];
rz(-0.78543109) q[3];
sx q[3];
rz(-0.75700179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66913644) q[0];
sx q[0];
rz(-1.0772912) q[0];
sx q[0];
rz(-2.6876167) q[0];
rz(-0.16481915) q[1];
sx q[1];
rz(-2.3655128) q[1];
sx q[1];
rz(-1.4705315) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31904083) q[0];
sx q[0];
rz(-1.3994354) q[0];
sx q[0];
rz(-1.1771855) q[0];
rz(-1.0376588) q[2];
sx q[2];
rz(-2.5488944) q[2];
sx q[2];
rz(-1.8225841) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.20023055) q[1];
sx q[1];
rz(-0.75298568) q[1];
sx q[1];
rz(2.689207) q[1];
rz(1.0461058) q[3];
sx q[3];
rz(-2.1614822) q[3];
sx q[3];
rz(-2.5928796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37956023) q[2];
sx q[2];
rz(-1.6635868) q[2];
sx q[2];
rz(-1.3564159) q[2];
rz(0.49946579) q[3];
sx q[3];
rz(-1.6127337) q[3];
sx q[3];
rz(-1.0511901) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22719638) q[0];
sx q[0];
rz(-2.2375186) q[0];
sx q[0];
rz(1.0585349) q[0];
rz(1.9805485) q[1];
sx q[1];
rz(-2.5807022) q[1];
sx q[1];
rz(-0.11722359) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0441705) q[0];
sx q[0];
rz(-1.6099842) q[0];
sx q[0];
rz(2.5990361) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2160596) q[2];
sx q[2];
rz(-1.4129854) q[2];
sx q[2];
rz(2.7881546) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99771755) q[1];
sx q[1];
rz(-1.8492325) q[1];
sx q[1];
rz(-0.77634676) q[1];
x q[2];
rz(-2.6113308) q[3];
sx q[3];
rz(-0.70928364) q[3];
sx q[3];
rz(1.5269296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.78493541) q[2];
sx q[2];
rz(-1.6931809) q[2];
sx q[2];
rz(-1.7735205) q[2];
rz(-1.28537) q[3];
sx q[3];
rz(-1.1077489) q[3];
sx q[3];
rz(-0.72803298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.042498978) q[0];
sx q[0];
rz(-2.2137401) q[0];
sx q[0];
rz(1.4121144) q[0];
rz(2.1389351) q[1];
sx q[1];
rz(-0.24191562) q[1];
sx q[1];
rz(-1.5826506) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2501734) q[0];
sx q[0];
rz(-2.1262263) q[0];
sx q[0];
rz(2.069418) q[0];
rz(-pi) q[1];
rz(-3.0465165) q[2];
sx q[2];
rz(-2.4209967) q[2];
sx q[2];
rz(-2.528185) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6465368) q[1];
sx q[1];
rz(-0.55083129) q[1];
sx q[1];
rz(1.5257833) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.094314625) q[3];
sx q[3];
rz(-1.2504761) q[3];
sx q[3];
rz(-1.8220779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5530508) q[2];
sx q[2];
rz(-1.7091227) q[2];
sx q[2];
rz(-2.7480965) q[2];
rz(-0.95997512) q[3];
sx q[3];
rz(-0.67592755) q[3];
sx q[3];
rz(-0.967832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(1.8680962) q[0];
sx q[0];
rz(-3.0143026) q[0];
sx q[0];
rz(0.18336503) q[0];
rz(-0.49194899) q[1];
sx q[1];
rz(-0.79194561) q[1];
sx q[1];
rz(1.9221745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3122092) q[0];
sx q[0];
rz(-1.6512733) q[0];
sx q[0];
rz(0.089437251) q[0];
x q[1];
rz(3.1111818) q[2];
sx q[2];
rz(-1.6492427) q[2];
sx q[2];
rz(-2.4837361) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9907836) q[1];
sx q[1];
rz(-1.1649302) q[1];
sx q[1];
rz(-1.1836786) q[1];
x q[2];
rz(-3.0571396) q[3];
sx q[3];
rz(-1.1005529) q[3];
sx q[3];
rz(3.0601644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.1313974) q[2];
sx q[2];
rz(-1.4364028) q[2];
sx q[2];
rz(0.08237002) q[2];
rz(-2.2149337) q[3];
sx q[3];
rz(-0.72946531) q[3];
sx q[3];
rz(1.6389219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53967404) q[0];
sx q[0];
rz(-0.32783666) q[0];
sx q[0];
rz(-0.71640054) q[0];
rz(-2.7377103) q[1];
sx q[1];
rz(-1.5733746) q[1];
sx q[1];
rz(1.4366038) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8536896) q[0];
sx q[0];
rz(-2.7349964) q[0];
sx q[0];
rz(-1.9415744) q[0];
rz(-pi) q[1];
x q[1];
rz(0.2242078) q[2];
sx q[2];
rz(-0.53935236) q[2];
sx q[2];
rz(-1.8856848) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0246504) q[1];
sx q[1];
rz(-2.7305909) q[1];
sx q[1];
rz(0.25115369) q[1];
rz(-pi) q[2];
rz(0.8088423) q[3];
sx q[3];
rz(-2.1385962) q[3];
sx q[3];
rz(-2.0140935) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7776103) q[2];
sx q[2];
rz(-2.1668375) q[2];
sx q[2];
rz(3.0762365) q[2];
rz(0.86723793) q[3];
sx q[3];
rz(-1.6403653) q[3];
sx q[3];
rz(1.8170099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48677483) q[0];
sx q[0];
rz(-1.7205394) q[0];
sx q[0];
rz(-2.4105893) q[0];
rz(-0.77516088) q[1];
sx q[1];
rz(-0.33439264) q[1];
sx q[1];
rz(0.41466546) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99546826) q[0];
sx q[0];
rz(-2.099461) q[0];
sx q[0];
rz(0.84923262) q[0];
x q[1];
rz(1.783466) q[2];
sx q[2];
rz(-1.050282) q[2];
sx q[2];
rz(1.5368324) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.4463628) q[1];
sx q[1];
rz(-1.8703504) q[1];
sx q[1];
rz(-1.3079554) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7348792) q[3];
sx q[3];
rz(-1.7030431) q[3];
sx q[3];
rz(-1.924072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.23184648) q[2];
sx q[2];
rz(-2.6555588) q[2];
sx q[2];
rz(0.092378423) q[2];
rz(2.3653024) q[3];
sx q[3];
rz(-1.4624566) q[3];
sx q[3];
rz(2.9379454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48731503) q[0];
sx q[0];
rz(-1.4947816) q[0];
sx q[0];
rz(-3.1238632) q[0];
rz(-1.0461294) q[1];
sx q[1];
rz(-1.7283864) q[1];
sx q[1];
rz(-1.0606934) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3002787) q[0];
sx q[0];
rz(-1.5848241) q[0];
sx q[0];
rz(0.068723516) q[0];
rz(-1.9477884) q[2];
sx q[2];
rz(-0.62219884) q[2];
sx q[2];
rz(-2.5925328) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.13421397) q[1];
sx q[1];
rz(-1.3985279) q[1];
sx q[1];
rz(-2.3399347) q[1];
rz(0.10194998) q[3];
sx q[3];
rz(-2.6155439) q[3];
sx q[3];
rz(2.7395484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.1549418) q[2];
sx q[2];
rz(-2.3255746) q[2];
sx q[2];
rz(-2.4972534) q[2];
rz(-2.2996357) q[3];
sx q[3];
rz(-1.5717477) q[3];
sx q[3];
rz(-0.90698609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7613206) q[0];
sx q[0];
rz(-2.705882) q[0];
sx q[0];
rz(2.6897588) q[0];
rz(-1.0093581) q[1];
sx q[1];
rz(-1.511907) q[1];
sx q[1];
rz(-1.5712646) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0382045) q[0];
sx q[0];
rz(-1.7409835) q[0];
sx q[0];
rz(1.5769995) q[0];
rz(-3.0250835) q[2];
sx q[2];
rz(-1.0986058) q[2];
sx q[2];
rz(2.0906868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6856076) q[1];
sx q[1];
rz(-2.7311014) q[1];
sx q[1];
rz(-1.2430771) q[1];
rz(-3.1327469) q[3];
sx q[3];
rz(-2.4191751) q[3];
sx q[3];
rz(-2.9904108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5181804) q[2];
sx q[2];
rz(-1.0717012) q[2];
sx q[2];
rz(-1.5706459) q[2];
rz(3.0359641) q[3];
sx q[3];
rz(-2.195334) q[3];
sx q[3];
rz(2.6251729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8138206) q[0];
sx q[0];
rz(-2.1229424) q[0];
sx q[0];
rz(0.13737296) q[0];
rz(2.9227921) q[1];
sx q[1];
rz(-1.7411502) q[1];
sx q[1];
rz(2.2314744) q[1];
rz(-1.995718) q[2];
sx q[2];
rz(-0.47987249) q[2];
sx q[2];
rz(2.1940827) q[2];
rz(-2.6811213) q[3];
sx q[3];
rz(-2.4302767) q[3];
sx q[3];
rz(-2.5828491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
