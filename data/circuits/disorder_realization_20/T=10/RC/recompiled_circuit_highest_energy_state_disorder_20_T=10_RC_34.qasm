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
rz(-2.1751997) q[0];
sx q[0];
rz(4.8973358) q[0];
sx q[0];
rz(10.020221) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(-0.59102494) q[1];
sx q[1];
rz(0.12330595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6791413) q[0];
sx q[0];
rz(-1.9872905) q[0];
sx q[0];
rz(0.97521675) q[0];
x q[1];
rz(0.78211981) q[2];
sx q[2];
rz(-2.189866) q[2];
sx q[2];
rz(1.9591046) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6001624) q[1];
sx q[1];
rz(-1.7847536) q[1];
sx q[1];
rz(-0.096264953) q[1];
rz(-pi) q[2];
rz(-1.418258) q[3];
sx q[3];
rz(-1.2087052) q[3];
sx q[3];
rz(0.83521508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3143602) q[2];
sx q[2];
rz(-1.102697) q[2];
sx q[2];
rz(-0.11191351) q[2];
rz(-1.7498451) q[3];
sx q[3];
rz(-1.1575969) q[3];
sx q[3];
rz(2.8351423) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4157442) q[0];
sx q[0];
rz(-0.84472504) q[0];
sx q[0];
rz(0.14990212) q[0];
rz(-2.8556178) q[1];
sx q[1];
rz(-1.7466702) q[1];
sx q[1];
rz(-1.1289977) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1492831) q[0];
sx q[0];
rz(-2.8516927) q[0];
sx q[0];
rz(1.6768181) q[0];
rz(-pi) q[1];
rz(-1.8104042) q[2];
sx q[2];
rz(-0.061366038) q[2];
sx q[2];
rz(2.1642016) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4399229) q[1];
sx q[1];
rz(-1.4050802) q[1];
sx q[1];
rz(2.7782337) q[1];
rz(-pi) q[2];
rz(2.1993447) q[3];
sx q[3];
rz(-1.0744922) q[3];
sx q[3];
rz(1.4013578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4631606) q[2];
sx q[2];
rz(-2.2191935) q[2];
sx q[2];
rz(1.9909667) q[2];
rz(1.1848508) q[3];
sx q[3];
rz(-1.7246282) q[3];
sx q[3];
rz(0.020603389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6612369) q[0];
sx q[0];
rz(-2.895597) q[0];
sx q[0];
rz(0.38159698) q[0];
rz(-1.2941788) q[1];
sx q[1];
rz(-1.1062063) q[1];
sx q[1];
rz(2.9433184) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7679032) q[0];
sx q[0];
rz(-1.7977909) q[0];
sx q[0];
rz(0.17902184) q[0];
x q[1];
rz(1.4691822) q[2];
sx q[2];
rz(-0.67199113) q[2];
sx q[2];
rz(-0.72173126) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5570017) q[1];
sx q[1];
rz(-1.3583359) q[1];
sx q[1];
rz(0.23834385) q[1];
rz(3.1214633) q[3];
sx q[3];
rz(-1.3960037) q[3];
sx q[3];
rz(-2.54769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.76564378) q[2];
sx q[2];
rz(-2.9889034) q[2];
sx q[2];
rz(0.51885968) q[2];
rz(-1.2505924) q[3];
sx q[3];
rz(-1.0873245) q[3];
sx q[3];
rz(0.58108228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4984703) q[0];
sx q[0];
rz(-0.20007087) q[0];
sx q[0];
rz(-2.6923687) q[0];
rz(1.6953702) q[1];
sx q[1];
rz(-2.4999373) q[1];
sx q[1];
rz(-2.5208688) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8410438) q[0];
sx q[0];
rz(-2.4660206) q[0];
sx q[0];
rz(-2.9094957) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63346699) q[2];
sx q[2];
rz(-1.6704428) q[2];
sx q[2];
rz(1.5758621) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6606969) q[1];
sx q[1];
rz(-2.6891853) q[1];
sx q[1];
rz(-0.47641944) q[1];
rz(-2.1093217) q[3];
sx q[3];
rz(-1.6163975) q[3];
sx q[3];
rz(-0.63889438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0351403) q[2];
sx q[2];
rz(-1.2713212) q[2];
sx q[2];
rz(2.980496) q[2];
rz(0.39786878) q[3];
sx q[3];
rz(-1.6631923) q[3];
sx q[3];
rz(-2.9919992) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75905269) q[0];
sx q[0];
rz(-1.7790786) q[0];
sx q[0];
rz(0.25704849) q[0];
rz(-0.88047475) q[1];
sx q[1];
rz(-0.6404666) q[1];
sx q[1];
rz(1.2535198) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.391025) q[0];
sx q[0];
rz(-1.5526958) q[0];
sx q[0];
rz(-1.6068293) q[0];
rz(0.23373105) q[2];
sx q[2];
rz(-0.82373754) q[2];
sx q[2];
rz(-1.1540138) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1977928) q[1];
sx q[1];
rz(-1.5256923) q[1];
sx q[1];
rz(2.3438441) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.569031) q[3];
sx q[3];
rz(-1.2236508) q[3];
sx q[3];
rz(-0.4345168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4132061) q[2];
sx q[2];
rz(-1.2089968) q[2];
sx q[2];
rz(-1.2665292) q[2];
rz(-0.62147102) q[3];
sx q[3];
rz(-2.5340243) q[3];
sx q[3];
rz(-2.6004041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40015873) q[0];
sx q[0];
rz(-0.83916894) q[0];
sx q[0];
rz(-0.025064502) q[0];
rz(1.9301682) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(0.2549583) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6699885) q[0];
sx q[0];
rz(-1.5657664) q[0];
sx q[0];
rz(0.049335376) q[0];
x q[1];
rz(1.0827052) q[2];
sx q[2];
rz(-2.1016663) q[2];
sx q[2];
rz(1.7386029) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7073696) q[1];
sx q[1];
rz(-1.0138144) q[1];
sx q[1];
rz(-2.2051478) q[1];
rz(-2.1433349) q[3];
sx q[3];
rz(-1.7504331) q[3];
sx q[3];
rz(1.4416173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.24594626) q[2];
sx q[2];
rz(-2.0993555) q[2];
sx q[2];
rz(0.45905217) q[2];
rz(1.4970655) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(-0.12115255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6727305) q[0];
sx q[0];
rz(-0.48968306) q[0];
sx q[0];
rz(-0.086061867) q[0];
rz(0.91880265) q[1];
sx q[1];
rz(-1.1163534) q[1];
sx q[1];
rz(-1.930621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0897601) q[0];
sx q[0];
rz(-0.38516949) q[0];
sx q[0];
rz(2.5776107) q[0];
x q[1];
rz(2.1384905) q[2];
sx q[2];
rz(-1.9493812) q[2];
sx q[2];
rz(2.1673983) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0667292) q[1];
sx q[1];
rz(-1.7742429) q[1];
sx q[1];
rz(-1.2438846) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3567621) q[3];
sx q[3];
rz(-2.6468224) q[3];
sx q[3];
rz(3.1113603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3226037) q[2];
sx q[2];
rz(-2.6456867) q[2];
sx q[2];
rz(-0.71713478) q[2];
rz(2.1142193) q[3];
sx q[3];
rz(-1.4077978) q[3];
sx q[3];
rz(-2.7023442) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1502007) q[0];
sx q[0];
rz(-2.1450277) q[0];
sx q[0];
rz(0.014658654) q[0];
rz(0.75153366) q[1];
sx q[1];
rz(-1.9207759) q[1];
sx q[1];
rz(1.470648) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6122307) q[0];
sx q[0];
rz(-2.9069773) q[0];
sx q[0];
rz(0.63572065) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4906497) q[2];
sx q[2];
rz(-1.8266457) q[2];
sx q[2];
rz(0.7776153) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.70794174) q[1];
sx q[1];
rz(-1.4580863) q[1];
sx q[1];
rz(1.9757084) q[1];
rz(-pi) q[2];
rz(1.5812381) q[3];
sx q[3];
rz(-2.7964948) q[3];
sx q[3];
rz(-2.0562248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.26059255) q[2];
sx q[2];
rz(-1.1129817) q[2];
sx q[2];
rz(-2.8680958) q[2];
rz(1.986844) q[3];
sx q[3];
rz(-0.65360779) q[3];
sx q[3];
rz(2.4912513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038789373) q[0];
sx q[0];
rz(-2.0162855) q[0];
sx q[0];
rz(-2.0661085) q[0];
rz(-3.0134046) q[1];
sx q[1];
rz(-2.2905541) q[1];
sx q[1];
rz(-0.34559616) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3227279) q[0];
sx q[0];
rz(-2.4562307) q[0];
sx q[0];
rz(-1.3166053) q[0];
x q[1];
rz(-2.1394452) q[2];
sx q[2];
rz(-1.2969742) q[2];
sx q[2];
rz(-1.3783) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.83452713) q[1];
sx q[1];
rz(-2.5967) q[1];
sx q[1];
rz(0.12172555) q[1];
rz(-pi) q[2];
rz(-0.59806268) q[3];
sx q[3];
rz(-1.7081714) q[3];
sx q[3];
rz(-0.029702317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0547611) q[2];
sx q[2];
rz(-2.0589477) q[2];
sx q[2];
rz(1.6205988) q[2];
rz(-1.7704376) q[3];
sx q[3];
rz(-0.75826472) q[3];
sx q[3];
rz(0.41485205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82776752) q[0];
sx q[0];
rz(-2.1791552) q[0];
sx q[0];
rz(0.23571043) q[0];
rz(2.56855) q[1];
sx q[1];
rz(-1.0083116) q[1];
sx q[1];
rz(-1.3689573) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.154207) q[0];
sx q[0];
rz(-1.9999749) q[0];
sx q[0];
rz(-2.3230419) q[0];
rz(-pi) q[1];
rz(-0.67566411) q[2];
sx q[2];
rz(-2.0620637) q[2];
sx q[2];
rz(1.5671687) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0693501) q[1];
sx q[1];
rz(-1.4003039) q[1];
sx q[1];
rz(-0.99872288) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4769796) q[3];
sx q[3];
rz(-1.2117077) q[3];
sx q[3];
rz(2.1916568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3706751) q[2];
sx q[2];
rz(-1.3266027) q[2];
sx q[2];
rz(3.0885922) q[2];
rz(-0.75183374) q[3];
sx q[3];
rz(-1.0337318) q[3];
sx q[3];
rz(-0.92541614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62405217) q[0];
sx q[0];
rz(-1.0160099) q[0];
sx q[0];
rz(2.1852063) q[0];
rz(1.3507631) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(-2.4443632) q[2];
sx q[2];
rz(-1.6710499) q[2];
sx q[2];
rz(1.2895126) q[2];
rz(0.57784537) q[3];
sx q[3];
rz(-2.2457849) q[3];
sx q[3];
rz(-0.59413322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
