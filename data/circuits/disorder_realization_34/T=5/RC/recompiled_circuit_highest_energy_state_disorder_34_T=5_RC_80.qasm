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
rz(-1.7582769) q[0];
sx q[0];
rz(4.5780616) q[0];
sx q[0];
rz(8.4599001) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(-0.42958346) q[1];
sx q[1];
rz(2.092195) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9193662) q[0];
sx q[0];
rz(-0.99577409) q[0];
sx q[0];
rz(2.7336043) q[0];
x q[1];
rz(1.8951178) q[2];
sx q[2];
rz(-1.4359546) q[2];
sx q[2];
rz(0.91712778) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8559554) q[1];
sx q[1];
rz(-2.1826943) q[1];
sx q[1];
rz(-0.4680856) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.5792867) q[3];
sx q[3];
rz(-1.6484954) q[3];
sx q[3];
rz(-2.169211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.45399484) q[2];
sx q[2];
rz(-1.5291841) q[2];
sx q[2];
rz(0.018608658) q[2];
rz(-2.5971557) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(0.64793599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7740771) q[0];
sx q[0];
rz(-1.6870966) q[0];
sx q[0];
rz(0.53502214) q[0];
rz(-0.53994838) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(1.7417057) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8606621) q[0];
sx q[0];
rz(-1.1854404) q[0];
sx q[0];
rz(-1.0866665) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97008791) q[2];
sx q[2];
rz(-2.0510489) q[2];
sx q[2];
rz(-0.18829543) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8230086) q[1];
sx q[1];
rz(-0.46097791) q[1];
sx q[1];
rz(-0.58580841) q[1];
x q[2];
rz(0.73604354) q[3];
sx q[3];
rz(-0.76220817) q[3];
sx q[3];
rz(-0.30022844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0672368) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(-0.92607099) q[2];
rz(0.98008424) q[3];
sx q[3];
rz(-2.1772549) q[3];
sx q[3];
rz(-0.66942352) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7922908) q[0];
sx q[0];
rz(-1.8632357) q[0];
sx q[0];
rz(-1.6287623) q[0];
rz(-0.28383645) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(1.1121174) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84217269) q[0];
sx q[0];
rz(-1.5779691) q[0];
sx q[0];
rz(-3.1345063) q[0];
rz(3.0653333) q[2];
sx q[2];
rz(-1.6789376) q[2];
sx q[2];
rz(1.0534444) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6472811) q[1];
sx q[1];
rz(-2.0485281) q[1];
sx q[1];
rz(-1.7812438) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96235449) q[3];
sx q[3];
rz(-0.43399226) q[3];
sx q[3];
rz(-0.8784465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6185559) q[2];
sx q[2];
rz(-1.3448389) q[2];
sx q[2];
rz(-2.0396566) q[2];
rz(0.88895041) q[3];
sx q[3];
rz(-2.6933653) q[3];
sx q[3];
rz(0.86047188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.76982826) q[0];
sx q[0];
rz(-2.9673321) q[0];
sx q[0];
rz(-0.68921047) q[0];
rz(0.31461942) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(-1.7291732) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1293761) q[0];
sx q[0];
rz(-2.3123992) q[0];
sx q[0];
rz(-1.1393121) q[0];
x q[1];
rz(0.41823776) q[2];
sx q[2];
rz(-0.23242885) q[2];
sx q[2];
rz(-1.7247891) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3312348) q[1];
sx q[1];
rz(-2.7572933) q[1];
sx q[1];
rz(-0.15366252) q[1];
x q[2];
rz(-1.0869157) q[3];
sx q[3];
rz(-1.7431419) q[3];
sx q[3];
rz(-0.92007557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.41669258) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(-2.5856384) q[2];
rz(-1.8703095) q[3];
sx q[3];
rz(-1.2621745) q[3];
sx q[3];
rz(0.2909734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3274662) q[0];
sx q[0];
rz(-0.1411345) q[0];
sx q[0];
rz(-0.59447527) q[0];
rz(2.5796083) q[1];
sx q[1];
rz(-2.2832506) q[1];
sx q[1];
rz(2.1655703) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41244477) q[0];
sx q[0];
rz(-1.2093822) q[0];
sx q[0];
rz(-1.1114208) q[0];
rz(-pi) q[1];
rz(-0.030841737) q[2];
sx q[2];
rz(-0.75428666) q[2];
sx q[2];
rz(3.0109143) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0208066) q[1];
sx q[1];
rz(-0.73985064) q[1];
sx q[1];
rz(2.6915361) q[1];
rz(-pi) q[2];
rz(-2.8401883) q[3];
sx q[3];
rz(-1.3632953) q[3];
sx q[3];
rz(2.1077771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.48656616) q[2];
sx q[2];
rz(-1.8883294) q[2];
sx q[2];
rz(-2.5308934) q[2];
rz(2.8357909) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(-1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3139451) q[0];
sx q[0];
rz(-0.018445404) q[0];
sx q[0];
rz(1.6754643) q[0];
rz(0.77955359) q[1];
sx q[1];
rz(-1.164271) q[1];
sx q[1];
rz(0.73878845) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95979106) q[0];
sx q[0];
rz(-2.7914146) q[0];
sx q[0];
rz(-2.2018593) q[0];
rz(-pi) q[1];
rz(1.5601842) q[2];
sx q[2];
rz(-1.877575) q[2];
sx q[2];
rz(-2.5576484) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2062981) q[1];
sx q[1];
rz(-1.2064465) q[1];
sx q[1];
rz(-0.97779001) q[1];
x q[2];
rz(0.8035369) q[3];
sx q[3];
rz(-0.49852405) q[3];
sx q[3];
rz(-2.5132708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5256727) q[2];
sx q[2];
rz(-1.9024666) q[2];
sx q[2];
rz(2.2789148) q[2];
rz(-2.681813) q[3];
sx q[3];
rz(-1.6254057) q[3];
sx q[3];
rz(1.9988683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2343242) q[0];
sx q[0];
rz(-2.7643804) q[0];
sx q[0];
rz(-2.1211076) q[0];
rz(3.0842969) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(-2.3540156) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6660992) q[0];
sx q[0];
rz(-1.8084053) q[0];
sx q[0];
rz(2.4852024) q[0];
x q[1];
rz(1.4992981) q[2];
sx q[2];
rz(-1.8282969) q[2];
sx q[2];
rz(-1.6704287) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5503834) q[1];
sx q[1];
rz(-2.0616643) q[1];
sx q[1];
rz(0.20259133) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.034463783) q[3];
sx q[3];
rz(-0.98317819) q[3];
sx q[3];
rz(1.8922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3202177) q[2];
sx q[2];
rz(-1.3221952) q[2];
sx q[2];
rz(0.58464948) q[2];
rz(2.4892877) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(-1.8023796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898191) q[0];
sx q[0];
rz(-1.3404055) q[0];
sx q[0];
rz(-2.8133494) q[0];
rz(2.7499061) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(1.4788871) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97349629) q[0];
sx q[0];
rz(-2.7910821) q[0];
sx q[0];
rz(0.71061937) q[0];
x q[1];
rz(1.0974036) q[2];
sx q[2];
rz(-2.0461296) q[2];
sx q[2];
rz(-1.9394099) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6819671) q[1];
sx q[1];
rz(-2.8468968) q[1];
sx q[1];
rz(1.9035643) q[1];
rz(-pi) q[2];
rz(-2.0746873) q[3];
sx q[3];
rz(-2.7253838) q[3];
sx q[3];
rz(-1.4598626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11999232) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(2.3528986) q[2];
rz(-0.24104077) q[3];
sx q[3];
rz(-0.92625109) q[3];
sx q[3];
rz(-0.81361667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64860827) q[0];
sx q[0];
rz(-2.6032175) q[0];
sx q[0];
rz(-0.96187821) q[0];
rz(2.9073763) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(-0.73807565) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80958145) q[0];
sx q[0];
rz(-1.1763078) q[0];
sx q[0];
rz(-1.3138735) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63491027) q[2];
sx q[2];
rz(-1.9378086) q[2];
sx q[2];
rz(-0.5762595) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16633655) q[1];
sx q[1];
rz(-1.8615906) q[1];
sx q[1];
rz(1.1329805) q[1];
x q[2];
rz(2.4752615) q[3];
sx q[3];
rz(-1.7937346) q[3];
sx q[3];
rz(1.8458888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6568079) q[2];
sx q[2];
rz(-2.1190376) q[2];
sx q[2];
rz(1.2602932) q[2];
rz(1.8079181) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(2.8792152) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(2.2743478) q[0];
rz(2.1278837) q[1];
sx q[1];
rz(-1.3288682) q[1];
sx q[1];
rz(2.0713461) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7263171) q[0];
sx q[0];
rz(-1.4196102) q[0];
sx q[0];
rz(0.12355208) q[0];
rz(-2.3142368) q[2];
sx q[2];
rz(-2.1340886) q[2];
sx q[2];
rz(-0.87633301) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.50430272) q[1];
sx q[1];
rz(-0.66537332) q[1];
sx q[1];
rz(0.65267434) q[1];
x q[2];
rz(1.6813047) q[3];
sx q[3];
rz(-1.7095057) q[3];
sx q[3];
rz(-1.8294301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9762207) q[2];
sx q[2];
rz(-0.82087159) q[2];
sx q[2];
rz(1.2316068) q[2];
rz(1.5008789) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.83723849) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(-0.41481836) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(-0.11484405) q[2];
sx q[2];
rz(-0.44841246) q[2];
sx q[2];
rz(2.8186225) q[2];
rz(2.1555156) q[3];
sx q[3];
rz(-2.232983) q[3];
sx q[3];
rz(1.9017526) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
