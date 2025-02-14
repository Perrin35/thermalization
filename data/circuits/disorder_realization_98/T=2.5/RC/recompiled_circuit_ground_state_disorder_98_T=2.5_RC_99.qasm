OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.15220517) q[0];
sx q[0];
rz(-0.21259354) q[0];
sx q[0];
rz(-2.4074182) q[0];
rz(-1.9250159) q[1];
sx q[1];
rz(2.7956378) q[1];
sx q[1];
rz(12.591127) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7914331) q[0];
sx q[0];
rz(-1.6140249) q[0];
sx q[0];
rz(-0.062383609) q[0];
rz(-2.0169746) q[2];
sx q[2];
rz(-2.8834448) q[2];
sx q[2];
rz(0.1218957) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7721482) q[1];
sx q[1];
rz(-1.8696864) q[1];
sx q[1];
rz(-1.6256096) q[1];
x q[2];
rz(1.3583899) q[3];
sx q[3];
rz(-1.6254566) q[3];
sx q[3];
rz(2.2509991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2572702) q[2];
sx q[2];
rz(-1.4929644) q[2];
sx q[2];
rz(2.8327668) q[2];
rz(-0.8257927) q[3];
sx q[3];
rz(-1.0079931) q[3];
sx q[3];
rz(-0.76396137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9708213) q[0];
sx q[0];
rz(-1.229137) q[0];
sx q[0];
rz(1.0276851) q[0];
rz(-0.81614256) q[1];
sx q[1];
rz(-1.1183389) q[1];
sx q[1];
rz(-2.3291086) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.621429) q[0];
sx q[0];
rz(-1.7296711) q[0];
sx q[0];
rz(-1.9245252) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.036244321) q[2];
sx q[2];
rz(-2.7680021) q[2];
sx q[2];
rz(-0.61818733) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.634944) q[1];
sx q[1];
rz(-0.5310943) q[1];
sx q[1];
rz(1.8455161) q[1];
rz(1.2151617) q[3];
sx q[3];
rz(-1.2421248) q[3];
sx q[3];
rz(-1.5394913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5673148) q[2];
sx q[2];
rz(-1.4894166) q[2];
sx q[2];
rz(-0.91892773) q[2];
rz(-0.53141665) q[3];
sx q[3];
rz(-1.0454949) q[3];
sx q[3];
rz(-0.04118583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70319217) q[0];
sx q[0];
rz(-2.0719318) q[0];
sx q[0];
rz(-2.0042787) q[0];
rz(-1.3905585) q[1];
sx q[1];
rz(-2.5847692) q[1];
sx q[1];
rz(2.4086319) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.304564) q[0];
sx q[0];
rz(-2.6246855) q[0];
sx q[0];
rz(-2.5190973) q[0];
rz(-2.3172313) q[2];
sx q[2];
rz(-2.3869053) q[2];
sx q[2];
rz(-3.1122308) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.93948352) q[1];
sx q[1];
rz(-2.0401478) q[1];
sx q[1];
rz(-2.6497627) q[1];
rz(-0.0065782733) q[3];
sx q[3];
rz(-1.1220782) q[3];
sx q[3];
rz(-2.6980163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.0040032337) q[2];
sx q[2];
rz(-1.7294451) q[2];
sx q[2];
rz(2.1027193) q[2];
rz(2.1393356) q[3];
sx q[3];
rz(-1.301731) q[3];
sx q[3];
rz(-2.8514298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20659474) q[0];
sx q[0];
rz(-2.0552141) q[0];
sx q[0];
rz(0.90718734) q[0];
rz(-2.9838003) q[1];
sx q[1];
rz(-1.0148427) q[1];
sx q[1];
rz(1.8603604) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84128252) q[0];
sx q[0];
rz(-2.1774946) q[0];
sx q[0];
rz(-0.9228031) q[0];
x q[1];
rz(2.7428076) q[2];
sx q[2];
rz(-2.246309) q[2];
sx q[2];
rz(0.28093064) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5425749) q[1];
sx q[1];
rz(-2.3532824) q[1];
sx q[1];
rz(-2.5224204) q[1];
x q[2];
rz(2.2415555) q[3];
sx q[3];
rz(-2.0931912) q[3];
sx q[3];
rz(-2.1417422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24718757) q[2];
sx q[2];
rz(-0.93623585) q[2];
sx q[2];
rz(-1.8518764) q[2];
rz(-1.3442518) q[3];
sx q[3];
rz(-1.4569837) q[3];
sx q[3];
rz(-2.4414506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6808788) q[0];
sx q[0];
rz(-0.49837708) q[0];
sx q[0];
rz(2.1233249) q[0];
rz(2.0385888) q[1];
sx q[1];
rz(-2.4296727) q[1];
sx q[1];
rz(-1.8011372) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0809652) q[0];
sx q[0];
rz(-0.61765352) q[0];
sx q[0];
rz(2.4897442) q[0];
rz(-pi) q[1];
rz(2.7619038) q[2];
sx q[2];
rz(-1.6616529) q[2];
sx q[2];
rz(1.0041217) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3065222) q[1];
sx q[1];
rz(-1.7002652) q[1];
sx q[1];
rz(-3.0135327) q[1];
x q[2];
rz(-3.0839447) q[3];
sx q[3];
rz(-2.1103922) q[3];
sx q[3];
rz(0.12069139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5995784) q[2];
sx q[2];
rz(-0.88853637) q[2];
sx q[2];
rz(2.1395903) q[2];
rz(2.4501948) q[3];
sx q[3];
rz(-1.9611497) q[3];
sx q[3];
rz(-1.5207759) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5387251) q[0];
sx q[0];
rz(-1.069101) q[0];
sx q[0];
rz(-0.62527239) q[0];
rz(-1.193115) q[1];
sx q[1];
rz(-0.68980491) q[1];
sx q[1];
rz(0.56732059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35463542) q[0];
sx q[0];
rz(-0.79878317) q[0];
sx q[0];
rz(-1.3243933) q[0];
rz(1.5789323) q[2];
sx q[2];
rz(-0.74441351) q[2];
sx q[2];
rz(1.2071963) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9822944) q[1];
sx q[1];
rz(-1.59677) q[1];
sx q[1];
rz(2.1893255) q[1];
rz(2.9653699) q[3];
sx q[3];
rz(-1.7473975) q[3];
sx q[3];
rz(-0.11208243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83711964) q[2];
sx q[2];
rz(-1.2834872) q[2];
sx q[2];
rz(1.5597957) q[2];
rz(-2.5290153) q[3];
sx q[3];
rz(-1.160683) q[3];
sx q[3];
rz(-0.15577236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0953858) q[0];
sx q[0];
rz(-1.1908197) q[0];
sx q[0];
rz(1.1268536) q[0];
rz(-1.2696179) q[1];
sx q[1];
rz(-1.9536628) q[1];
sx q[1];
rz(-0.52106214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4241735) q[0];
sx q[0];
rz(-1.9195286) q[0];
sx q[0];
rz(-2.2762389) q[0];
rz(-0.032835788) q[2];
sx q[2];
rz(-1.8996138) q[2];
sx q[2];
rz(1.1177899) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6566297) q[1];
sx q[1];
rz(-2.1466675) q[1];
sx q[1];
rz(2.5099975) q[1];
rz(1.8463797) q[3];
sx q[3];
rz(-2.0949938) q[3];
sx q[3];
rz(0.39488246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0236987) q[2];
sx q[2];
rz(-1.3641027) q[2];
sx q[2];
rz(1.2269616) q[2];
rz(1.6795233) q[3];
sx q[3];
rz(-1.6242124) q[3];
sx q[3];
rz(2.7000361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(3.0284477) q[0];
sx q[0];
rz(-2.1118836) q[0];
sx q[0];
rz(0.5471158) q[0];
rz(0.088134915) q[1];
sx q[1];
rz(-1.793975) q[1];
sx q[1];
rz(2.7190582) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28332253) q[0];
sx q[0];
rz(-1.8517324) q[0];
sx q[0];
rz(-3.1188117) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6190622) q[2];
sx q[2];
rz(-1.2471419) q[2];
sx q[2];
rz(0.98060545) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.48848885) q[1];
sx q[1];
rz(-0.088989921) q[1];
sx q[1];
rz(-2.2527534) q[1];
rz(-pi) q[2];
rz(-1.3874153) q[3];
sx q[3];
rz(-1.7088582) q[3];
sx q[3];
rz(-0.79425298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9162468) q[2];
sx q[2];
rz(-1.5044745) q[2];
sx q[2];
rz(2.8543191) q[2];
rz(0.54287994) q[3];
sx q[3];
rz(-0.77016872) q[3];
sx q[3];
rz(0.058102593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3995689) q[0];
sx q[0];
rz(-2.4574807) q[0];
sx q[0];
rz(-2.3642819) q[0];
rz(-0.9264535) q[1];
sx q[1];
rz(-1.1421721) q[1];
sx q[1];
rz(1.3949589) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0844432) q[0];
sx q[0];
rz(-1.2532338) q[0];
sx q[0];
rz(1.1614413) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.59916536) q[2];
sx q[2];
rz(-1.3627421) q[2];
sx q[2];
rz(0.64917246) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3923126) q[1];
sx q[1];
rz(-2.5651076) q[1];
sx q[1];
rz(-2.2063401) q[1];
x q[2];
rz(2.3353102) q[3];
sx q[3];
rz(-2.3438934) q[3];
sx q[3];
rz(-2.9258941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7159783) q[2];
sx q[2];
rz(-1.7691111) q[2];
sx q[2];
rz(2.051029) q[2];
rz(-0.96450949) q[3];
sx q[3];
rz(-2.1301853) q[3];
sx q[3];
rz(1.2151659) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(2.5665117) q[0];
sx q[0];
rz(-2.8031741) q[0];
sx q[0];
rz(-1.5100719) q[0];
rz(0.034320023) q[1];
sx q[1];
rz(-1.8020554) q[1];
sx q[1];
rz(-2.7095749) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98429798) q[0];
sx q[0];
rz(-2.108886) q[0];
sx q[0];
rz(-3.0731766) q[0];
rz(-pi) q[1];
rz(-1.7909348) q[2];
sx q[2];
rz(-2.7825522) q[2];
sx q[2];
rz(-2.9901148) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4544983) q[1];
sx q[1];
rz(-1.7155572) q[1];
sx q[1];
rz(1.471038) q[1];
x q[2];
rz(-0.98938359) q[3];
sx q[3];
rz(-1.793878) q[3];
sx q[3];
rz(-1.7599966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4711275) q[2];
sx q[2];
rz(-1.1521143) q[2];
sx q[2];
rz(-1.0738037) q[2];
rz(-0.18887575) q[3];
sx q[3];
rz(-0.60868588) q[3];
sx q[3];
rz(0.3824189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043924532) q[0];
sx q[0];
rz(-2.8207939) q[0];
sx q[0];
rz(2.4378142) q[0];
rz(1.7882998) q[1];
sx q[1];
rz(-2.4663993) q[1];
sx q[1];
rz(-0.83723062) q[1];
rz(2.4953669) q[2];
sx q[2];
rz(-0.85474174) q[2];
sx q[2];
rz(-1.3893736) q[2];
rz(2.682229) q[3];
sx q[3];
rz(-1.6323327) q[3];
sx q[3];
rz(-1.0428693) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
