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
rz(1.7088543) q[0];
sx q[0];
rz(-2.813485) q[0];
sx q[0];
rz(-0.83972591) q[0];
rz(1.4108763) q[1];
sx q[1];
rz(-1.5411935) q[1];
sx q[1];
rz(0.76959258) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5308204) q[0];
sx q[0];
rz(-1.4913173) q[0];
sx q[0];
rz(0.034850807) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5430449) q[2];
sx q[2];
rz(-0.74615462) q[2];
sx q[2];
rz(-2.1389824) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.85629) q[1];
sx q[1];
rz(-1.5701615) q[1];
sx q[1];
rz(1.2190422) q[1];
x q[2];
rz(-2.0887718) q[3];
sx q[3];
rz(-2.149579) q[3];
sx q[3];
rz(2.5797648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8771693) q[2];
sx q[2];
rz(-1.2494272) q[2];
sx q[2];
rz(1.1733615) q[2];
rz(2.7133283) q[3];
sx q[3];
rz(-0.56848017) q[3];
sx q[3];
rz(0.65304023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0083171) q[0];
sx q[0];
rz(-0.44047099) q[0];
sx q[0];
rz(1.0622729) q[0];
rz(1.5509037) q[1];
sx q[1];
rz(-0.56113243) q[1];
sx q[1];
rz(1.0947469) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4832981) q[0];
sx q[0];
rz(-1.6219369) q[0];
sx q[0];
rz(-1.3157033) q[0];
x q[1];
rz(1.0350448) q[2];
sx q[2];
rz(-2.3887815) q[2];
sx q[2];
rz(-2.8967146) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1909263) q[1];
sx q[1];
rz(-0.48787531) q[1];
sx q[1];
rz(0.64779727) q[1];
rz(-1.3571489) q[3];
sx q[3];
rz(-0.80677885) q[3];
sx q[3];
rz(2.2329604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67794472) q[2];
sx q[2];
rz(-2.8539168) q[2];
sx q[2];
rz(0.90052432) q[2];
rz(-1.0362222) q[3];
sx q[3];
rz(-0.13768727) q[3];
sx q[3];
rz(-3.1305967) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66505945) q[0];
sx q[0];
rz(-1.1256951) q[0];
sx q[0];
rz(1.4680468) q[0];
rz(2.2131069) q[1];
sx q[1];
rz(-2.1378345) q[1];
sx q[1];
rz(-1.4220062) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8967246) q[0];
sx q[0];
rz(-1.4804615) q[0];
sx q[0];
rz(-0.24855129) q[0];
x q[1];
rz(-2.3238238) q[2];
sx q[2];
rz(-1.0330794) q[2];
sx q[2];
rz(-1.0420711) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0735237) q[1];
sx q[1];
rz(-1.0367107) q[1];
sx q[1];
rz(1.992354) q[1];
x q[2];
rz(2.5676651) q[3];
sx q[3];
rz(-2.6477154) q[3];
sx q[3];
rz(-2.9536794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.15098393) q[2];
sx q[2];
rz(-1.7093806) q[2];
sx q[2];
rz(-0.81986156) q[2];
rz(0.12623434) q[3];
sx q[3];
rz(-1.9641967) q[3];
sx q[3];
rz(1.646515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6821297) q[0];
sx q[0];
rz(-1.7403025) q[0];
sx q[0];
rz(1.6023741) q[0];
rz(-1.0914717) q[1];
sx q[1];
rz(-2.3075054) q[1];
sx q[1];
rz(1.9713255) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051043432) q[0];
sx q[0];
rz(-1.9465595) q[0];
sx q[0];
rz(-2.7997478) q[0];
rz(1.865388) q[2];
sx q[2];
rz(-2.3995993) q[2];
sx q[2];
rz(0.8823673) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2728426) q[1];
sx q[1];
rz(-0.12050546) q[1];
sx q[1];
rz(0.7293479) q[1];
rz(-pi) q[2];
rz(0.50962944) q[3];
sx q[3];
rz(-0.42181236) q[3];
sx q[3];
rz(-3.075656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.9200661) q[2];
sx q[2];
rz(-0.92851323) q[2];
sx q[2];
rz(-0.083219223) q[2];
rz(0.65256882) q[3];
sx q[3];
rz(-0.021952732) q[3];
sx q[3];
rz(-0.86161247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6457152) q[0];
sx q[0];
rz(-0.28384122) q[0];
sx q[0];
rz(0.19677095) q[0];
rz(1.1605877) q[1];
sx q[1];
rz(-0.44191688) q[1];
sx q[1];
rz(1.388185) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3236448) q[0];
sx q[0];
rz(-0.16992386) q[0];
sx q[0];
rz(0.0090881149) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4193397) q[2];
sx q[2];
rz(-1.7709608) q[2];
sx q[2];
rz(1.6895001) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39972389) q[1];
sx q[1];
rz(-0.67662865) q[1];
sx q[1];
rz(-1.8477738) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84204414) q[3];
sx q[3];
rz(-1.8477412) q[3];
sx q[3];
rz(2.7634948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5768726) q[2];
sx q[2];
rz(-1.2083961) q[2];
sx q[2];
rz(1.8604856) q[2];
rz(2.7992904) q[3];
sx q[3];
rz(-0.050693158) q[3];
sx q[3];
rz(-2.2127693) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36011919) q[0];
sx q[0];
rz(-2.0442648) q[0];
sx q[0];
rz(-2.1088364) q[0];
rz(2.4646941) q[1];
sx q[1];
rz(-0.38336661) q[1];
sx q[1];
rz(-2.3597609) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9932258) q[0];
sx q[0];
rz(-1.9253823) q[0];
sx q[0];
rz(-0.605159) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1091613) q[2];
sx q[2];
rz(-1.5196745) q[2];
sx q[2];
rz(-1.0716719) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1429174) q[1];
sx q[1];
rz(-1.2641205) q[1];
sx q[1];
rz(-0.59345133) q[1];
x q[2];
rz(3.019612) q[3];
sx q[3];
rz(-2.8013419) q[3];
sx q[3];
rz(1.550479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8402164) q[2];
sx q[2];
rz(-0.8679114) q[2];
sx q[2];
rz(0.8737348) q[2];
rz(-0.78911632) q[3];
sx q[3];
rz(-1.5849761) q[3];
sx q[3];
rz(-0.48621392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9588722) q[0];
sx q[0];
rz(-3.0954269) q[0];
sx q[0];
rz(-3.0507372) q[0];
rz(0.60984045) q[1];
sx q[1];
rz(-1.5515386) q[1];
sx q[1];
rz(-0.59744936) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93865636) q[0];
sx q[0];
rz(-2.9466726) q[0];
sx q[0];
rz(2.418243) q[0];
x q[1];
rz(-1.1979073) q[2];
sx q[2];
rz(-1.2593566) q[2];
sx q[2];
rz(-0.58538891) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.3143464) q[1];
sx q[1];
rz(-1.1097849) q[1];
sx q[1];
rz(2.6516312) q[1];
rz(1.6555637) q[3];
sx q[3];
rz(-0.42094195) q[3];
sx q[3];
rz(1.3011025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.7647543) q[2];
sx q[2];
rz(-2.59616) q[2];
sx q[2];
rz(0.1405912) q[2];
rz(-1.8600672) q[3];
sx q[3];
rz(-1.3158653) q[3];
sx q[3];
rz(-2.6144821) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22290467) q[0];
sx q[0];
rz(-2.1253026) q[0];
sx q[0];
rz(-1.6131529) q[0];
rz(0.81360045) q[1];
sx q[1];
rz(-3.0366812) q[1];
sx q[1];
rz(-2.334107) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7225689) q[0];
sx q[0];
rz(-1.3815855) q[0];
sx q[0];
rz(1.36947) q[0];
rz(-pi) q[1];
rz(-3.0119704) q[2];
sx q[2];
rz(-2.0106914) q[2];
sx q[2];
rz(2.9578046) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2827816) q[1];
sx q[1];
rz(-0.75560299) q[1];
sx q[1];
rz(-2.1367461) q[1];
rz(-0.92370478) q[3];
sx q[3];
rz(-2.4974303) q[3];
sx q[3];
rz(1.8235064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7878824) q[2];
sx q[2];
rz(-1.7870125) q[2];
sx q[2];
rz(2.9969969) q[2];
rz(2.0374129) q[3];
sx q[3];
rz(-0.2140597) q[3];
sx q[3];
rz(1.0415227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5398194) q[0];
sx q[0];
rz(-1.1310534) q[0];
sx q[0];
rz(-0.55651504) q[0];
rz(0.76983184) q[1];
sx q[1];
rz(-0.12431215) q[1];
sx q[1];
rz(2.6803023) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249864) q[0];
sx q[0];
rz(-1.6184153) q[0];
sx q[0];
rz(2.141593) q[0];
rz(-pi) q[1];
x q[1];
rz(0.83314216) q[2];
sx q[2];
rz(-2.1066446) q[2];
sx q[2];
rz(-1.3539202) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0374679) q[1];
sx q[1];
rz(-0.79782432) q[1];
sx q[1];
rz(-2.7255157) q[1];
rz(2.5697124) q[3];
sx q[3];
rz(-1.2415621) q[3];
sx q[3];
rz(0.62599949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.43872675) q[2];
sx q[2];
rz(-2.8361969) q[2];
sx q[2];
rz(-2.3766282) q[2];
rz(-1.2139828) q[3];
sx q[3];
rz(-0.58911222) q[3];
sx q[3];
rz(-2.7629619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42534378) q[0];
sx q[0];
rz(-0.049976293) q[0];
sx q[0];
rz(0.36323994) q[0];
rz(-1.7323469) q[1];
sx q[1];
rz(-2.1556518) q[1];
sx q[1];
rz(2.9149616) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0186414) q[0];
sx q[0];
rz(-1.5737857) q[0];
sx q[0];
rz(-3.0063773) q[0];
x q[1];
rz(0.66618528) q[2];
sx q[2];
rz(-0.86951423) q[2];
sx q[2];
rz(1.2291069) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.94144635) q[1];
sx q[1];
rz(-2.2392803) q[1];
sx q[1];
rz(-2.723395) q[1];
rz(3.111242) q[3];
sx q[3];
rz(-2.8035946) q[3];
sx q[3];
rz(-0.81164593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5994485) q[2];
sx q[2];
rz(-0.52079529) q[2];
sx q[2];
rz(0.64860541) q[2];
rz(-0.058622807) q[3];
sx q[3];
rz(-2.5128745) q[3];
sx q[3];
rz(2.4962943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-0.36515737) q[0];
sx q[0];
rz(-1.2189652) q[0];
sx q[0];
rz(-0.49268876) q[0];
rz(2.0877214) q[1];
sx q[1];
rz(-0.55640472) q[1];
sx q[1];
rz(0.41592204) q[1];
rz(3.1171534) q[2];
sx q[2];
rz(-1.3346439) q[2];
sx q[2];
rz(1.1572184) q[2];
rz(-2.8063227) q[3];
sx q[3];
rz(-2.2380968) q[3];
sx q[3];
rz(2.5401881) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
