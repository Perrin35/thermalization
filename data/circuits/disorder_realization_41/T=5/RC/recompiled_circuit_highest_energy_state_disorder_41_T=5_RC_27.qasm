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
rz(-0.77223414) q[0];
sx q[0];
rz(-1.2152262) q[0];
sx q[0];
rz(-1.1507432) q[0];
rz(-0.28252217) q[1];
sx q[1];
rz(4.2673586) q[1];
sx q[1];
rz(10.537108) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38508666) q[0];
sx q[0];
rz(-1.6876564) q[0];
sx q[0];
rz(-2.8782842) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0289828) q[2];
sx q[2];
rz(-0.99163429) q[2];
sx q[2];
rz(-0.65161588) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9196786) q[1];
sx q[1];
rz(-2.3828607) q[1];
sx q[1];
rz(2.3054211) q[1];
x q[2];
rz(-2.2745773) q[3];
sx q[3];
rz(-2.2261282) q[3];
sx q[3];
rz(-2.1239779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34899601) q[2];
sx q[2];
rz(-2.0670321) q[2];
sx q[2];
rz(1.7462771) q[2];
rz(0.058517728) q[3];
sx q[3];
rz(-1.7250215) q[3];
sx q[3];
rz(-1.3099366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42404744) q[0];
sx q[0];
rz(-0.40322867) q[0];
sx q[0];
rz(0.44977093) q[0];
rz(2.6768118) q[1];
sx q[1];
rz(-1.3100781) q[1];
sx q[1];
rz(-0.25258499) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21097525) q[0];
sx q[0];
rz(-1.2196333) q[0];
sx q[0];
rz(3.1245078) q[0];
rz(-pi) q[1];
rz(2.4812883) q[2];
sx q[2];
rz(-1.2129158) q[2];
sx q[2];
rz(1.7972657) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.79073097) q[1];
sx q[1];
rz(-1.2892032) q[1];
sx q[1];
rz(3.1386203) q[1];
rz(-pi) q[2];
rz(2.1869867) q[3];
sx q[3];
rz(-0.86087117) q[3];
sx q[3];
rz(1.6740284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.8978591) q[2];
sx q[2];
rz(-2.2458138) q[2];
sx q[2];
rz(0.016156999) q[2];
rz(-0.92205087) q[3];
sx q[3];
rz(-2.1478896) q[3];
sx q[3];
rz(-3.0097358) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.048412662) q[0];
sx q[0];
rz(-1.6572297) q[0];
sx q[0];
rz(-0.68159252) q[0];
rz(2.5438578) q[1];
sx q[1];
rz(-1.455541) q[1];
sx q[1];
rz(-0.32424232) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77406835) q[0];
sx q[0];
rz(-1.9107959) q[0];
sx q[0];
rz(-0.46214477) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50679548) q[2];
sx q[2];
rz(-1.5118216) q[2];
sx q[2];
rz(0.61117327) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.26962599) q[1];
sx q[1];
rz(-1.5287279) q[1];
sx q[1];
rz(1.7394689) q[1];
rz(-pi) q[2];
rz(2.1799654) q[3];
sx q[3];
rz(-1.4318083) q[3];
sx q[3];
rz(-2.9442996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2250259) q[2];
sx q[2];
rz(-0.39174199) q[2];
sx q[2];
rz(-2.0908835) q[2];
rz(-2.4010036) q[3];
sx q[3];
rz(-1.2935484) q[3];
sx q[3];
rz(-0.061323969) q[3];
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
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4968313) q[0];
sx q[0];
rz(-2.3065541) q[0];
sx q[0];
rz(-2.187425) q[0];
rz(-2.0046115) q[1];
sx q[1];
rz(-1.6258207) q[1];
sx q[1];
rz(2.383393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9228464) q[0];
sx q[0];
rz(-0.2740261) q[0];
sx q[0];
rz(1.42204) q[0];
x q[1];
rz(1.3761282) q[2];
sx q[2];
rz(-2.3050781) q[2];
sx q[2];
rz(-1.0502953) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8710701) q[1];
sx q[1];
rz(-2.4409373) q[1];
sx q[1];
rz(1.4187705) q[1];
x q[2];
rz(-1.5605035) q[3];
sx q[3];
rz(-2.309628) q[3];
sx q[3];
rz(-2.9266299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0884462) q[2];
sx q[2];
rz(-2.5909178) q[2];
sx q[2];
rz(1.3147563) q[2];
rz(-0.30086073) q[3];
sx q[3];
rz(-1.8010537) q[3];
sx q[3];
rz(-0.68527591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.143173) q[0];
sx q[0];
rz(-1.5851861) q[0];
sx q[0];
rz(-1.6402624) q[0];
rz(2.7550664) q[1];
sx q[1];
rz(-2.0730348) q[1];
sx q[1];
rz(-1.5466669) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0879171) q[0];
sx q[0];
rz(-2.2639416) q[0];
sx q[0];
rz(-1.3455092) q[0];
rz(2.0208668) q[2];
sx q[2];
rz(-2.483027) q[2];
sx q[2];
rz(-2.3651802) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9068577) q[1];
sx q[1];
rz(-1.654792) q[1];
sx q[1];
rz(2.707721) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0647072) q[3];
sx q[3];
rz(-0.57514433) q[3];
sx q[3];
rz(-0.98371668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3950562) q[2];
sx q[2];
rz(-2.655513) q[2];
sx q[2];
rz(2.0716095) q[2];
rz(-0.88548958) q[3];
sx q[3];
rz(-2.077379) q[3];
sx q[3];
rz(-1.989367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92538658) q[0];
sx q[0];
rz(-0.33846551) q[0];
sx q[0];
rz(2.0881407) q[0];
rz(0.18928754) q[1];
sx q[1];
rz(-0.82866755) q[1];
sx q[1];
rz(-1.4922356) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5439592) q[0];
sx q[0];
rz(-1.2043556) q[0];
sx q[0];
rz(-2.9935915) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21250658) q[2];
sx q[2];
rz(-1.8706609) q[2];
sx q[2];
rz(2.1455255) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88613207) q[1];
sx q[1];
rz(-2.0719872) q[1];
sx q[1];
rz(1.4574128) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1023472) q[3];
sx q[3];
rz(-2.3330894) q[3];
sx q[3];
rz(-1.8448585) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3377127) q[2];
sx q[2];
rz(-2.1264075) q[2];
sx q[2];
rz(-2.0886776) q[2];
rz(-0.11689154) q[3];
sx q[3];
rz(-2.0785073) q[3];
sx q[3];
rz(-1.6803928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9601124) q[0];
sx q[0];
rz(-2.5616665) q[0];
sx q[0];
rz(2.5970698) q[0];
rz(-0.25360423) q[1];
sx q[1];
rz(-0.2519775) q[1];
sx q[1];
rz(0.25467083) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38886759) q[0];
sx q[0];
rz(-1.4988741) q[0];
sx q[0];
rz(-0.0062463721) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39310633) q[2];
sx q[2];
rz(-1.4908893) q[2];
sx q[2];
rz(0.42482947) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0510301) q[1];
sx q[1];
rz(-1.5080308) q[1];
sx q[1];
rz(2.323137) q[1];
x q[2];
rz(0.56470498) q[3];
sx q[3];
rz(-0.89313358) q[3];
sx q[3];
rz(1.1260179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2250681) q[2];
sx q[2];
rz(-2.2038286) q[2];
sx q[2];
rz(0.41898215) q[2];
rz(-0.032912832) q[3];
sx q[3];
rz(-1.907828) q[3];
sx q[3];
rz(-2.029443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(0.58881775) q[0];
sx q[0];
rz(-0.7426312) q[0];
sx q[0];
rz(1.3277998) q[0];
rz(-2.1090419) q[1];
sx q[1];
rz(-1.7951218) q[1];
sx q[1];
rz(2.5520777) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8701328) q[0];
sx q[0];
rz(-1.0268332) q[0];
sx q[0];
rz(-0.58596316) q[0];
rz(0.74331626) q[2];
sx q[2];
rz(-2.7545815) q[2];
sx q[2];
rz(-1.001337) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.81130109) q[1];
sx q[1];
rz(-1.9988235) q[1];
sx q[1];
rz(0.93702646) q[1];
x q[2];
rz(-1.0716075) q[3];
sx q[3];
rz(-1.2195865) q[3];
sx q[3];
rz(2.2413188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0699658) q[2];
sx q[2];
rz(-2.4667141) q[2];
sx q[2];
rz(0.91342941) q[2];
rz(-2.6895798) q[3];
sx q[3];
rz(-1.2830696) q[3];
sx q[3];
rz(1.7970596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87208676) q[0];
sx q[0];
rz(-1.0259314) q[0];
sx q[0];
rz(-1.6312697) q[0];
rz(-1.3144846) q[1];
sx q[1];
rz(-1.038082) q[1];
sx q[1];
rz(2.7297535) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7896693) q[0];
sx q[0];
rz(-2.4487307) q[0];
sx q[0];
rz(1.1293651) q[0];
rz(1.5225173) q[2];
sx q[2];
rz(-1.252151) q[2];
sx q[2];
rz(-1.4876721) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.85535586) q[1];
sx q[1];
rz(-2.0986631) q[1];
sx q[1];
rz(1.1638454) q[1];
rz(1.8985073) q[3];
sx q[3];
rz(-1.9372889) q[3];
sx q[3];
rz(1.9043363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0044378) q[2];
sx q[2];
rz(-1.4919446) q[2];
sx q[2];
rz(-0.23492661) q[2];
rz(-0.91819417) q[3];
sx q[3];
rz(-2.4166959) q[3];
sx q[3];
rz(-1.0880067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10201564) q[0];
sx q[0];
rz(-2.7871842) q[0];
sx q[0];
rz(1.2603731) q[0];
rz(1.6873987) q[1];
sx q[1];
rz(-1.6834384) q[1];
sx q[1];
rz(-1.6967324) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4871976) q[0];
sx q[0];
rz(-2.5154468) q[0];
sx q[0];
rz(-2.3930413) q[0];
rz(-pi) q[1];
rz(-2.1865322) q[2];
sx q[2];
rz(-0.79074016) q[2];
sx q[2];
rz(0.093971595) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9493172) q[1];
sx q[1];
rz(-2.2819464) q[1];
sx q[1];
rz(-1.7356731) q[1];
rz(0.71717307) q[3];
sx q[3];
rz(-0.48472084) q[3];
sx q[3];
rz(-2.6074099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5789648) q[2];
sx q[2];
rz(-0.29890385) q[2];
sx q[2];
rz(-3.01801) q[2];
rz(0.28892162) q[3];
sx q[3];
rz(-1.2762504) q[3];
sx q[3];
rz(2.6821274) q[3];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73744437) q[0];
sx q[0];
rz(-1.7973719) q[0];
sx q[0];
rz(1.6112882) q[0];
rz(0.96724802) q[1];
sx q[1];
rz(-2.1687242) q[1];
sx q[1];
rz(0.8868934) q[1];
rz(3.0616765) q[2];
sx q[2];
rz(-0.6007522) q[2];
sx q[2];
rz(-2.8795254) q[2];
rz(2.2933949) q[3];
sx q[3];
rz(-2.316939) q[3];
sx q[3];
rz(1.4311781) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
