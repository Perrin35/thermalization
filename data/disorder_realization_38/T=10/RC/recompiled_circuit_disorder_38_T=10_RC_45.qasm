OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71198553) q[0];
sx q[0];
rz(-2.7349732) q[0];
sx q[0];
rz(-0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7101947) q[0];
sx q[0];
rz(-0.22127998) q[0];
sx q[0];
rz(-1.6943323) q[0];
rz(-pi) q[1];
rz(-2.1044188) q[2];
sx q[2];
rz(-1.954477) q[2];
sx q[2];
rz(-1.2984315) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7632335) q[1];
sx q[1];
rz(-0.86508703) q[1];
sx q[1];
rz(-2.6325429) q[1];
rz(-pi) q[2];
rz(-0.083505587) q[3];
sx q[3];
rz(-2.180035) q[3];
sx q[3];
rz(-3.1374251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7258519) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(2.501781) q[2];
rz(0.85302991) q[3];
sx q[3];
rz(-2.5363688) q[3];
sx q[3];
rz(2.7602592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8752276) q[0];
sx q[0];
rz(-1.9748283) q[0];
sx q[0];
rz(0.27045989) q[0];
rz(-0.71331435) q[1];
sx q[1];
rz(-1.0353054) q[1];
sx q[1];
rz(-1.6289904) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4166116) q[0];
sx q[0];
rz(-0.98326251) q[0];
sx q[0];
rz(2.7532817) q[0];
rz(-pi) q[1];
rz(-1.5741882) q[2];
sx q[2];
rz(-1.8545824) q[2];
sx q[2];
rz(-1.1018745) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0257033) q[1];
sx q[1];
rz(-1.70417) q[1];
sx q[1];
rz(-0.89110903) q[1];
x q[2];
rz(-0.23962044) q[3];
sx q[3];
rz(-1.2926971) q[3];
sx q[3];
rz(-0.77277377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9397395) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.4927031) q[3];
sx q[3];
rz(-0.025618205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-2.0443679) q[0];
sx q[0];
rz(-0.29552466) q[0];
rz(2.9064536) q[1];
sx q[1];
rz(-1.7087015) q[1];
sx q[1];
rz(-0.74584109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.196516) q[0];
sx q[0];
rz(-1.3498422) q[0];
sx q[0];
rz(3.0483079) q[0];
rz(-0.52650555) q[2];
sx q[2];
rz(-1.7130934) q[2];
sx q[2];
rz(-1.7775747) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5902139) q[1];
sx q[1];
rz(-2.5775902) q[1];
sx q[1];
rz(-0.53932921) q[1];
x q[2];
rz(-1.3219222) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(-0.53019023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.088034257) q[2];
sx q[2];
rz(-0.025667889) q[2];
sx q[2];
rz(-2.4528465) q[2];
rz(0.050343242) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(-1.6200199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.922309) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(-0.10198378) q[0];
rz(-0.12022262) q[1];
sx q[1];
rz(-2.6133803) q[1];
sx q[1];
rz(2.8682958) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8768339) q[0];
sx q[0];
rz(-2.8059373) q[0];
sx q[0];
rz(-1.8434974) q[0];
rz(-1.0068514) q[2];
sx q[2];
rz(-2.1401005) q[2];
sx q[2];
rz(-0.068892613) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8371493) q[1];
sx q[1];
rz(-1.7177561) q[1];
sx q[1];
rz(1.8069581) q[1];
x q[2];
rz(2.4920739) q[3];
sx q[3];
rz(-1.4276541) q[3];
sx q[3];
rz(-0.79917819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.3499202) q[2];
sx q[2];
rz(-2.4035954) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(3.062011) q[3];
sx q[3];
rz(-1.9888398) q[3];
sx q[3];
rz(-0.3751522) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(-3.058847) q[0];
rz(-2.4619608) q[1];
sx q[1];
rz(-1.6487164) q[1];
sx q[1];
rz(-2.1544429) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50901978) q[0];
sx q[0];
rz(-2.1633254) q[0];
sx q[0];
rz(0.49125262) q[0];
rz(-pi) q[1];
rz(-3.1290595) q[2];
sx q[2];
rz(-2.4977411) q[2];
sx q[2];
rz(-2.8963793) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6189177) q[1];
sx q[1];
rz(-0.72026157) q[1];
sx q[1];
rz(2.656225) q[1];
x q[2];
rz(1.6721252) q[3];
sx q[3];
rz(-2.7469027) q[3];
sx q[3];
rz(0.24995382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8905028) q[2];
sx q[2];
rz(-0.7901929) q[2];
sx q[2];
rz(2.9303072) q[2];
rz(-2.7206897) q[3];
sx q[3];
rz(-2.588721) q[3];
sx q[3];
rz(0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.485065) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(-0.791839) q[0];
rz(0.99545288) q[1];
sx q[1];
rz(-2.1738926) q[1];
sx q[1];
rz(1.8621559) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6364481) q[0];
sx q[0];
rz(-2.8594058) q[0];
sx q[0];
rz(1.5178773) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40462599) q[2];
sx q[2];
rz(-0.64794108) q[2];
sx q[2];
rz(-2.4857156) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.93950677) q[1];
sx q[1];
rz(-2.0923951) q[1];
sx q[1];
rz(0.034394666) q[1];
rz(-0.1412973) q[3];
sx q[3];
rz(-1.391489) q[3];
sx q[3];
rz(0.80207223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8906158) q[2];
sx q[2];
rz(-1.0303409) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(1.4701014) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32271785) q[0];
sx q[0];
rz(-0.6495496) q[0];
sx q[0];
rz(0.055667002) q[0];
rz(-2.9260013) q[1];
sx q[1];
rz(-0.76342738) q[1];
sx q[1];
rz(3.1380222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5768893) q[0];
sx q[0];
rz(-1.2029552) q[0];
sx q[0];
rz(-0.44643114) q[0];
rz(-2.7262906) q[2];
sx q[2];
rz(-2.3447782) q[2];
sx q[2];
rz(1.0559527) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7230941) q[1];
sx q[1];
rz(-1.2730036) q[1];
sx q[1];
rz(0.94350068) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3725029) q[3];
sx q[3];
rz(-1.6453504) q[3];
sx q[3];
rz(1.9004746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5443762) q[2];
sx q[2];
rz(-2.1901013) q[2];
sx q[2];
rz(-0.92010951) q[2];
rz(-2.9428633) q[3];
sx q[3];
rz(-1.2343497) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51628095) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(1.4177119) q[0];
rz(-2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(-2.4628941) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5355276) q[0];
sx q[0];
rz(-1.3074271) q[0];
sx q[0];
rz(-2.8019964) q[0];
x q[1];
rz(1.0661725) q[2];
sx q[2];
rz(-0.80427158) q[2];
sx q[2];
rz(1.417516) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8933967) q[1];
sx q[1];
rz(-1.0475323) q[1];
sx q[1];
rz(-0.072695331) q[1];
rz(0.33903665) q[3];
sx q[3];
rz(-0.44796523) q[3];
sx q[3];
rz(-0.26089222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9644908) q[2];
sx q[2];
rz(-2.047838) q[2];
sx q[2];
rz(-0.91782451) q[2];
rz(1.9994036) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-0.93723047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98638242) q[0];
sx q[0];
rz(-2.4676403) q[0];
sx q[0];
rz(-0.25979364) q[0];
rz(2.4329176) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5567112) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(0.35993872) q[0];
x q[1];
rz(-1.9138463) q[2];
sx q[2];
rz(-0.68708778) q[2];
sx q[2];
rz(-0.81817852) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.867013) q[1];
sx q[1];
rz(-1.8362987) q[1];
sx q[1];
rz(0.32151476) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5832289) q[3];
sx q[3];
rz(-2.6287968) q[3];
sx q[3];
rz(-2.5496897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.32593411) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(-0.79088598) q[2];
rz(0.38665006) q[3];
sx q[3];
rz(-1.0145885) q[3];
sx q[3];
rz(-0.33716831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40795046) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-0.98544425) q[0];
rz(2.573029) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(2.7808166) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11689582) q[0];
sx q[0];
rz(-1.7381867) q[0];
sx q[0];
rz(-0.025339729) q[0];
rz(-pi) q[1];
rz(-3.0214494) q[2];
sx q[2];
rz(-1.3573682) q[2];
sx q[2];
rz(-0.10211589) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.68041486) q[1];
sx q[1];
rz(-2.3056681) q[1];
sx q[1];
rz(1.3780891) q[1];
x q[2];
rz(-1.7182699) q[3];
sx q[3];
rz(-1.4031938) q[3];
sx q[3];
rz(-2.9180805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0439904) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(-2.1981751) q[2];
rz(0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(1.7452314) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5744793) q[0];
sx q[0];
rz(-1.4470826) q[0];
sx q[0];
rz(-0.8599109) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(1.1311244) q[2];
sx q[2];
rz(-2.9154073) q[2];
sx q[2];
rz(0.67868457) q[2];
rz(-2.2684569) q[3];
sx q[3];
rz(-0.94948873) q[3];
sx q[3];
rz(-0.82617847) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
