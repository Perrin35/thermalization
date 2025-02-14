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
rz(-1.1308489) q[0];
sx q[0];
rz(4.9133416) q[0];
sx q[0];
rz(9.2815444) q[0];
rz(3.1205966) q[1];
sx q[1];
rz(-2.5591873) q[1];
sx q[1];
rz(2.6666759) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.582481) q[0];
sx q[0];
rz(-1.5149635) q[0];
sx q[0];
rz(0.076095993) q[0];
rz(-pi) q[1];
rz(-0.61491809) q[2];
sx q[2];
rz(-1.5919098) q[2];
sx q[2];
rz(2.9647765) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.21771297) q[1];
sx q[1];
rz(-0.7860187) q[1];
sx q[1];
rz(1.7255369) q[1];
rz(-2.3793061) q[3];
sx q[3];
rz(-0.62847465) q[3];
sx q[3];
rz(-2.3746128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1249866) q[2];
sx q[2];
rz(-1.6281444) q[2];
sx q[2];
rz(-2.5561257) q[2];
rz(2.3440907) q[3];
sx q[3];
rz(-1.5471285) q[3];
sx q[3];
rz(-0.40369478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1095235) q[0];
sx q[0];
rz(-1.7970947) q[0];
sx q[0];
rz(2.4865785) q[0];
rz(-0.0034927448) q[1];
sx q[1];
rz(-2.4066636) q[1];
sx q[1];
rz(-0.014911501) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92214092) q[0];
sx q[0];
rz(-0.08481124) q[0];
sx q[0];
rz(-2.831859) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2737421) q[2];
sx q[2];
rz(-1.2406435) q[2];
sx q[2];
rz(1.7148266) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2567581) q[1];
sx q[1];
rz(-1.2169098) q[1];
sx q[1];
rz(-2.3164151) q[1];
rz(-pi) q[2];
rz(0.80958812) q[3];
sx q[3];
rz(-1.6177655) q[3];
sx q[3];
rz(-1.498158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.94987503) q[2];
sx q[2];
rz(-1.8052552) q[2];
sx q[2];
rz(-2.8548262) q[2];
rz(0.030335434) q[3];
sx q[3];
rz(-1.5807296) q[3];
sx q[3];
rz(1.102977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2716118) q[0];
sx q[0];
rz(-2.5219707) q[0];
sx q[0];
rz(-1.9770589) q[0];
rz(-2.8097235) q[1];
sx q[1];
rz(-1.429456) q[1];
sx q[1];
rz(1.1480931) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5184421) q[0];
sx q[0];
rz(-1.7887576) q[0];
sx q[0];
rz(-0.017780546) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16638799) q[2];
sx q[2];
rz(-1.0999318) q[2];
sx q[2];
rz(-1.3883049) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5268574) q[1];
sx q[1];
rz(-1.2678483) q[1];
sx q[1];
rz(-2.9399859) q[1];
rz(1.7487995) q[3];
sx q[3];
rz(-1.5220259) q[3];
sx q[3];
rz(-2.3017987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0120734) q[2];
sx q[2];
rz(-1.4612863) q[2];
sx q[2];
rz(2.6626383) q[2];
rz(-0.78914133) q[3];
sx q[3];
rz(-2.456587) q[3];
sx q[3];
rz(-1.7710549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0399465) q[0];
sx q[0];
rz(-2.3532823) q[0];
sx q[0];
rz(2.3729861) q[0];
rz(0.45306122) q[1];
sx q[1];
rz(-1.0337579) q[1];
sx q[1];
rz(2.5624842) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3192037) q[0];
sx q[0];
rz(-1.7928078) q[0];
sx q[0];
rz(2.8846106) q[0];
rz(-1.0592346) q[2];
sx q[2];
rz(-1.2280012) q[2];
sx q[2];
rz(0.9139708) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.37984338) q[1];
sx q[1];
rz(-2.6452612) q[1];
sx q[1];
rz(1.6654254) q[1];
rz(-pi) q[2];
rz(-1.7594537) q[3];
sx q[3];
rz(-1.9629729) q[3];
sx q[3];
rz(-0.79526633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96827489) q[2];
sx q[2];
rz(-1.3710794) q[2];
sx q[2];
rz(-1.0260065) q[2];
rz(-2.6722028) q[3];
sx q[3];
rz(-1.0204851) q[3];
sx q[3];
rz(-0.49218407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.39407179) q[0];
sx q[0];
rz(-1.515027) q[0];
sx q[0];
rz(1.0006022) q[0];
rz(-1.8383149) q[1];
sx q[1];
rz(-2.3293827) q[1];
sx q[1];
rz(-2.2607048) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39932809) q[0];
sx q[0];
rz(-0.77631809) q[0];
sx q[0];
rz(0.72350435) q[0];
x q[1];
rz(2.8116436) q[2];
sx q[2];
rz(-2.4545728) q[2];
sx q[2];
rz(-0.52345146) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8288159) q[1];
sx q[1];
rz(-0.79598266) q[1];
sx q[1];
rz(0.46013855) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0468454) q[3];
sx q[3];
rz(-1.928471) q[3];
sx q[3];
rz(-0.053225191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96054968) q[2];
sx q[2];
rz(-1.1134104) q[2];
sx q[2];
rz(1.95365) q[2];
rz(-0.12399593) q[3];
sx q[3];
rz(-1.8143727) q[3];
sx q[3];
rz(-2.7441062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(0.7749629) q[0];
sx q[0];
rz(-2.9149945) q[0];
sx q[0];
rz(2.9887548) q[0];
rz(-0.14640181) q[1];
sx q[1];
rz(-1.791879) q[1];
sx q[1];
rz(2.4375367) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.31788) q[0];
sx q[0];
rz(-1.5410265) q[0];
sx q[0];
rz(-1.5533621) q[0];
rz(0.91135773) q[2];
sx q[2];
rz(-1.8375085) q[2];
sx q[2];
rz(3.0195723) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2347243) q[1];
sx q[1];
rz(-0.14130302) q[1];
sx q[1];
rz(1.0060746) q[1];
rz(-0.95100286) q[3];
sx q[3];
rz(-1.8547748) q[3];
sx q[3];
rz(-1.2172446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0600837) q[2];
sx q[2];
rz(-2.1174105) q[2];
sx q[2];
rz(2.7626959) q[2];
rz(-2.1156408) q[3];
sx q[3];
rz(-0.70464269) q[3];
sx q[3];
rz(-1.1161233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41554552) q[0];
sx q[0];
rz(-1.8141831) q[0];
sx q[0];
rz(2.5782247) q[0];
rz(-2.7401961) q[1];
sx q[1];
rz(-2.1386264) q[1];
sx q[1];
rz(-2.4664403) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0289265) q[0];
sx q[0];
rz(-3.108041) q[0];
sx q[0];
rz(-1.9153509) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60826081) q[2];
sx q[2];
rz(-1.741) q[2];
sx q[2];
rz(0.93675429) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.25738379) q[1];
sx q[1];
rz(-1.3335449) q[1];
sx q[1];
rz(-1.9006185) q[1];
rz(-pi) q[2];
x q[2];
rz(1.160333) q[3];
sx q[3];
rz(-0.77926765) q[3];
sx q[3];
rz(-1.1469649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.99951619) q[2];
sx q[2];
rz(-2.2579305) q[2];
sx q[2];
rz(-0.42259541) q[2];
rz(-2.8790867) q[3];
sx q[3];
rz(-2.0764949) q[3];
sx q[3];
rz(2.0601823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.066684) q[0];
sx q[0];
rz(-2.1796362) q[0];
sx q[0];
rz(-2.1426376) q[0];
rz(1.2133489) q[1];
sx q[1];
rz(-1.1937001) q[1];
sx q[1];
rz(-2.8003069) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0851819) q[0];
sx q[0];
rz(-0.23427948) q[0];
sx q[0];
rz(0.60472576) q[0];
rz(-2.2407124) q[2];
sx q[2];
rz(-1.7770801) q[2];
sx q[2];
rz(-1.2717977) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.44036699) q[1];
sx q[1];
rz(-1.4557252) q[1];
sx q[1];
rz(-0.23903317) q[1];
x q[2];
rz(-3.083087) q[3];
sx q[3];
rz(-1.2924177) q[3];
sx q[3];
rz(-2.6191421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21218941) q[2];
sx q[2];
rz(-1.5430278) q[2];
sx q[2];
rz(-0.25941485) q[2];
rz(0.2855531) q[3];
sx q[3];
rz(-2.4694337) q[3];
sx q[3];
rz(-1.0478421) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24006537) q[0];
sx q[0];
rz(-0.49235383) q[0];
sx q[0];
rz(2.0557851) q[0];
rz(-2.1366513) q[1];
sx q[1];
rz(-2.0045547) q[1];
sx q[1];
rz(2.5850632) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96592007) q[0];
sx q[0];
rz(-1.5300599) q[0];
sx q[0];
rz(1.7044742) q[0];
x q[1];
rz(-2.506571) q[2];
sx q[2];
rz(-2.1974034) q[2];
sx q[2];
rz(2.602586) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6145912) q[1];
sx q[1];
rz(-1.6629991) q[1];
sx q[1];
rz(-2.8648977) q[1];
rz(0.50332467) q[3];
sx q[3];
rz(-2.0317114) q[3];
sx q[3];
rz(0.028655076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0050547) q[2];
sx q[2];
rz(-1.7359066) q[2];
sx q[2];
rz(1.271099) q[2];
rz(-0.84730411) q[3];
sx q[3];
rz(-1.3264791) q[3];
sx q[3];
rz(2.8443851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.51410455) q[0];
sx q[0];
rz(-1.4493554) q[0];
sx q[0];
rz(-2.3741212) q[0];
rz(-2.1766359) q[1];
sx q[1];
rz(-1.8582397) q[1];
sx q[1];
rz(-2.8585785) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8022158) q[0];
sx q[0];
rz(-2.2664323) q[0];
sx q[0];
rz(-2.8921769) q[0];
rz(0.04472132) q[2];
sx q[2];
rz(-0.43202094) q[2];
sx q[2];
rz(0.17145874) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9099847) q[1];
sx q[1];
rz(-2.5043813) q[1];
sx q[1];
rz(2.5093394) q[1];
x q[2];
rz(1.0177294) q[3];
sx q[3];
rz(-0.8700287) q[3];
sx q[3];
rz(-1.7297945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.1538126) q[2];
sx q[2];
rz(-2.751838) q[2];
sx q[2];
rz(1.4022931) q[2];
rz(2.4757989) q[3];
sx q[3];
rz(-1.8699162) q[3];
sx q[3];
rz(1.8388892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.5182198) q[0];
sx q[0];
rz(-1.2517396) q[0];
sx q[0];
rz(-1.804833) q[0];
rz(1.5917336) q[1];
sx q[1];
rz(-1.5668329) q[1];
sx q[1];
rz(-0.11650539) q[1];
rz(-2.0143853) q[2];
sx q[2];
rz(-2.1090322) q[2];
sx q[2];
rz(2.3354989) q[2];
rz(2.3497818) q[3];
sx q[3];
rz(-1.7196426) q[3];
sx q[3];
rz(-0.28256217) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
