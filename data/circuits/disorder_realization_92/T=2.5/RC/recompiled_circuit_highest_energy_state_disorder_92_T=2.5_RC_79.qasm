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
rz(0.25843698) q[0];
sx q[0];
rz(-0.28764495) q[0];
sx q[0];
rz(-1.1076466) q[0];
rz(1.3232752) q[1];
sx q[1];
rz(3.4634436) q[1];
sx q[1];
rz(8.7965214) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2520066) q[0];
sx q[0];
rz(-1.4528251) q[0];
sx q[0];
rz(-1.3924166) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1789178) q[2];
sx q[2];
rz(-0.88681839) q[2];
sx q[2];
rz(2.1563403) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7767649) q[1];
sx q[1];
rz(-1.7994295) q[1];
sx q[1];
rz(-1.6079812) q[1];
rz(-pi) q[2];
rz(-2.2443549) q[3];
sx q[3];
rz(-2.4938886) q[3];
sx q[3];
rz(-0.11769262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.2884752) q[2];
sx q[2];
rz(-1.4318117) q[2];
sx q[2];
rz(-2.8004069) q[2];
rz(-0.8902542) q[3];
sx q[3];
rz(-1.2672) q[3];
sx q[3];
rz(2.6064579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.5219236) q[0];
sx q[0];
rz(-0.6837908) q[0];
sx q[0];
rz(-2.4271915) q[0];
rz(-1.5996541) q[1];
sx q[1];
rz(-2.6456412) q[1];
sx q[1];
rz(3.0468859) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15078397) q[0];
sx q[0];
rz(-1.2512733) q[0];
sx q[0];
rz(-0.099415778) q[0];
rz(-pi) q[1];
rz(-0.96927283) q[2];
sx q[2];
rz(-1.8516632) q[2];
sx q[2];
rz(-1.2796677) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4058806) q[1];
sx q[1];
rz(-0.42197126) q[1];
sx q[1];
rz(0.60776819) q[1];
rz(-1.4121145) q[3];
sx q[3];
rz(-2.074721) q[3];
sx q[3];
rz(-2.9014587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9793205) q[2];
sx q[2];
rz(-1.2459735) q[2];
sx q[2];
rz(-0.47404131) q[2];
rz(-2.3794409) q[3];
sx q[3];
rz(-0.63665974) q[3];
sx q[3];
rz(2.6704085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2786461) q[0];
sx q[0];
rz(-3.0248248) q[0];
sx q[0];
rz(2.2886724) q[0];
rz(2.9184753) q[1];
sx q[1];
rz(-2.6631963) q[1];
sx q[1];
rz(0.51582897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70948008) q[0];
sx q[0];
rz(-1.6211544) q[0];
sx q[0];
rz(-1.5169282) q[0];
rz(-pi) q[1];
rz(-1.91024) q[2];
sx q[2];
rz(-0.66563779) q[2];
sx q[2];
rz(-2.9817944) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.22166907) q[1];
sx q[1];
rz(-0.95429568) q[1];
sx q[1];
rz(-2.4024582) q[1];
rz(-2.7985057) q[3];
sx q[3];
rz(-0.90215397) q[3];
sx q[3];
rz(-0.17681387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8082661) q[2];
sx q[2];
rz(-0.9767248) q[2];
sx q[2];
rz(0.094956368) q[2];
rz(0.44148463) q[3];
sx q[3];
rz(-1.35651) q[3];
sx q[3];
rz(1.9879742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3917291) q[0];
sx q[0];
rz(-0.57498217) q[0];
sx q[0];
rz(-1.0561426) q[0];
rz(-0.9437584) q[1];
sx q[1];
rz(-0.25139233) q[1];
sx q[1];
rz(1.4518849) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7923342) q[0];
sx q[0];
rz(-0.99674388) q[0];
sx q[0];
rz(-2.3999155) q[0];
rz(-pi) q[1];
rz(-2.679616) q[2];
sx q[2];
rz(-1.2179795) q[2];
sx q[2];
rz(-1.3188254) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31409594) q[1];
sx q[1];
rz(-0.24330595) q[1];
sx q[1];
rz(-0.14863405) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11066505) q[3];
sx q[3];
rz(-1.5632716) q[3];
sx q[3];
rz(-1.8845817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4132495) q[2];
sx q[2];
rz(-0.80060935) q[2];
sx q[2];
rz(2.441414) q[2];
rz(-0.87107301) q[3];
sx q[3];
rz(-2.0756523) q[3];
sx q[3];
rz(-1.7843436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1313318) q[0];
sx q[0];
rz(-0.49322525) q[0];
sx q[0];
rz(2.6569271) q[0];
rz(-0.87801814) q[1];
sx q[1];
rz(-1.1596102) q[1];
sx q[1];
rz(0.39883089) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53610401) q[0];
sx q[0];
rz(-1.4358504) q[0];
sx q[0];
rz(0.55079726) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7473248) q[2];
sx q[2];
rz(-2.3986001) q[2];
sx q[2];
rz(0.34994469) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3276766) q[1];
sx q[1];
rz(-2.3688201) q[1];
sx q[1];
rz(-1.0303968) q[1];
rz(2.1966241) q[3];
sx q[3];
rz(-2.3381013) q[3];
sx q[3];
rz(-3.0253493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.47348076) q[2];
sx q[2];
rz(-1.5405737) q[2];
sx q[2];
rz(-0.87781805) q[2];
rz(-0.38464883) q[3];
sx q[3];
rz(-2.2967702) q[3];
sx q[3];
rz(-1.9627242) q[3];
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
rz(-0.28410742) q[0];
sx q[0];
rz(-2.6326038) q[0];
sx q[0];
rz(2.8711163) q[0];
rz(-2.8349304) q[1];
sx q[1];
rz(-0.51391927) q[1];
sx q[1];
rz(-0.57672966) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78806308) q[0];
sx q[0];
rz(-1.7559785) q[0];
sx q[0];
rz(3.0524521) q[0];
x q[1];
rz(-2.001695) q[2];
sx q[2];
rz(-2.4043519) q[2];
sx q[2];
rz(0.37079045) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9432183) q[1];
sx q[1];
rz(-0.70462185) q[1];
sx q[1];
rz(0.93438973) q[1];
rz(-0.33476013) q[3];
sx q[3];
rz(-1.5840414) q[3];
sx q[3];
rz(0.56559004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7755255) q[2];
sx q[2];
rz(-2.0719353) q[2];
sx q[2];
rz(0.97037399) q[2];
rz(2.7905285) q[3];
sx q[3];
rz(-0.23215663) q[3];
sx q[3];
rz(-2.8620913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.905726) q[0];
sx q[0];
rz(-0.38385639) q[0];
sx q[0];
rz(-0.50091499) q[0];
rz(2.4210335) q[1];
sx q[1];
rz(-2.2961398) q[1];
sx q[1];
rz(-0.93773425) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7739627) q[0];
sx q[0];
rz(-1.367462) q[0];
sx q[0];
rz(-0.042026566) q[0];
rz(-pi) q[1];
rz(0.47368424) q[2];
sx q[2];
rz(-1.2315183) q[2];
sx q[2];
rz(0.35655278) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5825469) q[1];
sx q[1];
rz(-1.9314879) q[1];
sx q[1];
rz(0.8949032) q[1];
rz(3.1295007) q[3];
sx q[3];
rz(-0.68688697) q[3];
sx q[3];
rz(1.6364607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0096036) q[2];
sx q[2];
rz(-2.9961573) q[2];
sx q[2];
rz(-0.84849882) q[2];
rz(-2.2961473) q[3];
sx q[3];
rz(-2.4852821) q[3];
sx q[3];
rz(-1.9158624) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41035143) q[0];
sx q[0];
rz(-2.1827965) q[0];
sx q[0];
rz(-2.248726) q[0];
rz(-3.0800173) q[1];
sx q[1];
rz(-0.67191809) q[1];
sx q[1];
rz(-2.9796013) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50349405) q[0];
sx q[0];
rz(-1.5862521) q[0];
sx q[0];
rz(-1.6729805) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4396877) q[2];
sx q[2];
rz(-1.6584975) q[2];
sx q[2];
rz(0.039393124) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29934853) q[1];
sx q[1];
rz(-0.33978251) q[1];
sx q[1];
rz(2.7187125) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6747159) q[3];
sx q[3];
rz(-1.4774408) q[3];
sx q[3];
rz(-2.4568077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.091247678) q[2];
sx q[2];
rz(-0.42751905) q[2];
sx q[2];
rz(-0.72489911) q[2];
rz(2.8643518) q[3];
sx q[3];
rz(-2.1767949) q[3];
sx q[3];
rz(-0.086517081) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8898833) q[0];
sx q[0];
rz(-0.6530264) q[0];
sx q[0];
rz(-3.1157893) q[0];
rz(1.7426527) q[1];
sx q[1];
rz(-1.502123) q[1];
sx q[1];
rz(0.10765156) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9190034) q[0];
sx q[0];
rz(-2.1371142) q[0];
sx q[0];
rz(-0.1583216) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4716267) q[2];
sx q[2];
rz(-1.9666858) q[2];
sx q[2];
rz(-2.5079923) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.15756179) q[1];
sx q[1];
rz(-1.6068629) q[1];
sx q[1];
rz(-0.47337516) q[1];
rz(-pi) q[2];
rz(-1.4571545) q[3];
sx q[3];
rz(-1.5885308) q[3];
sx q[3];
rz(-2.4271698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4699576) q[2];
sx q[2];
rz(-1.2592955) q[2];
sx q[2];
rz(1.0939481) q[2];
rz(-0.57540244) q[3];
sx q[3];
rz(-2.3016774) q[3];
sx q[3];
rz(2.0247139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53330082) q[0];
sx q[0];
rz(-2.5550483) q[0];
sx q[0];
rz(-0.71304148) q[0];
rz(-0.22480045) q[1];
sx q[1];
rz(-0.59200042) q[1];
sx q[1];
rz(2.0483268) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21205344) q[0];
sx q[0];
rz(-1.0997547) q[0];
sx q[0];
rz(0.77917288) q[0];
rz(-pi) q[1];
rz(-2.1987783) q[2];
sx q[2];
rz(-2.8066789) q[2];
sx q[2];
rz(-0.57491344) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7674637) q[1];
sx q[1];
rz(-2.0258365) q[1];
sx q[1];
rz(-1.4468852) q[1];
rz(-3.0244707) q[3];
sx q[3];
rz(-1.1599231) q[3];
sx q[3];
rz(2.9409941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1405868) q[2];
sx q[2];
rz(-2.4457377) q[2];
sx q[2];
rz(-2.801753) q[2];
rz(-0.042424399) q[3];
sx q[3];
rz(-1.8570615) q[3];
sx q[3];
rz(-2.5137918) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18702678) q[0];
sx q[0];
rz(-1.5196336) q[0];
sx q[0];
rz(-1.2484311) q[0];
rz(2.3595702) q[1];
sx q[1];
rz(-0.93688688) q[1];
sx q[1];
rz(-0.72601906) q[1];
rz(1.9044884) q[2];
sx q[2];
rz(-1.2174229) q[2];
sx q[2];
rz(0.65897043) q[2];
rz(1.1449161) q[3];
sx q[3];
rz(-1.799634) q[3];
sx q[3];
rz(-1.4992731) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
