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
rz(-1.7464632) q[0];
sx q[0];
rz(-0.75463086) q[0];
sx q[0];
rz(1.0838497) q[0];
rz(-0.87767449) q[1];
sx q[1];
rz(-1.2134774) q[1];
sx q[1];
rz(2.6649063) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86171976) q[0];
sx q[0];
rz(-1.828816) q[0];
sx q[0];
rz(1.9805542) q[0];
rz(-pi) q[1];
rz(2.347685) q[2];
sx q[2];
rz(-2.53144) q[2];
sx q[2];
rz(-2.8086503) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4922766) q[1];
sx q[1];
rz(-1.5536123) q[1];
sx q[1];
rz(-1.8284998) q[1];
rz(-2.6891119) q[3];
sx q[3];
rz(-2.6685331) q[3];
sx q[3];
rz(2.5291131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.013926355) q[2];
sx q[2];
rz(-1.4907336) q[2];
sx q[2];
rz(2.9909383) q[2];
rz(-1.6023747) q[3];
sx q[3];
rz(-2.7645002) q[3];
sx q[3];
rz(-0.25130513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0852614) q[0];
sx q[0];
rz(-1.3717317) q[0];
sx q[0];
rz(0.37013176) q[0];
rz(2.6689957) q[1];
sx q[1];
rz(-0.31817803) q[1];
sx q[1];
rz(-2.1841689) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9345624) q[0];
sx q[0];
rz(-1.0732722) q[0];
sx q[0];
rz(-0.091188641) q[0];
x q[1];
rz(-2.022012) q[2];
sx q[2];
rz(-1.7282082) q[2];
sx q[2];
rz(1.5296641) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.606876) q[1];
sx q[1];
rz(-0.83321111) q[1];
sx q[1];
rz(1.4556769) q[1];
rz(-pi) q[2];
rz(-2.4647275) q[3];
sx q[3];
rz(-1.3086623) q[3];
sx q[3];
rz(1.8755975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.0051673278) q[2];
sx q[2];
rz(-0.74446669) q[2];
sx q[2];
rz(2.1180507) q[2];
rz(1.7510022) q[3];
sx q[3];
rz(-1.3955045) q[3];
sx q[3];
rz(-1.3172147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35109529) q[0];
sx q[0];
rz(-2.1320765) q[0];
sx q[0];
rz(-0.56280953) q[0];
rz(-1.6273181) q[1];
sx q[1];
rz(-1.7104251) q[1];
sx q[1];
rz(-3.0624342) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.81632) q[0];
sx q[0];
rz(-2.5421341) q[0];
sx q[0];
rz(0.58167235) q[0];
rz(-pi) q[1];
rz(1.6116033) q[2];
sx q[2];
rz(-1.0804324) q[2];
sx q[2];
rz(0.23863579) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.2217949) q[1];
sx q[1];
rz(-0.88558965) q[1];
sx q[1];
rz(-2.5849839) q[1];
rz(-0.7003673) q[3];
sx q[3];
rz(-2.1383965) q[3];
sx q[3];
rz(-1.8710473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2341653) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(-1.3321715) q[2];
rz(-2.7857156) q[3];
sx q[3];
rz(-1.1938813) q[3];
sx q[3];
rz(0.98085421) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6818162) q[0];
sx q[0];
rz(-1.9720607) q[0];
sx q[0];
rz(2.0303149) q[0];
rz(0.92012826) q[1];
sx q[1];
rz(-1.5524813) q[1];
sx q[1];
rz(-2.8880602) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7074575) q[0];
sx q[0];
rz(-0.61096707) q[0];
sx q[0];
rz(-0.57662983) q[0];
rz(-pi) q[1];
rz(-2.5229613) q[2];
sx q[2];
rz(-2.9455373) q[2];
sx q[2];
rz(0.679099) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56139466) q[1];
sx q[1];
rz(-0.78402482) q[1];
sx q[1];
rz(-2.8785588) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0959804) q[3];
sx q[3];
rz(-0.56852698) q[3];
sx q[3];
rz(-0.15318682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8186875) q[2];
sx q[2];
rz(-0.56537586) q[2];
sx q[2];
rz(1.4595855) q[2];
rz(0.5736351) q[3];
sx q[3];
rz(-1.9394268) q[3];
sx q[3];
rz(-0.044376317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4742853) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(1.7942418) q[0];
rz(2.5509293) q[1];
sx q[1];
rz(-1.9213516) q[1];
sx q[1];
rz(1.4264872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9085981) q[0];
sx q[0];
rz(-0.54366606) q[0];
sx q[0];
rz(0.74169104) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9587742) q[2];
sx q[2];
rz(-0.96599865) q[2];
sx q[2];
rz(-1.4825578) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3326753) q[1];
sx q[1];
rz(-1.6170224) q[1];
sx q[1];
rz(2.8874365) q[1];
x q[2];
rz(-1.5065881) q[3];
sx q[3];
rz(-1.9252464) q[3];
sx q[3];
rz(1.3469157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.23919375) q[2];
sx q[2];
rz(-1.0603443) q[2];
sx q[2];
rz(-1.9690465) q[2];
rz(0.086056195) q[3];
sx q[3];
rz(-0.9525758) q[3];
sx q[3];
rz(0.15611592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0657144) q[0];
sx q[0];
rz(-1.5944163) q[0];
sx q[0];
rz(-0.86268798) q[0];
rz(2.8942096) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(-2.6695796) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87403546) q[0];
sx q[0];
rz(-0.78473259) q[0];
sx q[0];
rz(2.6536056) q[0];
rz(-2.0765993) q[2];
sx q[2];
rz(-1.4384631) q[2];
sx q[2];
rz(-1.1810609) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.029376349) q[1];
sx q[1];
rz(-0.35307717) q[1];
sx q[1];
rz(0.53332163) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9293044) q[3];
sx q[3];
rz(-0.95930305) q[3];
sx q[3];
rz(-1.4820198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2315959) q[2];
sx q[2];
rz(-0.96698499) q[2];
sx q[2];
rz(-1.3989353) q[2];
rz(1.6861964) q[3];
sx q[3];
rz(-2.3536436) q[3];
sx q[3];
rz(-2.5808891) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48651925) q[0];
sx q[0];
rz(-1.1460679) q[0];
sx q[0];
rz(-0.51862496) q[0];
rz(2.4229166) q[1];
sx q[1];
rz(-0.51281723) q[1];
sx q[1];
rz(1.3291043) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9921917) q[0];
sx q[0];
rz(-1.2103545) q[0];
sx q[0];
rz(-2.0130231) q[0];
rz(-1.3955922) q[2];
sx q[2];
rz(-0.28131286) q[2];
sx q[2];
rz(-0.11872053) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3490263) q[1];
sx q[1];
rz(-1.6036803) q[1];
sx q[1];
rz(2.3057641) q[1];
rz(-pi) q[2];
rz(1.5409327) q[3];
sx q[3];
rz(-0.58664188) q[3];
sx q[3];
rz(-1.3636953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3507877) q[2];
sx q[2];
rz(-2.4420276) q[2];
sx q[2];
rz(3.1034071) q[2];
rz(-2.3112678) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(3.0939046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6787978) q[0];
sx q[0];
rz(-2.9252453) q[0];
sx q[0];
rz(1.7713254) q[0];
rz(-2.7554152) q[1];
sx q[1];
rz(-1.736707) q[1];
sx q[1];
rz(2.4716299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7350175) q[0];
sx q[0];
rz(-1.1018714) q[0];
sx q[0];
rz(-1.0107299) q[0];
x q[1];
rz(-3.0490446) q[2];
sx q[2];
rz(-2.4227648) q[2];
sx q[2];
rz(-2.3859442) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48310623) q[1];
sx q[1];
rz(-1.4577333) q[1];
sx q[1];
rz(-2.1674564) q[1];
x q[2];
rz(-1.3988204) q[3];
sx q[3];
rz(-0.21504083) q[3];
sx q[3];
rz(0.57143927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2204444) q[2];
sx q[2];
rz(-2.8577652) q[2];
sx q[2];
rz(-2.0849126) q[2];
rz(-1.4165261) q[3];
sx q[3];
rz(-1.8959911) q[3];
sx q[3];
rz(1.2904803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.073762745) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(1.006806) q[0];
rz(-2.6489068) q[1];
sx q[1];
rz(-1.7308116) q[1];
sx q[1];
rz(-2.3115092) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.070979764) q[0];
sx q[0];
rz(-1.5375397) q[0];
sx q[0];
rz(-2.8625367) q[0];
x q[1];
rz(-2.8317881) q[2];
sx q[2];
rz(-2.5199157) q[2];
sx q[2];
rz(2.4182188) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.600297) q[1];
sx q[1];
rz(-1.4182555) q[1];
sx q[1];
rz(0.18220724) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1287635) q[3];
sx q[3];
rz(-1.1278102) q[3];
sx q[3];
rz(-1.6010546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3229708) q[2];
sx q[2];
rz(-1.3016204) q[2];
sx q[2];
rz(-2.0446365) q[2];
rz(1.8638301) q[3];
sx q[3];
rz(-0.68456972) q[3];
sx q[3];
rz(-0.6440312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85035664) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(-0.37503234) q[0];
rz(-3.0229783) q[1];
sx q[1];
rz(-1.7053441) q[1];
sx q[1];
rz(-0.92432252) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18436954) q[0];
sx q[0];
rz(-2.4878089) q[0];
sx q[0];
rz(-2.6517646) q[0];
rz(2.4690041) q[2];
sx q[2];
rz(-1.158357) q[2];
sx q[2];
rz(0.46063603) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0339151) q[1];
sx q[1];
rz(-2.7795791) q[1];
sx q[1];
rz(-3.0605143) q[1];
x q[2];
rz(1.9162991) q[3];
sx q[3];
rz(-1.9959108) q[3];
sx q[3];
rz(-0.88750741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.88860005) q[2];
sx q[2];
rz(-2.1414976) q[2];
sx q[2];
rz(2.2484153) q[2];
rz(-2.0231953) q[3];
sx q[3];
rz(-2.1233386) q[3];
sx q[3];
rz(-2.4848188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3746344) q[0];
sx q[0];
rz(-1.081291) q[0];
sx q[0];
rz(-1.3670856) q[0];
rz(-0.18192667) q[1];
sx q[1];
rz(-2.5839099) q[1];
sx q[1];
rz(-0.90669496) q[1];
rz(1.2557658) q[2];
sx q[2];
rz(-1.6581104) q[2];
sx q[2];
rz(2.8078512) q[2];
rz(2.6245194) q[3];
sx q[3];
rz(-1.2096249) q[3];
sx q[3];
rz(0.67740868) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
