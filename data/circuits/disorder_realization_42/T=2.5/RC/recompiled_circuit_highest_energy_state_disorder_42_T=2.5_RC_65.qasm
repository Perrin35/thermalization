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
rz(-1.090156) q[0];
sx q[0];
rz(-0.12026726) q[0];
sx q[0];
rz(-0.076844849) q[0];
rz(2.285217) q[1];
sx q[1];
rz(-1.0478579) q[1];
sx q[1];
rz(-0.82077789) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79362291) q[0];
sx q[0];
rz(-1.6467148) q[0];
sx q[0];
rz(2.6708219) q[0];
x q[1];
rz(1.3581543) q[2];
sx q[2];
rz(-2.4756788) q[2];
sx q[2];
rz(2.7261811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9498074) q[1];
sx q[1];
rz(-0.92647431) q[1];
sx q[1];
rz(-2.5271646) q[1];
x q[2];
rz(1.3935116) q[3];
sx q[3];
rz(-1.0330135) q[3];
sx q[3];
rz(0.017740109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1138175) q[2];
sx q[2];
rz(-2.5706988) q[2];
sx q[2];
rz(-1.5296193) q[2];
rz(0.87013236) q[3];
sx q[3];
rz(-1.2711997) q[3];
sx q[3];
rz(2.1421656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40365264) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(1.8988443) q[0];
rz(0.87951648) q[1];
sx q[1];
rz(-0.95006919) q[1];
sx q[1];
rz(1.3290728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91956988) q[0];
sx q[0];
rz(-2.0139512) q[0];
sx q[0];
rz(0.65624563) q[0];
x q[1];
rz(2.0602559) q[2];
sx q[2];
rz(-1.1266409) q[2];
sx q[2];
rz(1.7760488) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.14790598) q[1];
sx q[1];
rz(-2.4588486) q[1];
sx q[1];
rz(2.3629689) q[1];
rz(-pi) q[2];
rz(0.69345052) q[3];
sx q[3];
rz(-1.6957404) q[3];
sx q[3];
rz(-0.10378621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8326412) q[2];
sx q[2];
rz(-1.917104) q[2];
sx q[2];
rz(1.2844757) q[2];
rz(0.66713157) q[3];
sx q[3];
rz(-2.7693373) q[3];
sx q[3];
rz(-2.4767806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.852378) q[0];
sx q[0];
rz(-2.7279655) q[0];
sx q[0];
rz(-1.1997696) q[0];
rz(1.86357) q[1];
sx q[1];
rz(-0.27690241) q[1];
sx q[1];
rz(-0.75495458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4729378) q[0];
sx q[0];
rz(-1.1140214) q[0];
sx q[0];
rz(-3.1371389) q[0];
rz(-0.16904449) q[2];
sx q[2];
rz(-2.4158106) q[2];
sx q[2];
rz(0.58453416) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9782516) q[1];
sx q[1];
rz(-0.74865197) q[1];
sx q[1];
rz(1.5194375) q[1];
x q[2];
rz(-0.6471031) q[3];
sx q[3];
rz(-2.3386293) q[3];
sx q[3];
rz(0.76319198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7167012) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(-1.971604) q[2];
rz(-0.26015002) q[3];
sx q[3];
rz(-1.5258421) q[3];
sx q[3];
rz(-2.9366176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78885704) q[0];
sx q[0];
rz(-1.4997046) q[0];
sx q[0];
rz(0.66057551) q[0];
rz(-3.0109516) q[1];
sx q[1];
rz(-0.6159997) q[1];
sx q[1];
rz(2.7052243) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0562387) q[0];
sx q[0];
rz(-2.0699916) q[0];
sx q[0];
rz(-0.27277314) q[0];
rz(-pi) q[1];
rz(-0.5882259) q[2];
sx q[2];
rz(-0.67413051) q[2];
sx q[2];
rz(0.16954409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1088686) q[1];
sx q[1];
rz(-2.0744262) q[1];
sx q[1];
rz(-0.41028604) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.27159528) q[3];
sx q[3];
rz(-2.7885572) q[3];
sx q[3];
rz(-1.0853545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.067430647) q[2];
sx q[2];
rz(-2.6219411) q[2];
sx q[2];
rz(1.1379498) q[2];
rz(-2.2940995) q[3];
sx q[3];
rz(-1.7646101) q[3];
sx q[3];
rz(1.4165037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5402907) q[0];
sx q[0];
rz(-1.2492981) q[0];
sx q[0];
rz(0.20703319) q[0];
rz(1.8197458) q[1];
sx q[1];
rz(-1.9925947) q[1];
sx q[1];
rz(1.0386937) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2406851) q[0];
sx q[0];
rz(-1.0552013) q[0];
sx q[0];
rz(-2.2365966) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62049753) q[2];
sx q[2];
rz(-1.7417241) q[2];
sx q[2];
rz(2.2577499) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0952044) q[1];
sx q[1];
rz(-2.5051053) q[1];
sx q[1];
rz(-0.028179006) q[1];
rz(-pi) q[2];
rz(0.30666017) q[3];
sx q[3];
rz(-1.6741317) q[3];
sx q[3];
rz(-0.54319004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1553104) q[2];
sx q[2];
rz(-1.0460514) q[2];
sx q[2];
rz(2.5347064) q[2];
rz(0.33489975) q[3];
sx q[3];
rz(-1.5891985) q[3];
sx q[3];
rz(-1.5437532) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5034921) q[0];
sx q[0];
rz(-1.426921) q[0];
sx q[0];
rz(0.18216369) q[0];
rz(-2.4496574) q[1];
sx q[1];
rz(-1.7057995) q[1];
sx q[1];
rz(-1.8686132) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86827134) q[0];
sx q[0];
rz(-1.5813424) q[0];
sx q[0];
rz(-1.521973) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5820205) q[2];
sx q[2];
rz(-2.871964) q[2];
sx q[2];
rz(-2.0346682) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.5598076) q[1];
sx q[1];
rz(-0.81830561) q[1];
sx q[1];
rz(-0.35840987) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7837378) q[3];
sx q[3];
rz(-1.5113514) q[3];
sx q[3];
rz(2.418973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2481897) q[2];
sx q[2];
rz(-2.4290163) q[2];
sx q[2];
rz(2.2656061) q[2];
rz(-0.28410965) q[3];
sx q[3];
rz(-0.94485372) q[3];
sx q[3];
rz(2.7178154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1856325) q[0];
sx q[0];
rz(-0.9469339) q[0];
sx q[0];
rz(0.53674269) q[0];
rz(2.5770889) q[1];
sx q[1];
rz(-2.2269378) q[1];
sx q[1];
rz(-1.5435262) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48115981) q[0];
sx q[0];
rz(-0.33012046) q[0];
sx q[0];
rz(2.8274203) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6993611) q[2];
sx q[2];
rz(-0.45453192) q[2];
sx q[2];
rz(1.0584128) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7891414) q[1];
sx q[1];
rz(-2.0159715) q[1];
sx q[1];
rz(0.51207249) q[1];
rz(2.2143977) q[3];
sx q[3];
rz(-2.4193086) q[3];
sx q[3];
rz(2.3658906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87968612) q[2];
sx q[2];
rz(-1.8330816) q[2];
sx q[2];
rz(-1.7730664) q[2];
rz(2.8241099) q[3];
sx q[3];
rz(-1.8902794) q[3];
sx q[3];
rz(-2.4449463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19318652) q[0];
sx q[0];
rz(-0.53006154) q[0];
sx q[0];
rz(2.1004706) q[0];
rz(2.5913473) q[1];
sx q[1];
rz(-0.3717652) q[1];
sx q[1];
rz(2.9464029) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7025906) q[0];
sx q[0];
rz(-0.14678188) q[0];
sx q[0];
rz(0.37555666) q[0];
rz(-pi) q[1];
rz(0.11048079) q[2];
sx q[2];
rz(-2.049597) q[2];
sx q[2];
rz(0.88518054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4176538) q[1];
sx q[1];
rz(-2.854035) q[1];
sx q[1];
rz(-1.647526) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.45359277) q[3];
sx q[3];
rz(-2.2633584) q[3];
sx q[3];
rz(-2.3627687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13888415) q[2];
sx q[2];
rz(-2.0033074) q[2];
sx q[2];
rz(-1.2646593) q[2];
rz(1.5539315) q[3];
sx q[3];
rz(-1.2041661) q[3];
sx q[3];
rz(-1.1519661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.78466648) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(2.7333562) q[0];
rz(1.6881855) q[1];
sx q[1];
rz(-1.6994349) q[1];
sx q[1];
rz(-2.2834987) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6397809) q[0];
sx q[0];
rz(-0.94450356) q[0];
sx q[0];
rz(2.7485952) q[0];
rz(-pi) q[1];
rz(1.4420322) q[2];
sx q[2];
rz(-1.4742645) q[2];
sx q[2];
rz(-2.0842722) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0777867) q[1];
sx q[1];
rz(-1.9651455) q[1];
sx q[1];
rz(-1.1358244) q[1];
x q[2];
rz(1.4205877) q[3];
sx q[3];
rz(-0.69828639) q[3];
sx q[3];
rz(0.56415375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.475829) q[2];
sx q[2];
rz(-0.88226157) q[2];
sx q[2];
rz(-0.25064722) q[2];
rz(2.2541239) q[3];
sx q[3];
rz(-1.9991425) q[3];
sx q[3];
rz(2.7953776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2733521) q[0];
sx q[0];
rz(-2.1526985) q[0];
sx q[0];
rz(2.2762779) q[0];
rz(-2.6189651) q[1];
sx q[1];
rz(-0.80895439) q[1];
sx q[1];
rz(1.4563947) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0532074) q[0];
sx q[0];
rz(-1.5038948) q[0];
sx q[0];
rz(-1.7663203) q[0];
rz(-2.9695332) q[2];
sx q[2];
rz(-0.90751782) q[2];
sx q[2];
rz(-0.017939719) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8700614) q[1];
sx q[1];
rz(-1.9215596) q[1];
sx q[1];
rz(-2.6501772) q[1];
x q[2];
rz(1.9144451) q[3];
sx q[3];
rz(-2.0598186) q[3];
sx q[3];
rz(-2.1320313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7437848) q[2];
sx q[2];
rz(-0.91181552) q[2];
sx q[2];
rz(1.3249116) q[2];
rz(1.8600474) q[3];
sx q[3];
rz(-0.89373389) q[3];
sx q[3];
rz(0.76461422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.926173) q[0];
sx q[0];
rz(-1.8126491) q[0];
sx q[0];
rz(1.5480315) q[0];
rz(-2.678395) q[1];
sx q[1];
rz(-1.6013655) q[1];
sx q[1];
rz(-0.2688437) q[1];
rz(0.64434509) q[2];
sx q[2];
rz(-1.3186789) q[2];
sx q[2];
rz(-3.0659715) q[2];
rz(-0.67288053) q[3];
sx q[3];
rz(-1.4757446) q[3];
sx q[3];
rz(-0.12991005) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
