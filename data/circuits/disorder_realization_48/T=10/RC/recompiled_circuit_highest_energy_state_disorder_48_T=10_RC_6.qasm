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
rz(1.3951294) q[0];
sx q[0];
rz(3.8962235) q[0];
sx q[0];
rz(8.3409283) q[0];
rz(-0.87767449) q[1];
sx q[1];
rz(-1.2134774) q[1];
sx q[1];
rz(2.6649063) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81945588) q[0];
sx q[0];
rz(-1.9662204) q[0];
sx q[0];
rz(0.28015341) q[0];
x q[1];
rz(2.0332912) q[2];
sx q[2];
rz(-1.1574189) q[2];
sx q[2];
rz(-2.5819636) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0585441) q[1];
sx q[1];
rz(-1.8284608) q[1];
sx q[1];
rz(3.123822) q[1];
x q[2];
rz(0.45248078) q[3];
sx q[3];
rz(-0.47305952) q[3];
sx q[3];
rz(-2.5291131) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.013926355) q[2];
sx q[2];
rz(-1.4907336) q[2];
sx q[2];
rz(-2.9909383) q[2];
rz(-1.539218) q[3];
sx q[3];
rz(-0.37709245) q[3];
sx q[3];
rz(2.8902875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.056331228) q[0];
sx q[0];
rz(-1.3717317) q[0];
sx q[0];
rz(0.37013176) q[0];
rz(2.6689957) q[1];
sx q[1];
rz(-0.31817803) q[1];
sx q[1];
rz(0.95742375) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9345624) q[0];
sx q[0];
rz(-1.0732722) q[0];
sx q[0];
rz(-0.091188641) q[0];
rz(-pi) q[1];
rz(-1.2217058) q[2];
sx q[2];
rz(-2.6654976) q[2];
sx q[2];
rz(-0.27175909) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95847392) q[1];
sx q[1];
rz(-1.4856824) q[1];
sx q[1];
rz(0.74089153) q[1];
x q[2];
rz(-2.7368746) q[3];
sx q[3];
rz(-2.4232695) q[3];
sx q[3];
rz(2.5248034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1364253) q[2];
sx q[2];
rz(-2.397126) q[2];
sx q[2];
rz(-1.0235419) q[2];
rz(-1.7510022) q[3];
sx q[3];
rz(-1.3955045) q[3];
sx q[3];
rz(-1.8243779) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35109529) q[0];
sx q[0];
rz(-1.0095162) q[0];
sx q[0];
rz(-2.5787831) q[0];
rz(-1.5142745) q[1];
sx q[1];
rz(-1.4311675) q[1];
sx q[1];
rz(-3.0624342) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.81632) q[0];
sx q[0];
rz(-2.5421341) q[0];
sx q[0];
rz(2.5599203) q[0];
rz(-3.0653238) q[2];
sx q[2];
rz(-0.4919211) q[2];
sx q[2];
rz(2.9894376) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.2217949) q[1];
sx q[1];
rz(-0.88558965) q[1];
sx q[1];
rz(-0.55660875) q[1];
rz(2.2658562) q[3];
sx q[3];
rz(-2.1453224) q[3];
sx q[3];
rz(-0.72573001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9074273) q[2];
sx q[2];
rz(-1.5831466) q[2];
sx q[2];
rz(1.3321715) q[2];
rz(0.35587707) q[3];
sx q[3];
rz(-1.1938813) q[3];
sx q[3];
rz(-2.1607384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6818162) q[0];
sx q[0];
rz(-1.9720607) q[0];
sx q[0];
rz(-2.0303149) q[0];
rz(-2.2214644) q[1];
sx q[1];
rz(-1.5891113) q[1];
sx q[1];
rz(2.8880602) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64719114) q[0];
sx q[0];
rz(-1.8888942) q[0];
sx q[0];
rz(-2.6107016) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9811834) q[2];
sx q[2];
rz(-1.4575851) q[2];
sx q[2];
rz(2.8594174) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82099709) q[1];
sx q[1];
rz(-1.3861457) q[1];
sx q[1];
rz(-2.3750633) q[1];
x q[2];
rz(0.30999581) q[3];
sx q[3];
rz(-2.0553723) q[3];
sx q[3];
rz(2.692401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.8186875) q[2];
sx q[2];
rz(-0.56537586) q[2];
sx q[2];
rz(-1.4595855) q[2];
rz(2.5679576) q[3];
sx q[3];
rz(-1.2021659) q[3];
sx q[3];
rz(3.0972163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.4742853) q[0];
sx q[0];
rz(-1.0695589) q[0];
sx q[0];
rz(1.3473508) q[0];
rz(-0.59066331) q[1];
sx q[1];
rz(-1.2202411) q[1];
sx q[1];
rz(1.7151054) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9085981) q[0];
sx q[0];
rz(-0.54366606) q[0];
sx q[0];
rz(-2.3999016) q[0];
x q[1];
rz(-1.1828184) q[2];
sx q[2];
rz(-2.175594) q[2];
sx q[2];
rz(1.4825578) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25012452) q[1];
sx q[1];
rz(-1.824675) q[1];
sx q[1];
rz(-1.5230383) q[1];
rz(2.9699202) q[3];
sx q[3];
rz(-0.35997691) q[3];
sx q[3];
rz(-1.5300919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.23919375) q[2];
sx q[2];
rz(-2.0812483) q[2];
sx q[2];
rz(1.1725461) q[2];
rz(3.0555365) q[3];
sx q[3];
rz(-2.1890169) q[3];
sx q[3];
rz(0.15611592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0758783) q[0];
sx q[0];
rz(-1.5944163) q[0];
sx q[0];
rz(-0.86268798) q[0];
rz(-2.8942096) q[1];
sx q[1];
rz(-1.3037953) q[1];
sx q[1];
rz(2.6695796) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87403546) q[0];
sx q[0];
rz(-2.3568601) q[0];
sx q[0];
rz(2.6536056) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3026779) q[2];
sx q[2];
rz(-2.6202218) q[2];
sx q[2];
rz(0.15586317) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.029376349) q[1];
sx q[1];
rz(-0.35307717) q[1];
sx q[1];
rz(0.53332163) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9293044) q[3];
sx q[3];
rz(-0.95930305) q[3];
sx q[3];
rz(-1.6595728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2315959) q[2];
sx q[2];
rz(-2.1746077) q[2];
sx q[2];
rz(1.7426573) q[2];
rz(-1.6861964) q[3];
sx q[3];
rz(-0.78794909) q[3];
sx q[3];
rz(-2.5808891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48651925) q[0];
sx q[0];
rz(-1.9955248) q[0];
sx q[0];
rz(-2.6229677) q[0];
rz(0.71867603) q[1];
sx q[1];
rz(-2.6287754) q[1];
sx q[1];
rz(1.3291043) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14940093) q[0];
sx q[0];
rz(-1.9312381) q[0];
sx q[0];
rz(1.1285696) q[0];
rz(1.3955922) q[2];
sx q[2];
rz(-0.28131286) q[2];
sx q[2];
rz(0.11872053) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.79256637) q[1];
sx q[1];
rz(-1.6036803) q[1];
sx q[1];
rz(-2.3057641) q[1];
rz(-0.98436004) q[3];
sx q[3];
rz(-1.5873261) q[3];
sx q[3];
rz(-2.9096188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.79080498) q[2];
sx q[2];
rz(-2.4420276) q[2];
sx q[2];
rz(0.038185509) q[2];
rz(-2.3112678) q[3];
sx q[3];
rz(-1.3874976) q[3];
sx q[3];
rz(-0.047688095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6787978) q[0];
sx q[0];
rz(-0.21634732) q[0];
sx q[0];
rz(-1.3702673) q[0];
rz(-0.38617745) q[1];
sx q[1];
rz(-1.4048856) q[1];
sx q[1];
rz(2.4716299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7012335) q[0];
sx q[0];
rz(-1.0770849) q[0];
sx q[0];
rz(-0.53892737) q[0];
rz(-1.4901082) q[2];
sx q[2];
rz(-0.85570691) q[2];
sx q[2];
rz(-2.2632368) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.252287) q[1];
sx q[1];
rz(-0.60599594) q[1];
sx q[1];
rz(-1.3713981) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3588328) q[3];
sx q[3];
rz(-1.6073213) q[3];
sx q[3];
rz(1.167447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2204444) q[2];
sx q[2];
rz(-0.28382742) q[2];
sx q[2];
rz(1.05668) q[2];
rz(1.4165261) q[3];
sx q[3];
rz(-1.2456015) q[3];
sx q[3];
rz(-1.8511124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-3.0678299) q[0];
sx q[0];
rz(-2.1011598) q[0];
sx q[0];
rz(-2.1347866) q[0];
rz(0.49268588) q[1];
sx q[1];
rz(-1.410781) q[1];
sx q[1];
rz(2.3115092) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070979764) q[0];
sx q[0];
rz(-1.5375397) q[0];
sx q[0];
rz(2.8625367) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5987723) q[2];
sx q[2];
rz(-1.7493003) q[2];
sx q[2];
rz(1.1019966) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0574903) q[1];
sx q[1];
rz(-1.3907281) q[1];
sx q[1];
rz(-1.4157285) q[1];
rz(-pi) q[2];
rz(-0.012829145) q[3];
sx q[3];
rz(-2.0137824) q[3];
sx q[3];
rz(-1.540538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3229708) q[2];
sx q[2];
rz(-1.8399723) q[2];
sx q[2];
rz(2.0446365) q[2];
rz(1.2777626) q[3];
sx q[3];
rz(-0.68456972) q[3];
sx q[3];
rz(0.6440312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85035664) q[0];
sx q[0];
rz(-2.0391897) q[0];
sx q[0];
rz(0.37503234) q[0];
rz(3.0229783) q[1];
sx q[1];
rz(-1.4362486) q[1];
sx q[1];
rz(2.2172701) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40709201) q[0];
sx q[0];
rz(-1.0043) q[0];
sx q[0];
rz(1.2248216) q[0];
x q[1];
rz(2.0807939) q[2];
sx q[2];
rz(-0.96335232) q[2];
sx q[2];
rz(-0.80112544) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3872747) q[1];
sx q[1];
rz(-1.5994834) q[1];
sx q[1];
rz(-0.36092511) q[1];
rz(2.6931346) q[3];
sx q[3];
rz(-1.8844386) q[3];
sx q[3];
rz(0.83066108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.88860005) q[2];
sx q[2];
rz(-1.000095) q[2];
sx q[2];
rz(-0.89317733) q[2];
rz(1.1183974) q[3];
sx q[3];
rz(-1.018254) q[3];
sx q[3];
rz(2.4848188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-1.7669582) q[0];
sx q[0];
rz(-2.0603016) q[0];
sx q[0];
rz(1.7745071) q[0];
rz(-2.959666) q[1];
sx q[1];
rz(-0.55768273) q[1];
sx q[1];
rz(2.2348977) q[1];
rz(1.8858269) q[2];
sx q[2];
rz(-1.4834822) q[2];
sx q[2];
rz(-0.33374141) q[2];
rz(-2.4891067) q[3];
sx q[3];
rz(-2.5204044) q[3];
sx q[3];
rz(1.6922097) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
