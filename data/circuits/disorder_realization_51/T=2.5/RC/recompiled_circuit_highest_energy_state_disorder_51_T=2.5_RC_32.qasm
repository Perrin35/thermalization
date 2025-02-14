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
rz(-0.75495523) q[0];
sx q[0];
rz(3.8929953) q[0];
sx q[0];
rz(10.473517) q[0];
rz(-2.7222848) q[1];
sx q[1];
rz(4.5276548) q[1];
sx q[1];
rz(6.164896) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6615579) q[0];
sx q[0];
rz(-1.2999417) q[0];
sx q[0];
rz(-2.1476905) q[0];
rz(-pi) q[1];
rz(-1.5544791) q[2];
sx q[2];
rz(-1.6217578) q[2];
sx q[2];
rz(1.9555443) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.58913) q[1];
sx q[1];
rz(-1.5538061) q[1];
sx q[1];
rz(-3.1119027) q[1];
rz(-pi) q[2];
rz(1.4939601) q[3];
sx q[3];
rz(-1.8733403) q[3];
sx q[3];
rz(1.8267814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.92363858) q[2];
sx q[2];
rz(-0.0035692735) q[2];
sx q[2];
rz(0.16036073) q[2];
rz(1.1637566) q[3];
sx q[3];
rz(-1.0842423) q[3];
sx q[3];
rz(-2.3273996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1639975) q[0];
sx q[0];
rz(-1.4871335) q[0];
sx q[0];
rz(-0.98597041) q[0];
rz(1.5892971) q[1];
sx q[1];
rz(-0.28737107) q[1];
sx q[1];
rz(1.5602962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88402692) q[0];
sx q[0];
rz(-1.0959033) q[0];
sx q[0];
rz(2.7090461) q[0];
rz(3.1104187) q[2];
sx q[2];
rz(-0.59470648) q[2];
sx q[2];
rz(-3.1212357) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4680921) q[1];
sx q[1];
rz(-1.7335658) q[1];
sx q[1];
rz(-2.7767608) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.017707) q[3];
sx q[3];
rz(-2.200211) q[3];
sx q[3];
rz(2.7567425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.61957773) q[2];
sx q[2];
rz(-1.8176983) q[2];
sx q[2];
rz(-2.1066693) q[2];
rz(0.50152913) q[3];
sx q[3];
rz(-0.087567121) q[3];
sx q[3];
rz(-0.97626221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72306776) q[0];
sx q[0];
rz(-1.9880966) q[0];
sx q[0];
rz(-1.9368517) q[0];
rz(1.0459666) q[1];
sx q[1];
rz(-0.090066411) q[1];
sx q[1];
rz(0.14030309) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36049309) q[0];
sx q[0];
rz(-1.3083959) q[0];
sx q[0];
rz(2.187832) q[0];
x q[1];
rz(-2.7195752) q[2];
sx q[2];
rz(-2.2520503) q[2];
sx q[2];
rz(1.8651419) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.122815) q[1];
sx q[1];
rz(-0.63624708) q[1];
sx q[1];
rz(-2.754209) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9123593) q[3];
sx q[3];
rz(-1.7031125) q[3];
sx q[3];
rz(1.4190471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.77421618) q[2];
sx q[2];
rz(-2.1420631) q[2];
sx q[2];
rz(-0.2224758) q[2];
rz(0.065464822) q[3];
sx q[3];
rz(-1.2832578) q[3];
sx q[3];
rz(-0.84622598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(1.9533933) q[0];
sx q[0];
rz(-0.024217483) q[0];
sx q[0];
rz(0.54704332) q[0];
rz(0.25829265) q[1];
sx q[1];
rz(-0.02198418) q[1];
sx q[1];
rz(0.34150728) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9438081) q[0];
sx q[0];
rz(-3.1349234) q[0];
sx q[0];
rz(-1.7775373) q[0];
rz(-pi) q[1];
x q[1];
rz(0.90085588) q[2];
sx q[2];
rz(-2.6253999) q[2];
sx q[2];
rz(1.3951688) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.66879067) q[1];
sx q[1];
rz(-2.3040132) q[1];
sx q[1];
rz(2.0269143) q[1];
rz(-pi) q[2];
rz(2.0773402) q[3];
sx q[3];
rz(-1.9914802) q[3];
sx q[3];
rz(1.7947321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.42939886) q[2];
sx q[2];
rz(-1.8455576) q[2];
sx q[2];
rz(0.94924259) q[2];
rz(2.3927169) q[3];
sx q[3];
rz(-1.8604934) q[3];
sx q[3];
rz(-3.0794028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5267938) q[0];
sx q[0];
rz(-0.034788046) q[0];
sx q[0];
rz(-1.590796) q[0];
rz(1.7855478) q[1];
sx q[1];
rz(-0.0043914774) q[1];
sx q[1];
rz(3.0785676) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3857128) q[0];
sx q[0];
rz(-1.5405419) q[0];
sx q[0];
rz(-3.0755733) q[0];
rz(2.5777528) q[2];
sx q[2];
rz(-0.83602521) q[2];
sx q[2];
rz(-2.6631402) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.99044518) q[1];
sx q[1];
rz(-1.3626117) q[1];
sx q[1];
rz(3.1143536) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.78181569) q[3];
sx q[3];
rz(-1.2226579) q[3];
sx q[3];
rz(1.045351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.65681347) q[2];
sx q[2];
rz(-1.3098837) q[2];
sx q[2];
rz(-2.4978034) q[2];
rz(2.2667609) q[3];
sx q[3];
rz(-2.8372786) q[3];
sx q[3];
rz(2.3410102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0952045) q[0];
sx q[0];
rz(-3.0859741) q[0];
sx q[0];
rz(-0.355542) q[0];
rz(0.19861673) q[1];
sx q[1];
rz(-3.1348517) q[1];
sx q[1];
rz(0.14828646) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9614811) q[0];
sx q[0];
rz(-1.5683055) q[0];
sx q[0];
rz(-2.9585696) q[0];
rz(-pi) q[1];
rz(0.172074) q[2];
sx q[2];
rz(-0.41671696) q[2];
sx q[2];
rz(0.91815776) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7050138) q[1];
sx q[1];
rz(-1.8039898) q[1];
sx q[1];
rz(-1.4326653) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1790696) q[3];
sx q[3];
rz(-1.4279162) q[3];
sx q[3];
rz(2.8382728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6716914) q[2];
sx q[2];
rz(-2.9000059) q[2];
sx q[2];
rz(3.0174603) q[2];
rz(-2.5668674) q[3];
sx q[3];
rz(-2.997213) q[3];
sx q[3];
rz(-2.9706484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2588876) q[0];
sx q[0];
rz(-0.12508617) q[0];
sx q[0];
rz(-2.4001154) q[0];
rz(-2.8575836) q[1];
sx q[1];
rz(-3.1378855) q[1];
sx q[1];
rz(-0.31518087) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7717424) q[0];
sx q[0];
rz(-3.0707504) q[0];
sx q[0];
rz(1.1819345) q[0];
x q[1];
rz(2.0186606) q[2];
sx q[2];
rz(-1.383198) q[2];
sx q[2];
rz(0.59567398) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4969203) q[1];
sx q[1];
rz(-0.61750353) q[1];
sx q[1];
rz(2.8659077) q[1];
rz(1.3471425) q[3];
sx q[3];
rz(-1.5937433) q[3];
sx q[3];
rz(-0.060088559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.86889851) q[2];
sx q[2];
rz(-2.0468476) q[2];
sx q[2];
rz(-0.72186738) q[2];
rz(2.7591211) q[3];
sx q[3];
rz(-1.9964652) q[3];
sx q[3];
rz(-1.0684048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.566399) q[0];
sx q[0];
rz(-3.1167751) q[0];
sx q[0];
rz(-1.5750634) q[0];
rz(0.20340915) q[1];
sx q[1];
rz(-1.8433488) q[1];
sx q[1];
rz(2.4967616) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0087850182) q[0];
sx q[0];
rz(-1.5624579) q[0];
sx q[0];
rz(1.7783283) q[0];
rz(-pi) q[1];
rz(-1.5956085) q[2];
sx q[2];
rz(-0.96867079) q[2];
sx q[2];
rz(-0.40222049) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37721021) q[1];
sx q[1];
rz(-1.6444211) q[1];
sx q[1];
rz(2.1097357) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4025492) q[3];
sx q[3];
rz(-0.62646455) q[3];
sx q[3];
rz(-1.6221969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5975534) q[2];
sx q[2];
rz(-2.7880703) q[2];
sx q[2];
rz(2.7618347) q[2];
rz(-2.0960268) q[3];
sx q[3];
rz(-1.2328204) q[3];
sx q[3];
rz(-1.1988962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7757292) q[0];
sx q[0];
rz(-0.033670306) q[0];
sx q[0];
rz(1.3566383) q[0];
rz(2.7011073) q[1];
sx q[1];
rz(-1.0904652) q[1];
sx q[1];
rz(2.4408565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3627351) q[0];
sx q[0];
rz(-0.88934169) q[0];
sx q[0];
rz(-0.89618857) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9435025) q[2];
sx q[2];
rz(-0.89025324) q[2];
sx q[2];
rz(2.1777505) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.4266738) q[1];
sx q[1];
rz(-2.2853984) q[1];
sx q[1];
rz(1.8780519) q[1];
x q[2];
rz(1.5666991) q[3];
sx q[3];
rz(-2.0034241) q[3];
sx q[3];
rz(-3.0790902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.35357722) q[2];
sx q[2];
rz(-2.7699296) q[2];
sx q[2];
rz(1.8549982) q[2];
rz(-0.50518099) q[3];
sx q[3];
rz(-2.6954539) q[3];
sx q[3];
rz(-1.9193468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5196359) q[0];
sx q[0];
rz(-3.0918047) q[0];
sx q[0];
rz(-1.5420445) q[0];
rz(-0.75625769) q[1];
sx q[1];
rz(-0.007096346) q[1];
sx q[1];
rz(0.33682987) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9040457) q[0];
sx q[0];
rz(-1.8884475) q[0];
sx q[0];
rz(-1.5730412) q[0];
rz(-pi) q[1];
rz(-2.9354503) q[2];
sx q[2];
rz(-1.2954324) q[2];
sx q[2];
rz(-0.58616591) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4647563) q[1];
sx q[1];
rz(-0.14026961) q[1];
sx q[1];
rz(-0.91591473) q[1];
x q[2];
rz(-1.0098861) q[3];
sx q[3];
rz(-1.9282544) q[3];
sx q[3];
rz(-1.7021029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.51356641) q[2];
sx q[2];
rz(-2.1921373) q[2];
sx q[2];
rz(-0.27964082) q[2];
rz(2.5102992) q[3];
sx q[3];
rz(-2.203233) q[3];
sx q[3];
rz(0.62425557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30017988) q[0];
sx q[0];
rz(-1.5501839) q[0];
sx q[0];
rz(-1.3612904) q[0];
rz(-2.367876) q[1];
sx q[1];
rz(-0.63540375) q[1];
sx q[1];
rz(0.2159963) q[1];
rz(-2.28188) q[2];
sx q[2];
rz(-2.0794686) q[2];
sx q[2];
rz(2.5735264) q[2];
rz(-0.37626304) q[3];
sx q[3];
rz(-1.5496764) q[3];
sx q[3];
rz(-1.5839034) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
