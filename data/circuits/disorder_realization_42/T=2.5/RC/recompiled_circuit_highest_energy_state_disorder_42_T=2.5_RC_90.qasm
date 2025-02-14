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
rz(5.2353274) q[1];
sx q[1];
rz(8.6040001) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2164803) q[0];
sx q[0];
rz(-0.47639936) q[0];
sx q[0];
rz(0.16615482) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16430597) q[2];
sx q[2];
rz(-0.92245692) q[2];
sx q[2];
rz(0.68337464) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32626881) q[1];
sx q[1];
rz(-2.2827049) q[1];
sx q[1];
rz(-2.2253042) q[1];
x q[2];
rz(2.8540913) q[3];
sx q[3];
rz(-2.578082) q[3];
sx q[3];
rz(0.35421351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.0277752) q[2];
sx q[2];
rz(-2.5706988) q[2];
sx q[2];
rz(-1.6119733) q[2];
rz(0.87013236) q[3];
sx q[3];
rz(-1.2711997) q[3];
sx q[3];
rz(-0.99942708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.73794) q[0];
sx q[0];
rz(-1.0920748) q[0];
sx q[0];
rz(-1.8988443) q[0];
rz(-2.2620762) q[1];
sx q[1];
rz(-0.95006919) q[1];
sx q[1];
rz(1.3290728) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97016818) q[0];
sx q[0];
rz(-2.1546083) q[0];
sx q[0];
rz(-1.0310571) q[0];
rz(0.49449253) q[2];
sx q[2];
rz(-1.1323511) q[2];
sx q[2];
rz(-3.1218253) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.052292) q[1];
sx q[1];
rz(-2.0366021) q[1];
sx q[1];
rz(-1.0518846) q[1];
rz(0.19402282) q[3];
sx q[3];
rz(-0.7027773) q[3];
sx q[3];
rz(-1.3181835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.30895147) q[2];
sx q[2];
rz(-1.2244886) q[2];
sx q[2];
rz(1.857117) q[2];
rz(-0.66713157) q[3];
sx q[3];
rz(-2.7693373) q[3];
sx q[3];
rz(-0.66481203) q[3];
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
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852378) q[0];
sx q[0];
rz(-0.41362718) q[0];
sx q[0];
rz(1.941823) q[0];
rz(1.86357) q[1];
sx q[1];
rz(-2.8646902) q[1];
sx q[1];
rz(0.75495458) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65855712) q[0];
sx q[0];
rz(-2.6847976) q[0];
sx q[0];
rz(-1.5798588) q[0];
x q[1];
rz(0.71866106) q[2];
sx q[2];
rz(-1.6826944) q[2];
sx q[2];
rz(1.1132357) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1633411) q[1];
sx q[1];
rz(-0.74865197) q[1];
sx q[1];
rz(-1.5194375) q[1];
rz(-1.0126013) q[3];
sx q[3];
rz(-2.1821488) q[3];
sx q[3];
rz(-0.064289055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.42489147) q[2];
sx q[2];
rz(-0.5867914) q[2];
sx q[2];
rz(1.1699886) q[2];
rz(-0.26015002) q[3];
sx q[3];
rz(-1.5258421) q[3];
sx q[3];
rz(0.20497504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(0.78885704) q[0];
sx q[0];
rz(-1.4997046) q[0];
sx q[0];
rz(-0.66057551) q[0];
rz(3.0109516) q[1];
sx q[1];
rz(-2.525593) q[1];
sx q[1];
rz(2.7052243) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35231464) q[0];
sx q[0];
rz(-1.8095865) q[0];
sx q[0];
rz(1.055607) q[0];
rz(-pi) q[1];
rz(-0.58664586) q[2];
sx q[2];
rz(-1.2170976) q[2];
sx q[2];
rz(2.2207137) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.66884385) q[1];
sx q[1];
rz(-1.9276697) q[1];
sx q[1];
rz(-1.0297187) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4722669) q[3];
sx q[3];
rz(-1.2312345) q[3];
sx q[3];
rz(-1.3738541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.074162) q[2];
sx q[2];
rz(-0.51965153) q[2];
sx q[2];
rz(2.0036428) q[2];
rz(0.84749311) q[3];
sx q[3];
rz(-1.3769826) q[3];
sx q[3];
rz(-1.4165037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6013019) q[0];
sx q[0];
rz(-1.2492981) q[0];
sx q[0];
rz(2.9345595) q[0];
rz(-1.3218468) q[1];
sx q[1];
rz(-1.9925947) q[1];
sx q[1];
rz(-2.102899) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10932163) q[0];
sx q[0];
rz(-2.3242852) q[0];
sx q[0];
rz(0.82839806) q[0];
rz(0.28858273) q[2];
sx q[2];
rz(-2.5009857) q[2];
sx q[2];
rz(0.45329061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0952044) q[1];
sx q[1];
rz(-0.63648736) q[1];
sx q[1];
rz(-0.028179006) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8349325) q[3];
sx q[3];
rz(-1.467461) q[3];
sx q[3];
rz(-2.5984026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.98628226) q[2];
sx q[2];
rz(-2.0955413) q[2];
sx q[2];
rz(-2.5347064) q[2];
rz(-0.33489975) q[3];
sx q[3];
rz(-1.5891985) q[3];
sx q[3];
rz(1.5437532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5034921) q[0];
sx q[0];
rz(-1.426921) q[0];
sx q[0];
rz(-2.959429) q[0];
rz(-0.69193524) q[1];
sx q[1];
rz(-1.7057995) q[1];
sx q[1];
rz(1.8686132) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86827134) q[0];
sx q[0];
rz(-1.5813424) q[0];
sx q[0];
rz(1.521973) q[0];
rz(2.5820205) q[2];
sx q[2];
rz(-0.26962863) q[2];
sx q[2];
rz(2.0346682) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.581785) q[1];
sx q[1];
rz(-2.323287) q[1];
sx q[1];
rz(2.7831828) q[1];
x q[2];
rz(-1.6342513) q[3];
sx q[3];
rz(-1.2136019) q[3];
sx q[3];
rz(-2.315629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.89340297) q[2];
sx q[2];
rz(-2.4290163) q[2];
sx q[2];
rz(0.87598652) q[2];
rz(-0.28410965) q[3];
sx q[3];
rz(-2.1967389) q[3];
sx q[3];
rz(0.42377728) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1856325) q[0];
sx q[0];
rz(-0.9469339) q[0];
sx q[0];
rz(-2.60485) q[0];
rz(0.56450379) q[1];
sx q[1];
rz(-0.91465488) q[1];
sx q[1];
rz(1.5980665) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48115981) q[0];
sx q[0];
rz(-2.8114722) q[0];
sx q[0];
rz(-2.8274203) q[0];
x q[1];
rz(0.41588628) q[2];
sx q[2];
rz(-1.759811) q[2];
sx q[2];
rz(3.031446) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1607223) q[1];
sx q[1];
rz(-1.1127143) q[1];
sx q[1];
rz(-1.0700109) q[1];
x q[2];
rz(-2.1847637) q[3];
sx q[3];
rz(-1.162863) q[3];
sx q[3];
rz(1.3077426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2619065) q[2];
sx q[2];
rz(-1.3085111) q[2];
sx q[2];
rz(-1.7730664) q[2];
rz(-0.3174828) q[3];
sx q[3];
rz(-1.8902794) q[3];
sx q[3];
rz(0.69664636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
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
rz(-0.19318652) q[0];
sx q[0];
rz(-2.6115311) q[0];
sx q[0];
rz(-1.0411221) q[0];
rz(-2.5913473) q[1];
sx q[1];
rz(-2.7698275) q[1];
sx q[1];
rz(2.9464029) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0597417) q[0];
sx q[0];
rz(-1.4343111) q[0];
sx q[0];
rz(1.5166212) q[0];
rz(3.0311119) q[2];
sx q[2];
rz(-2.049597) q[2];
sx q[2];
rz(-0.88518054) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4176538) q[1];
sx q[1];
rz(-0.28755763) q[1];
sx q[1];
rz(1.4940667) q[1];
x q[2];
rz(2.3161668) q[3];
sx q[3];
rz(-1.9147827) q[3];
sx q[3];
rz(2.0478562) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13888415) q[2];
sx q[2];
rz(-1.1382853) q[2];
sx q[2];
rz(1.8769334) q[2];
rz(1.5539315) q[3];
sx q[3];
rz(-1.2041661) q[3];
sx q[3];
rz(1.9896265) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3569262) q[0];
sx q[0];
rz(-0.21360989) q[0];
sx q[0];
rz(-2.7333562) q[0];
rz(1.6881855) q[1];
sx q[1];
rz(-1.4421578) q[1];
sx q[1];
rz(2.2834987) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8342339) q[0];
sx q[0];
rz(-1.8862794) q[0];
sx q[0];
rz(-2.2351816) q[0];
rz(-pi) q[1];
rz(-3.0442601) q[2];
sx q[2];
rz(-1.6989577) q[2];
sx q[2];
rz(-2.6405957) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.063805997) q[1];
sx q[1];
rz(-1.9651455) q[1];
sx q[1];
rz(-1.1358244) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7210049) q[3];
sx q[3];
rz(-0.69828639) q[3];
sx q[3];
rz(2.5774389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.66576362) q[2];
sx q[2];
rz(-0.88226157) q[2];
sx q[2];
rz(2.8909454) q[2];
rz(0.88746873) q[3];
sx q[3];
rz(-1.1424501) q[3];
sx q[3];
rz(2.7953776) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86824054) q[0];
sx q[0];
rz(-2.1526985) q[0];
sx q[0];
rz(-0.86531472) q[0];
rz(0.52262753) q[1];
sx q[1];
rz(-2.3326383) q[1];
sx q[1];
rz(1.685198) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0883852) q[0];
sx q[0];
rz(-1.6376978) q[0];
sx q[0];
rz(1.3752723) q[0];
rz(-pi) q[1];
rz(-2.2413048) q[2];
sx q[2];
rz(-1.4354726) q[2];
sx q[2];
rz(-1.6953261) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6604546) q[1];
sx q[1];
rz(-2.0298971) q[1];
sx q[1];
rz(1.9641687) q[1];
x q[2];
rz(2.6271712) q[3];
sx q[3];
rz(-1.8728009) q[3];
sx q[3];
rz(-0.3946886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7437848) q[2];
sx q[2];
rz(-2.2297771) q[2];
sx q[2];
rz(-1.3249116) q[2];
rz(1.8600474) q[3];
sx q[3];
rz(-2.2478588) q[3];
sx q[3];
rz(2.3769784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(0.21541967) q[0];
sx q[0];
rz(-1.3289435) q[0];
sx q[0];
rz(-1.5935612) q[0];
rz(-0.46319766) q[1];
sx q[1];
rz(-1.5402272) q[1];
sx q[1];
rz(2.872749) q[1];
rz(0.40512591) q[2];
sx q[2];
rz(-0.68531698) q[2];
sx q[2];
rz(-1.1746049) q[2];
rz(0.67288053) q[3];
sx q[3];
rz(-1.665848) q[3];
sx q[3];
rz(3.0116826) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
