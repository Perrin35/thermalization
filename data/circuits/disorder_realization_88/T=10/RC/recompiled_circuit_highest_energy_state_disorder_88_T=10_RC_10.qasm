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
rz(2.9945381) q[0];
sx q[0];
rz(-1.7400063) q[0];
sx q[0];
rz(2.2214878) q[0];
rz(-0.080634557) q[1];
sx q[1];
rz(3.7205003) q[1];
sx q[1];
rz(10.401934) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.282519) q[0];
sx q[0];
rz(-0.4319829) q[0];
sx q[0];
rz(2.6175314) q[0];
rz(-pi) q[1];
rz(-1.5844272) q[2];
sx q[2];
rz(-1.7030099) q[2];
sx q[2];
rz(-1.4135828) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0339573) q[1];
sx q[1];
rz(-0.9743685) q[1];
sx q[1];
rz(1.7340842) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8771873) q[3];
sx q[3];
rz(-0.9054817) q[3];
sx q[3];
rz(0.09395919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5795634) q[2];
sx q[2];
rz(-1.9271489) q[2];
sx q[2];
rz(-2.5207632) q[2];
rz(0.51554716) q[3];
sx q[3];
rz(-1.4225057) q[3];
sx q[3];
rz(2.2990885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8949378) q[0];
sx q[0];
rz(-1.6064914) q[0];
sx q[0];
rz(2.5625693) q[0];
rz(2.31965) q[1];
sx q[1];
rz(-1.5162946) q[1];
sx q[1];
rz(0.45375219) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5229086) q[0];
sx q[0];
rz(-2.4755619) q[0];
sx q[0];
rz(-2.7722107) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46644174) q[2];
sx q[2];
rz(-1.4701519) q[2];
sx q[2];
rz(2.538344) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7612517) q[1];
sx q[1];
rz(-2.5793992) q[1];
sx q[1];
rz(-0.23558771) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1165127) q[3];
sx q[3];
rz(-1.1835872) q[3];
sx q[3];
rz(2.6301165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4072676) q[2];
sx q[2];
rz(-0.13886034) q[2];
sx q[2];
rz(0.56488758) q[2];
rz(2.7867553) q[3];
sx q[3];
rz(-0.9762888) q[3];
sx q[3];
rz(1.7179276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9722209) q[0];
sx q[0];
rz(-0.2114978) q[0];
sx q[0];
rz(-0.55150223) q[0];
rz(-0.36477271) q[1];
sx q[1];
rz(-0.33939895) q[1];
sx q[1];
rz(0.50484467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4992139) q[0];
sx q[0];
rz(-1.9758245) q[0];
sx q[0];
rz(-1.3971055) q[0];
rz(-pi) q[1];
rz(-0.21130685) q[2];
sx q[2];
rz(-2.720565) q[2];
sx q[2];
rz(-0.94079921) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.4144812) q[1];
sx q[1];
rz(-1.0950146) q[1];
sx q[1];
rz(-0.68051312) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48074333) q[3];
sx q[3];
rz(-1.9078443) q[3];
sx q[3];
rz(-2.2729006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7079033) q[2];
sx q[2];
rz(-1.4132376) q[2];
sx q[2];
rz(0.047253963) q[2];
rz(0.83682483) q[3];
sx q[3];
rz(-0.72332007) q[3];
sx q[3];
rz(1.0251454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.717201) q[0];
sx q[0];
rz(-1.0294788) q[0];
sx q[0];
rz(-1.3264054) q[0];
rz(-1.4855509) q[1];
sx q[1];
rz(-2.0573719) q[1];
sx q[1];
rz(-0.63492376) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6328652) q[0];
sx q[0];
rz(-1.5888402) q[0];
sx q[0];
rz(-3.1403072) q[0];
rz(-pi) q[1];
rz(-0.97135404) q[2];
sx q[2];
rz(-0.58341372) q[2];
sx q[2];
rz(0.82380481) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.078770854) q[1];
sx q[1];
rz(-2.7981542) q[1];
sx q[1];
rz(1.7867286) q[1];
rz(0.62795378) q[3];
sx q[3];
rz(-2.0139004) q[3];
sx q[3];
rz(1.7114044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9012458) q[2];
sx q[2];
rz(-2.2525658) q[2];
sx q[2];
rz(1.7690313) q[2];
rz(1.9339804) q[3];
sx q[3];
rz(-2.4977081) q[3];
sx q[3];
rz(2.8670368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32886252) q[0];
sx q[0];
rz(-0.48506081) q[0];
sx q[0];
rz(0.10511705) q[0];
rz(-0.33991995) q[1];
sx q[1];
rz(-0.9143908) q[1];
sx q[1];
rz(-2.0786659) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2539627) q[0];
sx q[0];
rz(-0.85940532) q[0];
sx q[0];
rz(1.5894058) q[0];
rz(-pi) q[1];
rz(2.8811997) q[2];
sx q[2];
rz(-0.77885926) q[2];
sx q[2];
rz(-3.0634524) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2731367) q[1];
sx q[1];
rz(-1.5850164) q[1];
sx q[1];
rz(-3.0830326) q[1];
rz(-pi) q[2];
rz(1.5664728) q[3];
sx q[3];
rz(-0.49884847) q[3];
sx q[3];
rz(-2.3259142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.062332705) q[2];
sx q[2];
rz(-1.7534813) q[2];
sx q[2];
rz(1.3366535) q[2];
rz(-0.054232728) q[3];
sx q[3];
rz(-1.9510061) q[3];
sx q[3];
rz(-2.5469053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9688251) q[0];
sx q[0];
rz(-0.69197881) q[0];
sx q[0];
rz(-1.2139976) q[0];
rz(-1.0460151) q[1];
sx q[1];
rz(-2.0966625) q[1];
sx q[1];
rz(-0.85375839) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7542242) q[0];
sx q[0];
rz(-1.7048536) q[0];
sx q[0];
rz(-0.06019528) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3860903) q[2];
sx q[2];
rz(-1.1248338) q[2];
sx q[2];
rz(-1.8309636) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56928192) q[1];
sx q[1];
rz(-1.0195059) q[1];
sx q[1];
rz(2.2612919) q[1];
rz(-pi) q[2];
rz(2.7925909) q[3];
sx q[3];
rz(-1.4274538) q[3];
sx q[3];
rz(0.81763148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(3.0115016) q[2];
sx q[2];
rz(-1.6513731) q[2];
sx q[2];
rz(2.3940274) q[2];
rz(-0.93303624) q[3];
sx q[3];
rz(-1.7739762) q[3];
sx q[3];
rz(0.30195495) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1099243) q[0];
sx q[0];
rz(-0.95369354) q[0];
sx q[0];
rz(-0.64144301) q[0];
rz(-0.96915069) q[1];
sx q[1];
rz(-1.0698003) q[1];
sx q[1];
rz(-1.1136805) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3308709) q[0];
sx q[0];
rz(-1.668005) q[0];
sx q[0];
rz(-1.9770798) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93514438) q[2];
sx q[2];
rz(-2.0144267) q[2];
sx q[2];
rz(1.1360687) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1195928) q[1];
sx q[1];
rz(-2.5869859) q[1];
sx q[1];
rz(1.9153992) q[1];
rz(-pi) q[2];
rz(-0.05984743) q[3];
sx q[3];
rz(-2.6204797) q[3];
sx q[3];
rz(1.9684362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.574719) q[2];
sx q[2];
rz(-0.69872624) q[2];
sx q[2];
rz(0.75801545) q[2];
rz(2.8356683) q[3];
sx q[3];
rz(-0.35861349) q[3];
sx q[3];
rz(2.6942286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-1.7253983) q[0];
sx q[0];
rz(-0.60878009) q[0];
sx q[0];
rz(2.388227) q[0];
rz(-2.2196409) q[1];
sx q[1];
rz(-1.4267068) q[1];
sx q[1];
rz(3.0070378) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.050722402) q[0];
sx q[0];
rz(-2.7742552) q[0];
sx q[0];
rz(1.814117) q[0];
rz(-0.33498165) q[2];
sx q[2];
rz(-2.4726395) q[2];
sx q[2];
rz(1.1990591) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0363732) q[1];
sx q[1];
rz(-1.4810307) q[1];
sx q[1];
rz(-1.5898934) q[1];
rz(-pi) q[2];
rz(-2.4681925) q[3];
sx q[3];
rz(-1.2030501) q[3];
sx q[3];
rz(2.0665702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3677463) q[2];
sx q[2];
rz(-1.4464804) q[2];
sx q[2];
rz(-1.9473677) q[2];
rz(-1.1456683) q[3];
sx q[3];
rz(-1.1402036) q[3];
sx q[3];
rz(2.0755419) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1016178) q[0];
sx q[0];
rz(-2.6823253) q[0];
sx q[0];
rz(-2.7591144) q[0];
rz(-2.3756012) q[1];
sx q[1];
rz(-1.4975558) q[1];
sx q[1];
rz(1.8006178) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1857698) q[0];
sx q[0];
rz(-2.9583724) q[0];
sx q[0];
rz(0.51143666) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.093042298) q[2];
sx q[2];
rz(-2.3822255) q[2];
sx q[2];
rz(2.0074646) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.36290088) q[1];
sx q[1];
rz(-2.3039989) q[1];
sx q[1];
rz(-1.2630839) q[1];
rz(-pi) q[2];
rz(-0.082499113) q[3];
sx q[3];
rz(-1.9256221) q[3];
sx q[3];
rz(0.88731836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0535023) q[2];
sx q[2];
rz(-2.221205) q[2];
sx q[2];
rz(-3.0255393) q[2];
rz(-1.2608438) q[3];
sx q[3];
rz(-1.1382269) q[3];
sx q[3];
rz(1.457823) q[3];
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
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5905404) q[0];
sx q[0];
rz(-0.11369471) q[0];
sx q[0];
rz(0.1846479) q[0];
rz(-2.2231936) q[1];
sx q[1];
rz(-1.5469488) q[1];
sx q[1];
rz(1.437423) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4686615) q[0];
sx q[0];
rz(-1.3672921) q[0];
sx q[0];
rz(-1.0887906) q[0];
rz(-3.060861) q[2];
sx q[2];
rz(-0.7237607) q[2];
sx q[2];
rz(-2.2249976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.559222) q[1];
sx q[1];
rz(-0.6383903) q[1];
sx q[1];
rz(-1.2450144) q[1];
x q[2];
rz(-0.40664169) q[3];
sx q[3];
rz(-0.71808091) q[3];
sx q[3];
rz(-2.6585916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.664428) q[2];
sx q[2];
rz(-1.2457341) q[2];
sx q[2];
rz(-0.63875443) q[2];
rz(-0.81513682) q[3];
sx q[3];
rz(-2.6705948) q[3];
sx q[3];
rz(-2.7208929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4169793) q[0];
sx q[0];
rz(-1.2495578) q[0];
sx q[0];
rz(0.15171262) q[0];
rz(-2.3607415) q[1];
sx q[1];
rz(-1.6897222) q[1];
sx q[1];
rz(2.2471468) q[1];
rz(1.2909605) q[2];
sx q[2];
rz(-1.4536413) q[2];
sx q[2];
rz(-2.8870256) q[2];
rz(-0.15624795) q[3];
sx q[3];
rz(-0.24689804) q[3];
sx q[3];
rz(-1.762184) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
