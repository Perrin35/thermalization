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
rz(0.084963381) q[0];
sx q[0];
rz(-2.8391916) q[0];
sx q[0];
rz(3.1095355) q[0];
rz(-0.0078460296) q[1];
sx q[1];
rz(-0.50385952) q[1];
sx q[1];
rz(2.388968) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0948715) q[0];
sx q[0];
rz(-1.3699475) q[0];
sx q[0];
rz(0.3263536) q[0];
rz(-pi) q[1];
x q[1];
rz(2.339509) q[2];
sx q[2];
rz(-0.49979106) q[2];
sx q[2];
rz(-1.0700723) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.17129414) q[1];
sx q[1];
rz(-0.42810218) q[1];
sx q[1];
rz(0.64117214) q[1];
rz(-pi) q[2];
rz(-1.2129477) q[3];
sx q[3];
rz(-1.3697624) q[3];
sx q[3];
rz(1.7261637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.96693119) q[2];
sx q[2];
rz(-1.175468) q[2];
sx q[2];
rz(2.3215129) q[2];
rz(-0.44102937) q[3];
sx q[3];
rz(-1.1038154) q[3];
sx q[3];
rz(-0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012861982) q[0];
sx q[0];
rz(-0.91330376) q[0];
sx q[0];
rz(-0.41020694) q[0];
rz(-2.7606616) q[1];
sx q[1];
rz(-1.610264) q[1];
sx q[1];
rz(-1.7832696) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5986431) q[0];
sx q[0];
rz(-0.99587593) q[0];
sx q[0];
rz(0.47358124) q[0];
x q[1];
rz(-1.0274506) q[2];
sx q[2];
rz(-1.1305692) q[2];
sx q[2];
rz(-1.5981975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1413026) q[1];
sx q[1];
rz(-1.4030255) q[1];
sx q[1];
rz(3.0649158) q[1];
rz(-1.3639327) q[3];
sx q[3];
rz(-0.65012041) q[3];
sx q[3];
rz(0.93322414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.35839781) q[2];
sx q[2];
rz(-2.7326549) q[2];
sx q[2];
rz(2.6386293) q[2];
rz(1.2991692) q[3];
sx q[3];
rz(-0.75853577) q[3];
sx q[3];
rz(0.23183307) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7314887) q[0];
sx q[0];
rz(-0.51173156) q[0];
sx q[0];
rz(2.1657535) q[0];
rz(0.93636912) q[1];
sx q[1];
rz(-1.5543289) q[1];
sx q[1];
rz(0.13793129) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5463211) q[0];
sx q[0];
rz(-2.0450651) q[0];
sx q[0];
rz(2.7300937) q[0];
rz(-pi) q[1];
rz(-2.8713869) q[2];
sx q[2];
rz(-2.7305354) q[2];
sx q[2];
rz(-0.11693987) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8616069) q[1];
sx q[1];
rz(-2.2642038) q[1];
sx q[1];
rz(0.1711636) q[1];
rz(-2.6757487) q[3];
sx q[3];
rz(-1.8154105) q[3];
sx q[3];
rz(-2.1810093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5522573) q[2];
sx q[2];
rz(-1.5312803) q[2];
sx q[2];
rz(1.662558) q[2];
rz(-2.3432483) q[3];
sx q[3];
rz(-2.2975497) q[3];
sx q[3];
rz(-3.121283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.289157) q[0];
sx q[0];
rz(-3.030179) q[0];
sx q[0];
rz(1.0668466) q[0];
rz(1.5486708) q[1];
sx q[1];
rz(-1.8916062) q[1];
sx q[1];
rz(2.7222395) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8035078) q[0];
sx q[0];
rz(-0.28059059) q[0];
sx q[0];
rz(2.1795033) q[0];
rz(-pi) q[1];
rz(-2.9871171) q[2];
sx q[2];
rz(-0.53336582) q[2];
sx q[2];
rz(1.5325002) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2372861) q[1];
sx q[1];
rz(-1.5066054) q[1];
sx q[1];
rz(-1.9285081) q[1];
rz(-pi) q[2];
rz(2.7953496) q[3];
sx q[3];
rz(-2.0308563) q[3];
sx q[3];
rz(1.8025043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7833917) q[2];
sx q[2];
rz(-2.6322067) q[2];
sx q[2];
rz(1.5436714) q[2];
rz(1.2265497) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(-1.8515057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.876038) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(-0.75620404) q[0];
rz(1.8887695) q[1];
sx q[1];
rz(-2.2105261) q[1];
sx q[1];
rz(-1.7255712) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6373581) q[0];
sx q[0];
rz(-1.7880482) q[0];
sx q[0];
rz(2.5385607) q[0];
x q[1];
rz(-1.4449233) q[2];
sx q[2];
rz(-1.9121661) q[2];
sx q[2];
rz(1.204513) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6292343) q[1];
sx q[1];
rz(-0.53881379) q[1];
sx q[1];
rz(-0.013756559) q[1];
rz(0.36325034) q[3];
sx q[3];
rz(-1.3923858) q[3];
sx q[3];
rz(1.6476064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8004134) q[2];
sx q[2];
rz(-2.4004553) q[2];
sx q[2];
rz(0.18079147) q[2];
rz(-1.6527269) q[3];
sx q[3];
rz(-1.7427665) q[3];
sx q[3];
rz(1.6256049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4329231) q[0];
sx q[0];
rz(-1.0754508) q[0];
sx q[0];
rz(0.80023009) q[0];
rz(-0.284614) q[1];
sx q[1];
rz(-2.0814643) q[1];
sx q[1];
rz(3.0175993) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2469047) q[0];
sx q[0];
rz(-3.0086811) q[0];
sx q[0];
rz(-1.7869496) q[0];
rz(0.059706177) q[2];
sx q[2];
rz(-1.8054031) q[2];
sx q[2];
rz(-1.4690746) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.968144) q[1];
sx q[1];
rz(-2.8255129) q[1];
sx q[1];
rz(-2.8760002) q[1];
rz(-0.37185566) q[3];
sx q[3];
rz(-1.0069435) q[3];
sx q[3];
rz(-1.2363889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.93512145) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(0.07621152) q[2];
rz(-1.291409) q[3];
sx q[3];
rz(-1.094787) q[3];
sx q[3];
rz(0.61856234) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9718219) q[0];
sx q[0];
rz(-1.5795647) q[0];
sx q[0];
rz(-2.3845657) q[0];
rz(-0.79611671) q[1];
sx q[1];
rz(-1.4069822) q[1];
sx q[1];
rz(-0.98181358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3552723) q[0];
sx q[0];
rz(-2.2844446) q[0];
sx q[0];
rz(-0.37843224) q[0];
rz(-pi) q[1];
rz(2.1481291) q[2];
sx q[2];
rz(-0.12778035) q[2];
sx q[2];
rz(-1.3410695) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0701778) q[1];
sx q[1];
rz(-1.9767673) q[1];
sx q[1];
rz(1.0557879) q[1];
rz(-1.8514093) q[3];
sx q[3];
rz(-1.4320489) q[3];
sx q[3];
rz(0.25158238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7299399) q[2];
sx q[2];
rz(-2.3306658) q[2];
sx q[2];
rz(0.5438424) q[2];
rz(3.1206711) q[3];
sx q[3];
rz(-1.7048416) q[3];
sx q[3];
rz(2.3232443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2127317) q[0];
sx q[0];
rz(-2.2305771) q[0];
sx q[0];
rz(-0.62335706) q[0];
rz(-2.8202672) q[1];
sx q[1];
rz(-2.5265381) q[1];
sx q[1];
rz(-2.8772433) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661115) q[0];
sx q[0];
rz(-0.74680416) q[0];
sx q[0];
rz(1.8854499) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2351722) q[2];
sx q[2];
rz(-1.5315637) q[2];
sx q[2];
rz(-2.7754663) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5874065) q[1];
sx q[1];
rz(-2.3465112) q[1];
sx q[1];
rz(1.1752179) q[1];
rz(-pi) q[2];
rz(-3.0086413) q[3];
sx q[3];
rz(-1.8859409) q[3];
sx q[3];
rz(-3.020435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.136772) q[2];
sx q[2];
rz(-2.4592168) q[2];
sx q[2];
rz(2.6668059) q[2];
rz(2.7759806) q[3];
sx q[3];
rz(-0.48706278) q[3];
sx q[3];
rz(-0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8356165) q[0];
sx q[0];
rz(-0.37537471) q[0];
sx q[0];
rz(-2.6557652) q[0];
rz(-2.0756508) q[1];
sx q[1];
rz(-1.424788) q[1];
sx q[1];
rz(-2.8071383) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7793559) q[0];
sx q[0];
rz(-1.4051172) q[0];
sx q[0];
rz(-0.30596531) q[0];
rz(-pi) q[1];
rz(0.29269258) q[2];
sx q[2];
rz(-0.45537696) q[2];
sx q[2];
rz(0.51560452) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.124467) q[1];
sx q[1];
rz(-1.4216058) q[1];
sx q[1];
rz(2.6890432) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0341088) q[3];
sx q[3];
rz(-2.3114738) q[3];
sx q[3];
rz(2.9024359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36755422) q[2];
sx q[2];
rz(-2.2944821) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(1.6537846) q[3];
sx q[3];
rz(-2.8130468) q[3];
sx q[3];
rz(-3.0706792) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3286572) q[0];
sx q[0];
rz(-2.2870977) q[0];
sx q[0];
rz(2.8712414) q[0];
rz(-1.5125754) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(-0.62634748) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8726845) q[0];
sx q[0];
rz(-1.88694) q[0];
sx q[0];
rz(0.85319467) q[0];
rz(-pi) q[1];
rz(-0.43105189) q[2];
sx q[2];
rz(-0.97857514) q[2];
sx q[2];
rz(0.86041245) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4463908) q[1];
sx q[1];
rz(-0.099041136) q[1];
sx q[1];
rz(2.7161069) q[1];
rz(-pi) q[2];
rz(-1.244943) q[3];
sx q[3];
rz(-1.9385466) q[3];
sx q[3];
rz(-0.25853423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.35499972) q[2];
sx q[2];
rz(-2.2248416) q[2];
sx q[2];
rz(1.5724486) q[2];
rz(-0.7507503) q[3];
sx q[3];
rz(-0.34199491) q[3];
sx q[3];
rz(-2.2641838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56562051) q[0];
sx q[0];
rz(-1.2979869) q[0];
sx q[0];
rz(-0.57814231) q[0];
rz(-1.7987953) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(-3.0192791) q[2];
sx q[2];
rz(-1.5707471) q[2];
sx q[2];
rz(-1.6483501) q[2];
rz(2.883705) q[3];
sx q[3];
rz(-1.5189239) q[3];
sx q[3];
rz(2.2699395) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
