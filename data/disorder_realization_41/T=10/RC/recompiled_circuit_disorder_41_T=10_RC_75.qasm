OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4364606) q[0];
sx q[0];
rz(-0.55186614) q[0];
sx q[0];
rz(-3.119757) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(4.6012576) q[1];
sx q[1];
rz(9.6396946) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41479933) q[0];
sx q[0];
rz(-2.7601295) q[0];
sx q[0];
rz(2.3261855) q[0];
x q[1];
rz(1.4161413) q[2];
sx q[2];
rz(-1.1064331) q[2];
sx q[2];
rz(-2.2048339) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.51337459) q[1];
sx q[1];
rz(-2.3623423) q[1];
sx q[1];
rz(-0.90374225) q[1];
rz(-pi) q[2];
x q[2];
rz(2.05902) q[3];
sx q[3];
rz(-0.91592741) q[3];
sx q[3];
rz(-2.0522576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4102143) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-2.0441701) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(2.2448418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39820406) q[0];
sx q[0];
rz(-2.2687015) q[0];
sx q[0];
rz(-0.25701216) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7361761) q[2];
sx q[2];
rz(-0.18341309) q[2];
sx q[2];
rz(-0.90454067) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.0838544) q[1];
sx q[1];
rz(-1.6946304) q[1];
sx q[1];
rz(2.3480575) q[1];
rz(-0.41665839) q[3];
sx q[3];
rz(-0.86371242) q[3];
sx q[3];
rz(1.7435031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(2.1014452) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74137694) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(1.0785412) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-2.7064586) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2745278) q[0];
sx q[0];
rz(-1.9221677) q[0];
sx q[0];
rz(1.685313) q[0];
rz(-pi) q[1];
rz(0.52559678) q[2];
sx q[2];
rz(-0.28738775) q[2];
sx q[2];
rz(-1.4488066) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.55089009) q[1];
sx q[1];
rz(-2.0520376) q[1];
sx q[1];
rz(-2.9571556) q[1];
x q[2];
rz(2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(-0.34624472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(0.30291525) q[2];
rz(1.3251925) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(-0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(-2.2241425) q[0];
rz(0.67287412) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-2.8767169) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0565856) q[0];
sx q[0];
rz(-1.6008953) q[0];
sx q[0];
rz(2.4823275) q[0];
x q[1];
rz(-0.27819602) q[2];
sx q[2];
rz(-1.2831266) q[2];
sx q[2];
rz(-0.96565914) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.781573) q[1];
sx q[1];
rz(-0.85687602) q[1];
sx q[1];
rz(-2.5931231) q[1];
x q[2];
rz(-2.360965) q[3];
sx q[3];
rz(-0.98140162) q[3];
sx q[3];
rz(-1.2300223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.778487) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(1.5765566) q[2];
rz(2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(-0.47873163) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(0.95265257) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8688696) q[0];
sx q[0];
rz(-2.5816744) q[0];
sx q[0];
rz(-2.9193004) q[0];
x q[1];
rz(-2.7258337) q[2];
sx q[2];
rz(-1.2649049) q[2];
sx q[2];
rz(1.3410459) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.18187411) q[1];
sx q[1];
rz(-1.3337143) q[1];
sx q[1];
rz(1.167776) q[1];
rz(-2.6635366) q[3];
sx q[3];
rz(-1.607967) q[3];
sx q[3];
rz(-0.93070785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.7197363) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(2.6110113) q[2];
rz(-1.7355708) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(0.83166844) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.574061) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.5166327) q[0];
rz(-1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(-0.17257246) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96777746) q[0];
sx q[0];
rz(-2.7249108) q[0];
sx q[0];
rz(-1.5370876) q[0];
rz(-0.2595915) q[2];
sx q[2];
rz(-2.0137557) q[2];
sx q[2];
rz(1.358658) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1291618) q[1];
sx q[1];
rz(-0.42814246) q[1];
sx q[1];
rz(1.8264324) q[1];
rz(1.3401838) q[3];
sx q[3];
rz(-0.95313822) q[3];
sx q[3];
rz(0.85268439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(-0.78222328) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(0.66147584) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(-0.75659928) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.341757) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(1.0210277) q[0];
x q[1];
rz(-2.7034764) q[2];
sx q[2];
rz(-2.6926059) q[2];
sx q[2];
rz(-0.24030906) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.36180624) q[1];
sx q[1];
rz(-0.27946073) q[1];
sx q[1];
rz(0.96868412) q[1];
rz(-pi) q[2];
rz(2.9499801) q[3];
sx q[3];
rz(-1.5125456) q[3];
sx q[3];
rz(-2.2239457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0044272) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(-0.19443092) q[2];
rz(-2.2284609) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(-0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(-0.98888046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7635927) q[0];
sx q[0];
rz(-2.6110296) q[0];
sx q[0];
rz(2.2647122) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9724389) q[2];
sx q[2];
rz(-3.0468468) q[2];
sx q[2];
rz(-0.46846889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5576396) q[1];
sx q[1];
rz(-1.9848616) q[1];
sx q[1];
rz(-0.072903452) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6259861) q[3];
sx q[3];
rz(-2.0022087) q[3];
sx q[3];
rz(-3.072217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(-0.40763339) q[2];
rz(-2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(-0.095104782) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(-0.87337714) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0489037) q[0];
sx q[0];
rz(-2.9498219) q[0];
sx q[0];
rz(-1.2094686) q[0];
x q[1];
rz(1.4344425) q[2];
sx q[2];
rz(-1.4246203) q[2];
sx q[2];
rz(1.2863976) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1210291) q[1];
sx q[1];
rz(-2.4596446) q[1];
sx q[1];
rz(2.6418583) q[1];
x q[2];
rz(2.14823) q[3];
sx q[3];
rz(-1.1953029) q[3];
sx q[3];
rz(-1.7253699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(2.4592887) q[2];
rz(-2.7673289) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-2.1886254) q[0];
rz(-2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.3964765) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29165927) q[0];
sx q[0];
rz(-0.33432654) q[0];
sx q[0];
rz(-2.8773984) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65812494) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(-2.4326774) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9979981) q[1];
sx q[1];
rz(-1.2600139) q[1];
sx q[1];
rz(2.8423611) q[1];
rz(0.26225984) q[3];
sx q[3];
rz(-1.5159303) q[3];
sx q[3];
rz(-2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(-2.9818025) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0970584) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(3.0083169) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(0.940154) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(0.084074323) q[3];
sx q[3];
rz(-2.6064059) q[3];
sx q[3];
rz(2.4168766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];