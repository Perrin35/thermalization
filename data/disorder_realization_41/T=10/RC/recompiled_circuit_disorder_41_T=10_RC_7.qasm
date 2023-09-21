OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.70513201) q[0];
sx q[0];
rz(-2.5897265) q[0];
sx q[0];
rz(3.119757) q[0];
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(-0.2149166) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41479933) q[0];
sx q[0];
rz(-2.7601295) q[0];
sx q[0];
rz(2.3261855) q[0];
rz(-pi) q[1];
rz(1.7254513) q[2];
sx q[2];
rz(-2.0351595) q[2];
sx q[2];
rz(-2.2048339) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.54675198) q[1];
sx q[1];
rz(-1.1210124) q[1];
sx q[1];
rz(-2.2307598) q[1];
rz(-pi) q[2];
rz(-2.05902) q[3];
sx q[3];
rz(-0.91592741) q[3];
sx q[3];
rz(-1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4102143) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(0.56420502) q[2];
rz(-1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0441701) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(-2.2136097) q[0];
rz(-1.9762951) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(2.2448418) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1363163) q[0];
sx q[0];
rz(-1.3747842) q[0];
sx q[0];
rz(-2.2851903) q[0];
rz(-pi) q[1];
x q[1];
rz(0.40541655) q[2];
sx q[2];
rz(-2.9581796) q[2];
sx q[2];
rz(-2.237052) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.36205081) q[1];
sx q[1];
rz(-2.3565787) q[1];
sx q[1];
rz(1.3951468) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3223022) q[3];
sx q[3];
rz(-1.2580401) q[3];
sx q[3];
rz(-2.6889338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26560489) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(1.0401475) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(-2.1133912) q[0];
rz(-2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-2.7064586) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9520156) q[0];
sx q[0];
rz(-0.36882419) q[0];
sx q[0];
rz(-2.8394305) q[0];
x q[1];
rz(-1.4235731) q[2];
sx q[2];
rz(-1.8185116) q[2];
sx q[2];
rz(1.148828) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.55089009) q[1];
sx q[1];
rz(-2.0520376) q[1];
sx q[1];
rz(2.9571556) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7828091) q[3];
sx q[3];
rz(-2.3453418) q[3];
sx q[3];
rz(-0.96804726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.53326398) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(0.30291525) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(2.2241425) q[0];
rz(-0.67287412) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(0.26487574) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0850071) q[0];
sx q[0];
rz(-1.5406973) q[0];
sx q[0];
rz(2.4823275) q[0];
x q[1];
rz(0.27819602) q[2];
sx q[2];
rz(-1.858466) q[2];
sx q[2];
rz(-0.96565914) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.59135624) q[1];
sx q[1];
rz(-1.1657506) q[1];
sx q[1];
rz(-2.3637799) q[1];
rz(-pi) q[2];
rz(0.75986741) q[3];
sx q[3];
rz(-0.93899512) q[3];
sx q[3];
rz(0.85216537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5650361) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(1.0331253) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(0.95265257) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27272308) q[0];
sx q[0];
rz(-0.55991828) q[0];
sx q[0];
rz(-2.9193004) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41575899) q[2];
sx q[2];
rz(-1.2649049) q[2];
sx q[2];
rz(-1.8005467) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.488727) q[1];
sx q[1];
rz(-1.9619202) q[1];
sx q[1];
rz(2.8847242) q[1];
rz(-pi) q[2];
rz(1.6126552) q[3];
sx q[3];
rz(-1.0930982) q[3];
sx q[3];
rz(-2.4822513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(2.6110113) q[2];
rz(1.4060219) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(-2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(1.5166327) q[0];
rz(-1.8364871) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(2.9690202) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1369551) q[0];
sx q[0];
rz(-1.9872268) q[0];
sx q[0];
rz(0.01491551) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8820011) q[2];
sx q[2];
rz(-2.0137557) q[2];
sx q[2];
rz(1.358658) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1291618) q[1];
sx q[1];
rz(-2.7134502) q[1];
sx q[1];
rz(-1.3151602) q[1];
rz(-pi) q[2];
rz(1.8014088) q[3];
sx q[3];
rz(-2.1884544) q[3];
sx q[3];
rz(0.85268439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(2.0149569) q[2];
rz(-0.78222328) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(-1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(0.66147584) q[0];
rz(0.96039564) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(2.3849934) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.341757) q[0];
sx q[0];
rz(-0.17827398) q[0];
sx q[0];
rz(-1.0210277) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3691749) q[2];
sx q[2];
rz(-1.1668418) q[2];
sx q[2];
rz(2.902365) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3488256) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(1.3385593) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6301304) q[3];
sx q[3];
rz(-1.3795128) q[3];
sx q[3];
rz(0.64185601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0044272) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(-2.9471617) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(0.98541361) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(-3.074926) q[0];
rz(-0.32456675) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(2.1527122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8152155) q[0];
sx q[0];
rz(-1.9003552) q[0];
sx q[0];
rz(-1.1471079) q[0];
x q[1];
rz(-0.093401508) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(-1.8708558) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.58395308) q[1];
sx q[1];
rz(-1.156731) q[1];
sx q[1];
rz(0.072903452) q[1];
rz(-0.51560651) q[3];
sx q[3];
rz(-1.139384) q[3];
sx q[3];
rz(3.072217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6797592) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(2.7339593) q[2];
rz(2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(-0.60920238) q[0];
rz(-3.0464879) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(0.87337714) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83335919) q[0];
sx q[0];
rz(-1.5033659) q[0];
sx q[0];
rz(1.7504577) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3960605) q[2];
sx q[2];
rz(-0.19956707) q[2];
sx q[2];
rz(2.0419288) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59233353) q[1];
sx q[1];
rz(-2.1570286) q[1];
sx q[1];
rz(-1.9418282) q[1];
rz(-0.94536762) q[3];
sx q[3];
rz(-2.4646467) q[3];
sx q[3];
rz(-0.35811801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(-0.68230391) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.7451161) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8499334) q[0];
sx q[0];
rz(-0.33432654) q[0];
sx q[0];
rz(2.8773984) q[0];
rz(-1.2730607) q[2];
sx q[2];
rz(-0.93431384) q[2];
sx q[2];
rz(1.0422848) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.64618387) q[1];
sx q[1];
rz(-0.42802654) q[1];
sx q[1];
rz(-0.82823786) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26225984) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(0.39615397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(-0.30188489) q[3];
sx q[3];
rz(-2.2146137) q[3];
sx q[3];
rz(2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(3.0970584) q[0];
sx q[0];
rz(-0.67561588) q[0];
sx q[0];
rz(-1.5560879) q[0];
rz(3.0083169) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-1.8272022) q[2];
sx q[2];
rz(-0.64654965) q[2];
sx q[2];
rz(1.8241573) q[2];
rz(-0.084074323) q[3];
sx q[3];
rz(-0.53518674) q[3];
sx q[3];
rz(-0.72471602) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
