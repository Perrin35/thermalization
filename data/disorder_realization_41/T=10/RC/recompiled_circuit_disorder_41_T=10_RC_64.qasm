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
rz(2.7472189) q[1];
sx q[1];
rz(-1.4596649) q[1];
sx q[1];
rz(-0.2149166) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2673187) q[0];
sx q[0];
rz(-1.8288757) q[0];
sx q[0];
rz(-1.8549071) q[0];
x q[1];
rz(-2.8432437) q[2];
sx q[2];
rz(-0.48765182) q[2];
sx q[2];
rz(1.2717441) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3494527) q[1];
sx q[1];
rz(-0.98590241) q[1];
sx q[1];
rz(2.5930415) q[1];
x q[2];
rz(-1.0825726) q[3];
sx q[3];
rz(-0.91592741) q[3];
sx q[3];
rz(1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(-2.5773876) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0974225) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(2.2136097) q[0];
rz(-1.1652975) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(2.2448418) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78643909) q[0];
sx q[0];
rz(-2.4053898) q[0];
sx q[0];
rz(1.2765221) q[0];
rz(-pi) q[1];
rz(-0.16883822) q[2];
sx q[2];
rz(-1.4988006) q[2];
sx q[2];
rz(-0.26693401) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.36205081) q[1];
sx q[1];
rz(-0.78501399) q[1];
sx q[1];
rz(-1.7464459) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7249343) q[3];
sx q[3];
rz(-2.2778802) q[3];
sx q[3];
rz(1.3980896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.26560489) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-1.0401475) q[2];
rz(-1.4552207) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(0.43513402) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9520156) q[0];
sx q[0];
rz(-0.36882419) q[0];
sx q[0];
rz(0.30216218) q[0];
x q[1];
rz(-2.6159959) q[2];
sx q[2];
rz(-0.28738775) q[2];
sx q[2];
rz(1.692786) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5907026) q[1];
sx q[1];
rz(-1.0895551) q[1];
sx q[1];
rz(-2.9571556) q[1];
rz(-pi) q[2];
rz(-0.76336236) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(-0.34624472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.53326398) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(0.30291525) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(-0.091025092) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19642297) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(-0.91745013) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-0.26487574) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53752758) q[0];
sx q[0];
rz(-0.91188216) q[0];
sx q[0];
rz(1.5327246) q[0];
rz(1.2722837) q[2];
sx q[2];
rz(-1.8372756) q[2];
sx q[2];
rz(-2.4556015) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5502364) q[1];
sx q[1];
rz(-1.9758421) q[1];
sx q[1];
rz(-0.77781271) q[1];
rz(-pi) q[2];
rz(-0.75986741) q[3];
sx q[3];
rz(-0.93899512) q[3];
sx q[3];
rz(-0.85216537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(-2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.251579) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(-1.0331253) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(0.95265257) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6081776) q[0];
sx q[0];
rz(-2.1153643) q[0];
sx q[0];
rz(-1.4334701) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6636159) q[2];
sx q[2];
rz(-0.51082078) q[2];
sx q[2];
rz(-2.7727327) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.18187411) q[1];
sx q[1];
rz(-1.8078783) q[1];
sx q[1];
rz(-1.167776) q[1];
x q[2];
rz(-1.6126552) q[3];
sx q[3];
rz(-2.0484945) q[3];
sx q[3];
rz(0.65934138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7197363) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(0.53058132) q[2];
rz(-1.7355708) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56753165) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(-1.5166327) q[0];
rz(-1.3051055) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(-2.9690202) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1738152) q[0];
sx q[0];
rz(-0.41668188) q[0];
sx q[0];
rz(-1.5370876) q[0];
rz(-pi) q[1];
rz(-1.1144981) q[2];
sx q[2];
rz(-1.8048394) q[2];
sx q[2];
rz(-0.32548387) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.3249358) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(-1.9865958) q[1];
x q[2];
rz(-1.8014088) q[3];
sx q[3];
rz(-2.1884544) q[3];
sx q[3];
rz(-0.85268439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.83795786) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(-2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(0.66147584) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(-0.75659928) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7717168) q[0];
sx q[0];
rz(-1.4780095) q[0];
sx q[0];
rz(-1.4183527) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7301894) q[2];
sx q[2];
rz(-1.7559933) q[2];
sx q[2];
rz(-1.7298557) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7927671) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(-1.8030333) q[1];
rz(2.8444418) q[3];
sx q[3];
rz(-0.20016709) q[3];
sx q[3];
rz(-0.9447007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-0.1903154) q[2];
sx q[2];
rz(0.19443092) q[2];
rz(0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(-0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(3.074926) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.5044183) q[1];
sx q[1];
rz(-0.98888046) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8152155) q[0];
sx q[0];
rz(-1.9003552) q[0];
sx q[0];
rz(-1.9944847) q[0];
rz(-pi) q[1];
rz(-3.0481911) q[2];
sx q[2];
rz(-1.5867234) q[2];
sx q[2];
rz(-1.8708558) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0162184) q[1];
sx q[1];
rz(-1.6375293) q[1];
sx q[1];
rz(1.1557505) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0841247) q[3];
sx q[3];
rz(-1.1063965) q[3];
sx q[3];
rz(-1.7341136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(0.76861012) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(-0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3826564) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(-3.0464879) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4163923) q[0];
sx q[0];
rz(-1.3915477) q[0];
sx q[0];
rz(-3.0730625) q[0];
x q[1];
rz(-2.3960605) q[2];
sx q[2];
rz(-2.9420256) q[2];
sx q[2];
rz(-2.0419288) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1904618) q[1];
sx q[1];
rz(-1.8776263) q[1];
sx q[1];
rz(0.61913403) q[1];
rz(-pi) q[2];
rz(2.14823) q[3];
sx q[3];
rz(-1.9462898) q[3];
sx q[3];
rz(-1.4162228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-0.78803524) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(-0.37426379) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(-0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436214) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(0.95296729) q[0];
rz(-2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.3964765) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29165927) q[0];
sx q[0];
rz(-2.8072661) q[0];
sx q[0];
rz(2.8773984) q[0];
x q[1];
rz(-0.37784414) q[2];
sx q[2];
rz(-2.4477738) q[2];
sx q[2];
rz(-0.56570429) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.64618387) q[1];
sx q[1];
rz(-0.42802654) q[1];
sx q[1];
rz(-0.82823786) q[1];
rz(1.6276007) q[3];
sx q[3];
rz(-1.3089404) q[3];
sx q[3];
rz(-1.9522304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(-0.15979016) q[2];
rz(0.30188489) q[3];
sx q[3];
rz(-2.2146137) q[3];
sx q[3];
rz(-2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044534279) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(-0.13327577) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(2.2014387) q[2];
sx q[2];
rz(-1.7241782) q[2];
sx q[2];
rz(0.45964514) q[2];
rz(2.6079569) q[3];
sx q[3];
rz(-1.5279557) q[3];
sx q[3];
rz(0.91844311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];