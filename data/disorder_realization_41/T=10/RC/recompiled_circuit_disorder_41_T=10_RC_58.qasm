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
rz(6.8350514) q[0];
sx q[0];
rz(9.4466136) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(4.6012576) q[1];
sx q[1];
rz(9.6396946) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7267933) q[0];
sx q[0];
rz(-2.7601295) q[0];
sx q[0];
rz(2.3261855) q[0];
rz(-pi) q[1];
rz(-0.29834892) q[2];
sx q[2];
rz(-0.48765182) q[2];
sx q[2];
rz(1.8698486) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.54675198) q[1];
sx q[1];
rz(-1.1210124) q[1];
sx q[1];
rz(-2.2307598) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.05902) q[3];
sx q[3];
rz(-0.91592741) q[3];
sx q[3];
rz(2.0522576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4102143) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(1.365186) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.9762951) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-2.2448418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78643909) q[0];
sx q[0];
rz(-0.73620287) q[0];
sx q[0];
rz(1.2765221) q[0];
rz(-pi) q[1];
rz(-1.6438269) q[2];
sx q[2];
rz(-1.739193) q[2];
sx q[2];
rz(1.3161236) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0838544) q[1];
sx q[1];
rz(-1.4469622) q[1];
sx q[1];
rz(-0.79353516) q[1];
rz(-2.7249343) q[3];
sx q[3];
rz(-0.86371242) q[3];
sx q[3];
rz(-1.7435031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(1.6863719) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(-2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(-1.0282015) q[0];
rz(-1.0785412) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(2.7064586) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66416392) q[0];
sx q[0];
rz(-1.4633044) q[0];
sx q[0];
rz(-2.7880923) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8912796) q[2];
sx q[2];
rz(-1.4280983) q[2];
sx q[2];
rz(-2.755969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1060461) q[1];
sx q[1];
rz(-1.4075081) q[1];
sx q[1];
rz(-2.0590904) q[1];
x q[2];
rz(-0.76336236) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9451697) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(2.2241425) q[0];
rz(-0.67287412) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(-2.8767169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6661975) q[0];
sx q[0];
rz(-0.65984939) q[0];
sx q[0];
rz(-0.049113627) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27819602) q[2];
sx q[2];
rz(-1.858466) q[2];
sx q[2];
rz(-2.1759335) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5323822) q[1];
sx q[1];
rz(-0.86984837) q[1];
sx q[1];
rz(-1.0290531) q[1];
rz(-pi) q[2];
rz(-0.75986741) q[3];
sx q[3];
rz(-0.93899512) q[3];
sx q[3];
rz(-0.85216537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(-1.5765566) q[2];
rz(-2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(-2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-0.95265257) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53341502) q[0];
sx q[0];
rz(-1.0262283) q[0];
sx q[0];
rz(-1.4334701) q[0];
rz(-pi) q[1];
rz(-1.2383934) q[2];
sx q[2];
rz(-1.1754416) q[2];
sx q[2];
rz(0.3619286) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.488727) q[1];
sx q[1];
rz(-1.1796724) q[1];
sx q[1];
rz(2.8847242) q[1];
x q[2];
rz(3.0609344) q[3];
sx q[3];
rz(-0.47938743) q[3];
sx q[3];
rz(-2.5731034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4218563) q[2];
sx q[2];
rz(-2.77878) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.5166327) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(2.9690202) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57219244) q[0];
sx q[0];
rz(-1.5844371) q[0];
sx q[0];
rz(1.1543247) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0270945) q[2];
sx q[2];
rz(-1.3367532) q[2];
sx q[2];
rz(2.8161088) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3249358) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(1.1549969) q[1];
rz(-pi) q[2];
rz(0.63038007) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(-2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(0.78222328) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(-2.4801168) q[0];
rz(-0.96039564) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(0.75659928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7717168) q[0];
sx q[0];
rz(-1.4780095) q[0];
sx q[0];
rz(-1.7232399) q[0];
rz(-pi) q[1];
rz(-2.7034764) q[2];
sx q[2];
rz(-0.44898673) q[2];
sx q[2];
rz(0.24030906) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3488256) q[1];
sx q[1];
rz(-1.4139237) q[1];
sx q[1];
rz(1.8030333) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9499801) q[3];
sx q[3];
rz(-1.5125456) q[3];
sx q[3];
rz(2.2239457) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0044272) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(-0.19443092) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(-0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(-2.1527122) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37799997) q[0];
sx q[0];
rz(-2.6110296) q[0];
sx q[0];
rz(-2.2647122) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.16915377) q[2];
sx q[2];
rz(-3.0468468) q[2];
sx q[2];
rz(0.46846889) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7372072) q[1];
sx q[1];
rz(-2.7215241) q[1];
sx q[1];
rz(1.7350446) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75108053) q[3];
sx q[3];
rz(-0.6595279) q[3];
sx q[3];
rz(2.2758323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4618335) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(2.7339593) q[2];
rz(0.76861012) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-2.2682155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3082335) q[0];
sx q[0];
rz(-1.6382268) q[0];
sx q[0];
rz(-1.7504577) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4344425) q[2];
sx q[2];
rz(-1.4246203) q[2];
sx q[2];
rz(-1.2863976) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5492591) q[1];
sx q[1];
rz(-0.98456406) q[1];
sx q[1];
rz(1.1997644) q[1];
rz(-pi) q[2];
x q[2];
rz(2.196225) q[3];
sx q[3];
rz(-0.67694596) q[3];
sx q[3];
rz(-2.7834746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-0.68230391) q[2];
rz(0.37426379) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436214) q[0];
sx q[0];
rz(-1.3598096) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.3964765) q[1];
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
rz(2.4834677) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(-2.4326774) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.1435946) q[1];
sx q[1];
rz(-1.8815787) q[1];
sx q[1];
rz(-0.29923156) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8793328) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(0.39615397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(2.9818025) q[2];
rz(2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.044534279) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(3.0083169) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-2.9524654) q[2];
sx q[2];
rz(-0.94869877) q[2];
sx q[2];
rz(-1.000065) q[2];
rz(3.0575183) q[3];
sx q[3];
rz(-0.53518674) q[3];
sx q[3];
rz(-0.72471602) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
