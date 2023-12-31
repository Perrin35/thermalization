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
rz(0.37786814) q[0];
sx q[0];
rz(-1.8452497) q[0];
sx q[0];
rz(2.8732357) q[0];
x q[1];
rz(2.8432437) q[2];
sx q[2];
rz(-2.6539408) q[2];
sx q[2];
rz(-1.8698486) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3494527) q[1];
sx q[1];
rz(-2.1556902) q[1];
sx q[1];
rz(2.5930415) q[1];
rz(-2.05902) q[3];
sx q[3];
rz(-2.2256652) q[3];
sx q[3];
rz(1.0893351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4102143) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(2.5773876) q[2];
rz(-1.365186) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(-1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0974225) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(0.89675084) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78643909) q[0];
sx q[0];
rz(-2.4053898) q[0];
sx q[0];
rz(-1.2765221) q[0];
rz(0.16883822) q[2];
sx q[2];
rz(-1.642792) q[2];
sx q[2];
rz(-0.26693401) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5337199) q[1];
sx q[1];
rz(-2.3405511) q[1];
sx q[1];
rz(0.17287066) q[1];
rz(2.3223022) q[3];
sx q[3];
rz(-1.8835526) q[3];
sx q[3];
rz(-0.4526588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.26560489) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(2.1014452) q[2];
rz(1.6863719) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(-2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74137694) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(2.1133912) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-0.43513402) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2745278) q[0];
sx q[0];
rz(-1.9221677) q[0];
sx q[0];
rz(-1.4562796) q[0];
rz(1.7180195) q[2];
sx q[2];
rz(-1.3230811) q[2];
sx q[2];
rz(1.9927646) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5907026) q[1];
sx q[1];
rz(-2.0520376) q[1];
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
rz(0.53326398) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(1.3251925) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(-0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19642297) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(-2.2241425) q[0];
rz(2.4687185) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(2.8767169) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6040651) q[0];
sx q[0];
rz(-0.91188216) q[0];
sx q[0];
rz(1.6088681) q[0];
x q[1];
rz(2.3189544) q[2];
sx q[2];
rz(-2.7441141) q[2];
sx q[2];
rz(0.17695225) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.59135624) q[1];
sx q[1];
rz(-1.9758421) q[1];
sx q[1];
rz(-0.77781271) q[1];
rz(-0.78062765) q[3];
sx q[3];
rz(-2.160191) q[3];
sx q[3];
rz(-1.2300223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.778487) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(2.1145084) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89001369) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(0.47873163) q[0];
rz(1.0331253) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-0.95265257) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53341502) q[0];
sx q[0];
rz(-1.0262283) q[0];
sx q[0];
rz(1.7081225) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4779768) q[2];
sx q[2];
rz(-0.51082078) q[2];
sx q[2];
rz(2.7727327) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.488727) q[1];
sx q[1];
rz(-1.9619202) q[1];
sx q[1];
rz(2.8847242) q[1];
x q[2];
rz(-1.5289375) q[3];
sx q[3];
rz(-1.0930982) q[3];
sx q[3];
rz(0.65934138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7197363) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(2.6110113) q[2];
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
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.6249599) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(0.17257246) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1738152) q[0];
sx q[0];
rz(-0.41668188) q[0];
sx q[0];
rz(1.5370876) q[0];
x q[1];
rz(-2.8820011) q[2];
sx q[2];
rz(-2.0137557) q[2];
sx q[2];
rz(1.7829347) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3249358) q[1];
sx q[1];
rz(-1.4656193) q[1];
sx q[1];
rz(-1.9865958) q[1];
rz(-pi) q[2];
rz(-0.63038007) q[3];
sx q[3];
rz(-1.3833589) q[3];
sx q[3];
rz(0.58296766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(-2.0149569) q[2];
rz(0.78222328) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(-2.4801168) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(2.3849934) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9264383) q[0];
sx q[0];
rz(-1.4190136) q[0];
sx q[0];
rz(3.0477235) q[0];
rz(-pi) q[1];
rz(2.7034764) q[2];
sx q[2];
rz(-2.6926059) q[2];
sx q[2];
rz(0.24030906) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7797864) q[1];
sx q[1];
rz(-0.27946073) q[1];
sx q[1];
rz(-2.1729085) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29715085) q[3];
sx q[3];
rz(-2.9414256) q[3];
sx q[3];
rz(-2.196892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0044272) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(-0.19443092) q[2];
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
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-0.98888046) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7522404) q[0];
sx q[0];
rz(-1.9703431) q[0];
sx q[0];
rz(-0.35895343) q[0];
x q[1];
rz(-3.0481911) q[2];
sx q[2];
rz(-1.5548692) q[2];
sx q[2];
rz(-1.2707368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7372072) q[1];
sx q[1];
rz(-0.42006856) q[1];
sx q[1];
rz(-1.406548) q[1];
x q[2];
rz(-2.6259861) q[3];
sx q[3];
rz(-1.139384) q[3];
sx q[3];
rz(0.069375667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6797592) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(-2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(-0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75893629) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(-2.5323903) q[0];
rz(-3.0464879) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4163923) q[0];
sx q[0];
rz(-1.7500449) q[0];
sx q[0];
rz(0.068530131) q[0];
x q[1];
rz(1.4344425) q[2];
sx q[2];
rz(-1.4246203) q[2];
sx q[2];
rz(1.2863976) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1904618) q[1];
sx q[1];
rz(-1.2639664) q[1];
sx q[1];
rz(-2.5224586) q[1];
rz(-pi) q[2];
rz(0.94536762) q[3];
sx q[3];
rz(-2.4646467) q[3];
sx q[3];
rz(0.35811801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(2.4592887) q[2];
rz(0.37426379) q[3];
sx q[3];
rz(-1.572861) q[3];
sx q[3];
rz(0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-2.1886254) q[0];
rz(2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.7451161) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5293225) q[0];
sx q[0];
rz(-1.6565874) q[0];
sx q[0];
rz(0.32353185) q[0];
rz(-pi) q[1];
rz(-0.65812494) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(0.70891526) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5212621) q[1];
sx q[1];
rz(-1.8552823) q[1];
sx q[1];
rz(1.2465338) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8793328) q[3];
sx q[3];
rz(-1.6256623) q[3];
sx q[3];
rz(-0.39615397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.8272022) q[2];
sx q[2];
rz(-0.64654965) q[2];
sx q[2];
rz(1.8241573) q[2];
rz(-3.0575183) q[3];
sx q[3];
rz(-2.6064059) q[3];
sx q[3];
rz(2.4168766) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
