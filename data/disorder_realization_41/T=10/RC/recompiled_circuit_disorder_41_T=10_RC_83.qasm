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
rz(2.9266761) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2673187) q[0];
sx q[0];
rz(-1.312717) q[0];
sx q[0];
rz(1.2866856) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4161413) q[2];
sx q[2];
rz(-1.1064331) q[2];
sx q[2];
rz(-2.2048339) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5948407) q[1];
sx q[1];
rz(-2.0205803) q[1];
sx q[1];
rz(-2.2307598) q[1];
rz(-2.4258852) q[3];
sx q[3];
rz(-1.189609) q[3];
sx q[3];
rz(-2.3472799) q[3];
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
rz(-1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.2692497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0441701) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(-1.1652975) q[1];
sx q[1];
rz(-1.5382643) q[1];
sx q[1];
rz(2.2448418) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3551536) q[0];
sx q[0];
rz(-0.73620287) q[0];
sx q[0];
rz(-1.8650706) q[0];
rz(-2.7361761) q[2];
sx q[2];
rz(-0.18341309) q[2];
sx q[2];
rz(-0.90454067) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5337199) q[1];
sx q[1];
rz(-2.3405511) q[1];
sx q[1];
rz(-2.968722) q[1];
rz(-pi) q[2];
rz(2.7249343) q[3];
sx q[3];
rz(-2.2778802) q[3];
sx q[3];
rz(-1.7435031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(-2.1014452) q[2];
rz(1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002157) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(-1.0282015) q[0];
rz(-1.0785412) q[1];
sx q[1];
rz(-0.56285793) q[1];
sx q[1];
rz(-2.7064586) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1895771) q[0];
sx q[0];
rz(-0.36882419) q[0];
sx q[0];
rz(2.8394305) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.52559678) q[2];
sx q[2];
rz(-0.28738775) q[2];
sx q[2];
rz(-1.692786) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1060461) q[1];
sx q[1];
rz(-1.4075081) q[1];
sx q[1];
rz(2.0590904) q[1];
rz(-pi) q[2];
rz(2.7828091) q[3];
sx q[3];
rz(-0.79625087) q[3];
sx q[3];
rz(0.96804726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.53326398) q[2];
sx q[2];
rz(-1.3003131) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(1.3251925) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(-3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
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
rz(2.2241425) q[0];
rz(-2.4687185) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-2.8767169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6040651) q[0];
sx q[0];
rz(-2.2297105) q[0];
sx q[0];
rz(-1.6088681) q[0];
x q[1];
rz(-1.869309) q[2];
sx q[2];
rz(-1.304317) q[2];
sx q[2];
rz(2.4556015) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5323822) q[1];
sx q[1];
rz(-0.86984837) q[1];
sx q[1];
rz(-1.0290531) q[1];
rz(-pi) q[2];
rz(2.360965) q[3];
sx q[3];
rz(-0.98140162) q[3];
sx q[3];
rz(-1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.778487) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(1.5650361) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89001369) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(-0.47873163) q[0];
rz(2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-2.1889401) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27272308) q[0];
sx q[0];
rz(-2.5816744) q[0];
sx q[0];
rz(0.22229226) q[0];
rz(-1.2383934) q[2];
sx q[2];
rz(-1.9661511) q[2];
sx q[2];
rz(-0.3619286) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2561803) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(-2.1229565) q[1];
x q[2];
rz(0.47805602) q[3];
sx q[3];
rz(-1.607967) q[3];
sx q[3];
rz(2.2108848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7197363) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(0.53058132) q[2];
rz(-1.7355708) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(0.83166844) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(1.5166327) q[0];
rz(-1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(2.9690202) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57219244) q[0];
sx q[0];
rz(-1.5571556) q[0];
sx q[0];
rz(1.9872679) q[0];
x q[1];
rz(2.0270945) q[2];
sx q[2];
rz(-1.3367532) q[2];
sx q[2];
rz(-2.8161088) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29218201) q[1];
sx q[1];
rz(-1.1574355) q[1];
sx q[1];
rz(3.0267015) q[1];
rz(-0.31130143) q[3];
sx q[3];
rz(-0.65400306) q[3];
sx q[3];
rz(-1.2378539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.83795786) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(-2.0149569) q[2];
rz(0.78222328) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.3379898) q[3];
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
rz(-0.28850266) q[0];
sx q[0];
rz(-0.30650109) q[0];
sx q[0];
rz(-2.4801168) q[0];
rz(0.96039564) q[1];
sx q[1];
rz(-1.7405225) q[1];
sx q[1];
rz(0.75659928) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9264383) q[0];
sx q[0];
rz(-1.4190136) q[0];
sx q[0];
rz(-3.0477235) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7034764) q[2];
sx q[2];
rz(-2.6926059) q[2];
sx q[2];
rz(-2.9012836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7797864) q[1];
sx q[1];
rz(-0.27946073) q[1];
sx q[1];
rz(-2.1729085) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19161253) q[3];
sx q[3];
rz(-1.629047) q[3];
sx q[3];
rz(-0.91764698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(2.9471617) q[2];
rz(-0.91313177) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(-2.156179) q[3];
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
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(0.32456675) q[1];
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
rz(2.7635927) q[0];
sx q[0];
rz(-2.6110296) q[0];
sx q[0];
rz(-0.87688045) q[0];
rz(-pi) q[1];
rz(-1.5867932) q[2];
sx q[2];
rz(-1.4774067) q[2];
sx q[2];
rz(0.29856759) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5576396) q[1];
sx q[1];
rz(-1.9848616) q[1];
sx q[1];
rz(-0.072903452) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0841247) q[3];
sx q[3];
rz(-1.1063965) q[3];
sx q[3];
rz(1.7341136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6797592) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(2.7339593) q[2];
rz(-0.76861012) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(0.9238981) q[3];
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
rz(-pi) q[0];
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
rz(3.0464879) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(2.2682155) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0489037) q[0];
sx q[0];
rz(-2.9498219) q[0];
sx q[0];
rz(-1.2094686) q[0];
x q[1];
rz(0.74553211) q[2];
sx q[2];
rz(-0.19956707) q[2];
sx q[2];
rz(-1.0996639) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.020563515) q[1];
sx q[1];
rz(-2.4596446) q[1];
sx q[1];
rz(2.6418583) q[1];
rz(-pi) q[2];
rz(2.14823) q[3];
sx q[3];
rz(-1.1953029) q[3];
sx q[3];
rz(1.4162228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-0.68230391) q[2];
rz(-2.7673289) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(2.9746829) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(2.3151746) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.7451161) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5293225) q[0];
sx q[0];
rz(-1.4850052) q[0];
sx q[0];
rz(0.32353185) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4834677) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(-0.70891526) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9979981) q[1];
sx q[1];
rz(-1.2600139) q[1];
sx q[1];
rz(2.8423611) q[1];
rz(0.20874899) q[3];
sx q[3];
rz(-2.8737846) q[3];
sx q[3];
rz(0.97313125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6754127) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-2.2146137) q[3];
sx q[3];
rz(-2.8543499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-2.9524654) q[2];
sx q[2];
rz(-0.94869877) q[2];
sx q[2];
rz(-1.000065) q[2];
rz(-1.6205447) q[3];
sx q[3];
rz(-2.1038901) q[3];
sx q[3];
rz(-0.62705561) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
