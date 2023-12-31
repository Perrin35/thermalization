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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2673187) q[0];
sx q[0];
rz(-1.8288757) q[0];
sx q[0];
rz(1.8549071) q[0];
rz(-2.8432437) q[2];
sx q[2];
rz(-0.48765182) q[2];
sx q[2];
rz(1.2717441) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.54675198) q[1];
sx q[1];
rz(-1.1210124) q[1];
sx q[1];
rz(2.2307598) q[1];
rz(-pi) q[2];
rz(2.4258852) q[3];
sx q[3];
rz(-1.189609) q[3];
sx q[3];
rz(2.3472799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4102143) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(-2.5773876) q[2];
rz(-1.365186) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974225) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-0.89675084) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78643909) q[0];
sx q[0];
rz(-0.73620287) q[0];
sx q[0];
rz(-1.8650706) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9727544) q[2];
sx q[2];
rz(-1.642792) q[2];
sx q[2];
rz(0.26693401) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.0577382) q[1];
sx q[1];
rz(-1.4469622) q[1];
sx q[1];
rz(0.79353516) q[1];
rz(0.41665839) q[3];
sx q[3];
rz(-0.86371242) q[3];
sx q[3];
rz(1.3980896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.26560489) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(1.0401475) q[2];
rz(-1.4552207) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-2.564036) q[0];
sx q[0];
rz(2.1133912) q[0];
rz(1.0785412) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-2.7064586) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9520156) q[0];
sx q[0];
rz(-2.7727685) q[0];
sx q[0];
rz(0.30216218) q[0];
x q[1];
rz(-2.8912796) q[2];
sx q[2];
rz(-1.7134943) q[2];
sx q[2];
rz(0.38562361) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.55089009) q[1];
sx q[1];
rz(-1.0895551) q[1];
sx q[1];
rz(2.9571556) q[1];
rz(-pi) q[2];
rz(-1.9153254) q[3];
sx q[3];
rz(-0.83762729) q[3];
sx q[3];
rz(1.6813577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(2.8386774) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
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
rz(0.26487574) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6661975) q[0];
sx q[0];
rz(-0.65984939) q[0];
sx q[0];
rz(3.092479) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2722837) q[2];
sx q[2];
rz(-1.304317) q[2];
sx q[2];
rz(0.68599115) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5323822) q[1];
sx q[1];
rz(-2.2717443) q[1];
sx q[1];
rz(-1.0290531) q[1];
x q[2];
rz(-2.360965) q[3];
sx q[3];
rz(-2.160191) q[3];
sx q[3];
rz(-1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.778487) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5650361) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
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
rz(-2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(2.1889401) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.032741) q[0];
sx q[0];
rz(-1.6881588) q[0];
sx q[0];
rz(-2.5928241) q[0];
x q[1];
rz(2.7258337) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(-1.8005467) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.88541234) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(1.0186362) q[1];
rz(3.0609344) q[3];
sx q[3];
rz(-2.6622052) q[3];
sx q[3];
rz(-0.56848923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7197363) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(-2.6110113) q[2];
rz(-1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.4968137) q[0];
sx q[0];
rz(-1.6249599) q[0];
rz(-1.3051055) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(2.9690202) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5694002) q[0];
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
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.8494106) q[1];
sx q[1];
rz(-1.1574355) q[1];
sx q[1];
rz(-0.11489111) q[1];
x q[2];
rz(0.63038007) q[3];
sx q[3];
rz(-1.3833589) q[3];
sx q[3];
rz(2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-1.1266358) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(-1.8036028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28850266) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(0.66147584) q[0];
rz(0.96039564) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(-0.75659928) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7717168) q[0];
sx q[0];
rz(-1.6635832) q[0];
sx q[0];
rz(-1.7232399) q[0];
rz(-pi) q[1];
x q[1];
rz(0.41140326) q[2];
sx q[2];
rz(-1.7559933) q[2];
sx q[2];
rz(-1.411737) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7927671) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(1.8030333) q[1];
x q[2];
rz(0.19161253) q[3];
sx q[3];
rz(-1.629047) q[3];
sx q[3];
rz(-2.2239457) q[3];
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
rz(0.91313177) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(0.066666691) q[0];
rz(0.32456675) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(-2.1527122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8152155) q[0];
sx q[0];
rz(-1.9003552) q[0];
sx q[0];
rz(1.1471079) q[0];
rz(2.9724389) q[2];
sx q[2];
rz(-0.094745853) q[2];
sx q[2];
rz(-0.46846889) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5576396) q[1];
sx q[1];
rz(-1.9848616) q[1];
sx q[1];
rz(-0.072903452) q[1];
rz(1.0841247) q[3];
sx q[3];
rz(-2.0351962) q[3];
sx q[3];
rz(-1.407479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6797592) q[2];
sx q[2];
rz(-0.89670783) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(-2.2176946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(-2.3826564) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(0.60920238) q[0];
rz(-3.0464879) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4163923) q[0];
sx q[0];
rz(-1.3915477) q[0];
sx q[0];
rz(-0.068530131) q[0];
rz(-pi) q[1];
rz(-0.74553211) q[2];
sx q[2];
rz(-2.9420256) q[2];
sx q[2];
rz(-1.0996639) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1904618) q[1];
sx q[1];
rz(-1.8776263) q[1];
sx q[1];
rz(-2.5224586) q[1];
rz(-pi) q[2];
rz(2.7018413) q[3];
sx q[3];
rz(-1.0381178) q[3];
sx q[3];
rz(2.7524878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1197027) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-2.1886254) q[0];
rz(-0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(1.7451161) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6122702) q[0];
sx q[0];
rz(-1.4850052) q[0];
sx q[0];
rz(2.8180608) q[0];
x q[1];
rz(1.868532) q[2];
sx q[2];
rz(-0.93431384) q[2];
sx q[2];
rz(1.0422848) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6203306) q[1];
sx q[1];
rz(-1.8552823) q[1];
sx q[1];
rz(1.2465338) q[1];
rz(1.5139919) q[3];
sx q[3];
rz(-1.8326522) q[3];
sx q[3];
rz(-1.9522304) q[3];
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
rz(0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
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
rz(3.0083169) q[1];
sx q[1];
rz(-1.517308) q[1];
sx q[1];
rz(3.0130253) q[1];
rz(-2.2014387) q[2];
sx q[2];
rz(-1.4174145) q[2];
sx q[2];
rz(-2.6819475) q[2];
rz(-2.6079569) q[3];
sx q[3];
rz(-1.613637) q[3];
sx q[3];
rz(-2.2231495) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
