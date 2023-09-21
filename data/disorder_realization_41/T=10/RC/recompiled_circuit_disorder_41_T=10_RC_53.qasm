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
rz(-0.39437374) q[1];
sx q[1];
rz(-1.6819277) q[1];
sx q[1];
rz(0.2149166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7637245) q[0];
sx q[0];
rz(-1.2963429) q[0];
sx q[0];
rz(2.8732357) q[0];
rz(2.6724042) q[2];
sx q[2];
rz(-1.4326296) q[2];
sx q[2];
rz(0.56433041) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.79214) q[1];
sx q[1];
rz(-0.98590241) q[1];
sx q[1];
rz(0.5485512) q[1];
rz(-pi) q[2];
rz(0.54833702) q[3];
sx q[3];
rz(-2.3468446) q[3];
sx q[3];
rz(-2.7693975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.73137838) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(-2.5773876) q[2];
rz(1.7764067) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.0441701) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(0.92798293) q[0];
rz(-1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(-2.2448418) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7433886) q[0];
sx q[0];
rz(-0.87289116) q[0];
sx q[0];
rz(2.8845805) q[0];
rz(-pi) q[1];
rz(-1.6438269) q[2];
sx q[2];
rz(-1.4023997) q[2];
sx q[2];
rz(1.8254691) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5337199) q[1];
sx q[1];
rz(-2.3405511) q[1];
sx q[1];
rz(-0.17287066) q[1];
rz(-pi) q[2];
rz(2.0131301) q[3];
sx q[3];
rz(-0.80207523) q[3];
sx q[3];
rz(2.341552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(2.1014452) q[2];
rz(1.4552207) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(-2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4002157) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(2.1133912) q[0];
rz(1.0785412) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-2.7064586) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66416392) q[0];
sx q[0];
rz(-1.6782883) q[0];
sx q[0];
rz(-2.7880923) q[0];
x q[1];
rz(-1.7180195) q[2];
sx q[2];
rz(-1.8185116) q[2];
sx q[2];
rz(-1.148828) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.5907026) q[1];
sx q[1];
rz(-1.0895551) q[1];
sx q[1];
rz(2.9571556) q[1];
rz(-pi) q[2];
rz(1.2262672) q[3];
sx q[3];
rz(-0.83762729) q[3];
sx q[3];
rz(-1.4602349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6083287) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(-0.30291525) q[2];
rz(1.8164002) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(0.091025092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19642297) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(2.2241425) q[0];
rz(0.67287412) q[1];
sx q[1];
rz(-2.0560975) q[1];
sx q[1];
rz(-2.8767169) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0565856) q[0];
sx q[0];
rz(-1.5406973) q[0];
sx q[0];
rz(2.4823275) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27819602) q[2];
sx q[2];
rz(-1.858466) q[2];
sx q[2];
rz(0.96565914) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.59135624) q[1];
sx q[1];
rz(-1.9758421) q[1];
sx q[1];
rz(-2.3637799) q[1];
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
rz(0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(1.5765566) q[2];
rz(-1.0270843) q[3];
sx q[3];
rz(-0.74147195) q[3];
sx q[3];
rz(-2.0402133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89001369) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(-0.47873163) q[0];
rz(-1.0331253) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(-2.1889401) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.032741) q[0];
sx q[0];
rz(-1.4534338) q[0];
sx q[0];
rz(-0.54876859) q[0];
x q[1];
rz(1.2383934) q[2];
sx q[2];
rz(-1.1754416) q[2];
sx q[2];
rz(2.7796641) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2561803) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(-1.0186362) q[1];
rz(-2.6635366) q[3];
sx q[3];
rz(-1.607967) q[3];
sx q[3];
rz(-0.93070785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7197363) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(-0.53058132) q[2];
rz(1.4060219) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(2.3099242) q[3];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(1.8364871) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(-2.9690202) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0046376) q[0];
sx q[0];
rz(-1.9872268) q[0];
sx q[0];
rz(-3.1266771) q[0];
rz(-1.1144981) q[2];
sx q[2];
rz(-1.3367532) q[2];
sx q[2];
rz(-2.8161088) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.1291618) q[1];
sx q[1];
rz(-2.7134502) q[1];
sx q[1];
rz(1.3151602) q[1];
rz(-pi) q[2];
rz(-2.5112126) q[3];
sx q[3];
rz(-1.3833589) q[3];
sx q[3];
rz(2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3036348) q[2];
sx q[2];
rz(-1.6875608) q[2];
sx q[2];
rz(2.0149569) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.8036028) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(-0.66147584) q[0];
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
rz(2.7998357) q[0];
sx q[0];
rz(-0.17827398) q[0];
sx q[0];
rz(-2.1205649) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7301894) q[2];
sx q[2];
rz(-1.7559933) q[2];
sx q[2];
rz(-1.7298557) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3488256) q[1];
sx q[1];
rz(-1.4139237) q[1];
sx q[1];
rz(1.8030333) q[1];
rz(-2.8444418) q[3];
sx q[3];
rz(-0.20016709) q[3];
sx q[3];
rz(0.9447007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(-2.9471617) q[2];
rz(2.2284609) q[3];
sx q[3];
rz(-1.7539932) q[3];
sx q[3];
rz(2.156179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59654355) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(3.074926) q[0];
rz(-0.32456675) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(-0.98888046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38935223) q[0];
sx q[0];
rz(-1.1712495) q[0];
sx q[0];
rz(-0.35895343) q[0];
x q[1];
rz(1.5547995) q[2];
sx q[2];
rz(-1.664186) q[2];
sx q[2];
rz(2.8430251) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0162184) q[1];
sx q[1];
rz(-1.5040633) q[1];
sx q[1];
rz(1.1557505) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51560651) q[3];
sx q[3];
rz(-1.139384) q[3];
sx q[3];
rz(-3.072217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(2.7339593) q[2];
rz(-2.3729825) q[3];
sx q[3];
rz(-1.3137484) q[3];
sx q[3];
rz(0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3826564) q[0];
sx q[0];
rz(-1.3019245) q[0];
sx q[0];
rz(-0.60920238) q[0];
rz(3.0464879) q[1];
sx q[1];
rz(-1.8895878) q[1];
sx q[1];
rz(-0.87337714) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72520032) q[0];
sx q[0];
rz(-1.7500449) q[0];
sx q[0];
rz(3.0730625) q[0];
rz(-pi) q[1];
x q[1];
rz(0.74553211) q[2];
sx q[2];
rz(-0.19956707) q[2];
sx q[2];
rz(-1.0996639) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.59233353) q[1];
sx q[1];
rz(-0.98456406) q[1];
sx q[1];
rz(1.1997644) q[1];
rz(-pi) q[2];
rz(0.94536762) q[3];
sx q[3];
rz(-0.67694596) q[3];
sx q[3];
rz(2.7834746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1197027) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(2.7673289) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-2.9746829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4436214) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-0.95296729) q[0];
rz(-2.3151746) q[1];
sx q[1];
rz(-2.4024139) q[1];
sx q[1];
rz(-1.3964765) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8499334) q[0];
sx q[0];
rz(-2.8072661) q[0];
sx q[0];
rz(2.8773984) q[0];
rz(-0.65812494) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(0.70891526) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.64618387) q[1];
sx q[1];
rz(-2.7135661) q[1];
sx q[1];
rz(-2.3133548) q[1];
rz(-pi) q[2];
rz(2.8793328) q[3];
sx q[3];
rz(-1.5159303) q[3];
sx q[3];
rz(2.7454387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.46618) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(0.30188489) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-1.3143905) q[2];
sx q[2];
rz(-2.495043) q[2];
sx q[2];
rz(-1.3174353) q[2];
rz(1.6205447) q[3];
sx q[3];
rz(-1.0377025) q[3];
sx q[3];
rz(2.514537) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];