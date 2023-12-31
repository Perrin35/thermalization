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
rz(-1.6819277) q[1];
sx q[1];
rz(0.2149166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2673187) q[0];
sx q[0];
rz(-1.8288757) q[0];
sx q[0];
rz(1.8549071) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6724042) q[2];
sx q[2];
rz(-1.4326296) q[2];
sx q[2];
rz(-0.56433041) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6282181) q[1];
sx q[1];
rz(-2.3623423) q[1];
sx q[1];
rz(-0.90374225) q[1];
x q[2];
rz(-0.7157075) q[3];
sx q[3];
rz(-1.189609) q[3];
sx q[3];
rz(-0.79431278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4102143) q[2];
sx q[2];
rz(-1.4593068) q[2];
sx q[2];
rz(-0.56420502) q[2];
rz(-1.365186) q[3];
sx q[3];
rz(-2.6919638) q[3];
sx q[3];
rz(-1.8723429) q[3];
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
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0974225) q[0];
sx q[0];
rz(-1.928227) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(2.2448418) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0052764) q[0];
sx q[0];
rz(-1.7668084) q[0];
sx q[0];
rz(0.85640237) q[0];
rz(-0.40541655) q[2];
sx q[2];
rz(-0.18341309) q[2];
sx q[2];
rz(-2.237052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7795418) q[1];
sx q[1];
rz(-2.3565787) q[1];
sx q[1];
rz(-1.3951468) q[1];
rz(2.3223022) q[3];
sx q[3];
rz(-1.2580401) q[3];
sx q[3];
rz(-2.6889338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8759878) q[2];
sx q[2];
rz(-0.53498712) q[2];
sx q[2];
rz(2.1014452) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.8486332) q[3];
sx q[3];
rz(0.35475981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74137694) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-0.43513402) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4774287) q[0];
sx q[0];
rz(-1.4633044) q[0];
sx q[0];
rz(2.7880923) q[0];
rz(-pi) q[1];
x q[1];
rz(0.25031309) q[2];
sx q[2];
rz(-1.4280983) q[2];
sx q[2];
rz(2.755969) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.16776925) q[1];
sx q[1];
rz(-2.6288189) q[1];
sx q[1];
rz(1.2330526) q[1];
x q[2];
rz(-2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(-2.7953479) q[3];
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
rz(-1.3251925) q[3];
sx q[3];
rz(-1.15851) q[3];
sx q[3];
rz(0.091025092) q[3];
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
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9451697) q[0];
sx q[0];
rz(-1.7315995) q[0];
sx q[0];
rz(-2.2241425) q[0];
rz(0.67287412) q[1];
sx q[1];
rz(-1.0854951) q[1];
sx q[1];
rz(-0.26487574) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53752758) q[0];
sx q[0];
rz(-2.2297105) q[0];
sx q[0];
rz(1.5327246) q[0];
rz(-pi) q[1];
rz(-1.2722837) q[2];
sx q[2];
rz(-1.304317) q[2];
sx q[2];
rz(0.68599115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.781573) q[1];
sx q[1];
rz(-2.2847166) q[1];
sx q[1];
rz(2.5931231) q[1];
rz(-0.78062765) q[3];
sx q[3];
rz(-0.98140162) q[3];
sx q[3];
rz(-1.9115703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.36310568) q[2];
sx q[2];
rz(-0.48831707) q[2];
sx q[2];
rz(1.5765566) q[2];
rz(1.0270843) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.251579) q[0];
sx q[0];
rz(-0.13680923) q[0];
sx q[0];
rz(-2.662861) q[0];
rz(-2.1084673) q[1];
sx q[1];
rz(-0.9712351) q[1];
sx q[1];
rz(2.1889401) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1088516) q[0];
sx q[0];
rz(-1.6881588) q[0];
sx q[0];
rz(0.54876859) q[0];
rz(-pi) q[1];
rz(-0.41575899) q[2];
sx q[2];
rz(-1.2649049) q[2];
sx q[2];
rz(1.8005467) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9597185) q[1];
sx q[1];
rz(-1.8078783) q[1];
sx q[1];
rz(1.167776) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6635366) q[3];
sx q[3];
rz(-1.5336256) q[3];
sx q[3];
rz(0.93070785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(0.53058132) q[2];
rz(-1.4060219) q[3];
sx q[3];
rz(-1.1281745) q[3];
sx q[3];
rz(2.3099242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.574061) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.6249599) q[0];
rz(1.3051055) q[1];
sx q[1];
rz(-1.3508947) q[1];
sx q[1];
rz(2.9690202) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96777746) q[0];
sx q[0];
rz(-2.7249108) q[0];
sx q[0];
rz(-1.5370876) q[0];
rz(-pi) q[1];
rz(0.2595915) q[2];
sx q[2];
rz(-1.1278369) q[2];
sx q[2];
rz(1.358658) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3249358) q[1];
sx q[1];
rz(-1.6759733) q[1];
sx q[1];
rz(-1.9865958) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31130143) q[3];
sx q[3];
rz(-2.4875896) q[3];
sx q[3];
rz(-1.9037387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83795786) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(1.1266358) q[2];
rz(2.3593694) q[3];
sx q[3];
rz(-1.9061079) q[3];
sx q[3];
rz(1.3379898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
sx q[0];
rz(-2.8350916) q[0];
sx q[0];
rz(2.4801168) q[0];
rz(-2.181197) q[1];
sx q[1];
rz(-1.4010701) q[1];
sx q[1];
rz(2.3849934) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7998357) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(-1.0210277) q[0];
rz(-pi) q[1];
rz(2.7034764) q[2];
sx q[2];
rz(-0.44898673) q[2];
sx q[2];
rz(2.9012836) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3488256) q[1];
sx q[1];
rz(-1.4139237) q[1];
sx q[1];
rz(1.8030333) q[1];
rz(-2.9499801) q[3];
sx q[3];
rz(-1.5125456) q[3];
sx q[3];
rz(-0.91764698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.1371655) q[2];
sx q[2];
rz(-2.9512773) q[2];
sx q[2];
rz(-0.19443092) q[2];
rz(2.2284609) q[3];
sx q[3];
rz(-1.3875995) q[3];
sx q[3];
rz(0.98541361) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59654355) q[0];
sx q[0];
rz(-2.5248435) q[0];
sx q[0];
rz(-0.066666691) q[0];
rz(-2.8170259) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(0.98888046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3263771) q[0];
sx q[0];
rz(-1.9003552) q[0];
sx q[0];
rz(1.9944847) q[0];
x q[1];
rz(-3.0481911) q[2];
sx q[2];
rz(-1.5867234) q[2];
sx q[2];
rz(-1.8708558) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0162184) q[1];
sx q[1];
rz(-1.5040633) q[1];
sx q[1];
rz(-1.1557505) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75108053) q[3];
sx q[3];
rz(-0.6595279) q[3];
sx q[3];
rz(-2.2758323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6797592) q[2];
sx q[2];
rz(-2.2448848) q[2];
sx q[2];
rz(0.40763339) q[2];
rz(2.3729825) q[3];
sx q[3];
rz(-1.8278443) q[3];
sx q[3];
rz(0.9238981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75893629) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(2.5323903) q[0];
rz(3.0464879) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(-2.2682155) q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(1.7071502) q[2];
sx q[2];
rz(-1.4246203) q[2];
sx q[2];
rz(-1.2863976) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1210291) q[1];
sx q[1];
rz(-0.68194807) q[1];
sx q[1];
rz(-2.6418583) q[1];
rz(-pi) q[2];
rz(2.196225) q[3];
sx q[3];
rz(-0.67694596) q[3];
sx q[3];
rz(-2.7834746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0218899) q[2];
sx q[2];
rz(-2.3535574) q[2];
sx q[2];
rz(-2.4592887) q[2];
rz(0.37426379) q[3];
sx q[3];
rz(-1.5687317) q[3];
sx q[3];
rz(-0.16690978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[2];
rz(-pi/2) q[2];
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
rz(-1.7451161) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8499334) q[0];
sx q[0];
rz(-0.33432654) q[0];
sx q[0];
rz(-2.8773984) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4834677) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(2.4326774) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9979981) q[1];
sx q[1];
rz(-1.8815787) q[1];
sx q[1];
rz(2.8423611) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-2.6754127) q[2];
sx q[2];
rz(-0.35623494) q[2];
sx q[2];
rz(0.15979016) q[2];
rz(-2.8397078) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(-0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(-3.0083169) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(1.3143905) q[2];
sx q[2];
rz(-0.64654965) q[2];
sx q[2];
rz(1.8241573) q[2];
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
