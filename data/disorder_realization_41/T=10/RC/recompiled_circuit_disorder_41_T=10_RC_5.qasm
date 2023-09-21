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
rz(-0.021835672) q[0];
rz(-0.39437374) q[1];
sx q[1];
rz(-1.6819277) q[1];
sx q[1];
rz(0.2149166) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37786814) q[0];
sx q[0];
rz(-1.2963429) q[0];
sx q[0];
rz(2.8732357) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6724042) q[2];
sx q[2];
rz(-1.708963) q[2];
sx q[2];
rz(-2.5772622) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6282181) q[1];
sx q[1];
rz(-0.77925032) q[1];
sx q[1];
rz(2.2378504) q[1];
rz(-pi) q[2];
x q[2];
rz(0.7157075) q[3];
sx q[3];
rz(-1.9519836) q[3];
sx q[3];
rz(-0.79431278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4102143) q[2];
sx q[2];
rz(-1.6822858) q[2];
sx q[2];
rz(2.5773876) q[2];
rz(-1.365186) q[3];
sx q[3];
rz(-0.44962883) q[3];
sx q[3];
rz(1.8723429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0441701) q[0];
sx q[0];
rz(-1.2133657) q[0];
sx q[0];
rz(-0.92798293) q[0];
rz(1.1652975) q[1];
sx q[1];
rz(-1.6033283) q[1];
sx q[1];
rz(2.2448418) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7433886) q[0];
sx q[0];
rz(-0.87289116) q[0];
sx q[0];
rz(-0.25701216) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7361761) q[2];
sx q[2];
rz(-0.18341309) q[2];
sx q[2];
rz(0.90454067) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7795418) q[1];
sx q[1];
rz(-0.78501399) q[1];
sx q[1];
rz(-1.3951468) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41665839) q[3];
sx q[3];
rz(-0.86371242) q[3];
sx q[3];
rz(-1.3980896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8759878) q[2];
sx q[2];
rz(-2.6066055) q[2];
sx q[2];
rz(-1.0401475) q[2];
rz(-1.6863719) q[3];
sx q[3];
rz(-1.2929595) q[3];
sx q[3];
rz(2.7868328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4002157) q[0];
sx q[0];
rz(-0.57755661) q[0];
sx q[0];
rz(1.0282015) q[0];
rz(2.0630515) q[1];
sx q[1];
rz(-2.5787347) q[1];
sx q[1];
rz(-0.43513402) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86706485) q[0];
sx q[0];
rz(-1.2194249) q[0];
sx q[0];
rz(-1.4562796) q[0];
x q[1];
rz(2.8912796) q[2];
sx q[2];
rz(-1.7134943) q[2];
sx q[2];
rz(-0.38562361) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.9738234) q[1];
sx q[1];
rz(-2.6288189) q[1];
sx q[1];
rz(-1.9085401) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3782303) q[3];
sx q[3];
rz(-1.3171139) q[3];
sx q[3];
rz(2.7953479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.53326398) q[2];
sx q[2];
rz(-1.8412795) q[2];
sx q[2];
rz(0.30291525) q[2];
rz(-1.8164002) q[3];
sx q[3];
rz(-1.9830827) q[3];
sx q[3];
rz(-3.0505676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9451697) q[0];
sx q[0];
rz(-1.4099932) q[0];
sx q[0];
rz(-2.2241425) q[0];
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
rz(-1.0565856) q[0];
sx q[0];
rz(-1.5406973) q[0];
sx q[0];
rz(2.4823275) q[0];
rz(-pi) q[1];
rz(-1.2722837) q[2];
sx q[2];
rz(-1.8372756) q[2];
sx q[2];
rz(-0.68599115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.59135624) q[1];
sx q[1];
rz(-1.9758421) q[1];
sx q[1];
rz(0.77781271) q[1];
x q[2];
rz(0.81569205) q[3];
sx q[3];
rz(-0.94592735) q[3];
sx q[3];
rz(2.9790785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.36310568) q[2];
sx q[2];
rz(-2.6532756) q[2];
sx q[2];
rz(-1.5650361) q[2];
rz(2.1145084) q[3];
sx q[3];
rz(-2.4001207) q[3];
sx q[3];
rz(-1.1013793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89001369) q[0];
sx q[0];
rz(-3.0047834) q[0];
sx q[0];
rz(-0.47873163) q[0];
rz(-1.0331253) q[1];
sx q[1];
rz(-2.1703576) q[1];
sx q[1];
rz(2.1889401) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6081776) q[0];
sx q[0];
rz(-1.0262283) q[0];
sx q[0];
rz(1.4334701) q[0];
rz(-pi) q[1];
rz(0.41575899) q[2];
sx q[2];
rz(-1.8766878) q[2];
sx q[2];
rz(1.8005467) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.88541234) q[1];
sx q[1];
rz(-2.6773239) q[1];
sx q[1];
rz(-2.1229565) q[1];
rz(-2.6635366) q[3];
sx q[3];
rz(-1.607967) q[3];
sx q[3];
rz(2.2108848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4218563) q[2];
sx q[2];
rz(-0.36281261) q[2];
sx q[2];
rz(0.53058132) q[2];
rz(-1.7355708) q[3];
sx q[3];
rz(-2.0134182) q[3];
sx q[3];
rz(-0.83166844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56753165) q[0];
sx q[0];
rz(-1.644779) q[0];
sx q[0];
rz(1.5166327) q[0];
rz(1.8364871) q[1];
sx q[1];
rz(-1.790698) q[1];
sx q[1];
rz(-0.17257246) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
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
rz(-1.8048394) q[2];
sx q[2];
rz(-0.32548387) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.012430819) q[1];
sx q[1];
rz(-0.42814246) q[1];
sx q[1];
rz(1.3151602) q[1];
rz(-2.5112126) q[3];
sx q[3];
rz(-1.7582338) q[3];
sx q[3];
rz(-2.558625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3036348) q[2];
sx q[2];
rz(-1.4540318) q[2];
sx q[2];
rz(-2.0149569) q[2];
rz(-2.3593694) q[3];
sx q[3];
rz(-1.2354847) q[3];
sx q[3];
rz(1.3379898) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.85309) q[0];
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
rz(2.7998357) q[0];
sx q[0];
rz(-2.9633187) q[0];
sx q[0];
rz(-1.0210277) q[0];
rz(-pi) q[1];
rz(-0.41140326) q[2];
sx q[2];
rz(-1.3855993) q[2];
sx q[2];
rz(-1.411737) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7927671) q[1];
sx q[1];
rz(-1.7276689) q[1];
sx q[1];
rz(1.8030333) q[1];
rz(-pi) q[2];
rz(-2.8444418) q[3];
sx q[3];
rz(-2.9414256) q[3];
sx q[3];
rz(2.196892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0044272) q[2];
sx q[2];
rz(-2.9512773) q[2];
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
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5450491) q[0];
sx q[0];
rz(-0.61674917) q[0];
sx q[0];
rz(-3.074926) q[0];
rz(2.8170259) q[1];
sx q[1];
rz(-1.6371744) q[1];
sx q[1];
rz(2.1527122) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8152155) q[0];
sx q[0];
rz(-1.9003552) q[0];
sx q[0];
rz(1.9944847) q[0];
x q[1];
rz(3.0481911) q[2];
sx q[2];
rz(-1.5867234) q[2];
sx q[2];
rz(1.8708558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7372072) q[1];
sx q[1];
rz(-2.7215241) q[1];
sx q[1];
rz(-1.7350446) q[1];
x q[2];
rz(2.3905121) q[3];
sx q[3];
rz(-2.4820648) q[3];
sx q[3];
rz(-0.86576033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.4618335) q[2];
sx q[2];
rz(-2.2448848) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(2.3826564) q[0];
sx q[0];
rz(-1.8396682) q[0];
sx q[0];
rz(-2.5323903) q[0];
rz(3.0464879) q[1];
sx q[1];
rz(-1.2520049) q[1];
sx q[1];
rz(0.87337714) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4163923) q[0];
sx q[0];
rz(-1.7500449) q[0];
sx q[0];
rz(-0.068530131) q[0];
x q[1];
rz(-2.9940669) q[2];
sx q[2];
rz(-1.4359056) q[2];
sx q[2];
rz(2.877176) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.59233353) q[1];
sx q[1];
rz(-2.1570286) q[1];
sx q[1];
rz(-1.9418282) q[1];
rz(-pi) q[2];
rz(-2.7018413) q[3];
sx q[3];
rz(-1.0381178) q[3];
sx q[3];
rz(0.38910481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(0.69797126) q[0];
sx q[0];
rz(-1.781783) q[0];
sx q[0];
rz(-2.1886254) q[0];
rz(0.8264181) q[1];
sx q[1];
rz(-0.73917878) q[1];
sx q[1];
rz(-1.7451161) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1288426) q[0];
sx q[0];
rz(-1.2484974) q[0];
sx q[0];
rz(1.4803356) q[0];
rz(-2.4834677) q[2];
sx q[2];
rz(-1.8089559) q[2];
sx q[2];
rz(2.4326774) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9979981) q[1];
sx q[1];
rz(-1.2600139) q[1];
sx q[1];
rz(-0.29923156) q[1];
rz(-pi) q[2];
rz(-2.9328437) q[3];
sx q[3];
rz(-0.26780805) q[3];
sx q[3];
rz(2.1684614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.46618) q[2];
sx q[2];
rz(-2.7853577) q[2];
sx q[2];
rz(2.9818025) q[2];
rz(-0.30188489) q[3];
sx q[3];
rz(-0.92697898) q[3];
sx q[3];
rz(0.2872428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0970584) q[0];
sx q[0];
rz(-2.4659768) q[0];
sx q[0];
rz(1.5855047) q[0];
rz(0.13327577) q[1];
sx q[1];
rz(-1.6242846) q[1];
sx q[1];
rz(-0.12856738) q[1];
rz(2.2014387) q[2];
sx q[2];
rz(-1.7241782) q[2];
sx q[2];
rz(0.45964514) q[2];
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
