OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.32831353) q[0];
sx q[0];
rz(5.1483122) q[0];
sx q[0];
rz(6.8490646) q[0];
rz(-2.8973051) q[1];
sx q[1];
rz(-1.7841508) q[1];
sx q[1];
rz(2.607333) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0097144011) q[0];
sx q[0];
rz(-1.3004954) q[0];
sx q[0];
rz(0.86608315) q[0];
rz(-pi) q[1];
rz(-2.7083394) q[2];
sx q[2];
rz(-0.64068551) q[2];
sx q[2];
rz(2.6034466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.104887) q[1];
sx q[1];
rz(-1.5944696) q[1];
sx q[1];
rz(2.8410683) q[1];
rz(-1.8194866) q[3];
sx q[3];
rz(-1.1258896) q[3];
sx q[3];
rz(1.1550732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.013022097) q[2];
sx q[2];
rz(-2.8910922) q[2];
sx q[2];
rz(0.24918431) q[2];
rz(1.7976044) q[3];
sx q[3];
rz(-1.598571) q[3];
sx q[3];
rz(0.070076076) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4942112) q[0];
sx q[0];
rz(-1.9152315) q[0];
sx q[0];
rz(0.98980728) q[0];
rz(0.79065943) q[1];
sx q[1];
rz(-2.374687) q[1];
sx q[1];
rz(-2.8299423) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3450736) q[0];
sx q[0];
rz(-2.2537044) q[0];
sx q[0];
rz(2.2497798) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0333854) q[2];
sx q[2];
rz(-1.756486) q[2];
sx q[2];
rz(0.78150392) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6271882) q[1];
sx q[1];
rz(-1.3986821) q[1];
sx q[1];
rz(-1.3283511) q[1];
rz(-pi) q[2];
rz(0.34087025) q[3];
sx q[3];
rz(-2.1580527) q[3];
sx q[3];
rz(1.0140918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1404169) q[2];
sx q[2];
rz(-1.9813462) q[2];
sx q[2];
rz(-2.1902093) q[2];
rz(-0.22826711) q[3];
sx q[3];
rz(-2.164866) q[3];
sx q[3];
rz(1.867928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839639) q[0];
sx q[0];
rz(-1.552859) q[0];
sx q[0];
rz(-0.11446318) q[0];
rz(-1.0022256) q[1];
sx q[1];
rz(-0.4267692) q[1];
sx q[1];
rz(0.1618298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1755668) q[0];
sx q[0];
rz(-1.9593666) q[0];
sx q[0];
rz(-2.2184664) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4828909) q[2];
sx q[2];
rz(-1.1832596) q[2];
sx q[2];
rz(2.0883002) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5787631) q[1];
sx q[1];
rz(-2.17318) q[1];
sx q[1];
rz(1.2788601) q[1];
rz(1.9818241) q[3];
sx q[3];
rz(-0.81361249) q[3];
sx q[3];
rz(0.66853722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.74730045) q[2];
sx q[2];
rz(-2.6960399) q[2];
sx q[2];
rz(3.0008345) q[2];
rz(-2.5898139) q[3];
sx q[3];
rz(-1.8675624) q[3];
sx q[3];
rz(2.3902635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.9363339) q[0];
sx q[0];
rz(-2.8755499) q[0];
sx q[0];
rz(0.67489135) q[0];
rz(0.15448013) q[1];
sx q[1];
rz(-1.7325502) q[1];
sx q[1];
rz(-2.6836269) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5602011) q[0];
sx q[0];
rz(-2.0580252) q[0];
sx q[0];
rz(-1.1283406) q[0];
rz(-pi) q[1];
rz(1.141233) q[2];
sx q[2];
rz(-1.2899961) q[2];
sx q[2];
rz(-2.8605638) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.030589015) q[1];
sx q[1];
rz(-1.2007502) q[1];
sx q[1];
rz(-1.3732598) q[1];
x q[2];
rz(2.9513957) q[3];
sx q[3];
rz(-0.51589291) q[3];
sx q[3];
rz(-0.35515337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0080879) q[2];
sx q[2];
rz(-2.4613481) q[2];
sx q[2];
rz(-0.31944719) q[2];
rz(1.8114629) q[3];
sx q[3];
rz(-2.1250686) q[3];
sx q[3];
rz(-2.5813848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9820246) q[0];
sx q[0];
rz(-1.1363131) q[0];
sx q[0];
rz(1.9200448) q[0];
rz(0.23280652) q[1];
sx q[1];
rz(-0.56929749) q[1];
sx q[1];
rz(-2.1585042) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1637835) q[0];
sx q[0];
rz(-2.0931232) q[0];
sx q[0];
rz(0.99951007) q[0];
rz(0.95011656) q[2];
sx q[2];
rz(-1.4779556) q[2];
sx q[2];
rz(2.8454121) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0016664) q[1];
sx q[1];
rz(-1.3845433) q[1];
sx q[1];
rz(1.7193394) q[1];
x q[2];
rz(1.9290673) q[3];
sx q[3];
rz(-2.6942109) q[3];
sx q[3];
rz(-3.0606396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7588707) q[2];
sx q[2];
rz(-1.657234) q[2];
sx q[2];
rz(-2.7663084) q[2];
rz(-1.075607) q[3];
sx q[3];
rz(-0.38638249) q[3];
sx q[3];
rz(-0.53494278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4009092) q[0];
sx q[0];
rz(-1.4367737) q[0];
sx q[0];
rz(-1.7042473) q[0];
rz(3.0628693) q[1];
sx q[1];
rz(-1.3767786) q[1];
sx q[1];
rz(-2.5934503) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.948471) q[0];
sx q[0];
rz(-1.7821494) q[0];
sx q[0];
rz(0.67074346) q[0];
rz(-0.70912045) q[2];
sx q[2];
rz(-0.68163423) q[2];
sx q[2];
rz(-0.97709419) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6316165) q[1];
sx q[1];
rz(-2.2215034) q[1];
sx q[1];
rz(-2.0664311) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75597922) q[3];
sx q[3];
rz(-2.0507617) q[3];
sx q[3];
rz(0.85473138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4569725) q[2];
sx q[2];
rz(-0.62070864) q[2];
sx q[2];
rz(2.4662628) q[2];
rz(1.5572549) q[3];
sx q[3];
rz(-1.8341589) q[3];
sx q[3];
rz(-2.5244782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5976582) q[0];
sx q[0];
rz(-1.1827396) q[0];
sx q[0];
rz(-2.6455998) q[0];
rz(-1.4542106) q[1];
sx q[1];
rz(-2.0861237) q[1];
sx q[1];
rz(1.1167663) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3163213) q[0];
sx q[0];
rz(-0.50646913) q[0];
sx q[0];
rz(-2.3536918) q[0];
x q[1];
rz(-2.871795) q[2];
sx q[2];
rz(-1.0493663) q[2];
sx q[2];
rz(-2.9228589) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3254814) q[1];
sx q[1];
rz(-1.2524091) q[1];
sx q[1];
rz(-0.35811425) q[1];
x q[2];
rz(-0.45601042) q[3];
sx q[3];
rz(-1.5152165) q[3];
sx q[3];
rz(-1.0694651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3108959) q[2];
sx q[2];
rz(-0.98429698) q[2];
sx q[2];
rz(-2.7867219) q[2];
rz(-0.76534671) q[3];
sx q[3];
rz(-2.2782875) q[3];
sx q[3];
rz(-0.089546831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.68882051) q[0];
sx q[0];
rz(-1.4790164) q[0];
sx q[0];
rz(2.3759957) q[0];
rz(-2.2701524) q[1];
sx q[1];
rz(-0.7656509) q[1];
sx q[1];
rz(1.8188459) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2521556) q[0];
sx q[0];
rz(-1.3103292) q[0];
sx q[0];
rz(1.5171264) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9051612) q[2];
sx q[2];
rz(-0.90518236) q[2];
sx q[2];
rz(-2.1588391) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2603399) q[1];
sx q[1];
rz(-1.5543206) q[1];
sx q[1];
rz(0.50838281) q[1];
rz(-pi) q[2];
rz(2.3115308) q[3];
sx q[3];
rz(-0.30127159) q[3];
sx q[3];
rz(-0.64727441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4837997) q[2];
sx q[2];
rz(-2.0988266) q[2];
sx q[2];
rz(-0.22171177) q[2];
rz(-2.0486369) q[3];
sx q[3];
rz(-1.4128128) q[3];
sx q[3];
rz(0.67813412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0677277) q[0];
sx q[0];
rz(-1.8949969) q[0];
sx q[0];
rz(-2.8644323) q[0];
rz(0.74835888) q[1];
sx q[1];
rz(-2.6942418) q[1];
sx q[1];
rz(1.1433196) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70049452) q[0];
sx q[0];
rz(-1.3702533) q[0];
sx q[0];
rz(1.8822549) q[0];
rz(0.21094811) q[2];
sx q[2];
rz(-0.8767414) q[2];
sx q[2];
rz(-1.2488493) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4797693) q[1];
sx q[1];
rz(-2.7854438) q[1];
sx q[1];
rz(-1.1060064) q[1];
rz(-2.7797749) q[3];
sx q[3];
rz(-1.731195) q[3];
sx q[3];
rz(-0.058506207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.9144168) q[2];
sx q[2];
rz(-1.02355) q[2];
sx q[2];
rz(3.0926256) q[2];
rz(-2.07552) q[3];
sx q[3];
rz(-0.62838307) q[3];
sx q[3];
rz(0.15570417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0186036) q[0];
sx q[0];
rz(-1.5276696) q[0];
sx q[0];
rz(-3.1220806) q[0];
rz(-0.77876577) q[1];
sx q[1];
rz(-2.4486783) q[1];
sx q[1];
rz(-1.5787517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1922958) q[0];
sx q[0];
rz(-1.6598633) q[0];
sx q[0];
rz(1.7013447) q[0];
x q[1];
rz(0.18264211) q[2];
sx q[2];
rz(-2.172875) q[2];
sx q[2];
rz(1.4878359) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.409118) q[1];
sx q[1];
rz(-0.98608398) q[1];
sx q[1];
rz(-0.52701108) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7354129) q[3];
sx q[3];
rz(-2.4009628) q[3];
sx q[3];
rz(-0.1333789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.94893) q[2];
sx q[2];
rz(-0.86343416) q[2];
sx q[2];
rz(0.16555244) q[2];
rz(0.60756573) q[3];
sx q[3];
rz(-1.2847885) q[3];
sx q[3];
rz(-2.6732388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3792569) q[0];
sx q[0];
rz(-2.6597334) q[0];
sx q[0];
rz(-0.045482176) q[0];
rz(-1.634585) q[1];
sx q[1];
rz(-1.4611117) q[1];
sx q[1];
rz(-1.5932105) q[1];
rz(-1.4277369) q[2];
sx q[2];
rz(-0.80052994) q[2];
sx q[2];
rz(1.1442727) q[2];
rz(-2.8589917) q[3];
sx q[3];
rz(-2.5167214) q[3];
sx q[3];
rz(-1.4759397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
