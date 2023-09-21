OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.83710837) q[0];
sx q[0];
rz(-1.4533071) q[0];
sx q[0];
rz(0.31153554) q[0];
rz(2.7040634) q[1];
sx q[1];
rz(-1.3181926) q[1];
sx q[1];
rz(-0.55895609) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5200978) q[0];
sx q[0];
rz(-1.1557475) q[0];
sx q[0];
rz(-0.15226224) q[0];
rz(-1.5833601) q[2];
sx q[2];
rz(-1.0135092) q[2];
sx q[2];
rz(-1.81665) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.641687) q[1];
sx q[1];
rz(-0.99398621) q[1];
sx q[1];
rz(2.4813502) q[1];
rz(-pi) q[2];
rz(0.43879126) q[3];
sx q[3];
rz(-1.3109129) q[3];
sx q[3];
rz(-0.70513844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.3922334) q[2];
sx q[2];
rz(-1.8584676) q[2];
sx q[2];
rz(2.5048845) q[2];
rz(2.2926245) q[3];
sx q[3];
rz(-0.62148062) q[3];
sx q[3];
rz(-0.23392114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37671509) q[0];
sx q[0];
rz(-2.8945518) q[0];
sx q[0];
rz(-2.9887181) q[0];
rz(-2.3846467) q[1];
sx q[1];
rz(-1.5870973) q[1];
sx q[1];
rz(-0.98639948) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7727535) q[0];
sx q[0];
rz(-1.5305688) q[0];
sx q[0];
rz(-3.0818135) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.7631626) q[2];
sx q[2];
rz(-1.1726511) q[2];
sx q[2];
rz(2.5415908) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1101148) q[1];
sx q[1];
rz(-1.8077393) q[1];
sx q[1];
rz(-2.7538339) q[1];
rz(-pi) q[2];
rz(-0.13568474) q[3];
sx q[3];
rz(-1.8968582) q[3];
sx q[3];
rz(0.56860926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5793005) q[2];
sx q[2];
rz(-1.9220756) q[2];
sx q[2];
rz(2.3584649) q[2];
rz(-3.1230208) q[3];
sx q[3];
rz(-1.637807) q[3];
sx q[3];
rz(-0.40772453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0531533) q[0];
sx q[0];
rz(-0.29215559) q[0];
sx q[0];
rz(2.1799178) q[0];
rz(-0.36034521) q[1];
sx q[1];
rz(-2.0397489) q[1];
sx q[1];
rz(-0.12869421) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6260687) q[0];
sx q[0];
rz(-1.6532073) q[0];
sx q[0];
rz(-3.0936196) q[0];
x q[1];
rz(1.0110858) q[2];
sx q[2];
rz(-0.84583827) q[2];
sx q[2];
rz(2.9618008) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2663181) q[1];
sx q[1];
rz(-2.6058063) q[1];
sx q[1];
rz(-2.5712625) q[1];
rz(-pi) q[2];
rz(-1.7101173) q[3];
sx q[3];
rz(-2.3686667) q[3];
sx q[3];
rz(-1.3611925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0814357) q[2];
sx q[2];
rz(-1.1971985) q[2];
sx q[2];
rz(1.241768) q[2];
rz(0.5870108) q[3];
sx q[3];
rz(-2.185052) q[3];
sx q[3];
rz(-0.97755066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.509165) q[0];
sx q[0];
rz(-2.2594663) q[0];
sx q[0];
rz(1.084491) q[0];
rz(1.658461) q[1];
sx q[1];
rz(-0.56743923) q[1];
sx q[1];
rz(0.09253563) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0412484) q[0];
sx q[0];
rz(-0.39533246) q[0];
sx q[0];
rz(-0.80907099) q[0];
x q[1];
rz(-2.7961568) q[2];
sx q[2];
rz(-1.1159117) q[2];
sx q[2];
rz(-2.5330184) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9738237) q[1];
sx q[1];
rz(-1.2281706) q[1];
sx q[1];
rz(1.6325566) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.16102287) q[3];
sx q[3];
rz(-0.92217126) q[3];
sx q[3];
rz(-0.30740689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3080421) q[2];
sx q[2];
rz(-1.4414859) q[2];
sx q[2];
rz(0.33205024) q[2];
rz(-2.0856693) q[3];
sx q[3];
rz(-0.27763405) q[3];
sx q[3];
rz(2.5312996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32245359) q[0];
sx q[0];
rz(-1.8503014) q[0];
sx q[0];
rz(0.21155587) q[0];
rz(1.8353204) q[1];
sx q[1];
rz(-1.897656) q[1];
sx q[1];
rz(2.4938915) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54007285) q[0];
sx q[0];
rz(-1.411502) q[0];
sx q[0];
rz(-1.457731) q[0];
rz(-pi) q[1];
rz(0.34824246) q[2];
sx q[2];
rz(-1.4255382) q[2];
sx q[2];
rz(0.20656221) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.60407818) q[1];
sx q[1];
rz(-1.068183) q[1];
sx q[1];
rz(-0.22488774) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1224498) q[3];
sx q[3];
rz(-1.5557516) q[3];
sx q[3];
rz(-2.6263833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.56090474) q[2];
sx q[2];
rz(-0.40955341) q[2];
sx q[2];
rz(2.4482751) q[2];
rz(2.4723315) q[3];
sx q[3];
rz(-1.4368493) q[3];
sx q[3];
rz(3.1366689) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3174021) q[0];
sx q[0];
rz(-0.22739246) q[0];
sx q[0];
rz(-1.2325226) q[0];
rz(-1.0725853) q[1];
sx q[1];
rz(-2.0697846) q[1];
sx q[1];
rz(-2.9673064) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6264544) q[0];
sx q[0];
rz(-0.8695375) q[0];
sx q[0];
rz(2.7324972) q[0];
rz(2.2767378) q[2];
sx q[2];
rz(-1.0360498) q[2];
sx q[2];
rz(2.6851482) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.6433405) q[1];
sx q[1];
rz(-0.87264112) q[1];
sx q[1];
rz(-2.6367285) q[1];
rz(1.918119) q[3];
sx q[3];
rz(-2.8109549) q[3];
sx q[3];
rz(-0.91352458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8217414) q[2];
sx q[2];
rz(-2.3345626) q[2];
sx q[2];
rz(2.9439587) q[2];
rz(2.8526784) q[3];
sx q[3];
rz(-2.2646326) q[3];
sx q[3];
rz(-1.4060085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.063868) q[0];
sx q[0];
rz(-2.6517695) q[0];
sx q[0];
rz(-0.20859627) q[0];
rz(-2.1754307) q[1];
sx q[1];
rz(-1.9816793) q[1];
sx q[1];
rz(1.6360412) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0783605) q[0];
sx q[0];
rz(-2.1713543) q[0];
sx q[0];
rz(-0.91852228) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43039544) q[2];
sx q[2];
rz(-1.5376523) q[2];
sx q[2];
rz(2.7660649) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5559616) q[1];
sx q[1];
rz(-1.0502083) q[1];
sx q[1];
rz(-0.52557892) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.286245) q[3];
sx q[3];
rz(-1.6702594) q[3];
sx q[3];
rz(3.0486097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.85764) q[2];
sx q[2];
rz(-0.74308926) q[2];
sx q[2];
rz(3.0440142) q[2];
rz(-1.7476667) q[3];
sx q[3];
rz(-1.7912309) q[3];
sx q[3];
rz(2.424749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7512648) q[0];
sx q[0];
rz(-1.8126235) q[0];
sx q[0];
rz(-0.51399291) q[0];
rz(0.12318525) q[1];
sx q[1];
rz(-0.24736483) q[1];
sx q[1];
rz(-2.2095912) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31873736) q[0];
sx q[0];
rz(-1.1294951) q[0];
sx q[0];
rz(0.1499099) q[0];
x q[1];
rz(1.9281689) q[2];
sx q[2];
rz(-2.4744611) q[2];
sx q[2];
rz(2.2156029) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.060170505) q[1];
sx q[1];
rz(-2.0875071) q[1];
sx q[1];
rz(-0.64232773) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85429116) q[3];
sx q[3];
rz(-1.4201418) q[3];
sx q[3];
rz(0.52779576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1121858) q[2];
sx q[2];
rz(-2.0307348) q[2];
sx q[2];
rz(2.712148) q[2];
rz(-1.2094234) q[3];
sx q[3];
rz(-2.7691787) q[3];
sx q[3];
rz(2.4485574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76686239) q[0];
sx q[0];
rz(-1.6748036) q[0];
sx q[0];
rz(-1.8027579) q[0];
rz(-2.4354637) q[1];
sx q[1];
rz(-1.2495722) q[1];
sx q[1];
rz(-0.95058092) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50981748) q[0];
sx q[0];
rz(-1.063856) q[0];
sx q[0];
rz(-0.99960534) q[0];
rz(-1.9865932) q[2];
sx q[2];
rz(-1.6396513) q[2];
sx q[2];
rz(-1.6378251) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.49950019) q[1];
sx q[1];
rz(-1.7576808) q[1];
sx q[1];
rz(-2.1128113) q[1];
rz(-pi) q[2];
rz(2.1080984) q[3];
sx q[3];
rz(-0.47035445) q[3];
sx q[3];
rz(0.62894097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4799698) q[2];
sx q[2];
rz(-1.6500436) q[2];
sx q[2];
rz(2.6573112) q[2];
rz(-0.92710036) q[3];
sx q[3];
rz(-1.8323332) q[3];
sx q[3];
rz(0.62121975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1442239) q[0];
sx q[0];
rz(-3.0529418) q[0];
sx q[0];
rz(-0.22928672) q[0];
rz(-2.7067822) q[1];
sx q[1];
rz(-1.2229342) q[1];
sx q[1];
rz(-0.71892175) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839448) q[0];
sx q[0];
rz(-0.63417182) q[0];
sx q[0];
rz(-1.7303403) q[0];
rz(-2.1585196) q[2];
sx q[2];
rz(-1.4721774) q[2];
sx q[2];
rz(-1.7416416) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.1246008) q[1];
sx q[1];
rz(-1.4233839) q[1];
sx q[1];
rz(2.7684545) q[1];
rz(-pi) q[2];
rz(0.50456725) q[3];
sx q[3];
rz(-2.1310398) q[3];
sx q[3];
rz(1.954078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3832613) q[2];
sx q[2];
rz(-1.9394082) q[2];
sx q[2];
rz(-0.74679217) q[2];
rz(0.87219277) q[3];
sx q[3];
rz(-0.82834297) q[3];
sx q[3];
rz(0.94223589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.223021) q[0];
sx q[0];
rz(-1.3970319) q[0];
sx q[0];
rz(-1.3402517) q[0];
rz(-2.7643798) q[1];
sx q[1];
rz(-1.6709534) q[1];
sx q[1];
rz(0.51660641) q[1];
rz(-0.51858356) q[2];
sx q[2];
rz(-1.2972144) q[2];
sx q[2];
rz(-1.5757061) q[2];
rz(0.71729284) q[3];
sx q[3];
rz(-1.4581231) q[3];
sx q[3];
rz(0.067618528) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
