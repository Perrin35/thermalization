OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(-3.0769899) q[0];
sx q[0];
rz(-0.021615418) q[0];
rz(2.1463483) q[1];
sx q[1];
rz(-1.8145476) q[1];
sx q[1];
rz(-1.8099161) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0639122) q[0];
sx q[0];
rz(-2.4644116) q[0];
sx q[0];
rz(2.1190686) q[0];
rz(-1.1515491) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(-2.575945) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.761258) q[1];
sx q[1];
rz(-2.3654656) q[1];
sx q[1];
rz(2.9556729) q[1];
rz(-pi) q[2];
rz(0.33820037) q[3];
sx q[3];
rz(-2.1669398) q[3];
sx q[3];
rz(-2.3232834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(-2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3246831) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(-1.0634134) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(-0.0016454776) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7749274) q[0];
sx q[0];
rz(-1.6158551) q[0];
sx q[0];
rz(0.010618322) q[0];
x q[1];
rz(-1.3130982) q[2];
sx q[2];
rz(-2.1971467) q[2];
sx q[2];
rz(1.4305654) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26047036) q[1];
sx q[1];
rz(-2.3350041) q[1];
sx q[1];
rz(1.9116198) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.315829) q[3];
sx q[3];
rz(-1.7353188) q[3];
sx q[3];
rz(-0.69873519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1559747) q[2];
sx q[2];
rz(-1.6321471) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(2.0478785) q[3];
sx q[3];
rz(-1.8023068) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22096069) q[0];
sx q[0];
rz(-1.4039803) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(-3.1365085) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(1.089383) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1339061) q[0];
sx q[0];
rz(-2.0984762) q[0];
sx q[0];
rz(2.6469127) q[0];
rz(-0.58615746) q[2];
sx q[2];
rz(-0.68550368) q[2];
sx q[2];
rz(3.0229178) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.032420302) q[1];
sx q[1];
rz(-1.65616) q[1];
sx q[1];
rz(-1.1847772) q[1];
rz(-pi) q[2];
rz(2.1316707) q[3];
sx q[3];
rz(-0.86285931) q[3];
sx q[3];
rz(-0.55299711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0638782) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(-0.88469488) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(2.4147721) q[0];
rz(-2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(0.35983905) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.737405) q[0];
sx q[0];
rz(-0.92659896) q[0];
sx q[0];
rz(-0.23097158) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9927822) q[2];
sx q[2];
rz(-2.4850922) q[2];
sx q[2];
rz(1.3615001) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6872014) q[1];
sx q[1];
rz(-1.193421) q[1];
sx q[1];
rz(2.6881933) q[1];
rz(1.9197649) q[3];
sx q[3];
rz(-1.8604606) q[3];
sx q[3];
rz(1.3834013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5741253) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(-2.3763669) q[2];
rz(-0.75677538) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9969479) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(-1.416052) q[0];
rz(-2.7744746) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(1.0353154) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028776289) q[0];
sx q[0];
rz(-2.007764) q[0];
sx q[0];
rz(-0.90941888) q[0];
rz(-pi) q[1];
rz(2.6501422) q[2];
sx q[2];
rz(-1.5361668) q[2];
sx q[2];
rz(2.2143242) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.1086515) q[1];
sx q[1];
rz(-1.6254394) q[1];
sx q[1];
rz(0.095468949) q[1];
x q[2];
rz(-2.8661795) q[3];
sx q[3];
rz(-1.3684891) q[3];
sx q[3];
rz(1.3261258) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3570024) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(2.3727097) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-2.3612645) q[3];
sx q[3];
rz(-1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.996421) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(2.9343658) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0475537) q[0];
sx q[0];
rz(-0.63699603) q[0];
sx q[0];
rz(2.9061222) q[0];
rz(-pi) q[1];
rz(1.1057165) q[2];
sx q[2];
rz(-1.0208703) q[2];
sx q[2];
rz(-1.20649) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8689649) q[1];
sx q[1];
rz(-2.7902863) q[1];
sx q[1];
rz(2.3428194) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3239922) q[3];
sx q[3];
rz(-1.6896276) q[3];
sx q[3];
rz(2.6314051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8032288) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(2.4198789) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(-2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24467829) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-2.4095643) q[0];
rz(3.1320944) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(2.8498555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4371944) q[0];
sx q[0];
rz(-1.6219553) q[0];
sx q[0];
rz(-0.41914661) q[0];
rz(3.0351699) q[2];
sx q[2];
rz(-1.2692045) q[2];
sx q[2];
rz(2.0472722) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.74624324) q[1];
sx q[1];
rz(-0.33848539) q[1];
sx q[1];
rz(2.6836718) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2613867) q[3];
sx q[3];
rz(-2.072812) q[3];
sx q[3];
rz(1.3537784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8774524) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(-0.83731246) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(0.062019197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432805) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(0.6638546) q[0];
rz(-0.10617667) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(2.1829139) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3631358) q[0];
sx q[0];
rz(-2.083337) q[0];
sx q[0];
rz(-0.97495671) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27997048) q[2];
sx q[2];
rz(-2.0915871) q[2];
sx q[2];
rz(0.94172804) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.9740323) q[1];
sx q[1];
rz(-1.6775727) q[1];
sx q[1];
rz(-1.0436996) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18387353) q[3];
sx q[3];
rz(-2.549578) q[3];
sx q[3];
rz(-1.5873991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.133193) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(-0.91910249) q[2];
rz(-1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(-1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(2.7046955) q[0];
rz(2.4412952) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(2.6760496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1818905) q[0];
sx q[0];
rz(-2.2756016) q[0];
sx q[0];
rz(1.4591135) q[0];
rz(-pi) q[1];
rz(1.2836254) q[2];
sx q[2];
rz(-0.42771491) q[2];
sx q[2];
rz(-2.5123793) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1411966) q[1];
sx q[1];
rz(-1.2277514) q[1];
sx q[1];
rz(2.6617906) q[1];
rz(-pi) q[2];
rz(1.9302619) q[3];
sx q[3];
rz(-1.3350447) q[3];
sx q[3];
rz(-2.204493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.4979866) q[2];
rz(2.9368029) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(-0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64514226) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(-1.1219332) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(2.5591992) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84401417) q[0];
sx q[0];
rz(-1.3676924) q[0];
sx q[0];
rz(-0.88589478) q[0];
x q[1];
rz(-2.9899644) q[2];
sx q[2];
rz(-1.0851589) q[2];
sx q[2];
rz(3.1079353) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0092587) q[1];
sx q[1];
rz(-1.7528105) q[1];
sx q[1];
rz(-0.11972129) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12606975) q[3];
sx q[3];
rz(-0.21890103) q[3];
sx q[3];
rz(0.88760469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(0.76254145) q[2];
rz(-1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5738752) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0653771) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(1.8021884) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(-2.8600678) q[2];
sx q[2];
rz(-0.83419656) q[2];
sx q[2];
rz(0.4783334) q[2];
rz(-0.4776095) q[3];
sx q[3];
rz(-2.0255247) q[3];
sx q[3];
rz(-0.075102641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];