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
rz(0.064602764) q[0];
sx q[0];
rz(6.3048007) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0927147) q[0];
sx q[0];
rz(-1.9034916) q[0];
sx q[0];
rz(-0.96941745) q[0];
x q[1];
rz(1.3433427) q[2];
sx q[2];
rz(-1.4706352) q[2];
sx q[2];
rz(-1.7286466) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6380438) q[1];
sx q[1];
rz(-2.330144) q[1];
sx q[1];
rz(1.391295) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33820037) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(-0.81830922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8109479) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3246831) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(2.0781793) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(3.1399472) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5981818) q[0];
sx q[0];
rz(-0.046292154) q[0];
sx q[0];
rz(-1.339519) q[0];
x q[1];
rz(-1.8284945) q[2];
sx q[2];
rz(-2.1971467) q[2];
sx q[2];
rz(-1.4305654) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8811223) q[1];
sx q[1];
rz(-2.3350041) q[1];
sx q[1];
rz(-1.9116198) q[1];
rz(-pi) q[2];
rz(-1.8257636) q[3];
sx q[3];
rz(-1.7353188) q[3];
sx q[3];
rz(-2.4428575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.985618) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(1.0937141) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(-2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(2.3143342) q[0];
rz(0.0050841252) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(-2.0522096) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1339061) q[0];
sx q[0];
rz(-1.0431164) q[0];
sx q[0];
rz(-2.6469127) q[0];
x q[1];
rz(2.5554352) q[2];
sx q[2];
rz(-2.456089) q[2];
sx q[2];
rz(0.11867487) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.032420302) q[1];
sx q[1];
rz(-1.65616) q[1];
sx q[1];
rz(1.9568155) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55604071) q[3];
sx q[3];
rz(-0.87198139) q[3];
sx q[3];
rz(1.8204821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0777145) q[2];
sx q[2];
rz(-2.2641247) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(1.9583154) q[3];
sx q[3];
rz(-1.3180472) q[3];
sx q[3];
rz(-1.6515091) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5723715) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-2.4147721) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(2.7817536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4041876) q[0];
sx q[0];
rz(-0.92659896) q[0];
sx q[0];
rz(-0.23097158) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.14881046) q[2];
sx q[2];
rz(-0.65650046) q[2];
sx q[2];
rz(-1.7800926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53107809) q[1];
sx q[1];
rz(-0.58137608) q[1];
sx q[1];
rz(0.73552144) q[1];
x q[2];
rz(2.8344645) q[3];
sx q[3];
rz(-1.2369452) q[3];
sx q[3];
rz(-3.0577554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5741253) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(2.3763669) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(-0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.1446447) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(1.7255406) q[0];
rz(2.7744746) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(-2.1062772) q[1];
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
x q[1];
rz(0.073280235) q[2];
sx q[2];
rz(-2.6490232) q[2];
sx q[2];
rz(0.57893334) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0329412) q[1];
sx q[1];
rz(-1.6254394) q[1];
sx q[1];
rz(3.0461237) q[1];
x q[2];
rz(-0.27541311) q[3];
sx q[3];
rz(-1.7731035) q[3];
sx q[3];
rz(-1.8154669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.7845903) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(-0.76888293) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.996421) q[0];
rz(-2.0369453) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-0.2072269) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80400318) q[0];
sx q[0];
rz(-0.9540671) q[0];
sx q[0];
rz(1.7417275) q[0];
x q[1];
rz(-0.63164288) q[2];
sx q[2];
rz(-0.70438671) q[2];
sx q[2];
rz(-2.7001675) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2726278) q[1];
sx q[1];
rz(-2.7902863) q[1];
sx q[1];
rz(2.3428194) q[1];
rz(-2.3239922) q[3];
sx q[3];
rz(-1.451965) q[3];
sx q[3];
rz(2.6314051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3383639) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(-0.72171372) q[2];
rz(1.8668113) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(-2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.24467829) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-0.73202837) q[0];
rz(3.1320944) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(0.2917372) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1563819) q[0];
sx q[0];
rz(-1.9893601) q[0];
sx q[0];
rz(-1.5147989) q[0];
rz(-pi) q[1];
rz(1.2417492) q[2];
sx q[2];
rz(-0.31927682) q[2];
sx q[2];
rz(0.74908756) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3953494) q[1];
sx q[1];
rz(-2.8031073) q[1];
sx q[1];
rz(-2.6836718) q[1];
rz(1.2613867) q[3];
sx q[3];
rz(-2.072812) q[3];
sx q[3];
rz(1.7878143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26414028) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(-1.9705747) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(2.4777381) q[0];
rz(3.035416) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(0.95867872) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7784568) q[0];
sx q[0];
rz(-1.0582557) q[0];
sx q[0];
rz(0.97495671) q[0];
rz(-0.27997048) q[2];
sx q[2];
rz(-1.0500056) q[2];
sx q[2];
rz(-2.1998646) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9194591) q[1];
sx q[1];
rz(-0.53680116) q[1];
sx q[1];
rz(-1.7807351) q[1];
rz(-2.9577191) q[3];
sx q[3];
rz(-0.5920147) q[3];
sx q[3];
rz(-1.5873991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.133193) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(-2.2224902) q[2];
rz(-1.3778) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(-1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8063426) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(0.43689716) q[0];
rz(-0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(2.6760496) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9597022) q[0];
sx q[0];
rz(-2.2756016) q[0];
sx q[0];
rz(1.6824791) q[0];
x q[1];
rz(0.12840694) q[2];
sx q[2];
rz(-1.9799211) q[2];
sx q[2];
rz(2.1985334) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.000396) q[1];
sx q[1];
rz(-1.2277514) q[1];
sx q[1];
rz(-0.47980206) q[1];
x q[2];
rz(-1.2113308) q[3];
sx q[3];
rz(-1.3350447) q[3];
sx q[3];
rz(0.93709968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-0.53933829) q[2];
sx q[2];
rz(1.643606) q[2];
rz(2.9368029) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(-0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4964504) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(-2.0196594) q[0];
rz(-0.75138584) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(-2.5591992) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2975785) q[0];
sx q[0];
rz(-1.7739002) q[0];
sx q[0];
rz(-2.2556979) q[0];
rz(-1.0803797) q[2];
sx q[2];
rz(-1.7047802) q[2];
sx q[2];
rz(-1.6756563) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45422428) q[1];
sx q[1];
rz(-2.9240989) q[1];
sx q[1];
rz(-0.99517676) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.598761) q[3];
sx q[3];
rz(-1.3536605) q[3];
sx q[3];
rz(-1.0167227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0104684) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(-0.76254145) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(-1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(1.8680686) q[2];
sx q[2];
rz(-2.3625629) q[2];
sx q[2];
rz(0.071803781) q[2];
rz(0.4776095) q[3];
sx q[3];
rz(-1.1160679) q[3];
sx q[3];
rz(3.06649) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];