OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4361753) q[0];
sx q[0];
rz(-0.56646148) q[0];
sx q[0];
rz(0.17106549) q[0];
rz(-1.6879727) q[1];
sx q[1];
rz(-2.7984518) q[1];
sx q[1];
rz(-1.8106102) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2301807) q[0];
sx q[0];
rz(-0.6517621) q[0];
sx q[0];
rz(-3.0551811) q[0];
rz(-3.1034307) q[2];
sx q[2];
rz(-2.6254203) q[2];
sx q[2];
rz(-1.6463726) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0806182) q[1];
sx q[1];
rz(-1.3197101) q[1];
sx q[1];
rz(0.017647839) q[1];
rz(-2.8475548) q[3];
sx q[3];
rz(-2.4774744) q[3];
sx q[3];
rz(-2.3703863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3657637) q[2];
sx q[2];
rz(-0.77314955) q[2];
sx q[2];
rz(-1.8120871) q[2];
rz(-1.4154411) q[3];
sx q[3];
rz(-1.5664145) q[3];
sx q[3];
rz(-2.9392488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1502894) q[0];
sx q[0];
rz(-0.54420272) q[0];
sx q[0];
rz(-1.6888899) q[0];
rz(-0.20092043) q[1];
sx q[1];
rz(-2.0566437) q[1];
sx q[1];
rz(2.8413049) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18314221) q[0];
sx q[0];
rz(-2.169325) q[0];
sx q[0];
rz(-1.7167702) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47313182) q[2];
sx q[2];
rz(-2.5169249) q[2];
sx q[2];
rz(1.2166789) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.049388) q[1];
sx q[1];
rz(-0.37282473) q[1];
sx q[1];
rz(2.2224109) q[1];
rz(1.3817915) q[3];
sx q[3];
rz(-1.731589) q[3];
sx q[3];
rz(2.8652428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.77256569) q[2];
sx q[2];
rz(-0.3521266) q[2];
sx q[2];
rz(1.0380113) q[2];
rz(-0.84233061) q[3];
sx q[3];
rz(-1.4146283) q[3];
sx q[3];
rz(-2.7712908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5791941) q[0];
sx q[0];
rz(-2.6694522) q[0];
sx q[0];
rz(0.38811362) q[0];
rz(-0.072470486) q[1];
sx q[1];
rz(-1.4265172) q[1];
sx q[1];
rz(-0.31633502) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64992031) q[0];
sx q[0];
rz(-2.8965817) q[0];
sx q[0];
rz(-1.5632167) q[0];
rz(0.57245589) q[2];
sx q[2];
rz(-1.2081895) q[2];
sx q[2];
rz(2.2628757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.12831941) q[1];
sx q[1];
rz(-2.0187906) q[1];
sx q[1];
rz(-2.0112579) q[1];
x q[2];
rz(1.9100788) q[3];
sx q[3];
rz(-1.8456036) q[3];
sx q[3];
rz(2.2846082) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0324273) q[2];
sx q[2];
rz(-2.9563603) q[2];
sx q[2];
rz(-2.9807828) q[2];
rz(0.12604776) q[3];
sx q[3];
rz(-1.7536609) q[3];
sx q[3];
rz(1.0236615) q[3];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2407103) q[0];
sx q[0];
rz(-2.4375589) q[0];
sx q[0];
rz(2.3098992) q[0];
rz(1.3446993) q[1];
sx q[1];
rz(-2.1422377) q[1];
sx q[1];
rz(1.5136738) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2154402) q[0];
sx q[0];
rz(-2.2427796) q[0];
sx q[0];
rz(1.5120904) q[0];
rz(-pi) q[1];
rz(-1.6985745) q[2];
sx q[2];
rz(-2.5737692) q[2];
sx q[2];
rz(0.77726269) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1308243) q[1];
sx q[1];
rz(-1.1225892) q[1];
sx q[1];
rz(-3.1229449) q[1];
rz(-0.45659153) q[3];
sx q[3];
rz(-0.78568229) q[3];
sx q[3];
rz(2.8031138) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4251129) q[2];
sx q[2];
rz(-1.0093062) q[2];
sx q[2];
rz(-1.9173737) q[2];
rz(-2.7246357) q[3];
sx q[3];
rz(-1.0427534) q[3];
sx q[3];
rz(-2.6127889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.058218) q[0];
sx q[0];
rz(-2.871802) q[0];
sx q[0];
rz(-1.303724) q[0];
rz(-2.6858792) q[1];
sx q[1];
rz(-2.8790751) q[1];
sx q[1];
rz(-0.051503332) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7464375) q[0];
sx q[0];
rz(-1.4711385) q[0];
sx q[0];
rz(-1.6874203) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9872143) q[2];
sx q[2];
rz(-1.8533684) q[2];
sx q[2];
rz(2.4992361) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0609378) q[1];
sx q[1];
rz(-0.47532156) q[1];
sx q[1];
rz(-1.9156384) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0052161) q[3];
sx q[3];
rz(-2.4032421) q[3];
sx q[3];
rz(2.7969558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.45433989) q[2];
sx q[2];
rz(-1.3663102) q[2];
sx q[2];
rz(-2.5435737) q[2];
rz(-2.5915742) q[3];
sx q[3];
rz(-2.9019182) q[3];
sx q[3];
rz(-1.1365183) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0705868) q[0];
sx q[0];
rz(-0.07659176) q[0];
sx q[0];
rz(0.20198527) q[0];
rz(-0.96549353) q[1];
sx q[1];
rz(-1.0243203) q[1];
sx q[1];
rz(-0.083267033) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.252713) q[0];
sx q[0];
rz(-1.7656592) q[0];
sx q[0];
rz(-1.7478419) q[0];
rz(1.8385356) q[2];
sx q[2];
rz(-1.8362852) q[2];
sx q[2];
rz(2.7051085) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6492918) q[1];
sx q[1];
rz(-1.9363083) q[1];
sx q[1];
rz(0.11458061) q[1];
rz(-0.45883026) q[3];
sx q[3];
rz(-2.5341431) q[3];
sx q[3];
rz(-1.6096887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0319556) q[2];
sx q[2];
rz(-1.457931) q[2];
sx q[2];
rz(1.3827682) q[2];
rz(-0.2494732) q[3];
sx q[3];
rz(-1.8981257) q[3];
sx q[3];
rz(-1.680254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7145342) q[0];
sx q[0];
rz(-2.336851) q[0];
sx q[0];
rz(0.89685857) q[0];
rz(-0.25009051) q[1];
sx q[1];
rz(-0.18083328) q[1];
sx q[1];
rz(-2.0239963) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2876802) q[0];
sx q[0];
rz(-2.2932862) q[0];
sx q[0];
rz(-1.7476837) q[0];
rz(-2.932981) q[2];
sx q[2];
rz(-2.2715855) q[2];
sx q[2];
rz(-2.029062) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.7270131) q[1];
sx q[1];
rz(-2.1257183) q[1];
sx q[1];
rz(-1.8068061) q[1];
rz(-pi) q[2];
rz(0.53447978) q[3];
sx q[3];
rz(-1.8090873) q[3];
sx q[3];
rz(1.0907382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0573132) q[2];
sx q[2];
rz(-1.3998312) q[2];
sx q[2];
rz(2.9220667) q[2];
rz(-0.38671842) q[3];
sx q[3];
rz(-0.76824776) q[3];
sx q[3];
rz(-2.1462671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0893843) q[0];
sx q[0];
rz(-1.5970255) q[0];
sx q[0];
rz(-0.21729939) q[0];
rz(-3.1106588) q[1];
sx q[1];
rz(-2.5035281) q[1];
sx q[1];
rz(2.4826179) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9650759) q[0];
sx q[0];
rz(-2.342431) q[0];
sx q[0];
rz(-1.6060711) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4023513) q[2];
sx q[2];
rz(-0.87202245) q[2];
sx q[2];
rz(-3.0657363) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.98098552) q[1];
sx q[1];
rz(-1.7874582) q[1];
sx q[1];
rz(-0.76088455) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6625255) q[3];
sx q[3];
rz(-0.64507161) q[3];
sx q[3];
rz(1.0969539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7028246) q[2];
sx q[2];
rz(-1.3805026) q[2];
sx q[2];
rz(0.32021114) q[2];
rz(-1.5252339) q[3];
sx q[3];
rz(-1.9459008) q[3];
sx q[3];
rz(-1.2493791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-2.3808688) q[0];
sx q[0];
rz(-2.9601233) q[0];
sx q[0];
rz(0.6828126) q[0];
rz(1.0702417) q[1];
sx q[1];
rz(-1.2425334) q[1];
sx q[1];
rz(-1.210093) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34459201) q[0];
sx q[0];
rz(-1.7779777) q[0];
sx q[0];
rz(1.209757) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2145043) q[2];
sx q[2];
rz(-0.63535832) q[2];
sx q[2];
rz(-2.2941342) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.71483892) q[1];
sx q[1];
rz(-2.4086082) q[1];
sx q[1];
rz(2.0183802) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61981598) q[3];
sx q[3];
rz(-1.2634988) q[3];
sx q[3];
rz(-2.4710771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5658297) q[2];
sx q[2];
rz(-2.1255707) q[2];
sx q[2];
rz(-2.5749717) q[2];
rz(2.2120655) q[3];
sx q[3];
rz(-0.98171392) q[3];
sx q[3];
rz(1.367759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.728445) q[0];
sx q[0];
rz(-2.8840273) q[0];
sx q[0];
rz(-1.43191) q[0];
rz(2.6152949) q[1];
sx q[1];
rz(-0.52733517) q[1];
sx q[1];
rz(-2.3419103) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0737338) q[0];
sx q[0];
rz(-1.4859383) q[0];
sx q[0];
rz(-1.2220864) q[0];
rz(-pi) q[1];
rz(0.46634679) q[2];
sx q[2];
rz(-2.9620167) q[2];
sx q[2];
rz(-1.3122383) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.666115) q[1];
sx q[1];
rz(-2.0787422) q[1];
sx q[1];
rz(2.75027) q[1];
rz(-pi) q[2];
rz(-1.4950072) q[3];
sx q[3];
rz(-1.7160551) q[3];
sx q[3];
rz(0.69243542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8563103) q[2];
sx q[2];
rz(-2.6976863) q[2];
sx q[2];
rz(0.69520673) q[2];
rz(-2.5837512) q[3];
sx q[3];
rz(-1.7772243) q[3];
sx q[3];
rz(0.69303304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0859062) q[0];
sx q[0];
rz(-1.9491371) q[0];
sx q[0];
rz(0.51167713) q[0];
rz(1.862539) q[1];
sx q[1];
rz(-2.3976354) q[1];
sx q[1];
rz(2.4639113) q[1];
rz(0.002481133) q[2];
sx q[2];
rz(-0.33591349) q[2];
sx q[2];
rz(-2.9813067) q[2];
rz(2.9813319) q[3];
sx q[3];
rz(-2.3370623) q[3];
sx q[3];
rz(-0.2094895) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
