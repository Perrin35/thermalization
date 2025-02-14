OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.81808972) q[0];
sx q[0];
rz(-0.40677318) q[0];
sx q[0];
rz(0.50954252) q[0];
rz(-0.31515631) q[1];
sx q[1];
rz(6.4767467) q[1];
sx q[1];
rz(11.560796) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1664274) q[0];
sx q[0];
rz(-1.7708798) q[0];
sx q[0];
rz(1.6471662) q[0];
x q[1];
rz(0.92410134) q[2];
sx q[2];
rz(-1.8515203) q[2];
sx q[2];
rz(-1.4626056) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.24304388) q[1];
sx q[1];
rz(-2.1412256) q[1];
sx q[1];
rz(-2.3977086) q[1];
rz(-0.96989743) q[3];
sx q[3];
rz(-0.88774046) q[3];
sx q[3];
rz(2.0252463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.84844184) q[2];
sx q[2];
rz(-1.2520496) q[2];
sx q[2];
rz(1.1085917) q[2];
rz(1.1340002) q[3];
sx q[3];
rz(-1.0568551) q[3];
sx q[3];
rz(3.0602684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14652458) q[0];
sx q[0];
rz(-1.1256555) q[0];
sx q[0];
rz(0.31196892) q[0];
rz(-0.23090714) q[1];
sx q[1];
rz(-1.1038019) q[1];
sx q[1];
rz(1.1280967) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9474079) q[0];
sx q[0];
rz(-2.0387406) q[0];
sx q[0];
rz(-2.9293961) q[0];
rz(-pi) q[1];
rz(1.1351734) q[2];
sx q[2];
rz(-0.95016236) q[2];
sx q[2];
rz(2.1114608) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1027185) q[1];
sx q[1];
rz(-2.6800214) q[1];
sx q[1];
rz(-0.25074236) q[1];
rz(-1.1328636) q[3];
sx q[3];
rz(-2.2859462) q[3];
sx q[3];
rz(-1.6555427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.4565304) q[2];
sx q[2];
rz(-1.5038467) q[2];
sx q[2];
rz(0.064229639) q[2];
rz(-1.8188933) q[3];
sx q[3];
rz(-2.5048246) q[3];
sx q[3];
rz(-2.2548811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5927758) q[0];
sx q[0];
rz(-0.23574695) q[0];
sx q[0];
rz(-1.892426) q[0];
rz(-0.33117548) q[1];
sx q[1];
rz(-2.0024029) q[1];
sx q[1];
rz(1.1963199) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0977444) q[0];
sx q[0];
rz(-3.1412438) q[0];
sx q[0];
rz(2.1972138) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5761137) q[2];
sx q[2];
rz(-1.3233159) q[2];
sx q[2];
rz(2.4385045) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0374245) q[1];
sx q[1];
rz(-1.0598039) q[1];
sx q[1];
rz(0.68242208) q[1];
rz(-pi) q[2];
rz(-2.7306741) q[3];
sx q[3];
rz(-2.344729) q[3];
sx q[3];
rz(1.381402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5195878) q[2];
sx q[2];
rz(-0.76498166) q[2];
sx q[2];
rz(2.2288442) q[2];
rz(1.7031472) q[3];
sx q[3];
rz(-1.6310952) q[3];
sx q[3];
rz(-1.335817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4790633) q[0];
sx q[0];
rz(-1.1071858) q[0];
sx q[0];
rz(2.4712439) q[0];
rz(2.4743075) q[1];
sx q[1];
rz(-2.1332462) q[1];
sx q[1];
rz(-2.537312) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33907783) q[0];
sx q[0];
rz(-0.83503631) q[0];
sx q[0];
rz(1.5123532) q[0];
rz(-pi) q[1];
rz(-1.7261581) q[2];
sx q[2];
rz(-1.1929907) q[2];
sx q[2];
rz(0.31766674) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5126219) q[1];
sx q[1];
rz(-1.781946) q[1];
sx q[1];
rz(-1.3107915) q[1];
rz(1.8205582) q[3];
sx q[3];
rz(-2.3951924) q[3];
sx q[3];
rz(-1.6370893) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0303354) q[2];
sx q[2];
rz(-1.5735441) q[2];
sx q[2];
rz(2.7643909) q[2];
rz(3.0180569) q[3];
sx q[3];
rz(-1.433452) q[3];
sx q[3];
rz(1.7326573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94443026) q[0];
sx q[0];
rz(-2.5640709) q[0];
sx q[0];
rz(0.29931983) q[0];
rz(-1.1209283) q[1];
sx q[1];
rz(-1.9363554) q[1];
sx q[1];
rz(-2.1790806) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.045557307) q[0];
sx q[0];
rz(-1.0518414) q[0];
sx q[0];
rz(-0.23764289) q[0];
x q[1];
rz(0.91470529) q[2];
sx q[2];
rz(-1.1924628) q[2];
sx q[2];
rz(-1.3303743) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.968041) q[1];
sx q[1];
rz(-0.7725726) q[1];
sx q[1];
rz(2.3501758) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.706418) q[3];
sx q[3];
rz(-2.4101541) q[3];
sx q[3];
rz(2.4789435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0943429) q[2];
sx q[2];
rz(-2.3361358) q[2];
sx q[2];
rz(2.1979525) q[2];
rz(-2.0761679) q[3];
sx q[3];
rz(-1.3506972) q[3];
sx q[3];
rz(-1.0984727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0672673) q[0];
sx q[0];
rz(-1.704498) q[0];
sx q[0];
rz(-0.0083010439) q[0];
rz(0.51757327) q[1];
sx q[1];
rz(-0.65474302) q[1];
sx q[1];
rz(-1.5932721) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28866523) q[0];
sx q[0];
rz(-1.9216683) q[0];
sx q[0];
rz(0.025625833) q[0];
rz(-pi) q[1];
rz(2.6363434) q[2];
sx q[2];
rz(-1.8191686) q[2];
sx q[2];
rz(-0.58282436) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8928194) q[1];
sx q[1];
rz(-1.237251) q[1];
sx q[1];
rz(1.1897586) q[1];
rz(-pi) q[2];
rz(-2.1861784) q[3];
sx q[3];
rz(-1.9166167) q[3];
sx q[3];
rz(2.7503848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3090618) q[2];
sx q[2];
rz(-1.4893724) q[2];
sx q[2];
rz(-1.3365411) q[2];
rz(3.0858223) q[3];
sx q[3];
rz(-0.99168188) q[3];
sx q[3];
rz(0.43558863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5830773) q[0];
sx q[0];
rz(-2.6641088) q[0];
sx q[0];
rz(-2.4537295) q[0];
rz(-1.6784003) q[1];
sx q[1];
rz(-0.72931591) q[1];
sx q[1];
rz(0.24737839) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8198422) q[0];
sx q[0];
rz(-0.30748707) q[0];
sx q[0];
rz(0.64096682) q[0];
rz(0.94893564) q[2];
sx q[2];
rz(-2.1315711) q[2];
sx q[2];
rz(-1.6552629) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2108442) q[1];
sx q[1];
rz(-1.9860876) q[1];
sx q[1];
rz(-2.2077435) q[1];
rz(1.803627) q[3];
sx q[3];
rz(-1.4484753) q[3];
sx q[3];
rz(-1.8450774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4443724) q[2];
sx q[2];
rz(-1.8070544) q[2];
sx q[2];
rz(0.42416254) q[2];
rz(-1.0423202) q[3];
sx q[3];
rz(-0.76516953) q[3];
sx q[3];
rz(0.55646363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1384077) q[0];
sx q[0];
rz(-0.98588949) q[0];
sx q[0];
rz(-2.5307122) q[0];
rz(2.8914087) q[1];
sx q[1];
rz(-1.7411722) q[1];
sx q[1];
rz(2.6649323) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8308423) q[0];
sx q[0];
rz(-2.0526969) q[0];
sx q[0];
rz(-3.0974814) q[0];
rz(-2.0221316) q[2];
sx q[2];
rz(-0.20344606) q[2];
sx q[2];
rz(-0.97285336) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6338773) q[1];
sx q[1];
rz(-0.65836473) q[1];
sx q[1];
rz(-0.69067278) q[1];
rz(-pi) q[2];
rz(2.1564756) q[3];
sx q[3];
rz(-1.2063615) q[3];
sx q[3];
rz(-1.6062615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99159795) q[2];
sx q[2];
rz(-1.0341045) q[2];
sx q[2];
rz(2.4617713) q[2];
rz(2.6796135) q[3];
sx q[3];
rz(-1.446412) q[3];
sx q[3];
rz(-1.5010887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81944549) q[0];
sx q[0];
rz(-1.3752022) q[0];
sx q[0];
rz(0.15400259) q[0];
rz(-0.18140659) q[1];
sx q[1];
rz(-2.0712974) q[1];
sx q[1];
rz(3.0214686) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4171158) q[0];
sx q[0];
rz(-1.6301148) q[0];
sx q[0];
rz(-1.1671221) q[0];
rz(-2.9445705) q[2];
sx q[2];
rz(-1.8294334) q[2];
sx q[2];
rz(-1.7816003) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57615818) q[1];
sx q[1];
rz(-1.4729958) q[1];
sx q[1];
rz(-1.7478155) q[1];
rz(1.4265911) q[3];
sx q[3];
rz(-2.0940082) q[3];
sx q[3];
rz(-1.6169523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4274365) q[2];
sx q[2];
rz(-1.3849881) q[2];
sx q[2];
rz(-2.6825421) q[2];
rz(2.6311724) q[3];
sx q[3];
rz(-0.84158689) q[3];
sx q[3];
rz(2.4418805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089461483) q[0];
sx q[0];
rz(-0.94877807) q[0];
sx q[0];
rz(0.38598886) q[0];
rz(1.5486859) q[1];
sx q[1];
rz(-0.49647757) q[1];
sx q[1];
rz(-1.5492424) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0141409) q[0];
sx q[0];
rz(-1.7397428) q[0];
sx q[0];
rz(-0.028860748) q[0];
rz(2.4108294) q[2];
sx q[2];
rz(-2.363702) q[2];
sx q[2];
rz(2.23627) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.71536556) q[1];
sx q[1];
rz(-1.5976397) q[1];
sx q[1];
rz(-1.4687085) q[1];
rz(-pi) q[2];
rz(0.066927197) q[3];
sx q[3];
rz(-1.4176344) q[3];
sx q[3];
rz(2.7471971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7768895) q[2];
sx q[2];
rz(-2.2025755) q[2];
sx q[2];
rz(-1.4466064) q[2];
rz(-1.0363091) q[3];
sx q[3];
rz(-0.88007897) q[3];
sx q[3];
rz(3.115263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25630367) q[0];
sx q[0];
rz(-0.97923179) q[0];
sx q[0];
rz(-1.6543065) q[0];
rz(2.9215095) q[1];
sx q[1];
rz(-1.1047803) q[1];
sx q[1];
rz(1.6688375) q[1];
rz(1.8840811) q[2];
sx q[2];
rz(-0.85713119) q[2];
sx q[2];
rz(-1.5766889) q[2];
rz(2.6406399) q[3];
sx q[3];
rz(-1.8617478) q[3];
sx q[3];
rz(-0.72248722) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
