OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.6686749) q[0];
sx q[0];
rz(-0.023107419) q[0];
sx q[0];
rz(0.90149108) q[0];
rz(1.3540406) q[1];
sx q[1];
rz(-1.5259589) q[1];
sx q[1];
rz(1.807133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62518277) q[0];
sx q[0];
rz(-1.3526655) q[0];
sx q[0];
rz(1.6263917) q[0];
rz(0.35144866) q[2];
sx q[2];
rz(-0.16021591) q[2];
sx q[2];
rz(-1.7376228) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6264947) q[1];
sx q[1];
rz(-0.52649311) q[1];
sx q[1];
rz(-1.4655617) q[1];
rz(-pi) q[2];
rz(2.9525312) q[3];
sx q[3];
rz(-2.8950518) q[3];
sx q[3];
rz(-2.2061359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2287075) q[2];
sx q[2];
rz(-0.097147377) q[2];
sx q[2];
rz(2.482282) q[2];
rz(-0.77624503) q[3];
sx q[3];
rz(-0.018748911) q[3];
sx q[3];
rz(-2.4565878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9126251) q[0];
sx q[0];
rz(-1.2094867) q[0];
sx q[0];
rz(-1.1932766) q[0];
rz(3.0972262) q[1];
sx q[1];
rz(-0.012244789) q[1];
sx q[1];
rz(-0.22656974) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3069585) q[0];
sx q[0];
rz(-3.1378085) q[0];
sx q[0];
rz(1.1192805) q[0];
rz(-pi) q[1];
rz(0.37906693) q[2];
sx q[2];
rz(-1.5549391) q[2];
sx q[2];
rz(-3.130558) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6084131) q[1];
sx q[1];
rz(-0.4371818) q[1];
sx q[1];
rz(2.9730148) q[1];
x q[2];
rz(0.84323378) q[3];
sx q[3];
rz(-1.7403462) q[3];
sx q[3];
rz(-0.42593004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6988397) q[2];
sx q[2];
rz(-1.5719465) q[2];
sx q[2];
rz(-1.6119831) q[2];
rz(0.93305856) q[3];
sx q[3];
rz(-1.6589087) q[3];
sx q[3];
rz(-2.9034767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9642692) q[0];
sx q[0];
rz(-0.015559109) q[0];
sx q[0];
rz(-0.200287) q[0];
rz(0.00042032584) q[1];
sx q[1];
rz(-2.2047408) q[1];
sx q[1];
rz(-0.013484152) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6157841) q[0];
sx q[0];
rz(-1.1466525) q[0];
sx q[0];
rz(-0.070930907) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5277943) q[2];
sx q[2];
rz(-1.5041877) q[2];
sx q[2];
rz(3.1284077) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.36253769) q[1];
sx q[1];
rz(-1.4351621) q[1];
sx q[1];
rz(-0.070222994) q[1];
x q[2];
rz(2.8369569) q[3];
sx q[3];
rz(-2.9713062) q[3];
sx q[3];
rz(-2.0326322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9424332) q[2];
sx q[2];
rz(-1.556267) q[2];
sx q[2];
rz(1.5419434) q[2];
rz(-1.0017627) q[3];
sx q[3];
rz(-0.21654138) q[3];
sx q[3];
rz(-2.969363) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7672985) q[0];
sx q[0];
rz(-2.9758487) q[0];
sx q[0];
rz(0.34749183) q[0];
rz(-0.63684288) q[1];
sx q[1];
rz(-0.0056191365) q[1];
sx q[1];
rz(1.1905131) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40909262) q[0];
sx q[0];
rz(-2.996759) q[0];
sx q[0];
rz(1.1454789) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5019102) q[2];
sx q[2];
rz(-1.578555) q[2];
sx q[2];
rz(-3.1338281) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3733565) q[1];
sx q[1];
rz(-2.3563085) q[1];
sx q[1];
rz(0.26765243) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6116294) q[3];
sx q[3];
rz(-0.41614446) q[3];
sx q[3];
rz(0.79485369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5845329) q[2];
sx q[2];
rz(-3.1019042) q[2];
sx q[2];
rz(1.7812799) q[2];
rz(1.4761866) q[3];
sx q[3];
rz(-1.5675661) q[3];
sx q[3];
rz(0.40569693) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0880459) q[0];
sx q[0];
rz(-2.4021554) q[0];
sx q[0];
rz(0.0060225688) q[0];
rz(-1.4206403) q[1];
sx q[1];
rz(-0.050844897) q[1];
sx q[1];
rz(-0.086070148) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47274855) q[0];
sx q[0];
rz(-1.5875582) q[0];
sx q[0];
rz(-3.0344323) q[0];
rz(-pi) q[1];
rz(-0.44186307) q[2];
sx q[2];
rz(-0.060906868) q[2];
sx q[2];
rz(1.4373844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9169315) q[1];
sx q[1];
rz(-1.5337394) q[1];
sx q[1];
rz(0.091450973) q[1];
rz(1.6313598) q[3];
sx q[3];
rz(-1.525536) q[3];
sx q[3];
rz(-2.0082983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.033279557) q[2];
sx q[2];
rz(-2.4187708) q[2];
sx q[2];
rz(-1.3839728) q[2];
rz(-2.7464187) q[3];
sx q[3];
rz(-3.0925909) q[3];
sx q[3];
rz(1.1672195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5351717) q[0];
sx q[0];
rz(-2.9223154) q[0];
sx q[0];
rz(2.1411335) q[0];
rz(-0.74042997) q[1];
sx q[1];
rz(-2.705997) q[1];
sx q[1];
rz(-2.7620517) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1943037) q[0];
sx q[0];
rz(-1.6662046) q[0];
sx q[0];
rz(1.545115) q[0];
rz(1.2984811) q[2];
sx q[2];
rz(-2.0305995) q[2];
sx q[2];
rz(-1.0155755) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4011456) q[1];
sx q[1];
rz(-1.5441431) q[1];
sx q[1];
rz(-1.0668287) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3336181) q[3];
sx q[3];
rz(-2.6814125) q[3];
sx q[3];
rz(-1.7250012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1257816) q[2];
sx q[2];
rz(-2.8534079) q[2];
sx q[2];
rz(0.077032653) q[2];
rz(-0.020126255) q[3];
sx q[3];
rz(-0.052611668) q[3];
sx q[3];
rz(0.7974112) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1064442) q[0];
sx q[0];
rz(-3.1290717) q[0];
sx q[0];
rz(1.7092108) q[0];
rz(-2.8445981) q[1];
sx q[1];
rz(-0.15428267) q[1];
sx q[1];
rz(-3.0800379) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8650353) q[0];
sx q[0];
rz(-1.1662959) q[0];
sx q[0];
rz(0.96488953) q[0];
rz(-pi) q[1];
rz(1.5764357) q[2];
sx q[2];
rz(-1.3562849) q[2];
sx q[2];
rz(2.0129865) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.48744943) q[1];
sx q[1];
rz(-1.8620849) q[1];
sx q[1];
rz(-1.6421516) q[1];
rz(3.0606817) q[3];
sx q[3];
rz(-0.43596163) q[3];
sx q[3];
rz(-1.6502768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84294549) q[2];
sx q[2];
rz(-2.8892398) q[2];
sx q[2];
rz(1.7227777) q[2];
rz(-1.8055387) q[3];
sx q[3];
rz(-0.027848363) q[3];
sx q[3];
rz(-1.4068039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0272738) q[0];
sx q[0];
rz(-0.71684664) q[0];
sx q[0];
rz(-1.5007716) q[0];
rz(-1.1227135) q[1];
sx q[1];
rz(-0.32674679) q[1];
sx q[1];
rz(1.8480802) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6155535) q[0];
sx q[0];
rz(-2.3991971) q[0];
sx q[0];
rz(-2.0406988) q[0];
x q[1];
rz(0.40385623) q[2];
sx q[2];
rz(-0.73550288) q[2];
sx q[2];
rz(2.1063358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.499524) q[1];
sx q[1];
rz(-1.4022938) q[1];
sx q[1];
rz(2.8649855) q[1];
rz(-0.82610749) q[3];
sx q[3];
rz(-2.2650121) q[3];
sx q[3];
rz(0.40761596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5571755) q[2];
sx q[2];
rz(-0.36089218) q[2];
sx q[2];
rz(1.3593675) q[2];
rz(-0.44698295) q[3];
sx q[3];
rz(-3.0990661) q[3];
sx q[3];
rz(0.16714787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8166703) q[0];
sx q[0];
rz(-1.336038) q[0];
sx q[0];
rz(-2.0409806) q[0];
rz(-2.0657516) q[1];
sx q[1];
rz(-2.4918719) q[1];
sx q[1];
rz(-2.4646387) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9788584) q[0];
sx q[0];
rz(-1.7679201) q[0];
sx q[0];
rz(0.39359025) q[0];
rz(2.3271884) q[2];
sx q[2];
rz(-2.9737817) q[2];
sx q[2];
rz(1.4175159) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9857603) q[1];
sx q[1];
rz(-1.5712156) q[1];
sx q[1];
rz(9.930519e-05) q[1];
x q[2];
rz(1.7221801) q[3];
sx q[3];
rz(-1.1202235) q[3];
sx q[3];
rz(2.7344658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0000275) q[2];
sx q[2];
rz(-0.0024777369) q[2];
sx q[2];
rz(0.72762093) q[2];
rz(1.242312) q[3];
sx q[3];
rz(-0.036402313) q[3];
sx q[3];
rz(-1.9696994) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90945554) q[0];
sx q[0];
rz(-1.0219034) q[0];
sx q[0];
rz(-0.68956462) q[0];
rz(1.6652971) q[1];
sx q[1];
rz(-2.8910525) q[1];
sx q[1];
rz(-2.9857059) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7321701) q[0];
sx q[0];
rz(-2.4722833) q[0];
sx q[0];
rz(0.21440345) q[0];
rz(-0.25145289) q[2];
sx q[2];
rz(-1.6864459) q[2];
sx q[2];
rz(-2.419099) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5248651) q[1];
sx q[1];
rz(-1.5749802) q[1];
sx q[1];
rz(-1.5692488) q[1];
rz(3.0468349) q[3];
sx q[3];
rz(-1.4156431) q[3];
sx q[3];
rz(-1.9431537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9547687) q[2];
sx q[2];
rz(-0.12207741) q[2];
sx q[2];
rz(2.1464777) q[2];
rz(-0.22594813) q[3];
sx q[3];
rz(-3.093231) q[3];
sx q[3];
rz(2.3154955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7786998) q[0];
sx q[0];
rz(-2.1669372) q[0];
sx q[0];
rz(-1.7397276) q[0];
rz(-1.7097991) q[1];
sx q[1];
rz(-1.2964389) q[1];
sx q[1];
rz(-2.5242205) q[1];
rz(-1.7079034) q[2];
sx q[2];
rz(-2.2513486) q[2];
sx q[2];
rz(0.90604102) q[2];
rz(1.2479242) q[3];
sx q[3];
rz(-0.40798305) q[3];
sx q[3];
rz(1.8508607) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
