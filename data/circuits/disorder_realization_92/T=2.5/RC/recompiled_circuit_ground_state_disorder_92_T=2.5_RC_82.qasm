OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2607245) q[0];
sx q[0];
rz(-0.27624929) q[0];
sx q[0];
rz(0.84375381) q[0];
rz(1.002797) q[1];
sx q[1];
rz(3.9763713) q[1];
sx q[1];
rz(3.6741771) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7363889) q[0];
sx q[0];
rz(-1.4215934) q[0];
sx q[0];
rz(2.2264541) q[0];
x q[1];
rz(-1.7108828) q[2];
sx q[2];
rz(-1.0647961) q[2];
sx q[2];
rz(1.9078474) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3644173) q[1];
sx q[1];
rz(-0.88349408) q[1];
sx q[1];
rz(-0.19608325) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37181969) q[3];
sx q[3];
rz(-1.216868) q[3];
sx q[3];
rz(1.33385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8325995) q[2];
sx q[2];
rz(-1.4048046) q[2];
sx q[2];
rz(0.12283202) q[2];
rz(0.041736688) q[3];
sx q[3];
rz(-1.5580956) q[3];
sx q[3];
rz(2.6532555) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1500583) q[0];
sx q[0];
rz(-2.8128862) q[0];
sx q[0];
rz(1.1154037) q[0];
rz(-2.5967122) q[1];
sx q[1];
rz(-1.3976401) q[1];
sx q[1];
rz(-0.56745183) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.713288) q[0];
sx q[0];
rz(-0.056763857) q[0];
sx q[0];
rz(-0.63637498) q[0];
rz(-pi) q[1];
rz(2.3653602) q[2];
sx q[2];
rz(-1.5934332) q[2];
sx q[2];
rz(0.55720383) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.27047711) q[1];
sx q[1];
rz(-2.1391627) q[1];
sx q[1];
rz(-1.4517484) q[1];
rz(1.2699158) q[3];
sx q[3];
rz(-2.0344276) q[3];
sx q[3];
rz(-2.5167572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51022092) q[2];
sx q[2];
rz(-0.12237445) q[2];
sx q[2];
rz(-0.1046293) q[2];
rz(-0.50466022) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(-2.4095355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1402682) q[0];
sx q[0];
rz(-2.9502385) q[0];
sx q[0];
rz(1.485317) q[0];
rz(-0.56400076) q[1];
sx q[1];
rz(-2.0878849) q[1];
sx q[1];
rz(-0.82411134) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7566273) q[0];
sx q[0];
rz(-0.55776309) q[0];
sx q[0];
rz(1.2560448) q[0];
rz(2.3468909) q[2];
sx q[2];
rz(-1.1247171) q[2];
sx q[2];
rz(0.54903713) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50269404) q[1];
sx q[1];
rz(-0.47595176) q[1];
sx q[1];
rz(2.6902728) q[1];
rz(-3.0715354) q[3];
sx q[3];
rz(-1.7825216) q[3];
sx q[3];
rz(0.61933653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.092502681) q[2];
sx q[2];
rz(-2.8503032) q[2];
sx q[2];
rz(-1.5352486) q[2];
rz(-0.90318471) q[3];
sx q[3];
rz(-1.3245557) q[3];
sx q[3];
rz(0.42455348) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20151888) q[0];
sx q[0];
rz(-2.1551082) q[0];
sx q[0];
rz(-2.7281813) q[0];
rz(-2.9460733) q[1];
sx q[1];
rz(-2.1045411) q[1];
sx q[1];
rz(-1.4541385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91823309) q[0];
sx q[0];
rz(-0.70230316) q[0];
sx q[0];
rz(1.7683517) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6624751) q[2];
sx q[2];
rz(-2.7450941) q[2];
sx q[2];
rz(-2.2355607) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.11024347) q[1];
sx q[1];
rz(-1.8520466) q[1];
sx q[1];
rz(-1.0737277) q[1];
x q[2];
rz(-0.24202427) q[3];
sx q[3];
rz(-0.70902642) q[3];
sx q[3];
rz(-2.5828862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2426408) q[2];
sx q[2];
rz(-0.61598888) q[2];
sx q[2];
rz(-1.8757437) q[2];
rz(-2.2856581) q[3];
sx q[3];
rz(-1.7549691) q[3];
sx q[3];
rz(-1.9267193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6619381) q[0];
sx q[0];
rz(-2.9582773) q[0];
sx q[0];
rz(1.7228175) q[0];
rz(-1.7105557) q[1];
sx q[1];
rz(-1.6173671) q[1];
sx q[1];
rz(2.4180744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2884525) q[0];
sx q[0];
rz(-0.61091629) q[0];
sx q[0];
rz(-2.8302961) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3143598) q[2];
sx q[2];
rz(-2.1783354) q[2];
sx q[2];
rz(-0.22996685) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12234593) q[1];
sx q[1];
rz(-0.86402383) q[1];
sx q[1];
rz(3.138423) q[1];
x q[2];
rz(0.35227065) q[3];
sx q[3];
rz(-1.5639502) q[3];
sx q[3];
rz(2.7628243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7001223) q[2];
sx q[2];
rz(-2.1258326) q[2];
sx q[2];
rz(-2.7777242) q[2];
rz(1.3077334) q[3];
sx q[3];
rz(-1.4738844) q[3];
sx q[3];
rz(-1.3522805) q[3];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42809197) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(0.8999024) q[0];
rz(2.0948386) q[1];
sx q[1];
rz(-0.37938198) q[1];
sx q[1];
rz(2.5894763) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22495843) q[0];
sx q[0];
rz(-0.33946589) q[0];
sx q[0];
rz(-0.084966226) q[0];
x q[1];
rz(-2.4399567) q[2];
sx q[2];
rz(-1.1404827) q[2];
sx q[2];
rz(-0.21077354) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0830517) q[1];
sx q[1];
rz(-1.6612435) q[1];
sx q[1];
rz(-1.7427216) q[1];
rz(-pi) q[2];
rz(-1.6523727) q[3];
sx q[3];
rz(-2.1595567) q[3];
sx q[3];
rz(-2.0643864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3209352) q[2];
sx q[2];
rz(-2.4652017) q[2];
sx q[2];
rz(-0.19860849) q[2];
rz(-2.0123539) q[3];
sx q[3];
rz(-1.6695453) q[3];
sx q[3];
rz(0.78013295) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8010913) q[0];
sx q[0];
rz(-1.8393562) q[0];
sx q[0];
rz(0.71402016) q[0];
rz(-1.9188312) q[1];
sx q[1];
rz(-2.265265) q[1];
sx q[1];
rz(0.27539918) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80197734) q[0];
sx q[0];
rz(-2.6875949) q[0];
sx q[0];
rz(-1.071644) q[0];
x q[1];
rz(-1.1366264) q[2];
sx q[2];
rz(-0.6587907) q[2];
sx q[2];
rz(-0.24497797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8620257) q[1];
sx q[1];
rz(-2.0525524) q[1];
sx q[1];
rz(2.606009) q[1];
x q[2];
rz(-0.72526284) q[3];
sx q[3];
rz(-1.5607087) q[3];
sx q[3];
rz(-0.46848224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3068646) q[2];
sx q[2];
rz(-1.4707668) q[2];
sx q[2];
rz(-1.6738221) q[2];
rz(1.62014) q[3];
sx q[3];
rz(-2.455267) q[3];
sx q[3];
rz(-0.93792382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9891605) q[0];
sx q[0];
rz(-2.2280362) q[0];
sx q[0];
rz(-2.1754225) q[0];
rz(0.78978157) q[1];
sx q[1];
rz(-2.8826394) q[1];
sx q[1];
rz(2.1678179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51289979) q[0];
sx q[0];
rz(-1.3256462) q[0];
sx q[0];
rz(-0.12314491) q[0];
rz(-1.9386693) q[2];
sx q[2];
rz(-1.972636) q[2];
sx q[2];
rz(1.3746868) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9653496) q[1];
sx q[1];
rz(-2.4479657) q[1];
sx q[1];
rz(2.0575869) q[1];
x q[2];
rz(-1.4548746) q[3];
sx q[3];
rz(-0.478607) q[3];
sx q[3];
rz(2.6461864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3070273) q[2];
sx q[2];
rz(-1.4296738) q[2];
sx q[2];
rz(-2.5359421) q[2];
rz(-2.0511131) q[3];
sx q[3];
rz(-2.8993789) q[3];
sx q[3];
rz(-3.0428913) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7379446) q[0];
sx q[0];
rz(-0.048796766) q[0];
sx q[0];
rz(-2.5700289) q[0];
rz(-0.38772186) q[1];
sx q[1];
rz(-2.3523836) q[1];
sx q[1];
rz(0.8943843) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8799202) q[0];
sx q[0];
rz(-2.399647) q[0];
sx q[0];
rz(-0.80961734) q[0];
rz(-pi) q[1];
rz(-0.11668514) q[2];
sx q[2];
rz(-1.0067847) q[2];
sx q[2];
rz(-0.91911585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.91796658) q[1];
sx q[1];
rz(-1.3296179) q[1];
sx q[1];
rz(1.9832703) q[1];
rz(-pi) q[2];
rz(-2.6329805) q[3];
sx q[3];
rz(-1.5023352) q[3];
sx q[3];
rz(0.48229846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9992708) q[2];
sx q[2];
rz(-0.88185507) q[2];
sx q[2];
rz(-0.87232653) q[2];
rz(3.0669323) q[3];
sx q[3];
rz(-1.9118237) q[3];
sx q[3];
rz(-2.5653896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9999303) q[0];
sx q[0];
rz(-2.7003728) q[0];
sx q[0];
rz(-0.73750752) q[0];
rz(-2.4420338) q[1];
sx q[1];
rz(-2.12205) q[1];
sx q[1];
rz(0.84651822) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93447157) q[0];
sx q[0];
rz(-1.9686598) q[0];
sx q[0];
rz(2.4475696) q[0];
rz(-pi) q[1];
x q[1];
rz(0.051324322) q[2];
sx q[2];
rz(-1.1451756) q[2];
sx q[2];
rz(1.0251306) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5246995) q[1];
sx q[1];
rz(-1.1438826) q[1];
sx q[1];
rz(-2.9906143) q[1];
x q[2];
rz(-0.41937866) q[3];
sx q[3];
rz(-2.7079795) q[3];
sx q[3];
rz(0.24672844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6503341) q[2];
sx q[2];
rz(-1.1639872) q[2];
sx q[2];
rz(-2.7872046) q[2];
rz(-0.77704159) q[3];
sx q[3];
rz(-0.6207501) q[3];
sx q[3];
rz(2.0508155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95018321) q[0];
sx q[0];
rz(-1.8671028) q[0];
sx q[0];
rz(-2.6381459) q[0];
rz(1.3632111) q[1];
sx q[1];
rz(-2.6111205) q[1];
sx q[1];
rz(3.0725239) q[1];
rz(-1.4466046) q[2];
sx q[2];
rz(-2.4249981) q[2];
sx q[2];
rz(-0.80719765) q[2];
rz(2.1997785) q[3];
sx q[3];
rz(-0.54472803) q[3];
sx q[3];
rz(1.7036078) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
