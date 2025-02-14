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
rz(0.077724783) q[0];
sx q[0];
rz(-2.1169777) q[0];
sx q[0];
rz(-1.2999363) q[0];
rz(-0.44252244) q[1];
sx q[1];
rz(-1.1392925) q[1];
sx q[1];
rz(-2.7055969) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.090442) q[0];
sx q[0];
rz(-1.6761177) q[0];
sx q[0];
rz(-1.1770583) q[0];
x q[1];
rz(2.4019064) q[2];
sx q[2];
rz(-1.4130744) q[2];
sx q[2];
rz(-0.24631234) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.30937815) q[1];
sx q[1];
rz(-1.16799) q[1];
sx q[1];
rz(1.554053) q[1];
x q[2];
rz(-1.3920062) q[3];
sx q[3];
rz(-0.62084508) q[3];
sx q[3];
rz(1.5736962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.062833) q[2];
sx q[2];
rz(-1.6597513) q[2];
sx q[2];
rz(0.01595846) q[2];
rz(1.4912841) q[3];
sx q[3];
rz(-1.242638) q[3];
sx q[3];
rz(2.7119467) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4780739) q[0];
sx q[0];
rz(-1.0915382) q[0];
sx q[0];
rz(1.7953405) q[0];
rz(-1.1675872) q[1];
sx q[1];
rz(-2.0474032) q[1];
sx q[1];
rz(1.8928554) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8909709) q[0];
sx q[0];
rz(-0.75031322) q[0];
sx q[0];
rz(-2.1560099) q[0];
rz(-pi) q[1];
rz(-2.4820868) q[2];
sx q[2];
rz(-2.9444866) q[2];
sx q[2];
rz(1.8243621) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.68945049) q[1];
sx q[1];
rz(-0.4045338) q[1];
sx q[1];
rz(-1.4305315) q[1];
rz(-pi) q[2];
rz(-2.1571659) q[3];
sx q[3];
rz(-1.8674486) q[3];
sx q[3];
rz(0.082060952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.1355797) q[2];
sx q[2];
rz(-0.69301444) q[2];
sx q[2];
rz(0.47323027) q[2];
rz(3.0114975) q[3];
sx q[3];
rz(-1.3869163) q[3];
sx q[3];
rz(3.1307898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3267645) q[0];
sx q[0];
rz(-0.15693754) q[0];
sx q[0];
rz(-0.21406315) q[0];
rz(1.3940943) q[1];
sx q[1];
rz(-0.97624818) q[1];
sx q[1];
rz(-0.086437978) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8562216) q[0];
sx q[0];
rz(-2.8204794) q[0];
sx q[0];
rz(1.1620528) q[0];
x q[1];
rz(2.5912924) q[2];
sx q[2];
rz(-0.96947602) q[2];
sx q[2];
rz(2.3148906) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.4519388) q[1];
sx q[1];
rz(-1.9898333) q[1];
sx q[1];
rz(1.7614014) q[1];
x q[2];
rz(2.2203234) q[3];
sx q[3];
rz(-1.1315529) q[3];
sx q[3];
rz(-2.4961703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.67132407) q[2];
sx q[2];
rz(-2.2094122) q[2];
sx q[2];
rz(2.0051125) q[2];
rz(1.350435) q[3];
sx q[3];
rz(-1.330749) q[3];
sx q[3];
rz(-2.7854846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7478624) q[0];
sx q[0];
rz(-0.50734729) q[0];
sx q[0];
rz(0.90079975) q[0];
rz(-1.9081217) q[1];
sx q[1];
rz(-2.2869459) q[1];
sx q[1];
rz(2.5305117) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38310941) q[0];
sx q[0];
rz(-1.0953259) q[0];
sx q[0];
rz(-0.0080541797) q[0];
rz(-pi) q[1];
rz(2.99154) q[2];
sx q[2];
rz(-2.9100579) q[2];
sx q[2];
rz(-2.6611212) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.722695) q[1];
sx q[1];
rz(-1.9729821) q[1];
sx q[1];
rz(-2.6821892) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7902725) q[3];
sx q[3];
rz(-1.6804983) q[3];
sx q[3];
rz(2.9370995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.63371199) q[2];
sx q[2];
rz(-0.66358006) q[2];
sx q[2];
rz(-3.0207685) q[2];
rz(0.027033022) q[3];
sx q[3];
rz(-0.20172541) q[3];
sx q[3];
rz(2.2656608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9487069) q[0];
sx q[0];
rz(-0.8256194) q[0];
sx q[0];
rz(-1.8008308) q[0];
rz(-1.6138529) q[1];
sx q[1];
rz(-1.8283045) q[1];
sx q[1];
rz(-2.5921879) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2717239) q[0];
sx q[0];
rz(-2.3921674) q[0];
sx q[0];
rz(-0.11336993) q[0];
rz(2.1020736) q[2];
sx q[2];
rz(-1.5434859) q[2];
sx q[2];
rz(2.4940048) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.57826248) q[1];
sx q[1];
rz(-1.5208929) q[1];
sx q[1];
rz(-2.0567377) q[1];
rz(-pi) q[2];
rz(-3.1079783) q[3];
sx q[3];
rz(-1.1096802) q[3];
sx q[3];
rz(-0.5881084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1250829) q[2];
sx q[2];
rz(-1.5089401) q[2];
sx q[2];
rz(-0.58270085) q[2];
rz(1.5199644) q[3];
sx q[3];
rz(-0.71526066) q[3];
sx q[3];
rz(1.3165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.80578605) q[0];
sx q[0];
rz(-2.057322) q[0];
sx q[0];
rz(-1.2544607) q[0];
rz(-1.9460024) q[1];
sx q[1];
rz(-1.4072199) q[1];
sx q[1];
rz(1.5171299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82255581) q[0];
sx q[0];
rz(-1.5040888) q[0];
sx q[0];
rz(-3.1298679) q[0];
rz(-pi) q[1];
rz(-2.1468494) q[2];
sx q[2];
rz(-0.85604224) q[2];
sx q[2];
rz(-1.7896259) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.62134777) q[1];
sx q[1];
rz(-2.1572182) q[1];
sx q[1];
rz(-1.919073) q[1];
rz(-pi) q[2];
rz(-0.020217309) q[3];
sx q[3];
rz(-2.6264694) q[3];
sx q[3];
rz(-2.8216854) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.967041) q[2];
sx q[2];
rz(-1.6074564) q[2];
sx q[2];
rz(-0.045844585) q[2];
rz(0.52715078) q[3];
sx q[3];
rz(-1.0322626) q[3];
sx q[3];
rz(-2.3003858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884278) q[0];
sx q[0];
rz(-0.66075745) q[0];
sx q[0];
rz(2.7556457) q[0];
rz(1.376232) q[1];
sx q[1];
rz(-2.0538581) q[1];
sx q[1];
rz(-1.8650581) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71300426) q[0];
sx q[0];
rz(-0.78996336) q[0];
sx q[0];
rz(-2.6053564) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0936095) q[2];
sx q[2];
rz(-1.9429824) q[2];
sx q[2];
rz(-1.3822777) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7915262) q[1];
sx q[1];
rz(-0.036553144) q[1];
sx q[1];
rz(1.790148) q[1];
x q[2];
rz(-1.0625528) q[3];
sx q[3];
rz(-0.073598737) q[3];
sx q[3];
rz(1.2307597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4253) q[2];
sx q[2];
rz(-1.425681) q[2];
sx q[2];
rz(-1.7745793) q[2];
rz(-3.0573209) q[3];
sx q[3];
rz(-2.027498) q[3];
sx q[3];
rz(-0.082898609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0609584) q[0];
sx q[0];
rz(-3.0301889) q[0];
sx q[0];
rz(1.2782619) q[0];
rz(-3.0807965) q[1];
sx q[1];
rz(-2.1863329) q[1];
sx q[1];
rz(2.0416868) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3344582) q[0];
sx q[0];
rz(-1.4239199) q[0];
sx q[0];
rz(-1.9554041) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79675891) q[2];
sx q[2];
rz(-2.1646059) q[2];
sx q[2];
rz(1.6006921) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.56406139) q[1];
sx q[1];
rz(-1.2857591) q[1];
sx q[1];
rz(0.37071812) q[1];
rz(-pi) q[2];
rz(2.4722443) q[3];
sx q[3];
rz(-0.7318157) q[3];
sx q[3];
rz(1.2199618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3595769) q[2];
sx q[2];
rz(-2.4961553) q[2];
sx q[2];
rz(0.47425708) q[2];
rz(2.5931902) q[3];
sx q[3];
rz(-0.67265284) q[3];
sx q[3];
rz(1.7614346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7062374) q[0];
sx q[0];
rz(-0.42689231) q[0];
sx q[0];
rz(-1.2109582) q[0];
rz(-1.8636761) q[1];
sx q[1];
rz(-1.5958818) q[1];
sx q[1];
rz(-0.011215297) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8497711) q[0];
sx q[0];
rz(-1.0106083) q[0];
sx q[0];
rz(-1.673242) q[0];
rz(-pi) q[1];
rz(-2.5915543) q[2];
sx q[2];
rz(-0.42114741) q[2];
sx q[2];
rz(2.5530346) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.01811951) q[1];
sx q[1];
rz(-1.6247066) q[1];
sx q[1];
rz(0.82236313) q[1];
x q[2];
rz(1.2645006) q[3];
sx q[3];
rz(-0.3031177) q[3];
sx q[3];
rz(-2.1949286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8212905) q[2];
sx q[2];
rz(-1.9052637) q[2];
sx q[2];
rz(-2.383929) q[2];
rz(-2.441794) q[3];
sx q[3];
rz(-2.564513) q[3];
sx q[3];
rz(-0.9052161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.5504172) q[0];
sx q[0];
rz(-1.2400405) q[0];
sx q[0];
rz(0.33429876) q[0];
rz(-2.7475157) q[1];
sx q[1];
rz(-0.95029345) q[1];
sx q[1];
rz(1.2324415) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98823791) q[0];
sx q[0];
rz(-1.352509) q[0];
sx q[0];
rz(1.2984896) q[0];
rz(-2.903995) q[2];
sx q[2];
rz(-0.71212775) q[2];
sx q[2];
rz(-0.99436307) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5345911) q[1];
sx q[1];
rz(-1.8465733) q[1];
sx q[1];
rz(-0.2291344) q[1];
rz(-pi) q[2];
rz(2.0471935) q[3];
sx q[3];
rz(-0.55022424) q[3];
sx q[3];
rz(-2.8617895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1831827) q[2];
sx q[2];
rz(-2.5739058) q[2];
sx q[2];
rz(-3.0779823) q[2];
rz(-2.375864) q[3];
sx q[3];
rz(-2.016341) q[3];
sx q[3];
rz(-0.061802797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80339377) q[0];
sx q[0];
rz(-2.3173208) q[0];
sx q[0];
rz(1.4650387) q[0];
rz(-1.6784531) q[1];
sx q[1];
rz(-1.2668162) q[1];
sx q[1];
rz(-0.75513671) q[1];
rz(2.4974291) q[2];
sx q[2];
rz(-1.436787) q[2];
sx q[2];
rz(0.53874642) q[2];
rz(1.75453) q[3];
sx q[3];
rz(-2.4020772) q[3];
sx q[3];
rz(-2.6030131) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
