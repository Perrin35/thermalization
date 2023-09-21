OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2774529) q[0];
sx q[0];
rz(-1.5885408) q[0];
sx q[0];
rz(1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(-2.526386) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078755137) q[0];
sx q[0];
rz(-2.3427999) q[0];
sx q[0];
rz(0.90057217) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35742128) q[2];
sx q[2];
rz(-1.7454141) q[2];
sx q[2];
rz(-2.0703966) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.33949172) q[1];
sx q[1];
rz(-1.1420982) q[1];
sx q[1];
rz(-2.2080253) q[1];
rz(-1.390236) q[3];
sx q[3];
rz(-1.872634) q[3];
sx q[3];
rz(-2.3328822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7011828) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(0.33828503) q[2];
rz(-1.4398549) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(2.2556944) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171339) q[0];
sx q[0];
rz(-0.71115029) q[0];
sx q[0];
rz(3.1112444) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-2.1538484) q[1];
sx q[1];
rz(1.617584) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47027662) q[0];
sx q[0];
rz(-0.19965262) q[0];
sx q[0];
rz(-1.5623564) q[0];
rz(-1.6127869) q[2];
sx q[2];
rz(-2.6868372) q[2];
sx q[2];
rz(3.0595879) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1542201) q[1];
sx q[1];
rz(-1.0013354) q[1];
sx q[1];
rz(0.082055883) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.094035427) q[3];
sx q[3];
rz(-1.7574851) q[3];
sx q[3];
rz(1.7969014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.38561884) q[2];
sx q[2];
rz(-0.9884584) q[2];
sx q[2];
rz(-1.1478109) q[2];
rz(1.3267481) q[3];
sx q[3];
rz(-1.8170522) q[3];
sx q[3];
rz(0.23708788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1266992) q[0];
sx q[0];
rz(-0.48148695) q[0];
sx q[0];
rz(0.31578627) q[0];
rz(-0.93859998) q[1];
sx q[1];
rz(-1.6789852) q[1];
sx q[1];
rz(-0.25207239) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.789327) q[0];
sx q[0];
rz(-2.0354712) q[0];
sx q[0];
rz(-2.4436823) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3825233) q[2];
sx q[2];
rz(-1.0319064) q[2];
sx q[2];
rz(-0.76314229) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.1356126) q[1];
sx q[1];
rz(-0.95893919) q[1];
sx q[1];
rz(-0.23000418) q[1];
x q[2];
rz(1.9275946) q[3];
sx q[3];
rz(-1.1988415) q[3];
sx q[3];
rz(1.7372204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1217653) q[2];
sx q[2];
rz(-2.6941507) q[2];
sx q[2];
rz(-3.1075409) q[2];
rz(0.017459067) q[3];
sx q[3];
rz(-1.358946) q[3];
sx q[3];
rz(1.0954558) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62717342) q[0];
sx q[0];
rz(-1.592941) q[0];
sx q[0];
rz(1.5267641) q[0];
rz(2.0544255) q[1];
sx q[1];
rz(-2.4612869) q[1];
sx q[1];
rz(2.4345051) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2061545) q[0];
sx q[0];
rz(-1.1797138) q[0];
sx q[0];
rz(2.6576256) q[0];
rz(-pi) q[1];
rz(-0.25755067) q[2];
sx q[2];
rz(-0.45986816) q[2];
sx q[2];
rz(-2.3515153) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69869631) q[1];
sx q[1];
rz(-0.92123182) q[1];
sx q[1];
rz(2.991308) q[1];
rz(-pi) q[2];
x q[2];
rz(0.21541689) q[3];
sx q[3];
rz(-1.8724752) q[3];
sx q[3];
rz(-0.33723436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2924071) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(0.46009955) q[2];
rz(1.397331) q[3];
sx q[3];
rz(-1.5887235) q[3];
sx q[3];
rz(2.8989255) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8206772) q[0];
sx q[0];
rz(-2.9292332) q[0];
sx q[0];
rz(-1.7472349) q[0];
rz(1.0955411) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-0.25462338) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4655612) q[0];
sx q[0];
rz(-1.6499632) q[0];
sx q[0];
rz(0.013750793) q[0];
x q[1];
rz(-0.014572797) q[2];
sx q[2];
rz(-2.1483148) q[2];
sx q[2];
rz(-2.2968959) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8498807) q[1];
sx q[1];
rz(-0.92392081) q[1];
sx q[1];
rz(1.8748267) q[1];
rz(1.0765431) q[3];
sx q[3];
rz(-0.77270618) q[3];
sx q[3];
rz(-2.5172174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1191117) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.7648034) q[2];
rz(-1.4962176) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(1.0866603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6901533) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(2.8421463) q[0];
rz(-2.1014138) q[1];
sx q[1];
rz(-1.6957915) q[1];
sx q[1];
rz(-0.20656955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5305938) q[0];
sx q[0];
rz(-1.8078513) q[0];
sx q[0];
rz(-2.3321652) q[0];
rz(0.95025392) q[2];
sx q[2];
rz(-2.5197919) q[2];
sx q[2];
rz(1.7813462) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.2628852) q[1];
sx q[1];
rz(-0.72808121) q[1];
sx q[1];
rz(1.8379777) q[1];
x q[2];
rz(-1.1849095) q[3];
sx q[3];
rz(-2.4486802) q[3];
sx q[3];
rz(-3.0544359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2999337) q[2];
sx q[2];
rz(-1.4175697) q[2];
sx q[2];
rz(-0.28277961) q[2];
rz(0.81280604) q[3];
sx q[3];
rz(-0.41906425) q[3];
sx q[3];
rz(-0.16684428) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92641002) q[0];
sx q[0];
rz(-2.0880501) q[0];
sx q[0];
rz(-2.7600631) q[0];
rz(-0.58386699) q[1];
sx q[1];
rz(-2.5983512) q[1];
sx q[1];
rz(-1.8136224) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0921811) q[0];
sx q[0];
rz(-1.6850527) q[0];
sx q[0];
rz(-1.1294424) q[0];
rz(-0.96372693) q[2];
sx q[2];
rz(-0.71787314) q[2];
sx q[2];
rz(1.3340064) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2139637) q[1];
sx q[1];
rz(-0.8968401) q[1];
sx q[1];
rz(2.9272635) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9973203) q[3];
sx q[3];
rz(-1.0970308) q[3];
sx q[3];
rz(-1.9667448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.61838377) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(1.7810812) q[2];
rz(1.4303738) q[3];
sx q[3];
rz(-1.0083219) q[3];
sx q[3];
rz(-2.2935304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(2.9597136) q[0];
rz(-2.6673642) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(0.95091933) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5247105) q[0];
sx q[0];
rz(-1.8572154) q[0];
sx q[0];
rz(-2.8911203) q[0];
rz(-1.3046706) q[2];
sx q[2];
rz(-1.5106491) q[2];
sx q[2];
rz(1.6595449) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2434395) q[1];
sx q[1];
rz(-0.21874084) q[1];
sx q[1];
rz(2.8380413) q[1];
rz(-pi) q[2];
rz(2.3202053) q[3];
sx q[3];
rz(-2.6139724) q[3];
sx q[3];
rz(-1.1452831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2237079) q[2];
sx q[2];
rz(-0.48030883) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(-2.5935796) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6178745) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.3762208) q[0];
rz(-0.41704047) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(2.4818647) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2840246) q[0];
sx q[0];
rz(-2.095247) q[0];
sx q[0];
rz(0.88733034) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3599239) q[2];
sx q[2];
rz(-1.4147007) q[2];
sx q[2];
rz(-0.26542703) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0162504) q[1];
sx q[1];
rz(-2.2012735) q[1];
sx q[1];
rz(2.7793105) q[1];
rz(-pi) q[2];
rz(-3.0845853) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(1.9407335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.518121) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(0.14885151) q[3];
sx q[3];
rz(-1.0326577) q[3];
sx q[3];
rz(-3.1159475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898107) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(-1.8359258) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.8635609) q[1];
sx q[1];
rz(1.0356888) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0337692) q[0];
sx q[0];
rz(-2.6080837) q[0];
sx q[0];
rz(-2.4643722) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5752441) q[2];
sx q[2];
rz(-2.1841335) q[2];
sx q[2];
rz(-1.7042421) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4020821) q[1];
sx q[1];
rz(-1.2716736) q[1];
sx q[1];
rz(1.3189391) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8412748) q[3];
sx q[3];
rz(-2.3016553) q[3];
sx q[3];
rz(2.4887816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.62853652) q[2];
sx q[2];
rz(-2.3302902) q[2];
sx q[2];
rz(1.9899842) q[2];
rz(-0.51268762) q[3];
sx q[3];
rz(-1.0980462) q[3];
sx q[3];
rz(-2.776896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37968996) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.6336541) q[1];
sx q[1];
rz(-0.58273756) q[1];
sx q[1];
rz(2.6599463) q[1];
rz(0.41082906) q[2];
sx q[2];
rz(-1.6098235) q[2];
sx q[2];
rz(3.1237684) q[2];
rz(2.5232265) q[3];
sx q[3];
rz(-1.5636087) q[3];
sx q[3];
rz(-2.9611361) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];