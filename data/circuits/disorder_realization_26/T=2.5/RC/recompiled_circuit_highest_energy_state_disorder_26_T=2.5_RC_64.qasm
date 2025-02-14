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
rz(-2.2505724) q[0];
sx q[0];
rz(-1.855259) q[0];
sx q[0];
rz(2.2556055) q[0];
rz(-2.6180144) q[1];
sx q[1];
rz(-0.55966592) q[1];
sx q[1];
rz(1.3203415) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26368982) q[0];
sx q[0];
rz(-0.56803507) q[0];
sx q[0];
rz(-1.9024416) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2314531) q[2];
sx q[2];
rz(-2.3187713) q[2];
sx q[2];
rz(2.0773599) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1883386) q[1];
sx q[1];
rz(-1.0085114) q[1];
sx q[1];
rz(0.072093318) q[1];
rz(2.3778649) q[3];
sx q[3];
rz(-0.9850174) q[3];
sx q[3];
rz(-1.8770742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.77871275) q[2];
sx q[2];
rz(-1.820463) q[2];
sx q[2];
rz(-2.2134181) q[2];
rz(1.5556395) q[3];
sx q[3];
rz(-1.0471683) q[3];
sx q[3];
rz(1.1699404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25598079) q[0];
sx q[0];
rz(-2.8035127) q[0];
sx q[0];
rz(2.0804491) q[0];
rz(-1.2127009) q[1];
sx q[1];
rz(-0.78293982) q[1];
sx q[1];
rz(1.3139668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3201225) q[0];
sx q[0];
rz(-2.263453) q[0];
sx q[0];
rz(-2.0351324) q[0];
rz(2.1232067) q[2];
sx q[2];
rz(-2.4168042) q[2];
sx q[2];
rz(1.9375436) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3281341) q[1];
sx q[1];
rz(-0.5603921) q[1];
sx q[1];
rz(-1.8916393) q[1];
rz(-2.6077929) q[3];
sx q[3];
rz(-0.91900037) q[3];
sx q[3];
rz(0.25215146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5706011) q[2];
sx q[2];
rz(-2.1676895) q[2];
sx q[2];
rz(2.2368597) q[2];
rz(1.5500801) q[3];
sx q[3];
rz(-1.855987) q[3];
sx q[3];
rz(1.2929644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.078481361) q[0];
sx q[0];
rz(-0.08992973) q[0];
sx q[0];
rz(-0.91823804) q[0];
rz(-0.69084424) q[1];
sx q[1];
rz(-1.3965239) q[1];
sx q[1];
rz(-2.3450559) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.349571) q[0];
sx q[0];
rz(-2.9137028) q[0];
sx q[0];
rz(-1.9877276) q[0];
x q[1];
rz(-1.5112707) q[2];
sx q[2];
rz(-1.5280485) q[2];
sx q[2];
rz(1.8190365) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.59350508) q[1];
sx q[1];
rz(-1.4626164) q[1];
sx q[1];
rz(0.87203474) q[1];
rz(-pi) q[2];
rz(2.2240766) q[3];
sx q[3];
rz(-1.2923354) q[3];
sx q[3];
rz(2.3115932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1306661) q[2];
sx q[2];
rz(-2.1911759) q[2];
sx q[2];
rz(1.2476904) q[2];
rz(-0.55825663) q[3];
sx q[3];
rz(-1.6853251) q[3];
sx q[3];
rz(-3.1202417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(2.438544) q[0];
sx q[0];
rz(-3.092364) q[0];
sx q[0];
rz(2.4141648) q[0];
rz(-0.1618596) q[1];
sx q[1];
rz(-1.9700123) q[1];
sx q[1];
rz(1.9836099) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3129565) q[0];
sx q[0];
rz(-1.6677756) q[0];
sx q[0];
rz(-3.0752379) q[0];
rz(-0.00096614758) q[2];
sx q[2];
rz(-1.4314993) q[2];
sx q[2];
rz(1.3749706) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.53044415) q[1];
sx q[1];
rz(-1.4475249) q[1];
sx q[1];
rz(-1.378343) q[1];
x q[2];
rz(1.653966) q[3];
sx q[3];
rz(-2.0953878) q[3];
sx q[3];
rz(-2.7980141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8689279) q[2];
sx q[2];
rz(-1.9395892) q[2];
sx q[2];
rz(1.7499917) q[2];
rz(-0.23022716) q[3];
sx q[3];
rz(-1.9672491) q[3];
sx q[3];
rz(-0.19190425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4153083) q[0];
sx q[0];
rz(-3.0400161) q[0];
sx q[0];
rz(-1.104785) q[0];
rz(0.11100189) q[1];
sx q[1];
rz(-1.1155201) q[1];
sx q[1];
rz(1.8288137) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5376624) q[0];
sx q[0];
rz(-0.27034187) q[0];
sx q[0];
rz(2.2946847) q[0];
rz(-pi) q[1];
rz(-0.13938015) q[2];
sx q[2];
rz(-2.5733893) q[2];
sx q[2];
rz(2.2903493) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5178419) q[1];
sx q[1];
rz(-1.1718796) q[1];
sx q[1];
rz(2.9241882) q[1];
rz(-pi) q[2];
rz(0.85185527) q[3];
sx q[3];
rz(-1.2491711) q[3];
sx q[3];
rz(2.3648398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.63867265) q[2];
sx q[2];
rz(-2.0827677) q[2];
sx q[2];
rz(2.1144833) q[2];
rz(0.98661679) q[3];
sx q[3];
rz(-2.8592181) q[3];
sx q[3];
rz(-1.7178887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33243772) q[0];
sx q[0];
rz(-1.9092535) q[0];
sx q[0];
rz(0.80068457) q[0];
rz(-2.370131) q[1];
sx q[1];
rz(-1.0779251) q[1];
sx q[1];
rz(-1.7652184) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3478394) q[0];
sx q[0];
rz(-0.79393605) q[0];
sx q[0];
rz(2.5695922) q[0];
rz(-0.034997392) q[2];
sx q[2];
rz(-2.0366324) q[2];
sx q[2];
rz(0.83110561) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4478057) q[1];
sx q[1];
rz(-2.4359833) q[1];
sx q[1];
rz(-0.020486995) q[1];
rz(-pi) q[2];
x q[2];
rz(0.82628754) q[3];
sx q[3];
rz(-1.0184231) q[3];
sx q[3];
rz(-0.9780405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.8419522) q[2];
sx q[2];
rz(-1.2066634) q[2];
sx q[2];
rz(0.017814962) q[2];
rz(2.8282015) q[3];
sx q[3];
rz(-1.0217977) q[3];
sx q[3];
rz(2.4221086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7844836) q[0];
sx q[0];
rz(-1.5743558) q[0];
sx q[0];
rz(-0.16726476) q[0];
rz(-2.3978865) q[1];
sx q[1];
rz(-0.88631648) q[1];
sx q[1];
rz(-3.0053265) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8945904) q[0];
sx q[0];
rz(-0.73406363) q[0];
sx q[0];
rz(1.1511001) q[0];
x q[1];
rz(-2.9685814) q[2];
sx q[2];
rz(-1.271651) q[2];
sx q[2];
rz(-1.1227705) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4140461) q[1];
sx q[1];
rz(-2.7767477) q[1];
sx q[1];
rz(2.9529497) q[1];
x q[2];
rz(1.5149917) q[3];
sx q[3];
rz(-1.2254834) q[3];
sx q[3];
rz(1.2962411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3027489) q[2];
sx q[2];
rz(-0.15041298) q[2];
sx q[2];
rz(2.0602843) q[2];
rz(2.9980764) q[3];
sx q[3];
rz(-1.2986526) q[3];
sx q[3];
rz(-1.6941841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3926587) q[0];
sx q[0];
rz(-1.3329788) q[0];
sx q[0];
rz(2.4984711) q[0];
rz(2.2195623) q[1];
sx q[1];
rz(-0.26607251) q[1];
sx q[1];
rz(1.5111074) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5095172) q[0];
sx q[0];
rz(-0.38783535) q[0];
sx q[0];
rz(-2.0616421) q[0];
rz(-pi) q[1];
rz(-3.0273075) q[2];
sx q[2];
rz(-1.292789) q[2];
sx q[2];
rz(-0.68113995) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.585941) q[1];
sx q[1];
rz(-1.562848) q[1];
sx q[1];
rz(1.7674602) q[1];
rz(-0.53427215) q[3];
sx q[3];
rz(-2.5639682) q[3];
sx q[3];
rz(2.3848611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1954605) q[2];
sx q[2];
rz(-1.5798502) q[2];
sx q[2];
rz(1.3004318) q[2];
rz(-0.91222936) q[3];
sx q[3];
rz(-1.2765086) q[3];
sx q[3];
rz(-0.31360489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6533971) q[0];
sx q[0];
rz(-2.6618239) q[0];
sx q[0];
rz(-2.524014) q[0];
rz(0.081534475) q[1];
sx q[1];
rz(-0.50059861) q[1];
sx q[1];
rz(0.68731442) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1253374) q[0];
sx q[0];
rz(-2.9000686) q[0];
sx q[0];
rz(-1.3708273) q[0];
rz(-pi) q[1];
rz(-0.388889) q[2];
sx q[2];
rz(-2.8538879) q[2];
sx q[2];
rz(0.012118169) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1256492) q[1];
sx q[1];
rz(-1.9852891) q[1];
sx q[1];
rz(1.8806367) q[1];
rz(3.0932759) q[3];
sx q[3];
rz(-2.1135277) q[3];
sx q[3];
rz(-0.84191546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11757892) q[2];
sx q[2];
rz(-1.1774096) q[2];
sx q[2];
rz(-3.0702316) q[2];
rz(1.9060382) q[3];
sx q[3];
rz(-2.8046298) q[3];
sx q[3];
rz(-2.5753042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5736893) q[0];
sx q[0];
rz(-0.25147831) q[0];
sx q[0];
rz(-3.0036744) q[0];
rz(0.57016405) q[1];
sx q[1];
rz(-1.0823931) q[1];
sx q[1];
rz(-3.1261442) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1874267) q[0];
sx q[0];
rz(-1.2492234) q[0];
sx q[0];
rz(-0.84521742) q[0];
rz(-pi) q[1];
rz(-0.76064588) q[2];
sx q[2];
rz(-1.9234675) q[2];
sx q[2];
rz(1.2404397) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9271157) q[1];
sx q[1];
rz(-1.0815911) q[1];
sx q[1];
rz(1.3838883) q[1];
rz(-pi) q[2];
rz(-3.0857419) q[3];
sx q[3];
rz(-1.7540364) q[3];
sx q[3];
rz(0.88411546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40921679) q[2];
sx q[2];
rz(-0.78745431) q[2];
sx q[2];
rz(0.43238762) q[2];
rz(2.2435097) q[3];
sx q[3];
rz(-0.87528527) q[3];
sx q[3];
rz(-1.557365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4846004) q[0];
sx q[0];
rz(-1.1285755) q[0];
sx q[0];
rz(-2.3791671) q[0];
rz(0.55105974) q[1];
sx q[1];
rz(-1.8012128) q[1];
sx q[1];
rz(2.6691379) q[1];
rz(3.0226784) q[2];
sx q[2];
rz(-2.7758895) q[2];
sx q[2];
rz(2.9056673) q[2];
rz(2.2021709) q[3];
sx q[3];
rz(-0.71790725) q[3];
sx q[3];
rz(0.090286615) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
