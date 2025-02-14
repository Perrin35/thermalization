OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7140952) q[0];
sx q[0];
rz(-2.5795955) q[0];
sx q[0];
rz(-0.23101097) q[0];
rz(-2.8958939) q[1];
sx q[1];
rz(-2.6872771) q[1];
sx q[1];
rz(1.8543724) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3481349) q[0];
sx q[0];
rz(-1.0890111) q[0];
sx q[0];
rz(2.8346377) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76164772) q[2];
sx q[2];
rz(-0.33312329) q[2];
sx q[2];
rz(-1.8322577) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.27632984) q[1];
sx q[1];
rz(-1.6693322) q[1];
sx q[1];
rz(-0.16698412) q[1];
rz(-pi) q[2];
rz(1.9735967) q[3];
sx q[3];
rz(-2.9514545) q[3];
sx q[3];
rz(-0.26241747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.40453688) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(-2.7810968) q[2];
rz(1.1245842) q[3];
sx q[3];
rz(-1.0485317) q[3];
sx q[3];
rz(1.2688961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(0.26038134) q[0];
sx q[0];
rz(-0.17558782) q[0];
sx q[0];
rz(-2.4847109) q[0];
rz(-2.8086713) q[1];
sx q[1];
rz(-2.0491144) q[1];
sx q[1];
rz(-2.2064256) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012197709) q[0];
sx q[0];
rz(-3.0327065) q[0];
sx q[0];
rz(-2.8588141) q[0];
rz(-pi) q[1];
rz(-2.8043973) q[2];
sx q[2];
rz(-2.8697578) q[2];
sx q[2];
rz(-1.4120917) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1166653) q[1];
sx q[1];
rz(-1.5313193) q[1];
sx q[1];
rz(-0.14703207) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6245429) q[3];
sx q[3];
rz(-0.40362261) q[3];
sx q[3];
rz(1.0583741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.3071345) q[2];
sx q[2];
rz(-1.5156526) q[2];
sx q[2];
rz(-1.9251941) q[2];
rz(-0.52792102) q[3];
sx q[3];
rz(-2.1025751) q[3];
sx q[3];
rz(1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5288178) q[0];
sx q[0];
rz(-0.82614326) q[0];
sx q[0];
rz(-2.5566027) q[0];
rz(3.0168369) q[1];
sx q[1];
rz(-2.545732) q[1];
sx q[1];
rz(1.6927208) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0199658) q[0];
sx q[0];
rz(-0.0040071132) q[0];
sx q[0];
rz(-0.83929707) q[0];
rz(-pi) q[1];
rz(0.18529202) q[2];
sx q[2];
rz(-3.0681562) q[2];
sx q[2];
rz(0.18723182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3347149) q[1];
sx q[1];
rz(-2.2141465) q[1];
sx q[1];
rz(-0.54912864) q[1];
x q[2];
rz(-0.65421974) q[3];
sx q[3];
rz(-0.4163792) q[3];
sx q[3];
rz(-0.49921303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2584194) q[2];
sx q[2];
rz(-1.0018188) q[2];
sx q[2];
rz(-1.4460571) q[2];
rz(0.7575194) q[3];
sx q[3];
rz(-1.5599374) q[3];
sx q[3];
rz(2.3022046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.3676753) q[0];
sx q[0];
rz(-2.2125419) q[0];
sx q[0];
rz(1.0876592) q[0];
rz(-0.37519535) q[1];
sx q[1];
rz(-1.1253858) q[1];
sx q[1];
rz(2.1813724) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4488569) q[0];
sx q[0];
rz(-1.5310643) q[0];
sx q[0];
rz(-1.8360774) q[0];
x q[1];
rz(1.8619991) q[2];
sx q[2];
rz(-0.10379496) q[2];
sx q[2];
rz(-1.4196996) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6672047) q[1];
sx q[1];
rz(-1.1423813) q[1];
sx q[1];
rz(-2.0640255) q[1];
rz(-0.26953617) q[3];
sx q[3];
rz(-1.4785267) q[3];
sx q[3];
rz(2.5470981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9356161) q[2];
sx q[2];
rz(-1.874186) q[2];
sx q[2];
rz(1.4975632) q[2];
rz(2.2802672) q[3];
sx q[3];
rz(-2.438811) q[3];
sx q[3];
rz(0.45708814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.977026) q[0];
sx q[0];
rz(-0.67104665) q[0];
sx q[0];
rz(0.60428756) q[0];
rz(1.9970278) q[1];
sx q[1];
rz(-1.6555758) q[1];
sx q[1];
rz(2.7191275) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0351663) q[0];
sx q[0];
rz(-1.4894333) q[0];
sx q[0];
rz(-1.401713) q[0];
rz(-pi) q[1];
rz(2.5786262) q[2];
sx q[2];
rz(-1.0721231) q[2];
sx q[2];
rz(2.6359205) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3389272) q[1];
sx q[1];
rz(-0.69310729) q[1];
sx q[1];
rz(-3.0285804) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.72679584) q[3];
sx q[3];
rz(-2.6495547) q[3];
sx q[3];
rz(0.43471042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3182688) q[2];
sx q[2];
rz(-0.6508998) q[2];
sx q[2];
rz(2.4853415) q[2];
rz(1.6107791) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(2.5608565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73606473) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(-0.62498012) q[0];
rz(0.75195733) q[1];
sx q[1];
rz(-2.0729005) q[1];
sx q[1];
rz(1.6065074) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4706375) q[0];
sx q[0];
rz(-0.29330253) q[0];
sx q[0];
rz(-2.9706012) q[0];
rz(-0.70930945) q[2];
sx q[2];
rz(-1.0934208) q[2];
sx q[2];
rz(-0.8698744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8934763) q[1];
sx q[1];
rz(-1.3424338) q[1];
sx q[1];
rz(0.62947692) q[1];
rz(-pi) q[2];
rz(-2.296359) q[3];
sx q[3];
rz(-0.74016011) q[3];
sx q[3];
rz(-2.6773767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9525166) q[2];
sx q[2];
rz(-1.7837046) q[2];
sx q[2];
rz(0.9291741) q[2];
rz(2.3164228) q[3];
sx q[3];
rz(-1.630183) q[3];
sx q[3];
rz(-2.0391803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55649844) q[0];
sx q[0];
rz(-0.80699054) q[0];
sx q[0];
rz(1.7671385) q[0];
rz(0.51042405) q[1];
sx q[1];
rz(-0.7904895) q[1];
sx q[1];
rz(-0.62320954) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5730815) q[0];
sx q[0];
rz(-0.36751908) q[0];
sx q[0];
rz(-1.1373684) q[0];
x q[1];
rz(2.6323967) q[2];
sx q[2];
rz(-1.3993345) q[2];
sx q[2];
rz(0.96656884) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12909266) q[1];
sx q[1];
rz(-2.9575037) q[1];
sx q[1];
rz(-2.2232008) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3680063) q[3];
sx q[3];
rz(-0.8335127) q[3];
sx q[3];
rz(2.0059507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0906543) q[2];
sx q[2];
rz(-2.3263558) q[2];
sx q[2];
rz(1.9742924) q[2];
rz(3.1189611) q[3];
sx q[3];
rz(-2.6204717) q[3];
sx q[3];
rz(-1.1289271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24564329) q[0];
sx q[0];
rz(-2.0063945) q[0];
sx q[0];
rz(2.1242712) q[0];
rz(1.7474878) q[1];
sx q[1];
rz(-0.21109763) q[1];
sx q[1];
rz(0.62754935) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16495569) q[0];
sx q[0];
rz(-2.4058488) q[0];
sx q[0];
rz(2.1829905) q[0];
rz(-1.6553788) q[2];
sx q[2];
rz(-0.89327565) q[2];
sx q[2];
rz(1.5280452) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.095762091) q[1];
sx q[1];
rz(-1.374829) q[1];
sx q[1];
rz(-0.24230559) q[1];
rz(-pi) q[2];
rz(-1.2305607) q[3];
sx q[3];
rz(-1.1200532) q[3];
sx q[3];
rz(-0.72117248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.62577406) q[2];
sx q[2];
rz(-0.63483441) q[2];
sx q[2];
rz(-0.41075692) q[2];
rz(3.0913894) q[3];
sx q[3];
rz(-1.8886731) q[3];
sx q[3];
rz(-3.0661809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5589767) q[0];
sx q[0];
rz(-2.2447383) q[0];
sx q[0];
rz(-0.26891747) q[0];
rz(2.8315663) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(3.0115829) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0144671) q[0];
sx q[0];
rz(-1.3499182) q[0];
sx q[0];
rz(-1.9015903) q[0];
x q[1];
rz(0.79370462) q[2];
sx q[2];
rz(-2.6598601) q[2];
sx q[2];
rz(-0.46221737) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2490319) q[1];
sx q[1];
rz(-1.4098123) q[1];
sx q[1];
rz(2.575483) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6428516) q[3];
sx q[3];
rz(-0.95512701) q[3];
sx q[3];
rz(-1.9174674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.93398634) q[2];
sx q[2];
rz(-0.76675582) q[2];
sx q[2];
rz(-3.0450191) q[2];
rz(-1.5038331) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(-2.3418929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0132975) q[0];
sx q[0];
rz(-2.9533563) q[0];
sx q[0];
rz(-0.53303322) q[0];
rz(-3.10532) q[1];
sx q[1];
rz(-2.3590922) q[1];
sx q[1];
rz(1.689555) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8743458) q[0];
sx q[0];
rz(-1.1935992) q[0];
sx q[0];
rz(2.1826571) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31923652) q[2];
sx q[2];
rz(-1.4672654) q[2];
sx q[2];
rz(1.5872623) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9033135) q[1];
sx q[1];
rz(-0.36064816) q[1];
sx q[1];
rz(0.084899501) q[1];
rz(-pi) q[2];
x q[2];
rz(2.96569) q[3];
sx q[3];
rz(-0.47187343) q[3];
sx q[3];
rz(2.7681818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.71296802) q[2];
sx q[2];
rz(-0.73178256) q[2];
sx q[2];
rz(0.0326322) q[2];
rz(0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(2.0613861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7964771) q[0];
sx q[0];
rz(-0.65237541) q[0];
sx q[0];
rz(2.8271578) q[0];
rz(-0.79700094) q[1];
sx q[1];
rz(-2.2140257) q[1];
sx q[1];
rz(0.066233403) q[1];
rz(-0.5354698) q[2];
sx q[2];
rz(-1.8269804) q[2];
sx q[2];
rz(1.2852558) q[2];
rz(0.12049051) q[3];
sx q[3];
rz(-2.3269666) q[3];
sx q[3];
rz(2.5757488) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
