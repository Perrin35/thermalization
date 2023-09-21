OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4296071) q[0];
sx q[0];
rz(-0.40661943) q[0];
sx q[0];
rz(0.24917319) q[0];
rz(3.0781526) q[1];
sx q[1];
rz(-0.97172207) q[1];
sx q[1];
rz(2.5914153) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.881641) q[0];
sx q[0];
rz(-1.5978442) q[0];
sx q[0];
rz(1.3511488) q[0];
rz(2.1044188) q[2];
sx q[2];
rz(-1.954477) q[2];
sx q[2];
rz(1.2984315) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6018511) q[1];
sx q[1];
rz(-1.1907693) q[1];
sx q[1];
rz(2.3439581) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.689765) q[3];
sx q[3];
rz(-2.527378) q[3];
sx q[3];
rz(0.14106942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.41574079) q[2];
sx q[2];
rz(-0.44885138) q[2];
sx q[2];
rz(-2.501781) q[2];
rz(2.2885627) q[3];
sx q[3];
rz(-0.60522389) q[3];
sx q[3];
rz(-0.38133347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2663651) q[0];
sx q[0];
rz(-1.1667644) q[0];
sx q[0];
rz(-0.27045989) q[0];
rz(-0.71331435) q[1];
sx q[1];
rz(-2.1062873) q[1];
sx q[1];
rz(-1.5126022) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7249811) q[0];
sx q[0];
rz(-2.1583301) q[0];
sx q[0];
rz(-2.7532817) q[0];
x q[1];
rz(0.28378758) q[2];
sx q[2];
rz(-1.5740526) q[2];
sx q[2];
rz(0.46797215) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.70817972) q[1];
sx q[1];
rz(-2.4509894) q[1];
sx q[1];
rz(-1.3604926) q[1];
x q[2];
rz(2.9019722) q[3];
sx q[3];
rz(-1.2926971) q[3];
sx q[3];
rz(-0.77277377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2018532) q[2];
sx q[2];
rz(-0.24213232) q[2];
sx q[2];
rz(0.80336037) q[2];
rz(-2.0837636) q[3];
sx q[3];
rz(-1.6488896) q[3];
sx q[3];
rz(-3.1159744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8818883) q[0];
sx q[0];
rz(-1.0972247) q[0];
sx q[0];
rz(-2.846068) q[0];
rz(0.23513901) q[1];
sx q[1];
rz(-1.4328911) q[1];
sx q[1];
rz(-0.74584109) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.196516) q[0];
sx q[0];
rz(-1.3498422) q[0];
sx q[0];
rz(0.093284746) q[0];
rz(-pi) q[1];
rz(-1.7350115) q[2];
sx q[2];
rz(-2.0914372) q[2];
sx q[2];
rz(-2.8525713) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.44887603) q[1];
sx q[1];
rz(-1.8489031) q[1];
sx q[1];
rz(0.49726185) q[1];
rz(1.8196705) q[3];
sx q[3];
rz(-1.1296774) q[3];
sx q[3];
rz(2.6114024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0535584) q[2];
sx q[2];
rz(-3.1159248) q[2];
sx q[2];
rz(-0.68874613) q[2];
rz(-3.0912494) q[3];
sx q[3];
rz(-0.91209948) q[3];
sx q[3];
rz(1.5215727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.922309) q[0];
sx q[0];
rz(-2.121121) q[0];
sx q[0];
rz(0.10198378) q[0];
rz(3.02137) q[1];
sx q[1];
rz(-0.52821237) q[1];
sx q[1];
rz(0.27329683) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057987) q[0];
sx q[0];
rz(-1.4819643) q[0];
sx q[0];
rz(1.8949132) q[0];
rz(-pi) q[1];
rz(-2.4934713) q[2];
sx q[2];
rz(-2.0378049) q[2];
sx q[2];
rz(-1.8304706) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3015775) q[1];
sx q[1];
rz(-1.8043648) q[1];
sx q[1];
rz(2.9905) q[1];
rz(-0.64951879) q[3];
sx q[3];
rz(-1.7139385) q[3];
sx q[3];
rz(-2.3424145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3499202) q[2];
sx q[2];
rz(-0.73799729) q[2];
sx q[2];
rz(2.6376574) q[2];
rz(0.079581633) q[3];
sx q[3];
rz(-1.1527529) q[3];
sx q[3];
rz(2.7664405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85916096) q[0];
sx q[0];
rz(-2.0895884) q[0];
sx q[0];
rz(-3.058847) q[0];
rz(2.4619608) q[1];
sx q[1];
rz(-1.4928763) q[1];
sx q[1];
rz(0.98714978) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3521096) q[0];
sx q[0];
rz(-1.9728567) q[0];
sx q[0];
rz(-2.2228918) q[0];
rz(-pi) q[1];
rz(-3.1290595) q[2];
sx q[2];
rz(-0.64385157) q[2];
sx q[2];
rz(2.8963793) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.0070237006) q[1];
sx q[1];
rz(-0.94788523) q[1];
sx q[1];
rz(1.182215) q[1];
x q[2];
rz(1.1779285) q[3];
sx q[3];
rz(-1.6097027) q[3];
sx q[3];
rz(1.2272569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2510898) q[2];
sx q[2];
rz(-2.3513998) q[2];
sx q[2];
rz(-2.9303072) q[2];
rz(2.7206897) q[3];
sx q[3];
rz(-0.55287164) q[3];
sx q[3];
rz(0.063407272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.485065) q[0];
sx q[0];
rz(-1.8873029) q[0];
sx q[0];
rz(-0.791839) q[0];
rz(2.1461398) q[1];
sx q[1];
rz(-0.96770006) q[1];
sx q[1];
rz(-1.2794367) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5813542) q[0];
sx q[0];
rz(-1.2890153) q[0];
sx q[0];
rz(3.1262585) q[0];
rz(-2.7369667) q[2];
sx q[2];
rz(-0.64794108) q[2];
sx q[2];
rz(0.65587703) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2710323) q[1];
sx q[1];
rz(-2.6189657) q[1];
sx q[1];
rz(1.6305627) q[1];
rz(-pi) q[2];
rz(3.0002954) q[3];
sx q[3];
rz(-1.7501037) q[3];
sx q[3];
rz(2.3395204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.25097686) q[2];
sx q[2];
rz(-2.1112517) q[2];
sx q[2];
rz(-2.3708169) q[2];
rz(1.6714913) q[3];
sx q[3];
rz(-2.7272868) q[3];
sx q[3];
rz(-0.18338403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8188748) q[0];
sx q[0];
rz(-2.4920431) q[0];
sx q[0];
rz(3.0859257) q[0];
rz(-2.9260013) q[1];
sx q[1];
rz(-2.3781653) q[1];
sx q[1];
rz(-3.1380222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5647033) q[0];
sx q[0];
rz(-1.2029552) q[0];
sx q[0];
rz(2.6951615) q[0];
rz(-pi) q[1];
rz(0.41530208) q[2];
sx q[2];
rz(-0.79681444) q[2];
sx q[2];
rz(2.08564) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7230941) q[1];
sx q[1];
rz(-1.2730036) q[1];
sx q[1];
rz(-2.198092) q[1];
rz(-pi) q[2];
rz(-3.0655541) q[3];
sx q[3];
rz(-1.7685316) q[3];
sx q[3];
rz(2.7969489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.59721649) q[2];
sx q[2];
rz(-0.95149136) q[2];
sx q[2];
rz(2.2214831) q[2];
rz(2.9428633) q[3];
sx q[3];
rz(-1.9072429) q[3];
sx q[3];
rz(0.41771093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6253117) q[0];
sx q[0];
rz(-1.5100864) q[0];
sx q[0];
rz(-1.4177119) q[0];
rz(2.7334546) q[1];
sx q[1];
rz(-1.0393655) q[1];
sx q[1];
rz(-0.67869854) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1269826) q[0];
sx q[0];
rz(-1.8982366) q[0];
sx q[0];
rz(-1.2922657) q[0];
rz(-pi) q[1];
rz(2.0754201) q[2];
sx q[2];
rz(-2.3373211) q[2];
sx q[2];
rz(-1.7240766) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7486836) q[1];
sx q[1];
rz(-2.6137685) q[1];
sx q[1];
rz(-1.6960359) q[1];
x q[2];
rz(0.42550605) q[3];
sx q[3];
rz(-1.7153499) q[3];
sx q[3];
rz(-1.0021462) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.17710182) q[2];
sx q[2];
rz(-1.0937546) q[2];
sx q[2];
rz(-2.2237681) q[2];
rz(-1.9994036) q[3];
sx q[3];
rz(-2.2873788) q[3];
sx q[3];
rz(-2.2043622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1552102) q[0];
sx q[0];
rz(-0.6739524) q[0];
sx q[0];
rz(-2.881799) q[0];
rz(-0.70867509) q[1];
sx q[1];
rz(-0.27888137) q[1];
sx q[1];
rz(-2.6146467) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5567112) q[0];
sx q[0];
rz(-1.7877794) q[0];
sx q[0];
rz(2.7816539) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.2692659) q[2];
sx q[2];
rz(-0.93062799) q[2];
sx q[2];
rz(0.38538853) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.37104169) q[1];
sx q[1];
rz(-0.41401225) q[1];
sx q[1];
rz(0.71055926) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.55836375) q[3];
sx q[3];
rz(-2.6287968) q[3];
sx q[3];
rz(-2.5496897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8156585) q[2];
sx q[2];
rz(-2.2951173) q[2];
sx q[2];
rz(0.79088598) q[2];
rz(-2.7549426) q[3];
sx q[3];
rz(-2.1270042) q[3];
sx q[3];
rz(-2.8044243) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7336422) q[0];
sx q[0];
rz(-0.17245094) q[0];
sx q[0];
rz(-0.98544425) q[0];
rz(-2.573029) q[1];
sx q[1];
rz(-2.0303346) q[1];
sx q[1];
rz(-2.7808166) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4496778) q[0];
sx q[0];
rz(-1.5957818) q[0];
sx q[0];
rz(-1.7382394) q[0];
rz(-1.785727) q[2];
sx q[2];
rz(-1.6882009) q[2];
sx q[2];
rz(-1.4431151) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.39721397) q[1];
sx q[1];
rz(-0.75512868) q[1];
sx q[1];
rz(-0.20882864) q[1];
rz(1.7182699) q[3];
sx q[3];
rz(-1.4031938) q[3];
sx q[3];
rz(-0.22351219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0439904) q[2];
sx q[2];
rz(-0.68949914) q[2];
sx q[2];
rz(-0.94341755) q[2];
rz(-0.57389456) q[3];
sx q[3];
rz(-0.4807764) q[3];
sx q[3];
rz(-1.7452314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5744793) q[0];
sx q[0];
rz(-1.6945101) q[0];
sx q[0];
rz(2.2816818) q[0];
rz(1.3600596) q[1];
sx q[1];
rz(-1.0018476) q[1];
sx q[1];
rz(-0.76029653) q[1];
rz(-1.1311244) q[2];
sx q[2];
rz(-0.22618539) q[2];
sx q[2];
rz(-2.4629081) q[2];
rz(-2.4102224) q[3];
sx q[3];
rz(-0.89805713) q[3];
sx q[3];
rz(-3.0039207) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];