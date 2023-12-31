OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2965887) q[0];
sx q[0];
rz(3.8656524) q[0];
sx q[0];
rz(11.081628) q[0];
rz(1.2031263) q[1];
sx q[1];
rz(-0.523518) q[1];
sx q[1];
rz(2.2533921) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48924016) q[0];
sx q[0];
rz(-1.7416735) q[0];
sx q[0];
rz(-1.5255552) q[0];
rz(-pi) q[1];
rz(-0.16562478) q[2];
sx q[2];
rz(-2.1016444) q[2];
sx q[2];
rz(0.80425516) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2759309) q[1];
sx q[1];
rz(-1.4038424) q[1];
sx q[1];
rz(1.5702412) q[1];
x q[2];
rz(-0.039460823) q[3];
sx q[3];
rz(-2.4426115) q[3];
sx q[3];
rz(2.3208502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.28329864) q[2];
sx q[2];
rz(-0.41559872) q[2];
sx q[2];
rz(2.0236012) q[2];
rz(2.9962712) q[3];
sx q[3];
rz(-1.558692) q[3];
sx q[3];
rz(3.0956691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5685101) q[0];
sx q[0];
rz(-1.8772323) q[0];
sx q[0];
rz(-0.65482393) q[0];
rz(1.2163935) q[1];
sx q[1];
rz(-1.1647859) q[1];
sx q[1];
rz(-0.12589802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1032216) q[0];
sx q[0];
rz(-1.9403337) q[0];
sx q[0];
rz(-2.080337) q[0];
x q[1];
rz(0.84463859) q[2];
sx q[2];
rz(-0.8234878) q[2];
sx q[2];
rz(2.5708432) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.241908) q[1];
sx q[1];
rz(-3.0113314) q[1];
sx q[1];
rz(-1.8036519) q[1];
rz(-pi) q[2];
rz(0.010766518) q[3];
sx q[3];
rz(-2.2942703) q[3];
sx q[3];
rz(0.41124287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.048432365) q[2];
sx q[2];
rz(-1.9530714) q[2];
sx q[2];
rz(-2.753567) q[2];
rz(-1.7175425) q[3];
sx q[3];
rz(-0.63801304) q[3];
sx q[3];
rz(0.58732906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5557142) q[0];
sx q[0];
rz(-2.3590187) q[0];
sx q[0];
rz(-3.0622603) q[0];
rz(-0.084005984) q[1];
sx q[1];
rz(-2.3386798) q[1];
sx q[1];
rz(-1.9817339) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40354363) q[0];
sx q[0];
rz(-1.6532716) q[0];
sx q[0];
rz(2.0820046) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94295393) q[2];
sx q[2];
rz(-1.1477594) q[2];
sx q[2];
rz(1.9726276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.68418903) q[1];
sx q[1];
rz(-2.0261129) q[1];
sx q[1];
rz(-1.1475569) q[1];
rz(-pi) q[2];
rz(1.0252762) q[3];
sx q[3];
rz(-1.6189515) q[3];
sx q[3];
rz(2.2930876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7362061) q[2];
sx q[2];
rz(-1.3233041) q[2];
sx q[2];
rz(-0.1082871) q[2];
rz(0.64374271) q[3];
sx q[3];
rz(-2.0635922) q[3];
sx q[3];
rz(0.61557499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.165034) q[0];
sx q[0];
rz(-1.3950011) q[0];
sx q[0];
rz(-1.6595586) q[0];
rz(0.7011134) q[1];
sx q[1];
rz(-1.2524403) q[1];
sx q[1];
rz(2.8569417) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92589256) q[0];
sx q[0];
rz(-1.4963576) q[0];
sx q[0];
rz(-2.0490993) q[0];
rz(-pi) q[1];
rz(1.8585763) q[2];
sx q[2];
rz(-1.3735526) q[2];
sx q[2];
rz(-1.6523199) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2440566) q[1];
sx q[1];
rz(-1.4923864) q[1];
sx q[1];
rz(-1.4439911) q[1];
rz(-2.2399726) q[3];
sx q[3];
rz(-1.9603143) q[3];
sx q[3];
rz(0.31182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.231679) q[2];
sx q[2];
rz(-0.96857962) q[2];
sx q[2];
rz(-2.5740734) q[2];
rz(2.7275758) q[3];
sx q[3];
rz(-1.1375789) q[3];
sx q[3];
rz(-1.1119941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.652997) q[0];
sx q[0];
rz(-2.1814006) q[0];
sx q[0];
rz(1.8150785) q[0];
rz(1.1524221) q[1];
sx q[1];
rz(-1.7633341) q[1];
sx q[1];
rz(-2.2036536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91468231) q[0];
sx q[0];
rz(-1.662928) q[0];
sx q[0];
rz(0.021962086) q[0];
rz(2.6229834) q[2];
sx q[2];
rz(-2.0687639) q[2];
sx q[2];
rz(1.9260315) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.282498) q[1];
sx q[1];
rz(-1.6051834) q[1];
sx q[1];
rz(-3.1153468) q[1];
rz(-0.33850833) q[3];
sx q[3];
rz(-0.37422985) q[3];
sx q[3];
rz(2.3238473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1753297) q[2];
sx q[2];
rz(-0.18925174) q[2];
sx q[2];
rz(-2.1002634) q[2];
rz(1.3025618) q[3];
sx q[3];
rz(-2.0077191) q[3];
sx q[3];
rz(-1.8235824) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43907169) q[0];
sx q[0];
rz(-1.9938001) q[0];
sx q[0];
rz(-0.55737108) q[0];
rz(2.5769261) q[1];
sx q[1];
rz(-2.4319885) q[1];
sx q[1];
rz(-0.55647892) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0628478) q[0];
sx q[0];
rz(-1.9468465) q[0];
sx q[0];
rz(-1.7929121) q[0];
rz(-pi) q[1];
rz(0.85530497) q[2];
sx q[2];
rz(-1.4391293) q[2];
sx q[2];
rz(1.5065187) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4887052) q[1];
sx q[1];
rz(-2.8294551) q[1];
sx q[1];
rz(2.1991792) q[1];
rz(1.6218833) q[3];
sx q[3];
rz(-0.79569492) q[3];
sx q[3];
rz(1.5339799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6094728) q[2];
sx q[2];
rz(-2.3807821) q[2];
sx q[2];
rz(1.9167985) q[2];
rz(-1.6312284) q[3];
sx q[3];
rz(-1.8211726) q[3];
sx q[3];
rz(1.5009376) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0790134) q[0];
sx q[0];
rz(-1.2591079) q[0];
sx q[0];
rz(3.0986837) q[0];
rz(0.91730109) q[1];
sx q[1];
rz(-2.5179472) q[1];
sx q[1];
rz(-0.50061217) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4780873) q[0];
sx q[0];
rz(-2.9604719) q[0];
sx q[0];
rz(1.3785133) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.074965) q[2];
sx q[2];
rz(-2.5410286) q[2];
sx q[2];
rz(1.9753089) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4088879) q[1];
sx q[1];
rz(-1.9703456) q[1];
sx q[1];
rz(-1.7016181) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41916267) q[3];
sx q[3];
rz(-2.3401642) q[3];
sx q[3];
rz(-2.6475346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5381955) q[2];
sx q[2];
rz(-2.0624702) q[2];
sx q[2];
rz(-0.67561692) q[2];
rz(-0.44089857) q[3];
sx q[3];
rz(-1.7023804) q[3];
sx q[3];
rz(-2.6202257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.065141) q[0];
sx q[0];
rz(-2.8171709) q[0];
sx q[0];
rz(-2.0741529) q[0];
rz(0.43287977) q[1];
sx q[1];
rz(-1.503711) q[1];
sx q[1];
rz(-1.4656461) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6976801) q[0];
sx q[0];
rz(-1.7051823) q[0];
sx q[0];
rz(1.105955) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6797811) q[2];
sx q[2];
rz(-1.7562859) q[2];
sx q[2];
rz(-2.2456004) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88361909) q[1];
sx q[1];
rz(-2.5273364) q[1];
sx q[1];
rz(2.1780464) q[1];
x q[2];
rz(-0.96887178) q[3];
sx q[3];
rz(-1.0191917) q[3];
sx q[3];
rz(-1.5920036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8026768) q[2];
sx q[2];
rz(-0.85614506) q[2];
sx q[2];
rz(-1.6652997) q[2];
rz(0.87351292) q[3];
sx q[3];
rz(-2.2647808) q[3];
sx q[3];
rz(-2.6749558) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8840238) q[0];
sx q[0];
rz(-2.5654061) q[0];
sx q[0];
rz(-2.4043758) q[0];
rz(-3.1228512) q[1];
sx q[1];
rz(-0.32967162) q[1];
sx q[1];
rz(0.92528701) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4243471) q[0];
sx q[0];
rz(-0.92110094) q[0];
sx q[0];
rz(0.93066494) q[0];
x q[1];
rz(0.66321744) q[2];
sx q[2];
rz(-1.6753917) q[2];
sx q[2];
rz(1.516972) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.195203) q[1];
sx q[1];
rz(-2.9159947) q[1];
sx q[1];
rz(2.5220847) q[1];
rz(0.85261811) q[3];
sx q[3];
rz(-1.2468306) q[3];
sx q[3];
rz(-0.94911239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.76901889) q[2];
sx q[2];
rz(-1.5344658) q[2];
sx q[2];
rz(1.0037237) q[2];
rz(0.090099661) q[3];
sx q[3];
rz(-3.1153479) q[3];
sx q[3];
rz(2.0285006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1019679) q[0];
sx q[0];
rz(-1.573338) q[0];
sx q[0];
rz(1.2596624) q[0];
rz(2.8758077) q[1];
sx q[1];
rz(-2.5140285) q[1];
sx q[1];
rz(-2.3840747) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3525317) q[0];
sx q[0];
rz(-1.4333945) q[0];
sx q[0];
rz(1.3760516) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1610356) q[2];
sx q[2];
rz(-1.7014116) q[2];
sx q[2];
rz(2.5622501) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1265035) q[1];
sx q[1];
rz(-1.1474097) q[1];
sx q[1];
rz(2.8979315) q[1];
rz(-2.431589) q[3];
sx q[3];
rz(-1.5922976) q[3];
sx q[3];
rz(-1.3083252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.99047986) q[2];
sx q[2];
rz(-1.3157142) q[2];
sx q[2];
rz(0.79375664) q[2];
rz(2.5027067) q[3];
sx q[3];
rz(-0.34404889) q[3];
sx q[3];
rz(1.0206153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6486075) q[0];
sx q[0];
rz(-0.47825559) q[0];
sx q[0];
rz(-0.91260845) q[0];
rz(1.6408625) q[1];
sx q[1];
rz(-2.224557) q[1];
sx q[1];
rz(1.7932737) q[1];
rz(-1.2407672) q[2];
sx q[2];
rz(-1.9874265) q[2];
sx q[2];
rz(-0.12513587) q[2];
rz(-2.0414447) q[3];
sx q[3];
rz(-0.63334076) q[3];
sx q[3];
rz(0.28950194) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
