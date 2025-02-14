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
rz(1.709197) q[0];
sx q[0];
rz(3.4442918) q[0];
sx q[0];
rz(11.533307) q[0];
rz(1.9167702) q[1];
sx q[1];
rz(-1.6515825) q[1];
sx q[1];
rz(0.19436793) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8914999) q[0];
sx q[0];
rz(-0.95026267) q[0];
sx q[0];
rz(1.7583855) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35440234) q[2];
sx q[2];
rz(-0.72847073) q[2];
sx q[2];
rz(3.1284863) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.74339515) q[1];
sx q[1];
rz(-1.6170073) q[1];
sx q[1];
rz(-0.034548918) q[1];
rz(-pi) q[2];
rz(-0.57036542) q[3];
sx q[3];
rz(-0.70934848) q[3];
sx q[3];
rz(1.5523497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4091461) q[2];
sx q[2];
rz(-1.1017841) q[2];
sx q[2];
rz(-2.5771602) q[2];
rz(-1.8730646) q[3];
sx q[3];
rz(-3.1334435) q[3];
sx q[3];
rz(-2.234999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.059435189) q[0];
sx q[0];
rz(-3.1340288) q[0];
sx q[0];
rz(-2.2285158) q[0];
rz(-0.68436855) q[1];
sx q[1];
rz(-0.00033907779) q[1];
sx q[1];
rz(0.93904644) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9371999) q[0];
sx q[0];
rz(-2.4449722) q[0];
sx q[0];
rz(-2.533769) q[0];
rz(-pi) q[1];
rz(-3.1153283) q[2];
sx q[2];
rz(-2.677478) q[2];
sx q[2];
rz(2.8796706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.26015684) q[1];
sx q[1];
rz(-1.9828567) q[1];
sx q[1];
rz(-2.7753633) q[1];
rz(-0.066870832) q[3];
sx q[3];
rz(-1.1648199) q[3];
sx q[3];
rz(1.2479046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.3188476) q[2];
sx q[2];
rz(-0.007195909) q[2];
sx q[2];
rz(-0.89246559) q[2];
rz(-1.5626296) q[3];
sx q[3];
rz(-0.024450863) q[3];
sx q[3];
rz(1.310937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40176323) q[0];
sx q[0];
rz(-2.6391397) q[0];
sx q[0];
rz(0.40973642) q[0];
rz(-3.1306664) q[1];
sx q[1];
rz(-0.22053638) q[1];
sx q[1];
rz(1.3440557) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9896802) q[0];
sx q[0];
rz(-1.5912959) q[0];
sx q[0];
rz(-1.0597888) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.063882752) q[2];
sx q[2];
rz(-1.4555571) q[2];
sx q[2];
rz(1.2196397) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.1045055) q[1];
sx q[1];
rz(-1.6356138) q[1];
sx q[1];
rz(1.3124036) q[1];
x q[2];
rz(-1.4555172) q[3];
sx q[3];
rz(-1.5805827) q[3];
sx q[3];
rz(-0.20231314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4695796) q[2];
sx q[2];
rz(-1.4521705) q[2];
sx q[2];
rz(-0.043188485) q[2];
rz(2.0016661) q[3];
sx q[3];
rz(-2.98525) q[3];
sx q[3];
rz(-1.9388916) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4129055) q[0];
sx q[0];
rz(-0.40189704) q[0];
sx q[0];
rz(-1.8034978) q[0];
rz(0.69522229) q[1];
sx q[1];
rz(-3.0137364) q[1];
sx q[1];
rz(-0.41740886) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0510088) q[0];
sx q[0];
rz(-0.90253297) q[0];
sx q[0];
rz(0.21940748) q[0];
rz(-pi) q[1];
x q[1];
rz(2.874006) q[2];
sx q[2];
rz(-0.94174615) q[2];
sx q[2];
rz(-2.3305394) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.7898832) q[1];
sx q[1];
rz(-1.7725866) q[1];
sx q[1];
rz(2.2191597) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8700852) q[3];
sx q[3];
rz(-2.6911754) q[3];
sx q[3];
rz(1.3189486) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7030316) q[2];
sx q[2];
rz(-0.87623864) q[2];
sx q[2];
rz(-1.38928) q[2];
rz(0.82052463) q[3];
sx q[3];
rz(-1.5503784) q[3];
sx q[3];
rz(-2.4337721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.163212) q[0];
sx q[0];
rz(-3.1040525) q[0];
sx q[0];
rz(-0.96330825) q[0];
rz(-2.9523201) q[1];
sx q[1];
rz(-3.1258588) q[1];
sx q[1];
rz(-2.9761369) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4625797) q[0];
sx q[0];
rz(-1.5272584) q[0];
sx q[0];
rz(-1.5971558) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2654183) q[2];
sx q[2];
rz(-1.3378222) q[2];
sx q[2];
rz(-2.7829952) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.4374809) q[1];
sx q[1];
rz(-2.3299814) q[1];
sx q[1];
rz(2.981217) q[1];
rz(2.1894375) q[3];
sx q[3];
rz(-0.51636558) q[3];
sx q[3];
rz(-1.0960033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8339707) q[2];
sx q[2];
rz(-0.51887363) q[2];
sx q[2];
rz(2.0818254) q[2];
rz(2.4506532) q[3];
sx q[3];
rz(-2.2773404) q[3];
sx q[3];
rz(1.2714269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6178013) q[0];
sx q[0];
rz(-0.093136223) q[0];
sx q[0];
rz(-1.6049438) q[0];
rz(-0.45283428) q[1];
sx q[1];
rz(-3.1335399) q[1];
sx q[1];
rz(-1.7013928) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33537597) q[0];
sx q[0];
rz(-1.7378238) q[0];
sx q[0];
rz(1.7352292) q[0];
rz(-pi) q[1];
rz(2.9199706) q[2];
sx q[2];
rz(-1.3563507) q[2];
sx q[2];
rz(2.1943486) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3611412) q[1];
sx q[1];
rz(-1.5781501) q[1];
sx q[1];
rz(2.9869534) q[1];
rz(1.4000638) q[3];
sx q[3];
rz(-1.6686107) q[3];
sx q[3];
rz(-2.3119237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77233934) q[2];
sx q[2];
rz(-2.1489096) q[2];
sx q[2];
rz(3.0625647) q[2];
rz(2.2760271) q[3];
sx q[3];
rz(-2.1871958) q[3];
sx q[3];
rz(-1.0237833) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2257776) q[0];
sx q[0];
rz(-3.1362035) q[0];
sx q[0];
rz(0.22501568) q[0];
rz(-0.30613884) q[1];
sx q[1];
rz(-3.1248416) q[1];
sx q[1];
rz(-0.87047815) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8057833) q[0];
sx q[0];
rz(-3.0558944) q[0];
sx q[0];
rz(0.58914574) q[0];
rz(-pi) q[1];
rz(2.5092431) q[2];
sx q[2];
rz(-2.3195757) q[2];
sx q[2];
rz(2.1123304) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.051389305) q[1];
sx q[1];
rz(-1.6312113) q[1];
sx q[1];
rz(2.3714972) q[1];
rz(-2.2629396) q[3];
sx q[3];
rz(-2.5402398) q[3];
sx q[3];
rz(-2.2215312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.47591448) q[2];
sx q[2];
rz(-3.0147538) q[2];
sx q[2];
rz(-1.0352943) q[2];
rz(-0.78694844) q[3];
sx q[3];
rz(-0.042777177) q[3];
sx q[3];
rz(-0.27534819) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.03054522) q[0];
sx q[0];
rz(-3.1304517) q[0];
sx q[0];
rz(-3.1159478) q[0];
rz(1.8772839) q[1];
sx q[1];
rz(-3.1176716) q[1];
sx q[1];
rz(0.69902507) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5584939) q[0];
sx q[0];
rz(-1.8374658) q[0];
sx q[0];
rz(1.3400274) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7412296) q[2];
sx q[2];
rz(-2.0336069) q[2];
sx q[2];
rz(1.5968666) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0720585) q[1];
sx q[1];
rz(-1.9542819) q[1];
sx q[1];
rz(-1.6504565) q[1];
x q[2];
rz(0.89939846) q[3];
sx q[3];
rz(-0.77435437) q[3];
sx q[3];
rz(-3.0029861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8890624) q[2];
sx q[2];
rz(-0.74873304) q[2];
sx q[2];
rz(-0.32386455) q[2];
rz(2.9845386) q[3];
sx q[3];
rz(-1.9040949) q[3];
sx q[3];
rz(0.15373716) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5667628) q[0];
sx q[0];
rz(-3.1269508) q[0];
sx q[0];
rz(2.5515442) q[0];
rz(0.7198965) q[1];
sx q[1];
rz(-3.0820334) q[1];
sx q[1];
rz(-0.86729008) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27408263) q[0];
sx q[0];
rz(-1.4658341) q[0];
sx q[0];
rz(2.5068953) q[0];
rz(-pi) q[1];
x q[1];
rz(0.9684674) q[2];
sx q[2];
rz(-1.3545389) q[2];
sx q[2];
rz(2.575084) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.54360753) q[1];
sx q[1];
rz(-1.6389009) q[1];
sx q[1];
rz(1.6833413) q[1];
rz(-1.7884364) q[3];
sx q[3];
rz(-3*pi/13) q[3];
sx q[3];
rz(0.21228895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1934293) q[2];
sx q[2];
rz(-2.2488504) q[2];
sx q[2];
rz(0.11290045) q[2];
rz(-2.0924977) q[3];
sx q[3];
rz(-1.5594522) q[3];
sx q[3];
rz(0.73672867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0650487) q[0];
sx q[0];
rz(-1.5600518) q[0];
sx q[0];
rz(-1.5034058) q[0];
rz(2.5050971) q[1];
sx q[1];
rz(-0.46000767) q[1];
sx q[1];
rz(-1.5696625) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8638679) q[0];
sx q[0];
rz(-1.6536599) q[0];
sx q[0];
rz(-3.0623097) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.0031329069) q[2];
sx q[2];
rz(-1.5701446) q[2];
sx q[2];
rz(-0.80564431) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7074938) q[1];
sx q[1];
rz(-1.5710305) q[1];
sx q[1];
rz(3.1386761) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51497634) q[3];
sx q[3];
rz(-1.4287717) q[3];
sx q[3];
rz(1.2339301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2703209) q[2];
sx q[2];
rz(-1.511829) q[2];
sx q[2];
rz(2.9809269) q[2];
rz(1.9130982) q[3];
sx q[3];
rz(-3.1077423) q[3];
sx q[3];
rz(0.20139774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0495618) q[0];
sx q[0];
rz(-1.1783538) q[0];
sx q[0];
rz(0.023068064) q[0];
rz(1.6231712) q[1];
sx q[1];
rz(-0.36133125) q[1];
sx q[1];
rz(0.28892118) q[1];
rz(-1.6466181) q[2];
sx q[2];
rz(-0.29360148) q[2];
sx q[2];
rz(0.32099024) q[2];
rz(-2.992441) q[3];
sx q[3];
rz(-2.0059398) q[3];
sx q[3];
rz(-2.5930349) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
