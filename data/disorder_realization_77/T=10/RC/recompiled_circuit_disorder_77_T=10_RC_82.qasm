OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(-0.84364426) q[0];
sx q[0];
rz(0.16790976) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(-0.056161031) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355609) q[0];
sx q[0];
rz(-2.6132085) q[0];
sx q[0];
rz(2.0023268) q[0];
rz(-0.92476966) q[2];
sx q[2];
rz(-1.7416818) q[2];
sx q[2];
rz(0.31121635) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.010477) q[1];
sx q[1];
rz(-1.8382204) q[1];
sx q[1];
rz(1.0130151) q[1];
rz(2.9763016) q[3];
sx q[3];
rz(-2.878302) q[3];
sx q[3];
rz(-0.96112448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.37796676) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(-1.1928605) q[3];
sx q[3];
rz(-1.9038707) q[3];
sx q[3];
rz(-2.7584934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(1.312785) q[0];
rz(2.9361172) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(1.1516494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5702471) q[0];
sx q[0];
rz(-0.76128188) q[0];
sx q[0];
rz(1.135457) q[0];
x q[1];
rz(1.0863016) q[2];
sx q[2];
rz(-1.0420348) q[2];
sx q[2];
rz(2.8690086) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.093106) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(2.6622245) q[1];
rz(-pi) q[2];
rz(0.63668164) q[3];
sx q[3];
rz(-0.59906206) q[3];
sx q[3];
rz(-2.3707795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1318704) q[2];
sx q[2];
rz(-1.4902318) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(2.7644073) q[3];
sx q[3];
rz(-0.42755869) q[3];
sx q[3];
rz(2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8310228) q[0];
sx q[0];
rz(-0.092386827) q[0];
sx q[0];
rz(-0.036852766) q[0];
rz(2.316078) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(0.056578606) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53783137) q[0];
sx q[0];
rz(-2.1268401) q[0];
sx q[0];
rz(2.6105196) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47358863) q[2];
sx q[2];
rz(-2.8550672) q[2];
sx q[2];
rz(-0.63602704) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1287071) q[1];
sx q[1];
rz(-1.7696847) q[1];
sx q[1];
rz(0.30602869) q[1];
rz(-0.64485456) q[3];
sx q[3];
rz(-1.8861024) q[3];
sx q[3];
rz(2.9063318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8686707) q[2];
sx q[2];
rz(-1.6438831) q[2];
sx q[2];
rz(2.2154714) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(1.042897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91519231) q[0];
sx q[0];
rz(-0.72015786) q[0];
sx q[0];
rz(-0.74209374) q[0];
rz(-1.1391976) q[1];
sx q[1];
rz(-2.6622055) q[1];
sx q[1];
rz(-2.6779968) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47953654) q[0];
sx q[0];
rz(-0.89131309) q[0];
sx q[0];
rz(2.7583073) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33275231) q[2];
sx q[2];
rz(-1.5126265) q[2];
sx q[2];
rz(-2.1259049) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61958085) q[1];
sx q[1];
rz(-1.7420235) q[1];
sx q[1];
rz(-0.79230688) q[1];
rz(-2.2673244) q[3];
sx q[3];
rz(-2.4618751) q[3];
sx q[3];
rz(3.0391039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7745557) q[2];
sx q[2];
rz(-2.98428) q[2];
sx q[2];
rz(3.0920933) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(0.29770011) q[0];
rz(2.659335) q[1];
sx q[1];
rz(-2.3869956) q[1];
sx q[1];
rz(-2.1972426) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34236318) q[0];
sx q[0];
rz(-2.9276597) q[0];
sx q[0];
rz(1.3571204) q[0];
rz(-pi) q[1];
rz(2.7371251) q[2];
sx q[2];
rz(-0.81768113) q[2];
sx q[2];
rz(-0.58194619) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9152865) q[1];
sx q[1];
rz(-1.5686791) q[1];
sx q[1];
rz(1.5096942) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.827042) q[3];
sx q[3];
rz(-1.1503997) q[3];
sx q[3];
rz(2.8375569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.258761) q[2];
sx q[2];
rz(-1.0927039) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(-0.0023068874) q[3];
sx q[3];
rz(-1.5294411) q[3];
sx q[3];
rz(-2.8172857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7047983) q[0];
sx q[0];
rz(-0.44610766) q[0];
sx q[0];
rz(-0.58445245) q[0];
rz(-2.2553518) q[1];
sx q[1];
rz(-2.5247572) q[1];
sx q[1];
rz(0.054919682) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45860896) q[0];
sx q[0];
rz(-0.28124547) q[0];
sx q[0];
rz(1.2693229) q[0];
x q[1];
rz(-1.5149649) q[2];
sx q[2];
rz(-2.8179114) q[2];
sx q[2];
rz(-1.2602381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9565441) q[1];
sx q[1];
rz(-1.946432) q[1];
sx q[1];
rz(2.5513785) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92298569) q[3];
sx q[3];
rz(-2.4499948) q[3];
sx q[3];
rz(1.637961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-0.12067623) q[2];
sx q[2];
rz(-1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(1.8090766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(-2.8570535) q[0];
rz(2.1971205) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(0.91032666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1032573) q[0];
sx q[0];
rz(-0.084780134) q[0];
sx q[0];
rz(-0.8707365) q[0];
x q[1];
rz(1.1548642) q[2];
sx q[2];
rz(-1.8850733) q[2];
sx q[2];
rz(-1.2736125) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6701723) q[1];
sx q[1];
rz(-2.6176665) q[1];
sx q[1];
rz(0.36797221) q[1];
rz(2.5038239) q[3];
sx q[3];
rz(-0.83003269) q[3];
sx q[3];
rz(-2.4142767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3900782) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(-2.2195623) q[2];
rz(-0.56728029) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(0.89408016) q[0];
sx q[0];
rz(-2.4649354) q[0];
sx q[0];
rz(-0.12938736) q[0];
rz(2.5091876) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(0.30050373) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.520641) q[0];
sx q[0];
rz(-1.2004939) q[0];
sx q[0];
rz(-2.305549) q[0];
x q[1];
rz(0.6408765) q[2];
sx q[2];
rz(-0.87101988) q[2];
sx q[2];
rz(-1.4027632) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0771675) q[1];
sx q[1];
rz(-2.7556813) q[1];
sx q[1];
rz(-0.47972958) q[1];
x q[2];
rz(0.31605966) q[3];
sx q[3];
rz(-0.83871597) q[3];
sx q[3];
rz(-0.72533208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58632103) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(-0.78197455) q[2];
rz(-0.55109763) q[3];
sx q[3];
rz(-1.3827773) q[3];
sx q[3];
rz(2.9836695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5683811) q[0];
sx q[0];
rz(-1.9848354) q[0];
sx q[0];
rz(-0.12776275) q[0];
rz(-2.5993775) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(-2.382747) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1369143) q[0];
sx q[0];
rz(-1.5877962) q[0];
sx q[0];
rz(1.1466115) q[0];
rz(-pi) q[1];
x q[1];
rz(0.27300948) q[2];
sx q[2];
rz(-0.62676478) q[2];
sx q[2];
rz(3.0221456) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85941852) q[1];
sx q[1];
rz(-0.71600435) q[1];
sx q[1];
rz(0.23846682) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5893448) q[3];
sx q[3];
rz(-0.66569257) q[3];
sx q[3];
rz(1.9513643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1252497) q[2];
sx q[2];
rz(-1.3602076) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(1.7193433) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(2.0786044) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816417) q[0];
sx q[0];
rz(-0.61976969) q[0];
sx q[0];
rz(-0.075335659) q[0];
rz(2.244859) q[1];
sx q[1];
rz(-1.1743841) q[1];
sx q[1];
rz(0.60992253) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7216435) q[0];
sx q[0];
rz(-2.3099265) q[0];
sx q[0];
rz(-1.9589817) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5209849) q[2];
sx q[2];
rz(-0.45244103) q[2];
sx q[2];
rz(-1.1021745) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8049106) q[1];
sx q[1];
rz(-1.8098127) q[1];
sx q[1];
rz(0.75093021) q[1];
rz(0.742357) q[3];
sx q[3];
rz(-1.9543813) q[3];
sx q[3];
rz(-2.50768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.909409) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(-2.4278736) q[2];
rz(2.7632726) q[3];
sx q[3];
rz(-0.49351966) q[3];
sx q[3];
rz(0.87987125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(2.4907885) q[1];
sx q[1];
rz(-1.4918068) q[1];
sx q[1];
rz(3.0202958) q[1];
rz(-3.1109839) q[2];
sx q[2];
rz(-1.7666713) q[2];
sx q[2];
rz(-0.91790744) q[2];
rz(1.1013423) q[3];
sx q[3];
rz(-0.51870844) q[3];
sx q[3];
rz(-1.4389256) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
