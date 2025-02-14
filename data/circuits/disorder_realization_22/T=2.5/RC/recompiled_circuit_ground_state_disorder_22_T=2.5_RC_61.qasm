OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.9353256) q[0];
sx q[0];
rz(-1.5201818) q[0];
sx q[0];
rz(-1.4186463) q[0];
rz(0.043579276) q[1];
sx q[1];
rz(-1.2830623) q[1];
sx q[1];
rz(-2.9820774) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1698604) q[0];
sx q[0];
rz(-1.5423449) q[0];
sx q[0];
rz(-1.736743) q[0];
rz(-pi) q[1];
rz(-2.7829917) q[2];
sx q[2];
rz(-2.498543) q[2];
sx q[2];
rz(-1.7687163) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7456513) q[1];
sx q[1];
rz(-2.6439754) q[1];
sx q[1];
rz(1.8630337) q[1];
rz(-pi) q[2];
rz(-2.2920424) q[3];
sx q[3];
rz(-1.2952943) q[3];
sx q[3];
rz(-1.7003789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.688545) q[2];
sx q[2];
rz(-1.1587554) q[2];
sx q[2];
rz(3.0513406) q[2];
rz(-0.074617535) q[3];
sx q[3];
rz(-1.6736232) q[3];
sx q[3];
rz(-0.41372764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8892882) q[0];
sx q[0];
rz(-1.1587208) q[0];
sx q[0];
rz(1.4537551) q[0];
rz(-2.3989035) q[1];
sx q[1];
rz(-1.7748666) q[1];
sx q[1];
rz(-2.3765366) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9738282) q[0];
sx q[0];
rz(-2.544738) q[0];
sx q[0];
rz(3.0227227) q[0];
rz(-0.68816784) q[2];
sx q[2];
rz(-0.98994561) q[2];
sx q[2];
rz(2.1877387) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7530427) q[1];
sx q[1];
rz(-1.2425798) q[1];
sx q[1];
rz(0.80372827) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5935947) q[3];
sx q[3];
rz(-1.5685076) q[3];
sx q[3];
rz(-1.8016004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9146156) q[2];
sx q[2];
rz(-2.1299629) q[2];
sx q[2];
rz(1.6870618) q[2];
rz(2.4948273) q[3];
sx q[3];
rz(-0.55964595) q[3];
sx q[3];
rz(1.2796992) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0890559) q[0];
sx q[0];
rz(-1.6707957) q[0];
sx q[0];
rz(-0.045150969) q[0];
rz(1.9106983) q[1];
sx q[1];
rz(-1.8321593) q[1];
sx q[1];
rz(-2.0848354) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.089356747) q[0];
sx q[0];
rz(-2.1187021) q[0];
sx q[0];
rz(-0.38378451) q[0];
rz(-0.62303876) q[2];
sx q[2];
rz(-2.3941561) q[2];
sx q[2];
rz(-1.7926271) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9581117) q[1];
sx q[1];
rz(-2.2327802) q[1];
sx q[1];
rz(-2.8173911) q[1];
rz(-1.8523434) q[3];
sx q[3];
rz(-0.70812884) q[3];
sx q[3];
rz(-0.81048465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7354273) q[2];
sx q[2];
rz(-0.96144095) q[2];
sx q[2];
rz(2.4889634) q[2];
rz(-1.2930219) q[3];
sx q[3];
rz(-0.76904622) q[3];
sx q[3];
rz(3.0434216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80766135) q[0];
sx q[0];
rz(-0.46285358) q[0];
sx q[0];
rz(1.7474784) q[0];
rz(-1.8380503) q[1];
sx q[1];
rz(-1.1640254) q[1];
sx q[1];
rz(1.0882783) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0387087) q[0];
sx q[0];
rz(-2.3293478) q[0];
sx q[0];
rz(-1.1884304) q[0];
rz(-pi) q[1];
rz(0.84962623) q[2];
sx q[2];
rz(-1.0996981) q[2];
sx q[2];
rz(-1.8249408) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.478178) q[1];
sx q[1];
rz(-0.82163787) q[1];
sx q[1];
rz(0.68382241) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.518098) q[3];
sx q[3];
rz(-2.0695588) q[3];
sx q[3];
rz(-0.76585211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3950222) q[2];
sx q[2];
rz(-1.3777379) q[2];
sx q[2];
rz(-2.6611888) q[2];
rz(-0.26614842) q[3];
sx q[3];
rz(-1.9726617) q[3];
sx q[3];
rz(2.2973785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0323623) q[0];
sx q[0];
rz(-2.5205595) q[0];
sx q[0];
rz(1.9367628) q[0];
rz(0.042639848) q[1];
sx q[1];
rz(-1.3748906) q[1];
sx q[1];
rz(2.1991594) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54357994) q[0];
sx q[0];
rz(-2.6930598) q[0];
sx q[0];
rz(-1.2537812) q[0];
rz(-0.82289121) q[2];
sx q[2];
rz(-2.3510755) q[2];
sx q[2];
rz(-0.86871494) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1548295) q[1];
sx q[1];
rz(-1.0752605) q[1];
sx q[1];
rz(-2.7096728) q[1];
x q[2];
rz(-0.14640866) q[3];
sx q[3];
rz(-0.64420986) q[3];
sx q[3];
rz(-0.25800426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.11613906) q[2];
sx q[2];
rz(-0.68486539) q[2];
sx q[2];
rz(1.2241414) q[2];
rz(0.72743607) q[3];
sx q[3];
rz(-1.6614611) q[3];
sx q[3];
rz(-3.091403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8023476) q[0];
sx q[0];
rz(-0.50528637) q[0];
sx q[0];
rz(2.8832866) q[0];
rz(-0.91336617) q[1];
sx q[1];
rz(-2.2759627) q[1];
sx q[1];
rz(-2.1736653) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0845522) q[0];
sx q[0];
rz(-1.0121317) q[0];
sx q[0];
rz(1.1112122) q[0];
rz(2.2931678) q[2];
sx q[2];
rz(-1.8471187) q[2];
sx q[2];
rz(-1.6538617) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.591394) q[1];
sx q[1];
rz(-1.0934783) q[1];
sx q[1];
rz(1.5309019) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0276035) q[3];
sx q[3];
rz(-0.99811998) q[3];
sx q[3];
rz(3.0610994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4919058) q[2];
sx q[2];
rz(-2.457149) q[2];
sx q[2];
rz(2.7223132) q[2];
rz(2.3573917) q[3];
sx q[3];
rz(-2.1214285) q[3];
sx q[3];
rz(-1.3721344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.075031) q[0];
sx q[0];
rz(-3.0176268) q[0];
sx q[0];
rz(-3.1228464) q[0];
rz(-0.092451267) q[1];
sx q[1];
rz(-1.8094614) q[1];
sx q[1];
rz(-0.094955347) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0571916) q[0];
sx q[0];
rz(-0.6633299) q[0];
sx q[0];
rz(-2.6003772) q[0];
rz(-pi) q[1];
rz(-1.1509535) q[2];
sx q[2];
rz(-2.125084) q[2];
sx q[2];
rz(-1.3406164) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2778451) q[1];
sx q[1];
rz(-1.8238457) q[1];
sx q[1];
rz(2.7051198) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9704814) q[3];
sx q[3];
rz(-1.7198472) q[3];
sx q[3];
rz(0.51250729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2350601) q[2];
sx q[2];
rz(-1.1842714) q[2];
sx q[2];
rz(0.85901421) q[2];
rz(0.61339316) q[3];
sx q[3];
rz(-0.20039138) q[3];
sx q[3];
rz(1.332351) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61619806) q[0];
sx q[0];
rz(-1.2735676) q[0];
sx q[0];
rz(-1.6599576) q[0];
rz(-2.0741277) q[1];
sx q[1];
rz(-0.72507247) q[1];
sx q[1];
rz(-1.3190837) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7049371) q[0];
sx q[0];
rz(-2.4709513) q[0];
sx q[0];
rz(1.4192307) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6254933) q[2];
sx q[2];
rz(-1.2428428) q[2];
sx q[2];
rz(2.536377) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6665277) q[1];
sx q[1];
rz(-1.804978) q[1];
sx q[1];
rz(0.98824309) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33101636) q[3];
sx q[3];
rz(-0.3754456) q[3];
sx q[3];
rz(-2.7956243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2655502) q[2];
sx q[2];
rz(-0.37028131) q[2];
sx q[2];
rz(-0.041157095) q[2];
rz(2.0187078) q[3];
sx q[3];
rz(-1.6582146) q[3];
sx q[3];
rz(2.6179204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6599643) q[0];
sx q[0];
rz(-2.1631503) q[0];
sx q[0];
rz(1.6690669) q[0];
rz(-0.72498471) q[1];
sx q[1];
rz(-1.9703194) q[1];
sx q[1];
rz(0.57458007) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0313505) q[0];
sx q[0];
rz(-1.5230122) q[0];
sx q[0];
rz(2.8035473) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8594744) q[2];
sx q[2];
rz(-1.1638255) q[2];
sx q[2];
rz(0.59772432) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7704627) q[1];
sx q[1];
rz(-1.5411301) q[1];
sx q[1];
rz(-0.71843036) q[1];
rz(2.2584142) q[3];
sx q[3];
rz(-1.4559345) q[3];
sx q[3];
rz(-2.3695947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.64728111) q[2];
sx q[2];
rz(-0.24792555) q[2];
sx q[2];
rz(-0.9606804) q[2];
rz(-2.658355) q[3];
sx q[3];
rz(-1.5944696) q[3];
sx q[3];
rz(1.5794992) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0110589) q[0];
sx q[0];
rz(-0.49376765) q[0];
sx q[0];
rz(0.44179398) q[0];
rz(-2.4590625) q[1];
sx q[1];
rz(-1.6923169) q[1];
sx q[1];
rz(2.0880879) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5614612) q[0];
sx q[0];
rz(-0.86505689) q[0];
sx q[0];
rz(2.8012026) q[0];
rz(-pi) q[1];
rz(1.7948661) q[2];
sx q[2];
rz(-2.7599979) q[2];
sx q[2];
rz(-0.14118186) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7451585) q[1];
sx q[1];
rz(-0.72275439) q[1];
sx q[1];
rz(2.0682343) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9656762) q[3];
sx q[3];
rz(-2.7421167) q[3];
sx q[3];
rz(-0.3829869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0423476) q[2];
sx q[2];
rz(-1.2688364) q[2];
sx q[2];
rz(2.919) q[2];
rz(-1.7706722) q[3];
sx q[3];
rz(-1.7226847) q[3];
sx q[3];
rz(-0.62209779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.901578) q[0];
sx q[0];
rz(-1.0116901) q[0];
sx q[0];
rz(-2.1160175) q[0];
rz(1.1657794) q[1];
sx q[1];
rz(-0.60973254) q[1];
sx q[1];
rz(2.9235074) q[1];
rz(-0.54033242) q[2];
sx q[2];
rz(-1.9987543) q[2];
sx q[2];
rz(1.0812154) q[2];
rz(-3.141116) q[3];
sx q[3];
rz(-1.8924185) q[3];
sx q[3];
rz(-2.3782829) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
