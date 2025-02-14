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
rz(3.7035898) q[0];
sx q[0];
rz(9.193767) q[0];
rz(0.24569874) q[1];
sx q[1];
rz(-0.45431554) q[1];
sx q[1];
rz(1.2872202) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79345771) q[0];
sx q[0];
rz(-2.0525816) q[0];
sx q[0];
rz(2.8346377) q[0];
rz(1.8051992) q[2];
sx q[2];
rz(-1.3318828) q[2];
sx q[2];
rz(2.6221681) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8652628) q[1];
sx q[1];
rz(-1.4722605) q[1];
sx q[1];
rz(-0.16698412) q[1];
x q[2];
rz(-1.3955529) q[3];
sx q[3];
rz(-1.6449494) q[3];
sx q[3];
rz(1.7046622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40453688) q[2];
sx q[2];
rz(-2.4344567) q[2];
sx q[2];
rz(0.36049584) q[2];
rz(-2.0170085) q[3];
sx q[3];
rz(-2.093061) q[3];
sx q[3];
rz(1.8726965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8812113) q[0];
sx q[0];
rz(-2.9660048) q[0];
sx q[0];
rz(-2.4847109) q[0];
rz(-0.33292133) q[1];
sx q[1];
rz(-1.0924783) q[1];
sx q[1];
rz(-2.2064256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29657179) q[0];
sx q[0];
rz(-1.6753418) q[0];
sx q[0];
rz(1.5403032) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33719535) q[2];
sx q[2];
rz(-0.27183485) q[2];
sx q[2];
rz(1.4120917) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6015687) q[1];
sx q[1];
rz(-1.717713) q[1];
sx q[1];
rz(1.5308892) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7788497) q[3];
sx q[3];
rz(-1.9192154) q[3];
sx q[3];
rz(-1.5294242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3071345) q[2];
sx q[2];
rz(-1.6259401) q[2];
sx q[2];
rz(1.2163986) q[2];
rz(2.6136716) q[3];
sx q[3];
rz(-2.1025751) q[3];
sx q[3];
rz(1.2549887) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5288178) q[0];
sx q[0];
rz(-2.3154494) q[0];
sx q[0];
rz(-0.58498996) q[0];
rz(3.0168369) q[1];
sx q[1];
rz(-2.545732) q[1];
sx q[1];
rz(1.6927208) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3901236) q[0];
sx q[0];
rz(-1.5678143) q[0];
sx q[0];
rz(-3.1389159) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5572433) q[2];
sx q[2];
rz(-1.4986191) q[2];
sx q[2];
rz(2.7685795) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.3347149) q[1];
sx q[1];
rz(-0.92744614) q[1];
sx q[1];
rz(-2.592464) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.65421974) q[3];
sx q[3];
rz(-2.7252135) q[3];
sx q[3];
rz(0.49921303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.8831732) q[2];
sx q[2];
rz(-2.1397739) q[2];
sx q[2];
rz(-1.6955356) q[2];
rz(-0.7575194) q[3];
sx q[3];
rz(-1.5816553) q[3];
sx q[3];
rz(2.3022046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3676753) q[0];
sx q[0];
rz(-0.92905074) q[0];
sx q[0];
rz(-1.0876592) q[0];
rz(-2.7663973) q[1];
sx q[1];
rz(-2.0162069) q[1];
sx q[1];
rz(-0.96022022) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86726928) q[0];
sx q[0];
rz(-1.835863) q[0];
sx q[0];
rz(0.041170711) q[0];
rz(-1.6702508) q[2];
sx q[2];
rz(-1.5410454) q[2];
sx q[2];
rz(0.1386252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.5615446) q[1];
sx q[1];
rz(-0.64148884) q[1];
sx q[1];
rz(-2.3381691) q[1];
rz(0.33444114) q[3];
sx q[3];
rz(-2.8570647) q[3];
sx q[3];
rz(-2.4872125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9356161) q[2];
sx q[2];
rz(-1.2674067) q[2];
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
x q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.977026) q[0];
sx q[0];
rz(-2.470546) q[0];
sx q[0];
rz(2.5373051) q[0];
rz(-1.9970278) q[1];
sx q[1];
rz(-1.4860169) q[1];
sx q[1];
rz(2.7191275) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6910962) q[0];
sx q[0];
rz(-1.7393149) q[0];
sx q[0];
rz(-0.082534747) q[0];
x q[1];
rz(2.3461012) q[2];
sx q[2];
rz(-2.4078712) q[2];
sx q[2];
rz(2.7249641) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.996811) q[1];
sx q[1];
rz(-1.6429122) q[1];
sx q[1];
rz(2.4516289) q[1];
rz(-pi) q[2];
rz(2.7606008) q[3];
sx q[3];
rz(-1.89011) q[3];
sx q[3];
rz(1.3407941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3182688) q[2];
sx q[2];
rz(-0.6508998) q[2];
sx q[2];
rz(-0.65625119) q[2];
rz(-1.5308135) q[3];
sx q[3];
rz(-1.2698413) q[3];
sx q[3];
rz(-0.58073616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73606473) q[0];
sx q[0];
rz(-1.0117714) q[0];
sx q[0];
rz(0.62498012) q[0];
rz(0.75195733) q[1];
sx q[1];
rz(-2.0729005) q[1];
sx q[1];
rz(1.6065074) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4055734) q[0];
sx q[0];
rz(-1.6200119) q[0];
sx q[0];
rz(0.28926138) q[0];
rz(-pi) q[1];
rz(-2.4322832) q[2];
sx q[2];
rz(-2.0481718) q[2];
sx q[2];
rz(-0.8698744) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8934763) q[1];
sx q[1];
rz(-1.3424338) q[1];
sx q[1];
rz(0.62947692) q[1];
rz(-pi) q[2];
x q[2];
rz(0.84523369) q[3];
sx q[3];
rz(-0.74016011) q[3];
sx q[3];
rz(-2.6773767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.1890761) q[2];
sx q[2];
rz(-1.357888) q[2];
sx q[2];
rz(0.9291741) q[2];
rz(2.3164228) q[3];
sx q[3];
rz(-1.5114096) q[3];
sx q[3];
rz(-1.1024124) q[3];
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
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55649844) q[0];
sx q[0];
rz(-2.3346021) q[0];
sx q[0];
rz(1.3744542) q[0];
rz(-0.51042405) q[1];
sx q[1];
rz(-2.3511032) q[1];
sx q[1];
rz(2.5183831) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1081559) q[0];
sx q[0];
rz(-1.2386453) q[0];
sx q[0];
rz(0.16031127) q[0];
rz(-pi) q[1];
rz(-0.50919598) q[2];
sx q[2];
rz(-1.3993345) q[2];
sx q[2];
rz(0.96656884) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0859061) q[1];
sx q[1];
rz(-1.4594363) q[1];
sx q[1];
rz(1.423905) q[1];
rz(-0.74759746) q[3];
sx q[3];
rz(-1.4211402) q[3];
sx q[3];
rz(-2.8437994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0906543) q[2];
sx q[2];
rz(-2.3263558) q[2];
sx q[2];
rz(1.1673002) q[2];
rz(-0.022631571) q[3];
sx q[3];
rz(-0.52112094) q[3];
sx q[3];
rz(-2.0126655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24564329) q[0];
sx q[0];
rz(-1.1351981) q[0];
sx q[0];
rz(-1.0173215) q[0];
rz(1.7474878) q[1];
sx q[1];
rz(-2.930495) q[1];
sx q[1];
rz(2.5140433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8857972) q[0];
sx q[0];
rz(-1.1748519) q[0];
sx q[0];
rz(-0.93314472) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4623221) q[2];
sx q[2];
rz(-1.5049266) q[2];
sx q[2];
rz(3.1312453) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.095762091) q[1];
sx q[1];
rz(-1.7667637) q[1];
sx q[1];
rz(-0.24230559) q[1];
x q[2];
rz(1.9110319) q[3];
sx q[3];
rz(-2.0215394) q[3];
sx q[3];
rz(-2.4204202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5158186) q[2];
sx q[2];
rz(-0.63483441) q[2];
sx q[2];
rz(-2.7308357) q[2];
rz(0.050203236) q[3];
sx q[3];
rz(-1.2529195) q[3];
sx q[3];
rz(0.075411782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.582616) q[0];
sx q[0];
rz(-0.89685431) q[0];
sx q[0];
rz(2.8726752) q[0];
rz(2.8315663) q[1];
sx q[1];
rz(-1.0870442) q[1];
sx q[1];
rz(-0.1300098) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0144671) q[0];
sx q[0];
rz(-1.3499182) q[0];
sx q[0];
rz(1.2400024) q[0];
x q[1];
rz(-1.9275877) q[2];
sx q[2];
rz(-1.2399106) q[2];
sx q[2];
rz(2.7499104) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.074689493) q[1];
sx q[1];
rz(-2.5554511) q[1];
sx q[1];
rz(-0.2939923) q[1];
x q[2];
rz(-0.89270182) q[3];
sx q[3];
rz(-1.1696276) q[3];
sx q[3];
rz(0.65142163) q[3];
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
rz(-1.6377595) q[3];
sx q[3];
rz(-1.8042754) q[3];
sx q[3];
rz(-0.79969978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0132975) q[0];
sx q[0];
rz(-0.18823637) q[0];
sx q[0];
rz(-2.6085594) q[0];
rz(-0.036272613) q[1];
sx q[1];
rz(-2.3590922) q[1];
sx q[1];
rz(-1.689555) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2672469) q[0];
sx q[0];
rz(-1.1935992) q[0];
sx q[0];
rz(-0.95893559) q[0];
rz(-pi) q[1];
rz(-1.6797941) q[2];
sx q[2];
rz(-1.8882635) q[2];
sx q[2];
rz(-3.1239117) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.8126077) q[1];
sx q[1];
rz(-1.2115062) q[1];
sx q[1];
rz(1.6027662) q[1];
x q[2];
rz(2.6759869) q[3];
sx q[3];
rz(-1.4911663) q[3];
sx q[3];
rz(-1.7871943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.71296802) q[2];
sx q[2];
rz(-2.4098101) q[2];
sx q[2];
rz(0.0326322) q[2];
rz(0.84515682) q[3];
sx q[3];
rz(-1.9149575) q[3];
sx q[3];
rz(-1.0802065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7964771) q[0];
sx q[0];
rz(-2.4892172) q[0];
sx q[0];
rz(-0.31443483) q[0];
rz(0.79700094) q[1];
sx q[1];
rz(-0.92756699) q[1];
sx q[1];
rz(-3.0753593) q[1];
rz(-2.6061229) q[2];
sx q[2];
rz(-1.3146123) q[2];
sx q[2];
rz(-1.8563369) q[2];
rz(1.6975523) q[3];
sx q[3];
rz(-0.76382617) q[3];
sx q[3];
rz(2.7504117) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
