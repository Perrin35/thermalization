OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.11290057) q[0];
sx q[0];
rz(-0.24554645) q[0];
sx q[0];
rz(3.088933) q[0];
rz(-2.0242937) q[1];
sx q[1];
rz(2.8007562) q[1];
sx q[1];
rz(10.513289) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3153305) q[0];
sx q[0];
rz(-1.1754987) q[0];
sx q[0];
rz(3.0733103) q[0];
rz(-pi) q[1];
rz(0.11012385) q[2];
sx q[2];
rz(-1.7118258) q[2];
sx q[2];
rz(1.5361496) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.12456914) q[1];
sx q[1];
rz(-2.2042511) q[1];
sx q[1];
rz(1.775746) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5616712) q[3];
sx q[3];
rz(-0.41968981) q[3];
sx q[3];
rz(1.7407103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6344305) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(2.5334899) q[2];
rz(2.0173343) q[3];
sx q[3];
rz(-1.2048771) q[3];
sx q[3];
rz(0.89157909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
rz(-2.7267777) q[0];
sx q[0];
rz(-1.5685273) q[0];
sx q[0];
rz(-2.538105) q[0];
rz(-0.51015774) q[1];
sx q[1];
rz(-0.77168232) q[1];
sx q[1];
rz(-1.0268432) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0870846) q[0];
sx q[0];
rz(-1.3985976) q[0];
sx q[0];
rz(1.4132176) q[0];
x q[1];
rz(-0.85058327) q[2];
sx q[2];
rz(-0.89909485) q[2];
sx q[2];
rz(-1.4631866) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4812117) q[1];
sx q[1];
rz(-0.91970316) q[1];
sx q[1];
rz(0.91459413) q[1];
rz(-pi) q[2];
rz(0.53431003) q[3];
sx q[3];
rz(-0.97917334) q[3];
sx q[3];
rz(-3.001377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5828731) q[2];
sx q[2];
rz(-1.1483973) q[2];
sx q[2];
rz(2.4951475) q[2];
rz(1.8630113) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(1.4518552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2090787) q[0];
sx q[0];
rz(-0.13020733) q[0];
sx q[0];
rz(1.5721488) q[0];
rz(-2.7945844) q[1];
sx q[1];
rz(-1.6667112) q[1];
sx q[1];
rz(-2.2775547) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5310609) q[0];
sx q[0];
rz(-2.5342434) q[0];
sx q[0];
rz(0.37215287) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7059495) q[2];
sx q[2];
rz(-1.3416939) q[2];
sx q[2];
rz(1.4080017) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7155418) q[1];
sx q[1];
rz(-2.155039) q[1];
sx q[1];
rz(0.76202315) q[1];
rz(-pi) q[2];
rz(1.044831) q[3];
sx q[3];
rz(-1.6623896) q[3];
sx q[3];
rz(2.4909702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0916834) q[2];
sx q[2];
rz(-0.20541643) q[2];
sx q[2];
rz(-1.4570215) q[2];
rz(1.0157887) q[3];
sx q[3];
rz(-0.53272811) q[3];
sx q[3];
rz(-2.8252025) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9006573) q[0];
sx q[0];
rz(-2.7641986) q[0];
sx q[0];
rz(1.578791) q[0];
rz(-1.0904301) q[1];
sx q[1];
rz(-2.5999887) q[1];
sx q[1];
rz(1.9900367) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3061515) q[0];
sx q[0];
rz(-2.5130898) q[0];
sx q[0];
rz(-2.7449904) q[0];
rz(-pi) q[1];
rz(-2.0456373) q[2];
sx q[2];
rz(-0.86855382) q[2];
sx q[2];
rz(-2.4202731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87428625) q[1];
sx q[1];
rz(-1.5360254) q[1];
sx q[1];
rz(-1.0465996) q[1];
rz(-pi) q[2];
rz(3.1138469) q[3];
sx q[3];
rz(-0.96409269) q[3];
sx q[3];
rz(-1.1535899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.84732032) q[2];
sx q[2];
rz(-0.5449833) q[2];
sx q[2];
rz(1.6418183) q[2];
rz(3.0070987) q[3];
sx q[3];
rz(-1.8650863) q[3];
sx q[3];
rz(-2.3597778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0624369) q[0];
sx q[0];
rz(-2.2298614) q[0];
sx q[0];
rz(0.4656747) q[0];
rz(1.8648719) q[1];
sx q[1];
rz(-1.0036422) q[1];
sx q[1];
rz(2.1086955) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5860234) q[0];
sx q[0];
rz(-2.9831671) q[0];
sx q[0];
rz(1.6371631) q[0];
x q[1];
rz(-3.0014702) q[2];
sx q[2];
rz(-1.6701785) q[2];
sx q[2];
rz(-1.4750208) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.031449854) q[1];
sx q[1];
rz(-1.8217161) q[1];
sx q[1];
rz(0.67278426) q[1];
x q[2];
rz(0.38893969) q[3];
sx q[3];
rz(-2.0802214) q[3];
sx q[3];
rz(-3.096599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.5501962) q[2];
sx q[2];
rz(-1.952221) q[2];
sx q[2];
rz(0.93185321) q[2];
rz(3.0469117) q[3];
sx q[3];
rz(-0.90016142) q[3];
sx q[3];
rz(1.8315106) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69022995) q[0];
sx q[0];
rz(-2.9424423) q[0];
sx q[0];
rz(0.1061826) q[0];
rz(2.5871318) q[1];
sx q[1];
rz(-0.78840557) q[1];
sx q[1];
rz(1.6000481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0672749) q[0];
sx q[0];
rz(-1.5948926) q[0];
sx q[0];
rz(0.39053183) q[0];
rz(1.1136069) q[2];
sx q[2];
rz(-1.3234252) q[2];
sx q[2];
rz(-2.8988225) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.2682802) q[1];
sx q[1];
rz(-2.2619216) q[1];
sx q[1];
rz(-1.4119488) q[1];
rz(-1.5149045) q[3];
sx q[3];
rz(-1.8168746) q[3];
sx q[3];
rz(2.6757307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.23325486) q[2];
sx q[2];
rz(-2.4007863) q[2];
sx q[2];
rz(-1.4976658) q[2];
rz(-2.6805367) q[3];
sx q[3];
rz(-2.4888829) q[3];
sx q[3];
rz(1.8259995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4696362) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(1.0787971) q[0];
rz(0.59393334) q[1];
sx q[1];
rz(-0.97336951) q[1];
sx q[1];
rz(2.914391) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35543252) q[0];
sx q[0];
rz(-0.18635145) q[0];
sx q[0];
rz(0.92599286) q[0];
rz(1.960177) q[2];
sx q[2];
rz(-2.7837278) q[2];
sx q[2];
rz(3.1111382) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6727184) q[1];
sx q[1];
rz(-1.0212082) q[1];
sx q[1];
rz(-2.7764715) q[1];
rz(1.9428857) q[3];
sx q[3];
rz(-1.477302) q[3];
sx q[3];
rz(-2.1402511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6195153) q[2];
sx q[2];
rz(-2.7232813) q[2];
sx q[2];
rz(2.412001) q[2];
rz(-1.6884165) q[3];
sx q[3];
rz(-1.6875234) q[3];
sx q[3];
rz(-0.61354536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5486117) q[0];
sx q[0];
rz(-2.225003) q[0];
sx q[0];
rz(2.7899637) q[0];
rz(0.31373203) q[1];
sx q[1];
rz(-1.2064563) q[1];
sx q[1];
rz(-1.3765913) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2257654) q[0];
sx q[0];
rz(-1.5236014) q[0];
sx q[0];
rz(2.7412358) q[0];
rz(0.11560924) q[2];
sx q[2];
rz(-1.7358276) q[2];
sx q[2];
rz(0.045265667) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.59960466) q[1];
sx q[1];
rz(-1.7270589) q[1];
sx q[1];
rz(0.57078741) q[1];
rz(-0.099011856) q[3];
sx q[3];
rz(-1.9388698) q[3];
sx q[3];
rz(-1.2457445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54479638) q[2];
sx q[2];
rz(-2.4185541) q[2];
sx q[2];
rz(1.5210305) q[2];
rz(0.14791402) q[3];
sx q[3];
rz(-1.3444129) q[3];
sx q[3];
rz(1.773905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3677597) q[0];
sx q[0];
rz(-1.8748883) q[0];
sx q[0];
rz(-1.2441147) q[0];
rz(-1.7077839) q[1];
sx q[1];
rz(-1.5807187) q[1];
sx q[1];
rz(-0.22020766) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51802902) q[0];
sx q[0];
rz(-1.9401778) q[0];
sx q[0];
rz(2.2165038) q[0];
x q[1];
rz(2.979109) q[2];
sx q[2];
rz(-1.1097317) q[2];
sx q[2];
rz(1.7562255) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23277321) q[1];
sx q[1];
rz(-0.72130433) q[1];
sx q[1];
rz(-2.2386293) q[1];
rz(-pi) q[2];
rz(0.76311402) q[3];
sx q[3];
rz(-0.13555758) q[3];
sx q[3];
rz(-0.28597304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9571017) q[2];
sx q[2];
rz(-1.1471006) q[2];
sx q[2];
rz(1.6415049) q[2];
rz(3.0575276) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(0.070092289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4641007) q[0];
sx q[0];
rz(-2.3832432) q[0];
sx q[0];
rz(2.7229934) q[0];
rz(0.94340008) q[1];
sx q[1];
rz(-1.6523596) q[1];
sx q[1];
rz(2.8387866) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7092512) q[0];
sx q[0];
rz(-2.925207) q[0];
sx q[0];
rz(-0.88309137) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3934028) q[2];
sx q[2];
rz(-2.8837969) q[2];
sx q[2];
rz(2.2141544) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9040519) q[1];
sx q[1];
rz(-1.6142323) q[1];
sx q[1];
rz(0.082708184) q[1];
x q[2];
rz(-1.915384) q[3];
sx q[3];
rz(-1.4610963) q[3];
sx q[3];
rz(-1.6747703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4960949) q[2];
sx q[2];
rz(-0.87254137) q[2];
sx q[2];
rz(-0.75054753) q[2];
rz(1.6802855) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(-2.6245978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0208329) q[0];
sx q[0];
rz(-1.5789541) q[0];
sx q[0];
rz(-1.5776237) q[0];
rz(-1.9137406) q[1];
sx q[1];
rz(-2.6422983) q[1];
sx q[1];
rz(-1.9849389) q[1];
rz(-0.56671178) q[2];
sx q[2];
rz(-1.4162345) q[2];
sx q[2];
rz(-0.92879374) q[2];
rz(1.0629366) q[3];
sx q[3];
rz(-2.776317) q[3];
sx q[3];
rz(-1.1767514) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
