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
rz(2.8960462) q[0];
sx q[0];
rz(9.4774376) q[0];
rz(-2.0242937) q[1];
sx q[1];
rz(2.8007562) q[1];
sx q[1];
rz(10.513289) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4910925) q[0];
sx q[0];
rz(-2.740747) q[0];
sx q[0];
rz(1.7328788) q[0];
rz(-pi) q[1];
rz(-2.229548) q[2];
sx q[2];
rz(-0.17870644) q[2];
sx q[2];
rz(-2.2026874) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6793078) q[1];
sx q[1];
rz(-2.4801697) q[1];
sx q[1];
rz(0.2703305) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9904706) q[3];
sx q[3];
rz(-1.5670781) q[3];
sx q[3];
rz(-2.9633455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5071621) q[2];
sx q[2];
rz(-1.3619225) q[2];
sx q[2];
rz(-2.5334899) q[2];
rz(-1.1242584) q[3];
sx q[3];
rz(-1.9367155) q[3];
sx q[3];
rz(-0.89157909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41481498) q[0];
sx q[0];
rz(-1.5685273) q[0];
sx q[0];
rz(-2.538105) q[0];
rz(-2.6314349) q[1];
sx q[1];
rz(-0.77168232) q[1];
sx q[1];
rz(-2.1147494) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3393766) q[0];
sx q[0];
rz(-0.23288865) q[0];
sx q[0];
rz(2.4075137) q[0];
x q[1];
rz(0.85058327) q[2];
sx q[2];
rz(-0.89909485) q[2];
sx q[2];
rz(-1.6784061) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4812117) q[1];
sx q[1];
rz(-0.91970316) q[1];
sx q[1];
rz(-2.2269985) q[1];
rz(-pi) q[2];
x q[2];
rz(0.90791865) q[3];
sx q[3];
rz(-2.0072122) q[3];
sx q[3];
rz(1.7493356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5828731) q[2];
sx q[2];
rz(-1.9931953) q[2];
sx q[2];
rz(2.4951475) q[2];
rz(-1.8630113) q[3];
sx q[3];
rz(-2.1918112) q[3];
sx q[3];
rz(-1.4518552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.2090787) q[0];
sx q[0];
rz(-3.0113853) q[0];
sx q[0];
rz(1.5694438) q[0];
rz(0.34700829) q[1];
sx q[1];
rz(-1.4748814) q[1];
sx q[1];
rz(2.2775547) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5310609) q[0];
sx q[0];
rz(-0.60734925) q[0];
sx q[0];
rz(-0.37215287) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7059495) q[2];
sx q[2];
rz(-1.7998988) q[2];
sx q[2];
rz(1.4080017) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5123295) q[1];
sx q[1];
rz(-2.1844668) q[1];
sx q[1];
rz(-0.83028173) q[1];
rz(1.044831) q[3];
sx q[3];
rz(-1.4792031) q[3];
sx q[3];
rz(0.65062246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.049909264) q[2];
sx q[2];
rz(-2.9361762) q[2];
sx q[2];
rz(1.6845711) q[2];
rz(-1.0157887) q[3];
sx q[3];
rz(-2.6088645) q[3];
sx q[3];
rz(-2.8252025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24093534) q[0];
sx q[0];
rz(-2.7641986) q[0];
sx q[0];
rz(-1.578791) q[0];
rz(-2.0511625) q[1];
sx q[1];
rz(-0.54160392) q[1];
sx q[1];
rz(-1.1515559) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3061515) q[0];
sx q[0];
rz(-2.5130898) q[0];
sx q[0];
rz(-0.39660227) q[0];
x q[1];
rz(0.76049034) q[2];
sx q[2];
rz(-1.2142688) q[2];
sx q[2];
rz(-0.52887756) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4249866) q[1];
sx q[1];
rz(-1.046949) q[1];
sx q[1];
rz(-3.1014342) q[1];
x q[2];
rz(1.5308446) q[3];
sx q[3];
rz(-0.60725799) q[3];
sx q[3];
rz(1.1049529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2942723) q[2];
sx q[2];
rz(-0.5449833) q[2];
sx q[2];
rz(1.6418183) q[2];
rz(-3.0070987) q[3];
sx q[3];
rz(-1.8650863) q[3];
sx q[3];
rz(2.3597778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0791557) q[0];
sx q[0];
rz(-0.9117313) q[0];
sx q[0];
rz(-0.4656747) q[0];
rz(-1.8648719) q[1];
sx q[1];
rz(-1.0036422) q[1];
sx q[1];
rz(1.0328971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55556923) q[0];
sx q[0];
rz(-0.15842552) q[0];
sx q[0];
rz(-1.6371631) q[0];
rz(-1.6711556) q[2];
sx q[2];
rz(-1.4313698) q[2];
sx q[2];
rz(-0.1097691) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.1101428) q[1];
sx q[1];
rz(-1.8217161) q[1];
sx q[1];
rz(-0.67278426) q[1];
rz(1.0276658) q[3];
sx q[3];
rz(-1.9082205) q[3];
sx q[3];
rz(1.8130482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.59139645) q[2];
sx q[2];
rz(-1.1893716) q[2];
sx q[2];
rz(-0.93185321) q[2];
rz(3.0469117) q[3];
sx q[3];
rz(-2.2414312) q[3];
sx q[3];
rz(1.3100821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69022995) q[0];
sx q[0];
rz(-0.19915038) q[0];
sx q[0];
rz(-0.1061826) q[0];
rz(-2.5871318) q[1];
sx q[1];
rz(-2.3531871) q[1];
sx q[1];
rz(-1.5415446) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6351955) q[0];
sx q[0];
rz(-1.180384) q[0];
sx q[0];
rz(1.5447389) q[0];
rz(1.0511398) q[2];
sx q[2];
rz(-0.51562947) q[2];
sx q[2];
rz(-1.7898498) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2682802) q[1];
sx q[1];
rz(-2.2619216) q[1];
sx q[1];
rz(-1.4119488) q[1];
rz(1.5149045) q[3];
sx q[3];
rz(-1.8168746) q[3];
sx q[3];
rz(-2.6757307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9083378) q[2];
sx q[2];
rz(-0.7408064) q[2];
sx q[2];
rz(1.6439269) q[2];
rz(2.6805367) q[3];
sx q[3];
rz(-0.65270972) q[3];
sx q[3];
rz(-1.3155931) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6719565) q[0];
sx q[0];
rz(-0.75103432) q[0];
sx q[0];
rz(-2.0627956) q[0];
rz(0.59393334) q[1];
sx q[1];
rz(-2.1682231) q[1];
sx q[1];
rz(0.22720164) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35543252) q[0];
sx q[0];
rz(-0.18635145) q[0];
sx q[0];
rz(2.2155998) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9038796) q[2];
sx q[2];
rz(-1.7041612) q[2];
sx q[2];
rz(-1.2343209) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2367475) q[1];
sx q[1];
rz(-1.2614025) q[1];
sx q[1];
rz(0.99037804) q[1];
x q[2];
rz(-3.041275) q[3];
sx q[3];
rz(-1.2004108) q[3];
sx q[3];
rz(-2.5357192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.6195153) q[2];
sx q[2];
rz(-0.41831133) q[2];
sx q[2];
rz(0.72959161) q[2];
rz(-1.4531762) q[3];
sx q[3];
rz(-1.6875234) q[3];
sx q[3];
rz(0.61354536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59298092) q[0];
sx q[0];
rz(-2.225003) q[0];
sx q[0];
rz(-0.35162893) q[0];
rz(-0.31373203) q[1];
sx q[1];
rz(-1.9351363) q[1];
sx q[1];
rz(-1.3765913) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5976395) q[0];
sx q[0];
rz(-0.40297976) q[0];
sx q[0];
rz(-0.12059327) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1765472) q[2];
sx q[2];
rz(-0.20119431) q[2];
sx q[2];
rz(-0.66058841) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.541988) q[1];
sx q[1];
rz(-1.7270589) q[1];
sx q[1];
rz(-2.5708052) q[1];
rz(-pi) q[2];
rz(-1.8217161) q[3];
sx q[3];
rz(-2.7610215) q[3];
sx q[3];
rz(-1.6264834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.54479638) q[2];
sx q[2];
rz(-0.72303855) q[2];
sx q[2];
rz(-1.5210305) q[2];
rz(-0.14791402) q[3];
sx q[3];
rz(-1.7971797) q[3];
sx q[3];
rz(-1.3676876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77383298) q[0];
sx q[0];
rz(-1.2667043) q[0];
sx q[0];
rz(1.897478) q[0];
rz(-1.7077839) q[1];
sx q[1];
rz(-1.5807187) q[1];
sx q[1];
rz(-0.22020766) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7871683) q[0];
sx q[0];
rz(-0.97499135) q[0];
sx q[0];
rz(0.45136772) q[0];
rz(-pi) q[1];
x q[1];
rz(2.037165) q[2];
sx q[2];
rz(-1.4254071) q[2];
sx q[2];
rz(0.25823247) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0427093) q[1];
sx q[1];
rz(-1.0257057) q[1];
sx q[1];
rz(-2.6429151) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0433398) q[3];
sx q[3];
rz(-1.4772526) q[3];
sx q[3];
rz(0.52631179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.1844909) q[2];
sx q[2];
rz(-1.1471006) q[2];
sx q[2];
rz(-1.6415049) q[2];
rz(-0.084065048) q[3];
sx q[3];
rz(-1.2596687) q[3];
sx q[3];
rz(-3.0715004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67749196) q[0];
sx q[0];
rz(-0.75834948) q[0];
sx q[0];
rz(0.41859928) q[0];
rz(0.94340008) q[1];
sx q[1];
rz(-1.6523596) q[1];
sx q[1];
rz(-0.30280608) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43234149) q[0];
sx q[0];
rz(-0.21638566) q[0];
sx q[0];
rz(-0.88309137) q[0];
rz(-pi) q[1];
rz(-2.3934028) q[2];
sx q[2];
rz(-2.8837969) q[2];
sx q[2];
rz(0.92743826) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8157685) q[1];
sx q[1];
rz(-0.093397141) q[1];
sx q[1];
rz(2.6572822) q[1];
rz(-pi) q[2];
rz(-1.2262086) q[3];
sx q[3];
rz(-1.6804964) q[3];
sx q[3];
rz(1.4668224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6454978) q[2];
sx q[2];
rz(-0.87254137) q[2];
sx q[2];
rz(2.3910451) q[2];
rz(-1.4613072) q[3];
sx q[3];
rz(-0.76992005) q[3];
sx q[3];
rz(-2.6245978) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0208329) q[0];
sx q[0];
rz(-1.5626386) q[0];
sx q[0];
rz(1.563969) q[0];
rz(1.2278521) q[1];
sx q[1];
rz(-2.6422983) q[1];
sx q[1];
rz(-1.9849389) q[1];
rz(-0.56671178) q[2];
sx q[2];
rz(-1.4162345) q[2];
sx q[2];
rz(-0.92879374) q[2];
rz(-2.078656) q[3];
sx q[3];
rz(-2.776317) q[3];
sx q[3];
rz(-1.1767514) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
