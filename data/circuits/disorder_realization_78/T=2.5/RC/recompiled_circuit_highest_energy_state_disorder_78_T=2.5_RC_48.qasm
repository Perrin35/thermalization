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
rz(0.084963381) q[0];
sx q[0];
rz(-2.8391916) q[0];
sx q[0];
rz(3.1095355) q[0];
rz(-0.0078460296) q[1];
sx q[1];
rz(-0.50385952) q[1];
sx q[1];
rz(-0.75262466) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1330452) q[0];
sx q[0];
rz(-2.7602782) q[0];
sx q[0];
rz(-0.56579907) q[0];
rz(-pi) q[1];
x q[1];
rz(2.339509) q[2];
sx q[2];
rz(-2.6418016) q[2];
sx q[2];
rz(1.0700723) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9702985) q[1];
sx q[1];
rz(-2.7134905) q[1];
sx q[1];
rz(-2.5004205) q[1];
x q[2];
rz(2.0977375) q[3];
sx q[3];
rz(-2.733272) q[3];
sx q[3];
rz(2.4957617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.96693119) q[2];
sx q[2];
rz(-1.9661247) q[2];
sx q[2];
rz(2.3215129) q[2];
rz(-2.7005633) q[3];
sx q[3];
rz(-1.1038154) q[3];
sx q[3];
rz(0.57344121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.012861982) q[0];
sx q[0];
rz(-2.2282889) q[0];
sx q[0];
rz(-2.7313857) q[0];
rz(-0.38093105) q[1];
sx q[1];
rz(-1.5313287) q[1];
sx q[1];
rz(1.358323) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5429496) q[0];
sx q[0];
rz(-2.1457167) q[0];
sx q[0];
rz(0.47358124) q[0];
rz(-1.0274506) q[2];
sx q[2];
rz(-2.0110235) q[2];
sx q[2];
rz(1.5981975) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.0002901) q[1];
sx q[1];
rz(-1.4030255) q[1];
sx q[1];
rz(3.0649158) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15492691) q[3];
sx q[3];
rz(-2.2048031) q[3];
sx q[3];
rz(0.67544322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7831948) q[2];
sx q[2];
rz(-0.40893778) q[2];
sx q[2];
rz(-2.6386293) q[2];
rz(-1.8424235) q[3];
sx q[3];
rz(-2.3830569) q[3];
sx q[3];
rz(2.9097596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4101039) q[0];
sx q[0];
rz(-2.6298611) q[0];
sx q[0];
rz(-0.9758392) q[0];
rz(0.93636912) q[1];
sx q[1];
rz(-1.5872637) q[1];
sx q[1];
rz(-0.13793129) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1722522) q[0];
sx q[0];
rz(-1.9346049) q[0];
sx q[0];
rz(1.0602289) q[0];
rz(0.27020578) q[2];
sx q[2];
rz(-0.41105726) q[2];
sx q[2];
rz(0.11693987) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.960818) q[1];
sx q[1];
rz(-1.7021693) q[1];
sx q[1];
rz(-2.2714494) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6343752) q[3];
sx q[3];
rz(-2.6196369) q[3];
sx q[3];
rz(-1.0591266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.58933538) q[2];
sx q[2];
rz(-1.6103123) q[2];
sx q[2];
rz(-1.662558) q[2];
rz(0.79834437) q[3];
sx q[3];
rz(-2.2975497) q[3];
sx q[3];
rz(-3.121283) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8524356) q[0];
sx q[0];
rz(-0.11141369) q[0];
sx q[0];
rz(-1.0668466) q[0];
rz(1.5486708) q[1];
sx q[1];
rz(-1.2499864) q[1];
sx q[1];
rz(-2.7222395) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3187554) q[0];
sx q[0];
rz(-1.4117804) q[0];
sx q[0];
rz(1.8029638) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.1544756) q[2];
sx q[2];
rz(-2.6082268) q[2];
sx q[2];
rz(1.5325002) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90430656) q[1];
sx q[1];
rz(-1.6349873) q[1];
sx q[1];
rz(1.9285081) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0859503) q[3];
sx q[3];
rz(-1.8797698) q[3];
sx q[3];
rz(3.0687208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.35820094) q[2];
sx q[2];
rz(-2.6322067) q[2];
sx q[2];
rz(1.5436714) q[2];
rz(-1.2265497) q[3];
sx q[3];
rz(-1.9251325) q[3];
sx q[3];
rz(1.8515057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2655547) q[0];
sx q[0];
rz(-2.4198678) q[0];
sx q[0];
rz(-0.75620404) q[0];
rz(1.2528231) q[1];
sx q[1];
rz(-0.93106657) q[1];
sx q[1];
rz(1.4160215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9276816) q[0];
sx q[0];
rz(-2.1577303) q[0];
sx q[0];
rz(1.8326493) q[0];
rz(0.3438832) q[2];
sx q[2];
rz(-1.4522219) q[2];
sx q[2];
rz(2.8176475) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6452613) q[1];
sx q[1];
rz(-1.0320391) q[1];
sx q[1];
rz(1.5625728) q[1];
x q[2];
rz(1.7613714) q[3];
sx q[3];
rz(-1.9280199) q[3];
sx q[3];
rz(0.009454184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.34117928) q[2];
sx q[2];
rz(-0.74113733) q[2];
sx q[2];
rz(0.18079147) q[2];
rz(1.4888658) q[3];
sx q[3];
rz(-1.3988262) q[3];
sx q[3];
rz(1.5159877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(-2.4329231) q[0];
sx q[0];
rz(-2.0661418) q[0];
sx q[0];
rz(-0.80023009) q[0];
rz(-0.284614) q[1];
sx q[1];
rz(-2.0814643) q[1];
sx q[1];
rz(-0.12399331) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.02888814) q[0];
sx q[0];
rz(-1.7005973) q[0];
sx q[0];
rz(3.1129254) q[0];
rz(3.0818865) q[2];
sx q[2];
rz(-1.8054031) q[2];
sx q[2];
rz(1.4690746) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2468894) q[1];
sx q[1];
rz(-1.2661625) q[1];
sx q[1];
rz(-1.6564293) q[1];
x q[2];
rz(0.37185566) q[3];
sx q[3];
rz(-2.1346492) q[3];
sx q[3];
rz(1.9052037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.93512145) q[2];
sx q[2];
rz(-1.4149041) q[2];
sx q[2];
rz(-3.0653811) q[2];
rz(1.8501836) q[3];
sx q[3];
rz(-2.0468057) q[3];
sx q[3];
rz(2.5230303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9718219) q[0];
sx q[0];
rz(-1.562028) q[0];
sx q[0];
rz(-2.3845657) q[0];
rz(2.3454759) q[1];
sx q[1];
rz(-1.4069822) q[1];
sx q[1];
rz(-0.98181358) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7863203) q[0];
sx q[0];
rz(-2.2844446) q[0];
sx q[0];
rz(2.7631604) q[0];
x q[1];
rz(-1.4635529) q[2];
sx q[2];
rz(-1.5011884) q[2];
sx q[2];
rz(2.7977242) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.71926861) q[1];
sx q[1];
rz(-2.0403721) q[1];
sx q[1];
rz(2.6828241) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1037681) q[3];
sx q[3];
rz(-0.31222725) q[3];
sx q[3];
rz(-1.3749773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.41165274) q[2];
sx q[2];
rz(-2.3306658) q[2];
sx q[2];
rz(-2.5977503) q[2];
rz(0.020921556) q[3];
sx q[3];
rz(-1.436751) q[3];
sx q[3];
rz(2.3232443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92886096) q[0];
sx q[0];
rz(-0.91101557) q[0];
sx q[0];
rz(-0.62335706) q[0];
rz(2.8202672) q[1];
sx q[1];
rz(-0.61505452) q[1];
sx q[1];
rz(-2.8772433) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8835444) q[0];
sx q[0];
rz(-0.86853456) q[0];
sx q[0];
rz(2.8625898) q[0];
rz(-1.4521763) q[2];
sx q[2];
rz(-2.8037694) q[2];
sx q[2];
rz(-2.0488957) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26786727) q[1];
sx q[1];
rz(-1.2920989) q[1];
sx q[1];
rz(-0.81588094) q[1];
x q[2];
rz(0.13295138) q[3];
sx q[3];
rz(-1.8859409) q[3];
sx q[3];
rz(0.12115762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.136772) q[2];
sx q[2];
rz(-0.68237582) q[2];
sx q[2];
rz(0.47478673) q[2];
rz(-0.36561203) q[3];
sx q[3];
rz(-2.6545299) q[3];
sx q[3];
rz(0.18317187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8356165) q[0];
sx q[0];
rz(-2.7662179) q[0];
sx q[0];
rz(-2.6557652) q[0];
rz(-1.0659418) q[1];
sx q[1];
rz(-1.7168047) q[1];
sx q[1];
rz(-2.8071383) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4518648) q[0];
sx q[0];
rz(-2.7948927) q[0];
sx q[0];
rz(0.50677104) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7031519) q[2];
sx q[2];
rz(-1.4435569) q[2];
sx q[2];
rz(1.3194989) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0171256) q[1];
sx q[1];
rz(-1.7199868) q[1];
sx q[1];
rz(-2.6890432) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6869557) q[3];
sx q[3];
rz(-2.2918923) q[3];
sx q[3];
rz(0.39814147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7740384) q[2];
sx q[2];
rz(-0.84711051) q[2];
sx q[2];
rz(-0.96662194) q[2];
rz(1.6537846) q[3];
sx q[3];
rz(-0.32854587) q[3];
sx q[3];
rz(-0.070913471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8129355) q[0];
sx q[0];
rz(-0.85449496) q[0];
sx q[0];
rz(0.27035126) q[0];
rz(-1.6290172) q[1];
sx q[1];
rz(-0.83643475) q[1];
sx q[1];
rz(0.62634748) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26890818) q[0];
sx q[0];
rz(-1.88694) q[0];
sx q[0];
rz(-2.288398) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93339351) q[2];
sx q[2];
rz(-1.2167756) q[2];
sx q[2];
rz(2.6824981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.69520187) q[1];
sx q[1];
rz(-3.0425515) q[1];
sx q[1];
rz(0.42548577) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8966497) q[3];
sx q[3];
rz(-1.2030461) q[3];
sx q[3];
rz(2.8830584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7865929) q[2];
sx q[2];
rz(-2.2248416) q[2];
sx q[2];
rz(1.5691441) q[2];
rz(2.3908424) q[3];
sx q[3];
rz(-2.7995977) q[3];
sx q[3];
rz(-0.87740889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56562051) q[0];
sx q[0];
rz(-1.2979869) q[0];
sx q[0];
rz(-0.57814231) q[0];
rz(1.3427973) q[1];
sx q[1];
rz(-2.8248351) q[1];
sx q[1];
rz(0.14541365) q[1];
rz(1.5708459) q[2];
sx q[2];
rz(-1.6931099) q[2];
sx q[2];
rz(-0.077559774) q[2];
rz(-2.883705) q[3];
sx q[3];
rz(-1.6226688) q[3];
sx q[3];
rz(-0.87165312) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
