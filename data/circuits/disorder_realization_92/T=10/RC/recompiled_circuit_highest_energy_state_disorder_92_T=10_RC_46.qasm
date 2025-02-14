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
rz(-3.0109316) q[0];
sx q[0];
rz(-1.8008512) q[0];
sx q[0];
rz(1.8774348) q[0];
rz(0.90509993) q[1];
sx q[1];
rz(-1.0882508) q[1];
sx q[1];
rz(-3.1286214) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91468477) q[0];
sx q[0];
rz(-2.6400262) q[0];
sx q[0];
rz(-2.8109994) q[0];
x q[1];
rz(1.6799203) q[2];
sx q[2];
rz(-1.732601) q[2];
sx q[2];
rz(-1.2798123) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5789507) q[1];
sx q[1];
rz(-1.2516252) q[1];
sx q[1];
rz(0.36960543) q[1];
x q[2];
rz(-2.5586224) q[3];
sx q[3];
rz(-1.7503925) q[3];
sx q[3];
rz(-2.666134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.26244792) q[2];
sx q[2];
rz(-1.6877354) q[2];
sx q[2];
rz(-0.095414735) q[2];
rz(2.6573507) q[3];
sx q[3];
rz(-2.3384194) q[3];
sx q[3];
rz(-1.6835083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5937623) q[0];
sx q[0];
rz(-1.7050803) q[0];
sx q[0];
rz(2.4875212) q[0];
rz(-1.3520799) q[1];
sx q[1];
rz(-2.4857931) q[1];
sx q[1];
rz(-0.10993122) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1760565) q[0];
sx q[0];
rz(-1.5829464) q[0];
sx q[0];
rz(1.5909082) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6835126) q[2];
sx q[2];
rz(-0.9914248) q[2];
sx q[2];
rz(1.9922076) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.054903237) q[1];
sx q[1];
rz(-2.6536021) q[1];
sx q[1];
rz(-2.7665374) q[1];
x q[2];
rz(-1.8255005) q[3];
sx q[3];
rz(-1.1858018) q[3];
sx q[3];
rz(-2.8196872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8249417) q[2];
sx q[2];
rz(-1.2085088) q[2];
sx q[2];
rz(1.6744772) q[2];
rz(-1.0351099) q[3];
sx q[3];
rz(-0.78101522) q[3];
sx q[3];
rz(-1.8234113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3159981) q[0];
sx q[0];
rz(-2.9075629) q[0];
sx q[0];
rz(-0.79063928) q[0];
rz(2.3979893) q[1];
sx q[1];
rz(-1.5433106) q[1];
sx q[1];
rz(2.5620983) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52344874) q[0];
sx q[0];
rz(-1.0001527) q[0];
sx q[0];
rz(-0.6499001) q[0];
rz(2.8786009) q[2];
sx q[2];
rz(-1.8557669) q[2];
sx q[2];
rz(-2.6927352) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.0708376) q[1];
sx q[1];
rz(-2.8808314) q[1];
sx q[1];
rz(-1.7160796) q[1];
rz(-1.0813339) q[3];
sx q[3];
rz(-1.0940942) q[3];
sx q[3];
rz(0.34927327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4641331) q[2];
sx q[2];
rz(-0.83659283) q[2];
sx q[2];
rz(-0.38132384) q[2];
rz(-2.2903806) q[3];
sx q[3];
rz(-2.2150453) q[3];
sx q[3];
rz(0.73494953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5320936) q[0];
sx q[0];
rz(-2.5272326) q[0];
sx q[0];
rz(-2.7771948) q[0];
rz(-2.9413307) q[1];
sx q[1];
rz(-1.7297144) q[1];
sx q[1];
rz(-3.0416378) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9635506) q[0];
sx q[0];
rz(-0.077712312) q[0];
sx q[0];
rz(1.8366209) q[0];
x q[1];
rz(-0.1047524) q[2];
sx q[2];
rz(-1.9974553) q[2];
sx q[2];
rz(-1.1327281) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.050592) q[1];
sx q[1];
rz(-0.09517955) q[1];
sx q[1];
rz(-2.9042713) q[1];
x q[2];
rz(2.3175008) q[3];
sx q[3];
rz(-0.76585117) q[3];
sx q[3];
rz(-1.6834109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0692856) q[2];
sx q[2];
rz(-2.0733209) q[2];
sx q[2];
rz(-0.11492534) q[2];
rz(-2.0011486) q[3];
sx q[3];
rz(-1.8585669) q[3];
sx q[3];
rz(2.1663402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2045778) q[0];
sx q[0];
rz(-2.1890722) q[0];
sx q[0];
rz(-0.17247795) q[0];
rz(-1.0109488) q[1];
sx q[1];
rz(-0.5609678) q[1];
sx q[1];
rz(-0.15484658) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5342543) q[0];
sx q[0];
rz(-1.3528498) q[0];
sx q[0];
rz(1.2675257) q[0];
rz(-pi) q[1];
rz(-0.24598083) q[2];
sx q[2];
rz(-0.94106758) q[2];
sx q[2];
rz(-0.069815947) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.46609572) q[1];
sx q[1];
rz(-1.4267529) q[1];
sx q[1];
rz(-2.7269468) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17821594) q[3];
sx q[3];
rz(-1.4608188) q[3];
sx q[3];
rz(-2.4770155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.054472063) q[2];
sx q[2];
rz(-2.8643769) q[2];
sx q[2];
rz(0.37230125) q[2];
rz(-0.39792684) q[3];
sx q[3];
rz(-1.1436661) q[3];
sx q[3];
rz(-0.00042375617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18846866) q[0];
sx q[0];
rz(-0.57250452) q[0];
sx q[0];
rz(-1.5484126) q[0];
rz(-0.36987034) q[1];
sx q[1];
rz(-1.7431755) q[1];
sx q[1];
rz(-0.63873783) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4700482) q[0];
sx q[0];
rz(-2.8570768) q[0];
sx q[0];
rz(-2.1485062) q[0];
rz(-2.0962786) q[2];
sx q[2];
rz(-1.7470226) q[2];
sx q[2];
rz(0.8187364) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.35013851) q[1];
sx q[1];
rz(-1.119507) q[1];
sx q[1];
rz(-0.43124388) q[1];
x q[2];
rz(-0.5587479) q[3];
sx q[3];
rz(-2.2674446) q[3];
sx q[3];
rz(0.071867094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0686191) q[2];
sx q[2];
rz(-2.0899453) q[2];
sx q[2];
rz(-0.2529141) q[2];
rz(-0.84826338) q[3];
sx q[3];
rz(-1.6631283) q[3];
sx q[3];
rz(-3.0566791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10097583) q[0];
sx q[0];
rz(-0.210013) q[0];
sx q[0];
rz(-0.14895359) q[0];
rz(1.3111929) q[1];
sx q[1];
rz(-1.4102178) q[1];
sx q[1];
rz(3.0874918) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.220517) q[0];
sx q[0];
rz(-2.069325) q[0];
sx q[0];
rz(-1.4709298) q[0];
rz(-pi) q[1];
rz(-2.3391244) q[2];
sx q[2];
rz(-1.5395428) q[2];
sx q[2];
rz(-0.6337589) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74443597) q[1];
sx q[1];
rz(-2.2658093) q[1];
sx q[1];
rz(-1.7139462) q[1];
x q[2];
rz(1.0955174) q[3];
sx q[3];
rz(-1.5500808) q[3];
sx q[3];
rz(-1.4521862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.80186239) q[2];
sx q[2];
rz(-1.170155) q[2];
sx q[2];
rz(2.1745963) q[2];
rz(2.4407834) q[3];
sx q[3];
rz(-2.1613224) q[3];
sx q[3];
rz(0.35082671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(3.0130149) q[0];
sx q[0];
rz(-1.1739434) q[0];
sx q[0];
rz(-1.5933734) q[0];
rz(2.7633527) q[1];
sx q[1];
rz(-2.389237) q[1];
sx q[1];
rz(0.38984782) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0336817) q[0];
sx q[0];
rz(-2.0706688) q[0];
sx q[0];
rz(1.3477911) q[0];
rz(-pi) q[1];
rz(-1.5178096) q[2];
sx q[2];
rz(-0.30445004) q[2];
sx q[2];
rz(1.9805769) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88943849) q[1];
sx q[1];
rz(-2.5233626) q[1];
sx q[1];
rz(-0.63207027) q[1];
rz(1.3482565) q[3];
sx q[3];
rz(-1.8197818) q[3];
sx q[3];
rz(-2.1771262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.0491911) q[2];
sx q[2];
rz(-2.0202049) q[2];
sx q[2];
rz(1.8828877) q[2];
rz(-2.8973268) q[3];
sx q[3];
rz(-1.3552856) q[3];
sx q[3];
rz(2.145483) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32023892) q[0];
sx q[0];
rz(-1.9122253) q[0];
sx q[0];
rz(-1.7852596) q[0];
rz(-2.9755196) q[1];
sx q[1];
rz(-2.1898654) q[1];
sx q[1];
rz(0.43620268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22101519) q[0];
sx q[0];
rz(-1.8816784) q[0];
sx q[0];
rz(0.47500821) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8955375) q[2];
sx q[2];
rz(-0.51273275) q[2];
sx q[2];
rz(0.89287478) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2959152) q[1];
sx q[1];
rz(-1.3406173) q[1];
sx q[1];
rz(-2.3451817) q[1];
rz(-2.9968676) q[3];
sx q[3];
rz(-0.44071482) q[3];
sx q[3];
rz(-0.056469744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0210375) q[2];
sx q[2];
rz(-1.1641116) q[2];
sx q[2];
rz(-2.5980921) q[2];
rz(-3.0920658) q[3];
sx q[3];
rz(-1.4855569) q[3];
sx q[3];
rz(1.4878976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5307584) q[0];
sx q[0];
rz(-3.0093091) q[0];
sx q[0];
rz(-0.079205967) q[0];
rz(2.7414956) q[1];
sx q[1];
rz(-2.2875417) q[1];
sx q[1];
rz(-2.1150463) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0745747) q[0];
sx q[0];
rz(-1.8047143) q[0];
sx q[0];
rz(-0.83252711) q[0];
x q[1];
rz(2.8978149) q[2];
sx q[2];
rz(-0.37636435) q[2];
sx q[2];
rz(1.7610904) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.77791926) q[1];
sx q[1];
rz(-1.8936186) q[1];
sx q[1];
rz(-1.670091) q[1];
rz(-pi) q[2];
rz(1.6502569) q[3];
sx q[3];
rz(-1.2209519) q[3];
sx q[3];
rz(-2.3235278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2058699) q[2];
sx q[2];
rz(-1.3047855) q[2];
sx q[2];
rz(1.6892461) q[2];
rz(-0.82990372) q[3];
sx q[3];
rz(-0.87703505) q[3];
sx q[3];
rz(-0.95480603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6185388) q[0];
sx q[0];
rz(-1.1445615) q[0];
sx q[0];
rz(-1.6339697) q[0];
rz(-1.8745096) q[1];
sx q[1];
rz(-2.0590084) q[1];
sx q[1];
rz(0.72437292) q[1];
rz(-2.6020423) q[2];
sx q[2];
rz(-1.0536516) q[2];
sx q[2];
rz(-0.45894844) q[2];
rz(-2.6250056) q[3];
sx q[3];
rz(-0.62415515) q[3];
sx q[3];
rz(-2.2976807) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
