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
rz(2.7569438) q[0];
sx q[0];
rz(-0.83847133) q[0];
sx q[0];
rz(-0.60929259) q[0];
rz(-2.4203909) q[1];
sx q[1];
rz(-2.2130122) q[1];
sx q[1];
rz(2.0452926) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4247835) q[0];
sx q[0];
rz(-1.4583298) q[0];
sx q[0];
rz(-2.2676668) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.9303665) q[2];
sx q[2];
rz(-1.1784679) q[2];
sx q[2];
rz(-2.0517672) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9723477) q[1];
sx q[1];
rz(-1.2903288) q[1];
sx q[1];
rz(0.53963668) q[1];
rz(-pi) q[2];
rz(-0.15518409) q[3];
sx q[3];
rz(-1.6801356) q[3];
sx q[3];
rz(-0.21462277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2234552) q[2];
sx q[2];
rz(-1.3104985) q[2];
sx q[2];
rz(-1.1892595) q[2];
rz(-2.7506645) q[3];
sx q[3];
rz(-1.9783741) q[3];
sx q[3];
rz(1.1322359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44553462) q[0];
sx q[0];
rz(-1.7867418) q[0];
sx q[0];
rz(-1.498244) q[0];
rz(2.4636726) q[1];
sx q[1];
rz(-1.8410212) q[1];
sx q[1];
rz(0.71151412) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16879119) q[0];
sx q[0];
rz(-2.0341221) q[0];
sx q[0];
rz(0.40131779) q[0];
rz(1.8078126) q[2];
sx q[2];
rz(-2.0227602) q[2];
sx q[2];
rz(-2.1353561) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.9715648) q[1];
sx q[1];
rz(-1.4907279) q[1];
sx q[1];
rz(-0.34645924) q[1];
rz(-pi) q[2];
rz(1.6145124) q[3];
sx q[3];
rz(-0.62407848) q[3];
sx q[3];
rz(1.8246375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70637643) q[2];
sx q[2];
rz(-1.7855568) q[2];
sx q[2];
rz(1.0630652) q[2];
rz(-1.4937909) q[3];
sx q[3];
rz(-0.46896514) q[3];
sx q[3];
rz(-0.74023214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4110334) q[0];
sx q[0];
rz(-2.7506802) q[0];
sx q[0];
rz(0.60182369) q[0];
rz(0.14566323) q[1];
sx q[1];
rz(-2.115695) q[1];
sx q[1];
rz(0.62785968) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3981741) q[0];
sx q[0];
rz(-1.6891915) q[0];
sx q[0];
rz(1.5009686) q[0];
rz(-pi) q[1];
rz(1.026903) q[2];
sx q[2];
rz(-2.5182407) q[2];
sx q[2];
rz(2.198213) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.8257432) q[1];
sx q[1];
rz(-1.9239195) q[1];
sx q[1];
rz(-0.24314116) q[1];
rz(0.95471621) q[3];
sx q[3];
rz(-2.4297415) q[3];
sx q[3];
rz(-2.4095132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.38498983) q[2];
sx q[2];
rz(-1.0641229) q[2];
sx q[2];
rz(0.21558726) q[2];
rz(-2.0032517) q[3];
sx q[3];
rz(-1.8214858) q[3];
sx q[3];
rz(-2.5707572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(0.1348006) q[0];
sx q[0];
rz(-1.726806) q[0];
sx q[0];
rz(-0.11882812) q[0];
rz(1.6670594) q[1];
sx q[1];
rz(-0.88913616) q[1];
sx q[1];
rz(2.3337591) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0519179) q[0];
sx q[0];
rz(-1.0052983) q[0];
sx q[0];
rz(-0.90948186) q[0];
rz(-pi) q[1];
rz(1.4363244) q[2];
sx q[2];
rz(-0.32609144) q[2];
sx q[2];
rz(-1.4876613) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1203907) q[1];
sx q[1];
rz(-0.63285108) q[1];
sx q[1];
rz(-1.8561884) q[1];
x q[2];
rz(2.4633154) q[3];
sx q[3];
rz(-1.5587574) q[3];
sx q[3];
rz(-1.2348242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.64041758) q[2];
sx q[2];
rz(-1.7265604) q[2];
sx q[2];
rz(2.625722) q[2];
rz(3.0473895) q[3];
sx q[3];
rz(-0.7754063) q[3];
sx q[3];
rz(1.5532106) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3696988) q[0];
sx q[0];
rz(-2.0409245) q[0];
sx q[0];
rz(-1.1871673) q[0];
rz(0.59133235) q[1];
sx q[1];
rz(-1.2966803) q[1];
sx q[1];
rz(-0.82685131) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88604414) q[0];
sx q[0];
rz(-0.40103087) q[0];
sx q[0];
rz(0.31440763) q[0];
x q[1];
rz(-1.0127147) q[2];
sx q[2];
rz(-1.6396171) q[2];
sx q[2];
rz(2.0134913) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8137774) q[1];
sx q[1];
rz(-2.2499488) q[1];
sx q[1];
rz(0.16736729) q[1];
x q[2];
rz(-0.50179568) q[3];
sx q[3];
rz(-1.0579946) q[3];
sx q[3];
rz(-0.38450188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7190711) q[2];
sx q[2];
rz(-2.8870388) q[2];
sx q[2];
rz(2.1264326) q[2];
rz(2.2319131) q[3];
sx q[3];
rz(-1.9010952) q[3];
sx q[3];
rz(0.66143405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61843094) q[0];
sx q[0];
rz(-1.9310512) q[0];
sx q[0];
rz(2.8316408) q[0];
rz(0.018772086) q[1];
sx q[1];
rz(-1.4614146) q[1];
sx q[1];
rz(-0.087336691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76742889) q[0];
sx q[0];
rz(-1.5383676) q[0];
sx q[0];
rz(-0.57291605) q[0];
x q[1];
rz(-2.2703253) q[2];
sx q[2];
rz(-2.0972898) q[2];
sx q[2];
rz(-0.73824182) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.046459196) q[1];
sx q[1];
rz(-1.1756056) q[1];
sx q[1];
rz(-1.0254775) q[1];
rz(-pi) q[2];
x q[2];
rz(2.099192) q[3];
sx q[3];
rz(-1.434552) q[3];
sx q[3];
rz(2.6353177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.35149082) q[2];
sx q[2];
rz(-2.6678706) q[2];
sx q[2];
rz(2.6215485) q[2];
rz(0.13488787) q[3];
sx q[3];
rz(-1.3210195) q[3];
sx q[3];
rz(-1.7322056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0670369) q[0];
sx q[0];
rz(-1.073607) q[0];
sx q[0];
rz(2.7556162) q[0];
rz(1.0999673) q[1];
sx q[1];
rz(-2.4620582) q[1];
sx q[1];
rz(-2.3540672) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5533977) q[0];
sx q[0];
rz(-2.6752691) q[0];
sx q[0];
rz(-0.88769261) q[0];
x q[1];
rz(-0.50561302) q[2];
sx q[2];
rz(-2.4596301) q[2];
sx q[2];
rz(-1.4346892) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0037076) q[1];
sx q[1];
rz(-0.862993) q[1];
sx q[1];
rz(-2.6801609) q[1];
rz(-pi) q[2];
x q[2];
rz(0.16997108) q[3];
sx q[3];
rz(-1.8026226) q[3];
sx q[3];
rz(-1.9915723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0218574) q[2];
sx q[2];
rz(-1.3081552) q[2];
sx q[2];
rz(-1.5137399) q[2];
rz(-1.2210023) q[3];
sx q[3];
rz(-1.1648213) q[3];
sx q[3];
rz(0.47202078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30962238) q[0];
sx q[0];
rz(-1.4291052) q[0];
sx q[0];
rz(2.9085462) q[0];
rz(-1.4876935) q[1];
sx q[1];
rz(-2.1079886) q[1];
sx q[1];
rz(-1.0071365) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24503532) q[0];
sx q[0];
rz(-0.63525891) q[0];
sx q[0];
rz(-0.59478514) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8677668) q[2];
sx q[2];
rz(-2.179702) q[2];
sx q[2];
rz(2.0021653) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1888652) q[1];
sx q[1];
rz(-2.7161274) q[1];
sx q[1];
rz(-1.9317606) q[1];
rz(-0.24407401) q[3];
sx q[3];
rz(-1.5904038) q[3];
sx q[3];
rz(-1.5527035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2795589) q[2];
sx q[2];
rz(-2.0290012) q[2];
sx q[2];
rz(-0.88279185) q[2];
rz(-2.8629996) q[3];
sx q[3];
rz(-1.9735347) q[3];
sx q[3];
rz(0.45894233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.19926628) q[0];
sx q[0];
rz(-0.47573221) q[0];
sx q[0];
rz(1.5717773) q[0];
rz(2.3528631) q[1];
sx q[1];
rz(-2.1756344) q[1];
sx q[1];
rz(-2.4615361) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2325033) q[0];
sx q[0];
rz(-1.6304558) q[0];
sx q[0];
rz(-0.76086126) q[0];
x q[1];
rz(-3.0228258) q[2];
sx q[2];
rz(-1.6965908) q[2];
sx q[2];
rz(2.072352) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9100807) q[1];
sx q[1];
rz(-0.68889131) q[1];
sx q[1];
rz(2.5206294) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2979638) q[3];
sx q[3];
rz(-1.1427726) q[3];
sx q[3];
rz(0.53800636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.98836977) q[2];
sx q[2];
rz(-0.92066568) q[2];
sx q[2];
rz(1.9752768) q[2];
rz(-1.204528) q[3];
sx q[3];
rz(-1.322999) q[3];
sx q[3];
rz(0.62674633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.523943) q[0];
sx q[0];
rz(-0.88215041) q[0];
sx q[0];
rz(-2.7045265) q[0];
rz(-2.3433459) q[1];
sx q[1];
rz(-0.50004807) q[1];
sx q[1];
rz(-0.22013586) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6740538) q[0];
sx q[0];
rz(-2.3706145) q[0];
sx q[0];
rz(-1.5037554) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9371168) q[2];
sx q[2];
rz(-1.2618974) q[2];
sx q[2];
rz(2.1155865) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4365873) q[1];
sx q[1];
rz(-1.2975177) q[1];
sx q[1];
rz(-0.80555861) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6159358) q[3];
sx q[3];
rz(-0.76971005) q[3];
sx q[3];
rz(1.5509645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2685214) q[2];
sx q[2];
rz(-1.9067418) q[2];
sx q[2];
rz(0.29042563) q[2];
rz(2.5051266) q[3];
sx q[3];
rz(-1.3822184) q[3];
sx q[3];
rz(-1.9071473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3270522) q[0];
sx q[0];
rz(-2.4335813) q[0];
sx q[0];
rz(-2.0306564) q[0];
rz(0.064432714) q[1];
sx q[1];
rz(-1.4705407) q[1];
sx q[1];
rz(2.4992117) q[1];
rz(0.9759554) q[2];
sx q[2];
rz(-1.5447164) q[2];
sx q[2];
rz(-1.0093052) q[2];
rz(-0.53917428) q[3];
sx q[3];
rz(-1.8001582) q[3];
sx q[3];
rz(-3.0348626) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
