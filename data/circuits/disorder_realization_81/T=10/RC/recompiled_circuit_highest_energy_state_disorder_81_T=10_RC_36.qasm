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
rz(-3.0523678) q[0];
sx q[0];
rz(-2.3951946) q[0];
sx q[0];
rz(2.449692) q[0];
rz(3.6045868) q[1];
sx q[1];
rz(4.7582518) q[1];
sx q[1];
rz(7.313348) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6917186) q[0];
sx q[0];
rz(-1.0152467) q[0];
sx q[0];
rz(-1.8434974) q[0];
rz(2.9582623) q[2];
sx q[2];
rz(-1.1784484) q[2];
sx q[2];
rz(3.0073056) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.2166442) q[1];
sx q[1];
rz(-0.34662592) q[1];
sx q[1];
rz(-1.9199172) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9943732) q[3];
sx q[3];
rz(-0.69425827) q[3];
sx q[3];
rz(-2.8822299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2618711) q[2];
sx q[2];
rz(-0.77360952) q[2];
sx q[2];
rz(-0.36641463) q[2];
rz(-1.0877747) q[3];
sx q[3];
rz(-0.78947869) q[3];
sx q[3];
rz(-2.4521949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.825603) q[0];
sx q[0];
rz(-3.0593384) q[0];
sx q[0];
rz(0.91973037) q[0];
rz(2.3986744) q[1];
sx q[1];
rz(-1.2580322) q[1];
sx q[1];
rz(-1.119335) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21788439) q[0];
sx q[0];
rz(-2.9842434) q[0];
sx q[0];
rz(-1.9718534) q[0];
rz(0.2074457) q[2];
sx q[2];
rz(-2.4118057) q[2];
sx q[2];
rz(-1.6285673) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.38779681) q[1];
sx q[1];
rz(-0.8555853) q[1];
sx q[1];
rz(1.254735) q[1];
rz(-2.2833586) q[3];
sx q[3];
rz(-0.7131812) q[3];
sx q[3];
rz(-2.3767917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6433581) q[2];
sx q[2];
rz(-1.1946119) q[2];
sx q[2];
rz(0.052113459) q[2];
rz(2.0335967) q[3];
sx q[3];
rz(-2.2823157) q[3];
sx q[3];
rz(-3.0767379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4614918) q[0];
sx q[0];
rz(-1.9100186) q[0];
sx q[0];
rz(1.289806) q[0];
rz(2.3884804) q[1];
sx q[1];
rz(-1.6878004) q[1];
sx q[1];
rz(-2.6444816) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6996993) q[0];
sx q[0];
rz(-1.178368) q[0];
sx q[0];
rz(0.59458574) q[0];
rz(-pi) q[1];
rz(0.67694734) q[2];
sx q[2];
rz(-1.0407037) q[2];
sx q[2];
rz(2.9924336) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0935082) q[1];
sx q[1];
rz(-1.1968379) q[1];
sx q[1];
rz(-1.611534) q[1];
rz(1.6670973) q[3];
sx q[3];
rz(-0.9203921) q[3];
sx q[3];
rz(1.4092829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.46811732) q[2];
sx q[2];
rz(-2.868729) q[2];
sx q[2];
rz(-2.4172778) q[2];
rz(-0.22504462) q[3];
sx q[3];
rz(-1.2757755) q[3];
sx q[3];
rz(-0.79506522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1555136) q[0];
sx q[0];
rz(-0.88548311) q[0];
sx q[0];
rz(-2.8610863) q[0];
rz(1.6579423) q[1];
sx q[1];
rz(-1.9397441) q[1];
sx q[1];
rz(2.6703506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0126033) q[0];
sx q[0];
rz(-1.5935197) q[0];
sx q[0];
rz(-1.9381592) q[0];
rz(3.0755694) q[2];
sx q[2];
rz(-1.4581265) q[2];
sx q[2];
rz(2.6875935) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1859101) q[1];
sx q[1];
rz(-1.5985399) q[1];
sx q[1];
rz(1.300059) q[1];
x q[2];
rz(-2.1625278) q[3];
sx q[3];
rz(-2.2869898) q[3];
sx q[3];
rz(-1.5997052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.19807854) q[2];
sx q[2];
rz(-0.55800617) q[2];
sx q[2];
rz(-2.1264709) q[2];
rz(0.73457581) q[3];
sx q[3];
rz(-1.3646804) q[3];
sx q[3];
rz(-2.6034897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9839142) q[0];
sx q[0];
rz(-1.2175125) q[0];
sx q[0];
rz(0.96018803) q[0];
rz(-3.0958815) q[1];
sx q[1];
rz(-2.263133) q[1];
sx q[1];
rz(-3.1083621) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3094745) q[0];
sx q[0];
rz(-2.4178979) q[0];
sx q[0];
rz(0.81731082) q[0];
x q[1];
rz(1.303238) q[2];
sx q[2];
rz(-1.5853527) q[2];
sx q[2];
rz(-2.6637258) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7864125) q[1];
sx q[1];
rz(-0.35874736) q[1];
sx q[1];
rz(1.9500109) q[1];
x q[2];
rz(1.1636249) q[3];
sx q[3];
rz(-1.1615175) q[3];
sx q[3];
rz(-0.57151505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1632605) q[2];
sx q[2];
rz(-0.99908081) q[2];
sx q[2];
rz(0.31326374) q[2];
rz(0.35798171) q[3];
sx q[3];
rz(-1.0397747) q[3];
sx q[3];
rz(0.09859214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1545496) q[0];
sx q[0];
rz(-1.510226) q[0];
sx q[0];
rz(-2.6556067) q[0];
rz(1.1827353) q[1];
sx q[1];
rz(-2.4458838) q[1];
sx q[1];
rz(2.8114496) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3957452) q[0];
sx q[0];
rz(-1.8680649) q[0];
sx q[0];
rz(0.031121522) q[0];
x q[1];
rz(-0.92440549) q[2];
sx q[2];
rz(-0.87714783) q[2];
sx q[2];
rz(-2.2104757) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9818519) q[1];
sx q[1];
rz(-1.0969437) q[1];
sx q[1];
rz(-2.5706446) q[1];
x q[2];
rz(2.7948912) q[3];
sx q[3];
rz(-1.3606405) q[3];
sx q[3];
rz(-2.8488248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6351472) q[2];
sx q[2];
rz(-0.61656419) q[2];
sx q[2];
rz(-0.35663024) q[2];
rz(0.59456524) q[3];
sx q[3];
rz(-1.4831355) q[3];
sx q[3];
rz(-2.4026292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4643788) q[0];
sx q[0];
rz(-2.3742299) q[0];
sx q[0];
rz(-0.58694029) q[0];
rz(-0.37239536) q[1];
sx q[1];
rz(-1.7814691) q[1];
sx q[1];
rz(2.8940103) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42668396) q[0];
sx q[0];
rz(-1.539402) q[0];
sx q[0];
rz(2.035729) q[0];
rz(2.8147698) q[2];
sx q[2];
rz(-2.4832209) q[2];
sx q[2];
rz(-3.1382552) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.5953427) q[1];
sx q[1];
rz(-0.74186013) q[1];
sx q[1];
rz(-1.0651903) q[1];
rz(-pi) q[2];
rz(-2.5334397) q[3];
sx q[3];
rz(-2.7233015) q[3];
sx q[3];
rz(-0.88514144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8783012) q[2];
sx q[2];
rz(-1.9204488) q[2];
sx q[2];
rz(2.9325874) q[2];
rz(-2.5448997) q[3];
sx q[3];
rz(-0.34691072) q[3];
sx q[3];
rz(1.6056812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40808943) q[0];
sx q[0];
rz(-0.96036378) q[0];
sx q[0];
rz(0.39304131) q[0];
rz(-2.4709002) q[1];
sx q[1];
rz(-1.1436983) q[1];
sx q[1];
rz(1.6938946) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.873998) q[0];
sx q[0];
rz(-2.0646668) q[0];
sx q[0];
rz(-0.14125342) q[0];
rz(-pi) q[1];
rz(2.1850519) q[2];
sx q[2];
rz(-1.8280085) q[2];
sx q[2];
rz(-1.1957912) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5911342) q[1];
sx q[1];
rz(-2.0789642) q[1];
sx q[1];
rz(2.152446) q[1];
rz(-pi) q[2];
rz(-0.46341571) q[3];
sx q[3];
rz(-0.8532195) q[3];
sx q[3];
rz(-1.0427102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4286917) q[2];
sx q[2];
rz(-0.8608326) q[2];
sx q[2];
rz(-3.0248896) q[2];
rz(1.0137001) q[3];
sx q[3];
rz(-1.9096392) q[3];
sx q[3];
rz(1.3056171) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882465) q[0];
sx q[0];
rz(-1.4660864) q[0];
sx q[0];
rz(-1.4353132) q[0];
rz(-2.1460136) q[1];
sx q[1];
rz(-1.3187871) q[1];
sx q[1];
rz(2.6011655) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30975809) q[0];
sx q[0];
rz(-2.0481526) q[0];
sx q[0];
rz(-1.274513) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2733354) q[2];
sx q[2];
rz(-2.3405055) q[2];
sx q[2];
rz(-2.4641387) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1241972) q[1];
sx q[1];
rz(-2.3216343) q[1];
sx q[1];
rz(2.6142121) q[1];
rz(1.1630658) q[3];
sx q[3];
rz(-0.98568788) q[3];
sx q[3];
rz(-1.2297022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2702177) q[2];
sx q[2];
rz(-1.8367218) q[2];
sx q[2];
rz(-2.3649141) q[2];
rz(-3.1089879) q[3];
sx q[3];
rz(-1.5332581) q[3];
sx q[3];
rz(1.5672654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8107574) q[0];
sx q[0];
rz(-2.7726655) q[0];
sx q[0];
rz(0.43357968) q[0];
rz(-0.66917229) q[1];
sx q[1];
rz(-1.9586261) q[1];
sx q[1];
rz(0.42632595) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.688072) q[0];
sx q[0];
rz(-1.768508) q[0];
sx q[0];
rz(-1.3526249) q[0];
rz(-pi) q[1];
rz(0.74684116) q[2];
sx q[2];
rz(-1.8364808) q[2];
sx q[2];
rz(-2.8248252) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.9712228) q[1];
sx q[1];
rz(-1.0863546) q[1];
sx q[1];
rz(-1.96934) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4678585) q[3];
sx q[3];
rz(-0.70286432) q[3];
sx q[3];
rz(-1.7029312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3907889) q[2];
sx q[2];
rz(-0.28826737) q[2];
sx q[2];
rz(-1.1061579) q[2];
rz(0.11354167) q[3];
sx q[3];
rz(-0.93324408) q[3];
sx q[3];
rz(-3.09789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9857367) q[0];
sx q[0];
rz(-0.899213) q[0];
sx q[0];
rz(-0.56094299) q[0];
rz(2.7583495) q[1];
sx q[1];
rz(-2.0189197) q[1];
sx q[1];
rz(-0.22519208) q[1];
rz(-1.8763992) q[2];
sx q[2];
rz(-1.0355179) q[2];
sx q[2];
rz(-2.1260362) q[2];
rz(2.0416971) q[3];
sx q[3];
rz(-2.3318137) q[3];
sx q[3];
rz(-1.6398738) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
