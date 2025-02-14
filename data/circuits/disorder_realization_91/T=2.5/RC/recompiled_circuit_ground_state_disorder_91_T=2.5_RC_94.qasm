OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8102326) q[0];
sx q[0];
rz(4.3963764) q[0];
sx q[0];
rz(8.7788361) q[0];
rz(1.6826001) q[1];
sx q[1];
rz(3.269722) q[1];
sx q[1];
rz(10.092957) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39820078) q[0];
sx q[0];
rz(-1.8132134) q[0];
sx q[0];
rz(-1.4546118) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8142848) q[2];
sx q[2];
rz(-2.1919554) q[2];
sx q[2];
rz(0.46681625) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1937374) q[1];
sx q[1];
rz(-1.9661709) q[1];
sx q[1];
rz(2.1389524) q[1];
rz(-pi) q[2];
rz(2.6553538) q[3];
sx q[3];
rz(-0.91043962) q[3];
sx q[3];
rz(2.5442991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51840034) q[2];
sx q[2];
rz(-2.4861768) q[2];
sx q[2];
rz(-1.0830967) q[2];
rz(2.8862503) q[3];
sx q[3];
rz(-0.89699236) q[3];
sx q[3];
rz(-1.9820836) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1622247) q[0];
sx q[0];
rz(-0.14142445) q[0];
sx q[0];
rz(-0.033578385) q[0];
rz(2.0032739) q[1];
sx q[1];
rz(-0.99457026) q[1];
sx q[1];
rz(0.23695645) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36719272) q[0];
sx q[0];
rz(-2.8731308) q[0];
sx q[0];
rz(2.5949097) q[0];
x q[1];
rz(-2.9541624) q[2];
sx q[2];
rz(-1.1016365) q[2];
sx q[2];
rz(-0.7721061) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1019056) q[1];
sx q[1];
rz(-0.19409212) q[1];
sx q[1];
rz(-1.7029352) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44741607) q[3];
sx q[3];
rz(-2.1037115) q[3];
sx q[3];
rz(-0.79870236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.529155) q[2];
sx q[2];
rz(-2.0945022) q[2];
sx q[2];
rz(0.1203514) q[2];
rz(-0.58593166) q[3];
sx q[3];
rz(-1.7610995) q[3];
sx q[3];
rz(-1.2683292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.878207) q[0];
sx q[0];
rz(-2.1757941) q[0];
sx q[0];
rz(0.98814386) q[0];
rz(-1.490961) q[1];
sx q[1];
rz(-2.4246876) q[1];
sx q[1];
rz(1.5536701) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33086209) q[0];
sx q[0];
rz(-1.5867044) q[0];
sx q[0];
rz(1.5696947) q[0];
rz(-pi) q[1];
rz(-1.6397255) q[2];
sx q[2];
rz(-1.4367099) q[2];
sx q[2];
rz(-0.20251911) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8804259) q[1];
sx q[1];
rz(-2.8787887) q[1];
sx q[1];
rz(1.6815503) q[1];
x q[2];
rz(1.9918898) q[3];
sx q[3];
rz(-2.097887) q[3];
sx q[3];
rz(2.2452109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3042018) q[2];
sx q[2];
rz(-1.7981671) q[2];
sx q[2];
rz(3.0710132) q[2];
rz(-0.48921674) q[3];
sx q[3];
rz(-2.1507806) q[3];
sx q[3];
rz(2.9941471) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9551142) q[0];
sx q[0];
rz(-0.69545737) q[0];
sx q[0];
rz(1.7406933) q[0];
rz(-0.1700302) q[1];
sx q[1];
rz(-1.4100807) q[1];
sx q[1];
rz(1.9084575) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2182541) q[0];
sx q[0];
rz(-1.6955351) q[0];
sx q[0];
rz(1.1472923) q[0];
rz(-pi) q[1];
x q[1];
rz(0.71738404) q[2];
sx q[2];
rz(-1.1108494) q[2];
sx q[2];
rz(-1.1954952) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4608232) q[1];
sx q[1];
rz(-1.5857547) q[1];
sx q[1];
rz(-1.2911002) q[1];
x q[2];
rz(-1.5828269) q[3];
sx q[3];
rz(-1.1063442) q[3];
sx q[3];
rz(2.8248276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.79973811) q[2];
sx q[2];
rz(-1.4350472) q[2];
sx q[2];
rz(2.5119761) q[2];
rz(3.0047505) q[3];
sx q[3];
rz(-0.62961951) q[3];
sx q[3];
rz(-1.7262044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8585696) q[0];
sx q[0];
rz(-0.50823277) q[0];
sx q[0];
rz(0.10230219) q[0];
rz(1.6434068) q[1];
sx q[1];
rz(-1.841265) q[1];
sx q[1];
rz(2.7844875) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6320366) q[0];
sx q[0];
rz(-0.3118383) q[0];
sx q[0];
rz(-2.3760892) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5979366) q[2];
sx q[2];
rz(-1.4476336) q[2];
sx q[2];
rz(0.20448286) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.50430027) q[1];
sx q[1];
rz(-1.2724814) q[1];
sx q[1];
rz(0.53108414) q[1];
rz(-pi) q[2];
rz(-2.6980347) q[3];
sx q[3];
rz(-2.0152115) q[3];
sx q[3];
rz(-2.9068974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6685278) q[2];
sx q[2];
rz(-1.7230956) q[2];
sx q[2];
rz(-1.8190039) q[2];
rz(1.4332917) q[3];
sx q[3];
rz(-1.8247484) q[3];
sx q[3];
rz(-0.1639666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.598031) q[0];
sx q[0];
rz(-2.141641) q[0];
sx q[0];
rz(1.7033956) q[0];
rz(1.9012798) q[1];
sx q[1];
rz(-0.65483171) q[1];
sx q[1];
rz(-1.7887438) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2738631) q[0];
sx q[0];
rz(-0.41711787) q[0];
sx q[0];
rz(1.3365082) q[0];
rz(-pi) q[1];
rz(-0.5261658) q[2];
sx q[2];
rz(-1.7169098) q[2];
sx q[2];
rz(-2.7046159) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6412615) q[1];
sx q[1];
rz(-1.2029543) q[1];
sx q[1];
rz(-1.0235051) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7463673) q[3];
sx q[3];
rz(-1.4661643) q[3];
sx q[3];
rz(-2.6992309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.214434) q[2];
sx q[2];
rz(-0.36403251) q[2];
sx q[2];
rz(0.39043179) q[2];
rz(-1.8292142) q[3];
sx q[3];
rz(-0.81109154) q[3];
sx q[3];
rz(-0.81982476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.696233) q[0];
sx q[0];
rz(-0.92388988) q[0];
sx q[0];
rz(-1.4129289) q[0];
rz(1.6784809) q[1];
sx q[1];
rz(-1.585958) q[1];
sx q[1];
rz(-2.5020592) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7752658) q[0];
sx q[0];
rz(-1.9444939) q[0];
sx q[0];
rz(-2.213272) q[0];
rz(-1.495218) q[2];
sx q[2];
rz(-1.7210782) q[2];
sx q[2];
rz(1.4270475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9078341) q[1];
sx q[1];
rz(-0.69926622) q[1];
sx q[1];
rz(-2.722867) q[1];
x q[2];
rz(-3.022081) q[3];
sx q[3];
rz(-1.9627213) q[3];
sx q[3];
rz(0.78531314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21706906) q[2];
sx q[2];
rz(-1.3288493) q[2];
sx q[2];
rz(0.45219839) q[2];
rz(-2.9186115) q[3];
sx q[3];
rz(-0.28135869) q[3];
sx q[3];
rz(-0.15534672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13101354) q[0];
sx q[0];
rz(-2.2192945) q[0];
sx q[0];
rz(2.9796694) q[0];
rz(0.34795347) q[1];
sx q[1];
rz(-1.8674928) q[1];
sx q[1];
rz(-2.9071992) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1944626) q[0];
sx q[0];
rz(-0.85887733) q[0];
sx q[0];
rz(-2.3361337) q[0];
x q[1];
rz(1.8816119) q[2];
sx q[2];
rz(-2.2665885) q[2];
sx q[2];
rz(-1.1281079) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.80658) q[1];
sx q[1];
rz(-1.3836821) q[1];
sx q[1];
rz(-2.2172865) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9107844) q[3];
sx q[3];
rz(-0.2838906) q[3];
sx q[3];
rz(-2.6687255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.28635412) q[2];
sx q[2];
rz(-2.0704634) q[2];
sx q[2];
rz(0.62412778) q[2];
rz(0.74328077) q[3];
sx q[3];
rz(-2.1478839) q[3];
sx q[3];
rz(-0.68964094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2766992) q[0];
sx q[0];
rz(-1.6962637) q[0];
sx q[0];
rz(1.0869166) q[0];
rz(-1.2326321) q[1];
sx q[1];
rz(-2.0338567) q[1];
sx q[1];
rz(-1.3023652) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24827458) q[0];
sx q[0];
rz(-2.9030307) q[0];
sx q[0];
rz(0.49805157) q[0];
rz(2.0696569) q[2];
sx q[2];
rz(-2.4833224) q[2];
sx q[2];
rz(-1.9287623) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.86587807) q[1];
sx q[1];
rz(-0.68531407) q[1];
sx q[1];
rz(0.65071836) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.62303657) q[3];
sx q[3];
rz(-1.5861694) q[3];
sx q[3];
rz(-2.8492209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6752601) q[2];
sx q[2];
rz(-1.8856498) q[2];
sx q[2];
rz(-2.7272398) q[2];
rz(-0.77053344) q[3];
sx q[3];
rz(-1.7194175) q[3];
sx q[3];
rz(0.93366247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3122124) q[0];
sx q[0];
rz(-0.82013622) q[0];
sx q[0];
rz(2.6729551) q[0];
rz(1.8679484) q[1];
sx q[1];
rz(-1.8654774) q[1];
sx q[1];
rz(0.31148568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0429223) q[0];
sx q[0];
rz(-0.85349529) q[0];
sx q[0];
rz(3.0649661) q[0];
x q[1];
rz(2.5045372) q[2];
sx q[2];
rz(-0.14547507) q[2];
sx q[2];
rz(1.5381952) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.905924) q[1];
sx q[1];
rz(-2.0579268) q[1];
sx q[1];
rz(-3.0887927) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52759513) q[3];
sx q[3];
rz(-2.0012534) q[3];
sx q[3];
rz(-0.93387077) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.0353388) q[2];
sx q[2];
rz(-1.9474578) q[2];
sx q[2];
rz(1.7485471) q[2];
rz(-2.2140908) q[3];
sx q[3];
rz(-1.5774957) q[3];
sx q[3];
rz(0.86021304) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2911693) q[0];
sx q[0];
rz(-2.8130154) q[0];
sx q[0];
rz(-1.6541506) q[0];
rz(-2.9604079) q[1];
sx q[1];
rz(-0.94534992) q[1];
sx q[1];
rz(2.2015991) q[1];
rz(2.1045024) q[2];
sx q[2];
rz(-1.7052393) q[2];
sx q[2];
rz(1.3503804) q[2];
rz(-2.5355382) q[3];
sx q[3];
rz(-2.1447976) q[3];
sx q[3];
rz(-1.6404649) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
