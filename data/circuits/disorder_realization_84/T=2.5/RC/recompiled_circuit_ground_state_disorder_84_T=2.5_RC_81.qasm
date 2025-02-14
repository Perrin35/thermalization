OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.328188) q[0];
sx q[0];
rz(-2.5321811) q[0];
sx q[0];
rz(-3.1031026) q[0];
rz(2.6961532) q[1];
sx q[1];
rz(-0.69094509) q[1];
sx q[1];
rz(3.0145187) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84665438) q[0];
sx q[0];
rz(-1.9401258) q[0];
sx q[0];
rz(-2.3890004) q[0];
rz(1.5760001) q[2];
sx q[2];
rz(-1.5706816) q[2];
sx q[2];
rz(-0.031129908) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.91158479) q[1];
sx q[1];
rz(-0.42427166) q[1];
sx q[1];
rz(-2.2390635) q[1];
rz(-pi) q[2];
rz(2.6576275) q[3];
sx q[3];
rz(-2.1089156) q[3];
sx q[3];
rz(0.78236587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8700063) q[2];
sx q[2];
rz(-2.8903457) q[2];
sx q[2];
rz(-2.0269488) q[2];
rz(-0.33846578) q[3];
sx q[3];
rz(-1.7287799) q[3];
sx q[3];
rz(-0.72845212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2289497) q[0];
sx q[0];
rz(-2.8680608) q[0];
sx q[0];
rz(-1.1852513) q[0];
rz(1.7456906) q[1];
sx q[1];
rz(-1.6604559) q[1];
sx q[1];
rz(3.0286068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11155726) q[0];
sx q[0];
rz(-1.5148692) q[0];
sx q[0];
rz(1.5779172) q[0];
x q[1];
rz(1.127203) q[2];
sx q[2];
rz(-1.796486) q[2];
sx q[2];
rz(-2.0114102) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.005882) q[1];
sx q[1];
rz(-0.23488472) q[1];
sx q[1];
rz(-1.2967111) q[1];
x q[2];
rz(-2.7725622) q[3];
sx q[3];
rz(-1.1347767) q[3];
sx q[3];
rz(-2.8216803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.2286171) q[2];
sx q[2];
rz(-1.7824219) q[2];
sx q[2];
rz(1.011298) q[2];
rz(3.0806115) q[3];
sx q[3];
rz(-2.0296378) q[3];
sx q[3];
rz(-2.8673867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7762452) q[0];
sx q[0];
rz(-1.341935) q[0];
sx q[0];
rz(1.9976529) q[0];
rz(-1.9260319) q[1];
sx q[1];
rz(-2.2823915) q[1];
sx q[1];
rz(-2.5124195) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1224642) q[0];
sx q[0];
rz(-0.97342426) q[0];
sx q[0];
rz(3.1015784) q[0];
x q[1];
rz(-2.0780656) q[2];
sx q[2];
rz(-1.4794297) q[2];
sx q[2];
rz(-0.12742119) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.53051004) q[1];
sx q[1];
rz(-1.9323009) q[1];
sx q[1];
rz(1.481249) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8820417) q[3];
sx q[3];
rz(-1.5207612) q[3];
sx q[3];
rz(-0.52122859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4048142) q[2];
sx q[2];
rz(-1.7958769) q[2];
sx q[2];
rz(0.2573615) q[2];
rz(-1.8836053) q[3];
sx q[3];
rz(-1.6519929) q[3];
sx q[3];
rz(0.45132288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8266206) q[0];
sx q[0];
rz(-0.74084145) q[0];
sx q[0];
rz(-1.6175057) q[0];
rz(-1.357088) q[1];
sx q[1];
rz(-0.70248228) q[1];
sx q[1];
rz(-1.2290899) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8600677) q[0];
sx q[0];
rz(-1.901214) q[0];
sx q[0];
rz(1.2072674) q[0];
rz(0.76009373) q[2];
sx q[2];
rz(-2.1017535) q[2];
sx q[2];
rz(-2.107055) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.22013979) q[1];
sx q[1];
rz(-1.8760257) q[1];
sx q[1];
rz(-1.6294999) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6544615) q[3];
sx q[3];
rz(-2.122195) q[3];
sx q[3];
rz(-2.7745807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7288397) q[2];
sx q[2];
rz(-1.3430026) q[2];
sx q[2];
rz(-0.084224852) q[2];
rz(1.6526875) q[3];
sx q[3];
rz(-0.050857734) q[3];
sx q[3];
rz(-0.97676718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.65234891) q[0];
sx q[0];
rz(-0.69559613) q[0];
sx q[0];
rz(-3.1074281) q[0];
rz(1.4802406) q[1];
sx q[1];
rz(-0.40373293) q[1];
sx q[1];
rz(1.9306978) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4896442) q[0];
sx q[0];
rz(-1.9495954) q[0];
sx q[0];
rz(-0.31812152) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53721526) q[2];
sx q[2];
rz(-2.5491733) q[2];
sx q[2];
rz(-2.7835961) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38773266) q[1];
sx q[1];
rz(-1.4974623) q[1];
sx q[1];
rz(-2.0078986) q[1];
x q[2];
rz(1.3493531) q[3];
sx q[3];
rz(-1.4998271) q[3];
sx q[3];
rz(-2.8896795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5898798) q[2];
sx q[2];
rz(-2.5248542) q[2];
sx q[2];
rz(1.4025154) q[2];
rz(-1.2276522) q[3];
sx q[3];
rz(-2.472885) q[3];
sx q[3];
rz(-0.68784586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27125204) q[0];
sx q[0];
rz(-1.4372062) q[0];
sx q[0];
rz(1.2649076) q[0];
rz(1.854031) q[1];
sx q[1];
rz(-0.93390211) q[1];
sx q[1];
rz(0.63281473) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70775565) q[0];
sx q[0];
rz(-0.3908866) q[0];
sx q[0];
rz(-0.99796064) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1870474) q[2];
sx q[2];
rz(-0.9034268) q[2];
sx q[2];
rz(-1.2000893) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8532475) q[1];
sx q[1];
rz(-1.4215905) q[1];
sx q[1];
rz(0.38686727) q[1];
x q[2];
rz(-0.37042494) q[3];
sx q[3];
rz(-1.0492477) q[3];
sx q[3];
rz(-2.5481147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5536993) q[2];
sx q[2];
rz(-2.9554695) q[2];
sx q[2];
rz(-2.5659989) q[2];
rz(1.1147739) q[3];
sx q[3];
rz(-1.9284748) q[3];
sx q[3];
rz(2.7093844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1681528) q[0];
sx q[0];
rz(-1.7519209) q[0];
sx q[0];
rz(2.6218276) q[0];
rz(1.4555367) q[1];
sx q[1];
rz(-0.74570233) q[1];
sx q[1];
rz(-2.0533452) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1092508) q[0];
sx q[0];
rz(-2.2025488) q[0];
sx q[0];
rz(2.0103309) q[0];
rz(-pi) q[1];
x q[1];
rz(2.112816) q[2];
sx q[2];
rz(-1.7732829) q[2];
sx q[2];
rz(0.40371343) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.714915) q[1];
sx q[1];
rz(-1.2132056) q[1];
sx q[1];
rz(2.5044051) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3081) q[3];
sx q[3];
rz(-0.98219959) q[3];
sx q[3];
rz(-1.780542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.18152848) q[2];
sx q[2];
rz(-2.0025415) q[2];
sx q[2];
rz(-3.0354434) q[2];
rz(1.6541121) q[3];
sx q[3];
rz(-2.5663576) q[3];
sx q[3];
rz(2.9851448) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0625668) q[0];
sx q[0];
rz(-1.2747958) q[0];
sx q[0];
rz(-1.6612735) q[0];
rz(1.5777499) q[1];
sx q[1];
rz(-0.5636951) q[1];
sx q[1];
rz(3.1196583) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3907174) q[0];
sx q[0];
rz(-0.97033935) q[0];
sx q[0];
rz(-0.12156528) q[0];
x q[1];
rz(1.1925132) q[2];
sx q[2];
rz(-1.6334051) q[2];
sx q[2];
rz(-0.68990842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1749461) q[1];
sx q[1];
rz(-1.7632329) q[1];
sx q[1];
rz(2.220146) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.094385967) q[3];
sx q[3];
rz(-1.1991074) q[3];
sx q[3];
rz(-0.92640141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1070626) q[2];
sx q[2];
rz(-1.7835534) q[2];
sx q[2];
rz(-3.0800842) q[2];
rz(0.48418489) q[3];
sx q[3];
rz(-1.8248841) q[3];
sx q[3];
rz(0.22783247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.43559647) q[0];
sx q[0];
rz(-1.0974925) q[0];
sx q[0];
rz(-0.47501269) q[0];
rz(-1.8571521) q[1];
sx q[1];
rz(-0.78452763) q[1];
sx q[1];
rz(2.6433943) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45979701) q[0];
sx q[0];
rz(-1.4180935) q[0];
sx q[0];
rz(0.55241779) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.118843) q[2];
sx q[2];
rz(-2.2265847) q[2];
sx q[2];
rz(-1.5178258) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.47583844) q[1];
sx q[1];
rz(-1.7301754) q[1];
sx q[1];
rz(-1.5965553) q[1];
rz(-pi) q[2];
rz(-0.9356168) q[3];
sx q[3];
rz(-0.86305076) q[3];
sx q[3];
rz(1.3570116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66607928) q[2];
sx q[2];
rz(-2.3011415) q[2];
sx q[2];
rz(-1.6125512) q[2];
rz(-2.9764002) q[3];
sx q[3];
rz(-1.2673204) q[3];
sx q[3];
rz(2.7388549) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2861479) q[0];
sx q[0];
rz(-2.5832472) q[0];
sx q[0];
rz(2.1530491) q[0];
rz(0.63829154) q[1];
sx q[1];
rz(-2.0604362) q[1];
sx q[1];
rz(-3.1054896) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7571018) q[0];
sx q[0];
rz(-0.57390139) q[0];
sx q[0];
rz(-1.6792273) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1616012) q[2];
sx q[2];
rz(-2.0124751) q[2];
sx q[2];
rz(0.74143386) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.83602025) q[1];
sx q[1];
rz(-1.0936833) q[1];
sx q[1];
rz(0.95997196) q[1];
rz(-2.9265981) q[3];
sx q[3];
rz(-1.3500596) q[3];
sx q[3];
rz(-0.97209111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.710076) q[2];
sx q[2];
rz(-0.71147951) q[2];
sx q[2];
rz(2.9284488) q[2];
rz(1.3171116) q[3];
sx q[3];
rz(-2.9340543) q[3];
sx q[3];
rz(-0.801314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8217736) q[0];
sx q[0];
rz(-1.736301) q[0];
sx q[0];
rz(2.2176493) q[0];
rz(-2.3453103) q[1];
sx q[1];
rz(-2.2932107) q[1];
sx q[1];
rz(-0.36416818) q[1];
rz(2.4305565) q[2];
sx q[2];
rz(-1.0461251) q[2];
sx q[2];
rz(0.82991215) q[2];
rz(-2.2880461) q[3];
sx q[3];
rz(-0.86409909) q[3];
sx q[3];
rz(1.7649337) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
