OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.1420105) q[0];
sx q[0];
rz(-2.1394696) q[0];
sx q[0];
rz(-2.2440417) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(2.0770567) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9165186) q[0];
sx q[0];
rz(-1.5009891) q[0];
sx q[0];
rz(-0.9401456) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2034982) q[2];
sx q[2];
rz(-1.860306) q[2];
sx q[2];
rz(-3.0245568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.41649109) q[1];
sx q[1];
rz(-0.93826586) q[1];
sx q[1];
rz(1.7822669) q[1];
x q[2];
rz(1.7552056) q[3];
sx q[3];
rz(-0.82114906) q[3];
sx q[3];
rz(-2.8781995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66427461) q[2];
sx q[2];
rz(-0.6330108) q[2];
sx q[2];
rz(0.73195362) q[2];
rz(2.1814363) q[3];
sx q[3];
rz(-2.319016) q[3];
sx q[3];
rz(-1.4320954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(-0.91645855) q[0];
rz(0.48049277) q[1];
sx q[1];
rz(-0.57467159) q[1];
sx q[1];
rz(0.8786456) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7322313) q[0];
sx q[0];
rz(-1.1890113) q[0];
sx q[0];
rz(-2.2344927) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7777113) q[2];
sx q[2];
rz(-2.4907787) q[2];
sx q[2];
rz(1.770307) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.87989992) q[1];
sx q[1];
rz(-0.42527929) q[1];
sx q[1];
rz(0.17875032) q[1];
rz(-1.3602123) q[3];
sx q[3];
rz(-2.8254291) q[3];
sx q[3];
rz(0.27526835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8615222) q[2];
sx q[2];
rz(-1.4081988) q[2];
sx q[2];
rz(0.64727616) q[2];
rz(0.17368008) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(-2.98996) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(-0.24599427) q[0];
rz(1.4100769) q[1];
sx q[1];
rz(-1.1743816) q[1];
sx q[1];
rz(-2.0203967) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3504346) q[0];
sx q[0];
rz(-2.3179623) q[0];
sx q[0];
rz(-0.50823786) q[0];
x q[1];
rz(-3.0776575) q[2];
sx q[2];
rz(-1.1786412) q[2];
sx q[2];
rz(0.25142297) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.725863) q[1];
sx q[1];
rz(-1.2408857) q[1];
sx q[1];
rz(1.2582474) q[1];
rz(-2.6477473) q[3];
sx q[3];
rz(-1.6742799) q[3];
sx q[3];
rz(0.34085694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(2.1901954) q[2];
rz(0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51405108) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(0.78805584) q[0];
rz(2.0846941) q[1];
sx q[1];
rz(-1.9582656) q[1];
sx q[1];
rz(-3.025211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43198904) q[0];
sx q[0];
rz(-2.5076712) q[0];
sx q[0];
rz(-0.043179913) q[0];
rz(-pi) q[1];
rz(-2.7599081) q[2];
sx q[2];
rz(-1.0978062) q[2];
sx q[2];
rz(-0.7009398) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9434005) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(-2.0342159) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6223525) q[3];
sx q[3];
rz(-1.0815074) q[3];
sx q[3];
rz(1.4435022) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5300166) q[2];
sx q[2];
rz(-1.2449896) q[2];
sx q[2];
rz(-2.8895565) q[2];
rz(0.37825545) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(2.7799515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56548059) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(2.6089923) q[0];
rz(1.700092) q[1];
sx q[1];
rz(-2.7756727) q[1];
sx q[1];
rz(-1.1486357) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0557077) q[0];
sx q[0];
rz(-0.62725337) q[0];
sx q[0];
rz(-1.4447312) q[0];
rz(-pi) q[1];
rz(2.2448036) q[2];
sx q[2];
rz(-1.2090346) q[2];
sx q[2];
rz(-1.0587495) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7733113) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(-1.7721304) q[1];
x q[2];
rz(0.72539056) q[3];
sx q[3];
rz(-1.2687506) q[3];
sx q[3];
rz(0.52291742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0107515) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(-2.1095236) q[2];
rz(2.4268835) q[3];
sx q[3];
rz(-1.8604449) q[3];
sx q[3];
rz(-3.1177974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.5500568) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(0.39598879) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.4424125) q[1];
sx q[1];
rz(-0.2063624) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43564046) q[0];
sx q[0];
rz(-0.75169509) q[0];
sx q[0];
rz(-0.12775001) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.80870734) q[2];
sx q[2];
rz(-0.84918298) q[2];
sx q[2];
rz(-1.2869814) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.98112647) q[1];
sx q[1];
rz(-1.6674575) q[1];
sx q[1];
rz(-1.6756945) q[1];
x q[2];
rz(2.0993125) q[3];
sx q[3];
rz(-2.7665666) q[3];
sx q[3];
rz(0.11881766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.6017194) q[2];
sx q[2];
rz(-2.0809934) q[2];
sx q[2];
rz(0.27077857) q[2];
rz(-0.74832908) q[3];
sx q[3];
rz(-1.3724519) q[3];
sx q[3];
rz(2.0628827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2449743) q[0];
sx q[0];
rz(-1.1029607) q[0];
sx q[0];
rz(0.57624972) q[0];
rz(-0.17310625) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(-1.2111838) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9528708) q[0];
sx q[0];
rz(-2.0824008) q[0];
sx q[0];
rz(0.94922025) q[0];
rz(-pi) q[1];
rz(-1.7089825) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(-1.9463584) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8853828) q[1];
sx q[1];
rz(-1.7269644) q[1];
sx q[1];
rz(-0.068710879) q[1];
x q[2];
rz(1.1189506) q[3];
sx q[3];
rz(-0.65330905) q[3];
sx q[3];
rz(-2.2484231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(-3.1170735) q[2];
rz(-2.426614) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(-1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0054758469) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(-2.1210282) q[0];
rz(-0.15469805) q[1];
sx q[1];
rz(-1.6212515) q[1];
sx q[1];
rz(1.221009) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1655501) q[0];
sx q[0];
rz(-1.8414458) q[0];
sx q[0];
rz(-0.793215) q[0];
x q[1];
rz(-0.74322015) q[2];
sx q[2];
rz(-0.23822242) q[2];
sx q[2];
rz(-1.5889578) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3386821) q[1];
sx q[1];
rz(-1.7065587) q[1];
sx q[1];
rz(-1.6822097) q[1];
rz(2.8052748) q[3];
sx q[3];
rz(-2.291403) q[3];
sx q[3];
rz(-1.99828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8528379) q[2];
sx q[2];
rz(-1.6483665) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(2.5382036) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(1.6220629) q[3];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.110638) q[0];
sx q[0];
rz(-1.5663261) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(0.6048454) q[1];
sx q[1];
rz(-2.4046661) q[1];
sx q[1];
rz(0.12577122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14995689) q[0];
sx q[0];
rz(-1.3949252) q[0];
sx q[0];
rz(-1.7624904) q[0];
rz(-pi) q[1];
rz(1.3734829) q[2];
sx q[2];
rz(-1.243405) q[2];
sx q[2];
rz(1.2158074) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8217433) q[1];
sx q[1];
rz(-0.84801596) q[1];
sx q[1];
rz(-0.38139947) q[1];
x q[2];
rz(-1.7726692) q[3];
sx q[3];
rz(-1.2807506) q[3];
sx q[3];
rz(-1.9162852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13828364) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(2.4323145) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-2.2642093) q[3];
sx q[3];
rz(-0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.1518635) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(1.9412769) q[0];
rz(0.26750803) q[1];
sx q[1];
rz(-0.8539044) q[1];
sx q[1];
rz(1.1402003) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70779078) q[0];
sx q[0];
rz(-1.3085438) q[0];
sx q[0];
rz(3.0387525) q[0];
x q[1];
rz(1.1085547) q[2];
sx q[2];
rz(-1.8046364) q[2];
sx q[2];
rz(-3.0548981) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.4789341) q[1];
sx q[1];
rz(-1.7591898) q[1];
sx q[1];
rz(1.4807832) q[1];
rz(1.8885918) q[3];
sx q[3];
rz(-0.80298775) q[3];
sx q[3];
rz(3.0190937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7662979) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-0.89861384) q[2];
rz(0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.7858508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2469149) q[0];
sx q[0];
rz(-1.7264195) q[0];
sx q[0];
rz(1.3608426) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(2.649879) q[2];
sx q[2];
rz(-1.1469054) q[2];
sx q[2];
rz(2.7804874) q[2];
rz(-0.58581523) q[3];
sx q[3];
rz(-1.9184434) q[3];
sx q[3];
rz(-0.8003269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];