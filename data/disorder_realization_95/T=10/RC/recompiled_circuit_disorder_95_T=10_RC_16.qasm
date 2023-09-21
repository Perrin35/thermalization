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
rz(0.89755091) q[0];
rz(-0.23437962) q[1];
sx q[1];
rz(-0.27581629) q[1];
sx q[1];
rz(-1.0645359) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7449887) q[0];
sx q[0];
rz(-2.19967) q[0];
sx q[0];
rz(0.086358503) q[0];
x q[1];
rz(2.2034982) q[2];
sx q[2];
rz(-1.860306) q[2];
sx q[2];
rz(-3.0245568) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.41649109) q[1];
sx q[1];
rz(-2.2033268) q[1];
sx q[1];
rz(1.3593258) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1944794) q[3];
sx q[3];
rz(-0.76768657) q[3];
sx q[3];
rz(0.0038113468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.66427461) q[2];
sx q[2];
rz(-2.5085818) q[2];
sx q[2];
rz(-2.409639) q[2];
rz(0.96015635) q[3];
sx q[3];
rz(-0.82257661) q[3];
sx q[3];
rz(1.7094973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(2.9663548) q[0];
sx q[0];
rz(-0.8834928) q[0];
sx q[0];
rz(0.91645855) q[0];
rz(2.6610999) q[1];
sx q[1];
rz(-2.5669211) q[1];
sx q[1];
rz(-2.2629471) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40936138) q[0];
sx q[0];
rz(-1.9525813) q[0];
sx q[0];
rz(-2.2344927) q[0];
x q[1];
rz(1.8354561) q[2];
sx q[2];
rz(-2.172643) q[2];
sx q[2];
rz(-2.2167609) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0658873) q[1];
sx q[1];
rz(-1.9888708) q[1];
sx q[1];
rz(1.4904406) q[1];
rz(-pi) q[2];
rz(3.0733172) q[3];
sx q[3];
rz(-1.8797415) q[3];
sx q[3];
rz(0.49648778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2800704) q[2];
sx q[2];
rz(-1.7333938) q[2];
sx q[2];
rz(-2.4943165) q[2];
rz(0.17368008) q[3];
sx q[3];
rz(-2.0208385) q[3];
sx q[3];
rz(0.15163264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52755255) q[0];
sx q[0];
rz(-0.63967597) q[0];
sx q[0];
rz(-0.24599427) q[0];
rz(1.7315158) q[1];
sx q[1];
rz(-1.9672111) q[1];
sx q[1];
rz(-2.0203967) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58223984) q[0];
sx q[0];
rz(-1.9358578) q[0];
sx q[0];
rz(2.3854726) q[0];
rz(-1.4175225) q[2];
sx q[2];
rz(-2.7445265) q[2];
sx q[2];
rz(0.085445554) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0908302) q[1];
sx q[1];
rz(-1.8659667) q[1];
sx q[1];
rz(0.34543085) q[1];
x q[2];
rz(-1.6882012) q[3];
sx q[3];
rz(-2.061764) q[3];
sx q[3];
rz(-1.1743869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.40144172) q[2];
sx q[2];
rz(-1.4462877) q[2];
sx q[2];
rz(-2.1901954) q[2];
rz(0.65008632) q[3];
sx q[3];
rz(-1.8875467) q[3];
sx q[3];
rz(2.8459809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6275416) q[0];
sx q[0];
rz(-0.61215949) q[0];
sx q[0];
rz(-0.78805584) q[0];
rz(1.0568985) q[1];
sx q[1];
rz(-1.1833271) q[1];
sx q[1];
rz(0.11638164) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7096036) q[0];
sx q[0];
rz(-0.63392144) q[0];
sx q[0];
rz(-0.043179913) q[0];
x q[1];
rz(-1.0668683) q[2];
sx q[2];
rz(-1.908784) q[2];
sx q[2];
rz(-2.4525814) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1981922) q[1];
sx q[1];
rz(-2.0640089) q[1];
sx q[1];
rz(2.0342159) q[1];
x q[2];
rz(-0.096480358) q[3];
sx q[3];
rz(-2.6498142) q[3];
sx q[3];
rz(-1.5887367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6115761) q[2];
sx q[2];
rz(-1.896603) q[2];
sx q[2];
rz(-0.25203618) q[2];
rz(-2.7633372) q[3];
sx q[3];
rz(-0.16246048) q[3];
sx q[3];
rz(-0.3616412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5761121) q[0];
sx q[0];
rz(-1.4030554) q[0];
sx q[0];
rz(-0.53260032) q[0];
rz(-1.700092) q[1];
sx q[1];
rz(-0.36591995) q[1];
sx q[1];
rz(1.9929569) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.085885) q[0];
sx q[0];
rz(-0.62725337) q[0];
sx q[0];
rz(1.6968615) q[0];
x q[1];
rz(-2.2448036) q[2];
sx q[2];
rz(-1.9325581) q[2];
sx q[2];
rz(2.0828431) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3682813) q[1];
sx q[1];
rz(-1.2174264) q[1];
sx q[1];
rz(-1.3694622) q[1];
rz(-pi) q[2];
rz(1.9653737) q[3];
sx q[3];
rz(-0.88486457) q[3];
sx q[3];
rz(1.3057614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.1308412) q[2];
sx q[2];
rz(-0.66117078) q[2];
sx q[2];
rz(1.032069) q[2];
rz(-2.4268835) q[3];
sx q[3];
rz(-1.2811477) q[3];
sx q[3];
rz(0.023795279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
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
rz(-2.5500568) q[0];
sx q[0];
rz(-1.93601) q[0];
sx q[0];
rz(-2.7456039) q[0];
rz(-1.6962601) q[1];
sx q[1];
rz(-1.6991801) q[1];
sx q[1];
rz(-2.9352303) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2287184) q[0];
sx q[0];
rz(-1.4836856) q[0];
sx q[0];
rz(-2.3939783) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3328853) q[2];
sx q[2];
rz(-0.84918298) q[2];
sx q[2];
rz(-1.8546113) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.57950912) q[1];
sx q[1];
rz(-1.4663896) q[1];
sx q[1];
rz(-0.097192055) q[1];
rz(-1.2431074) q[3];
sx q[3];
rz(-1.7565691) q[3];
sx q[3];
rz(2.1873308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8966184) q[0];
sx q[0];
rz(-2.0386319) q[0];
sx q[0];
rz(-0.57624972) q[0];
rz(-2.9684864) q[1];
sx q[1];
rz(-2.3997967) q[1];
sx q[1];
rz(-1.9304088) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21748397) q[0];
sx q[0];
rz(-2.3586914) q[0];
sx q[0];
rz(-2.3379675) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4326101) q[2];
sx q[2];
rz(-2.0565363) q[2];
sx q[2];
rz(1.1952343) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8163029) q[1];
sx q[1];
rz(-1.6386697) q[1];
sx q[1];
rz(-1.4142649) q[1];
rz(-pi) q[2];
rz(-2.022642) q[3];
sx q[3];
rz(-2.4882836) q[3];
sx q[3];
rz(-0.89316955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8048191) q[2];
sx q[2];
rz(-0.65827289) q[2];
sx q[2];
rz(-3.1170735) q[2];
rz(2.426614) q[3];
sx q[3];
rz(-1.67778) q[3];
sx q[3];
rz(1.5589176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1361168) q[0];
sx q[0];
rz(-1.550721) q[0];
sx q[0];
rz(2.1210282) q[0];
rz(2.9868946) q[1];
sx q[1];
rz(-1.5203412) q[1];
sx q[1];
rz(1.9205836) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8120136) q[0];
sx q[0];
rz(-2.3276969) q[0];
sx q[0];
rz(-1.194186) q[0];
x q[1];
rz(2.9646655) q[2];
sx q[2];
rz(-1.7311586) q[2];
sx q[2];
rz(-0.71080506) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.78302661) q[1];
sx q[1];
rz(-1.4604124) q[1];
sx q[1];
rz(-0.13659887) q[1];
rz(-0.33631781) q[3];
sx q[3];
rz(-2.291403) q[3];
sx q[3];
rz(1.1433126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8528379) q[2];
sx q[2];
rz(-1.4932262) q[2];
sx q[2];
rz(-0.70518804) q[2];
rz(-0.60338902) q[3];
sx q[3];
rz(-0.861895) q[3];
sx q[3];
rz(-1.5195297) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.110638) q[0];
sx q[0];
rz(-1.5752666) q[0];
sx q[0];
rz(0.13701339) q[0];
rz(-0.6048454) q[1];
sx q[1];
rz(-0.73692656) q[1];
sx q[1];
rz(0.12577122) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14995689) q[0];
sx q[0];
rz(-1.7466674) q[0];
sx q[0];
rz(-1.3791023) q[0];
rz(-pi) q[1];
rz(2.618082) q[2];
sx q[2];
rz(-0.38041174) q[2];
sx q[2];
rz(1.3695804) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.51024918) q[1];
sx q[1];
rz(-1.2878839) q[1];
sx q[1];
rz(-2.3307073) q[1];
rz(2.845876) q[3];
sx q[3];
rz(-1.7641281) q[3];
sx q[3];
rz(2.8545692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.003309) q[2];
sx q[2];
rz(-0.97390276) q[2];
sx q[2];
rz(-0.70927817) q[2];
rz(-2.5214031) q[3];
sx q[3];
rz(-0.87738335) q[3];
sx q[3];
rz(0.15670776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9897292) q[0];
sx q[0];
rz(-2.3487838) q[0];
sx q[0];
rz(-1.9412769) q[0];
rz(-0.26750803) q[1];
sx q[1];
rz(-2.2876883) q[1];
sx q[1];
rz(-2.0013924) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2518371) q[0];
sx q[0];
rz(-1.6701084) q[0];
sx q[0];
rz(-1.3072144) q[0];
rz(-pi) q[1];
rz(-2.033038) q[2];
sx q[2];
rz(-1.8046364) q[2];
sx q[2];
rz(0.086694593) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4789341) q[1];
sx q[1];
rz(-1.3824029) q[1];
sx q[1];
rz(1.4807832) q[1];
rz(-2.8285698) q[3];
sx q[3];
rz(-0.81837294) q[3];
sx q[3];
rz(-2.5767874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7662979) q[2];
sx q[2];
rz(-1.6833498) q[2];
sx q[2];
rz(-2.2429788) q[2];
rz(-0.0071772655) q[3];
sx q[3];
rz(-0.74331784) q[3];
sx q[3];
rz(-1.3557419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89467775) q[0];
sx q[0];
rz(-1.4151731) q[0];
sx q[0];
rz(-1.7807501) q[0];
rz(-2.6208411) q[1];
sx q[1];
rz(-1.3856577) q[1];
sx q[1];
rz(1.7210977) q[1];
rz(0.49171369) q[2];
sx q[2];
rz(-1.9946873) q[2];
sx q[2];
rz(-0.36110525) q[2];
rz(-2.5557774) q[3];
sx q[3];
rz(-1.2231493) q[3];
sx q[3];
rz(2.3412658) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
