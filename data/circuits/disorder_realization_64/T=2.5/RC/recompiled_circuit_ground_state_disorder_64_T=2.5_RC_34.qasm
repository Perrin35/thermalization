OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.46848133) q[0];
sx q[0];
rz(3.0335479) q[0];
sx q[0];
rz(11.57668) q[0];
rz(1.5653079) q[1];
sx q[1];
rz(-1.5578527) q[1];
sx q[1];
rz(-1.4563814) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0385168) q[0];
sx q[0];
rz(-1.4355735) q[0];
sx q[0];
rz(1.200202) q[0];
x q[1];
rz(-1.5338495) q[2];
sx q[2];
rz(-1.7040898) q[2];
sx q[2];
rz(0.61090602) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.072468199) q[1];
sx q[1];
rz(-2.1241423) q[1];
sx q[1];
rz(2.8072559) q[1];
x q[2];
rz(0.043660284) q[3];
sx q[3];
rz(-2.147445) q[3];
sx q[3];
rz(0.94983627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9911389) q[2];
sx q[2];
rz(-0.011971124) q[2];
sx q[2];
rz(1.0452622) q[2];
rz(0.996905) q[3];
sx q[3];
rz(-0.0054587047) q[3];
sx q[3];
rz(-1.3158984) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5528706) q[0];
sx q[0];
rz(-1.2405688) q[0];
sx q[0];
rz(1.3528104) q[0];
rz(-0.040933985) q[1];
sx q[1];
rz(-1.9238238) q[1];
sx q[1];
rz(-1.5997684) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41087655) q[0];
sx q[0];
rz(-1.7109032) q[0];
sx q[0];
rz(0.1430169) q[0];
x q[1];
rz(3.1253042) q[2];
sx q[2];
rz(-1.5495346) q[2];
sx q[2];
rz(-2.4580815) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.99166283) q[1];
sx q[1];
rz(-0.59715358) q[1];
sx q[1];
rz(-1.5157022) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1349234) q[3];
sx q[3];
rz(-0.65186497) q[3];
sx q[3];
rz(1.8396371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.505595) q[2];
sx q[2];
rz(-3.1090241) q[2];
sx q[2];
rz(0.51844281) q[2];
rz(2.017766) q[3];
sx q[3];
rz(-2.3321407) q[3];
sx q[3];
rz(-2.7037485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464624) q[0];
sx q[0];
rz(-0.050124425) q[0];
sx q[0];
rz(-1.6019524) q[0];
rz(0.72499544) q[1];
sx q[1];
rz(-0.031818964) q[1];
sx q[1];
rz(-0.67951387) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5711229) q[0];
sx q[0];
rz(-0.75665604) q[0];
sx q[0];
rz(2.4837912) q[0];
rz(-pi) q[1];
rz(2.9271487) q[2];
sx q[2];
rz(-1.5661998) q[2];
sx q[2];
rz(1.4397804) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.589349) q[1];
sx q[1];
rz(-1.8830944) q[1];
sx q[1];
rz(2.8683788) q[1];
x q[2];
rz(-1.1876538) q[3];
sx q[3];
rz(-2.5184298) q[3];
sx q[3];
rz(-0.23322695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9069549) q[2];
sx q[2];
rz(-0.24181557) q[2];
sx q[2];
rz(0.50092906) q[2];
rz(0.45589724) q[3];
sx q[3];
rz(-3.1152476) q[3];
sx q[3];
rz(-0.19550368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86922115) q[0];
sx q[0];
rz(-3.0932194) q[0];
sx q[0];
rz(-0.79917556) q[0];
rz(-2.8436106) q[1];
sx q[1];
rz(-2.8831392) q[1];
sx q[1];
rz(2.2395649) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.395754) q[0];
sx q[0];
rz(-0.82358783) q[0];
sx q[0];
rz(1.7584778) q[0];
rz(-2.600014) q[2];
sx q[2];
rz(-1.7006093) q[2];
sx q[2];
rz(-2.1101348) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8405172) q[1];
sx q[1];
rz(-1.7184192) q[1];
sx q[1];
rz(-0.0096489659) q[1];
rz(-pi) q[2];
rz(-0.055156339) q[3];
sx q[3];
rz(-1.6481019) q[3];
sx q[3];
rz(-2.5173924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4681089) q[2];
sx q[2];
rz(-3.1099042) q[2];
sx q[2];
rz(1.6710949) q[2];
rz(2.9521613) q[3];
sx q[3];
rz(-0.048308689) q[3];
sx q[3];
rz(2.3725177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22537941) q[0];
sx q[0];
rz(-3.0140641) q[0];
sx q[0];
rz(-0.39176971) q[0];
rz(-1.0852934) q[1];
sx q[1];
rz(-0.009805209) q[1];
sx q[1];
rz(-0.36878961) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1802496) q[0];
sx q[0];
rz(-1.215544) q[0];
sx q[0];
rz(2.536098) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1283705) q[2];
sx q[2];
rz(-0.48071024) q[2];
sx q[2];
rz(1.7465357) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49807355) q[1];
sx q[1];
rz(-0.0066702492) q[1];
sx q[1];
rz(1.4590864) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4318352) q[3];
sx q[3];
rz(-1.3922979) q[3];
sx q[3];
rz(0.75542007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.58086777) q[2];
sx q[2];
rz(-0.11595011) q[2];
sx q[2];
rz(1.7725393) q[2];
rz(-3.0154058) q[3];
sx q[3];
rz(-0.43712619) q[3];
sx q[3];
rz(2.0807467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5996025) q[0];
sx q[0];
rz(-0.76102155) q[0];
sx q[0];
rz(1.5916995) q[0];
rz(2.1688993) q[1];
sx q[1];
rz(-0.21316554) q[1];
sx q[1];
rz(-1.0290283) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81188503) q[0];
sx q[0];
rz(-0.68674478) q[0];
sx q[0];
rz(-2.6558557) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9036361) q[2];
sx q[2];
rz(-1.4334442) q[2];
sx q[2];
rz(-2.1417625) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29063836) q[1];
sx q[1];
rz(-0.096157638) q[1];
sx q[1];
rz(-1.5753217) q[1];
rz(-0.036110445) q[3];
sx q[3];
rz(-1.7381769) q[3];
sx q[3];
rz(-3.0073364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.782393) q[2];
sx q[2];
rz(-1.426037) q[2];
sx q[2];
rz(-2.7101809) q[2];
rz(-0.99724489) q[3];
sx q[3];
rz(-3.1044208) q[3];
sx q[3];
rz(0.55716151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.7998841) q[0];
sx q[0];
rz(-0.48483098) q[0];
sx q[0];
rz(-2.0836015) q[0];
rz(0.82408389) q[1];
sx q[1];
rz(-3.1415756) q[1];
sx q[1];
rz(-0.81738671) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5800487) q[0];
sx q[0];
rz(-1.9684682) q[0];
sx q[0];
rz(-0.22821594) q[0];
rz(-pi) q[1];
rz(-1.5777052) q[2];
sx q[2];
rz(-1.5890317) q[2];
sx q[2];
rz(-1.4438786) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.89725607) q[1];
sx q[1];
rz(-0.24459363) q[1];
sx q[1];
rz(1.7331428) q[1];
x q[2];
rz(0.55900662) q[3];
sx q[3];
rz(-0.83448258) q[3];
sx q[3];
rz(-2.112039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.99007964) q[2];
sx q[2];
rz(-1.238287) q[2];
sx q[2];
rz(1.5129169) q[2];
rz(2.7401636) q[3];
sx q[3];
rz(-0.028086834) q[3];
sx q[3];
rz(-0.95429558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7734739) q[0];
sx q[0];
rz(-3.0768657) q[0];
sx q[0];
rz(-1.3637967) q[0];
rz(-3.0996481) q[1];
sx q[1];
rz(-0.13101235) q[1];
sx q[1];
rz(-2.1379474) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.100228) q[0];
sx q[0];
rz(-1.4778839) q[0];
sx q[0];
rz(2.5071457) q[0];
rz(-pi) q[1];
rz(3.1334468) q[2];
sx q[2];
rz(-1.9081433) q[2];
sx q[2];
rz(-1.0376736) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0692301) q[1];
sx q[1];
rz(-1.8025928) q[1];
sx q[1];
rz(0.20834558) q[1];
rz(1.8742626) q[3];
sx q[3];
rz(-0.93765536) q[3];
sx q[3];
rz(-1.9078524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7046788) q[2];
sx q[2];
rz(-3.0824326) q[2];
sx q[2];
rz(1.7436279) q[2];
rz(0.49020234) q[3];
sx q[3];
rz(-0.043488113) q[3];
sx q[3];
rz(1.9628261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.040084664) q[0];
sx q[0];
rz(-0.12797102) q[0];
sx q[0];
rz(-0.21139938) q[0];
rz(-1.118411) q[1];
sx q[1];
rz(-0.011592955) q[1];
sx q[1];
rz(-1.6308019) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3113041) q[0];
sx q[0];
rz(-1.2823866) q[0];
sx q[0];
rz(-2.8187241) q[0];
x q[1];
rz(0.57447471) q[2];
sx q[2];
rz(-2.4316638) q[2];
sx q[2];
rz(-2.1394875) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.80434752) q[1];
sx q[1];
rz(-1.7253863) q[1];
sx q[1];
rz(-1.5881601) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6602732) q[3];
sx q[3];
rz(-1.4050639) q[3];
sx q[3];
rz(0.16752088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7819034) q[2];
sx q[2];
rz(-3.1236881) q[2];
sx q[2];
rz(2.8622799) q[2];
rz(0.8638047) q[3];
sx q[3];
rz(-3.1371208) q[3];
sx q[3];
rz(-2.4590676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.913468) q[0];
sx q[0];
rz(-2.0371912) q[0];
sx q[0];
rz(1.8260691) q[0];
rz(2.3663991) q[1];
sx q[1];
rz(-2.9029791) q[1];
sx q[1];
rz(-1.7528037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0553) q[0];
sx q[0];
rz(-2.820083) q[0];
sx q[0];
rz(-2.6525524) q[0];
rz(-pi) q[1];
rz(-1.3451683) q[2];
sx q[2];
rz(-1.1755845) q[2];
sx q[2];
rz(-1.164142) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5701616) q[1];
sx q[1];
rz(-1.5699016) q[1];
sx q[1];
rz(-1.5714297) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3540415) q[3];
sx q[3];
rz(-2.5803704) q[3];
sx q[3];
rz(2.7360327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3751601) q[2];
sx q[2];
rz(-3.1258686) q[2];
sx q[2];
rz(1.6895705) q[2];
rz(2.8826513) q[3];
sx q[3];
rz(-2.9539234) q[3];
sx q[3];
rz(1.9848721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6257085) q[0];
sx q[0];
rz(-0.71737552) q[0];
sx q[0];
rz(1.4364568) q[0];
rz(1.6547849) q[1];
sx q[1];
rz(-2.868352) q[1];
sx q[1];
rz(0.19788338) q[1];
rz(-1.5639772) q[2];
sx q[2];
rz(-1.7710254) q[2];
sx q[2];
rz(0.22631021) q[2];
rz(-3.1155125) q[3];
sx q[3];
rz(-1.9175538) q[3];
sx q[3];
rz(-3.1160311) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
