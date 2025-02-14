OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.34831369) q[0];
sx q[0];
rz(8.0170595) q[0];
sx q[0];
rz(9.4693139) q[0];
rz(2.0308004) q[1];
sx q[1];
rz(-1.9171311) q[1];
sx q[1];
rz(0.78723025) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6985748) q[0];
sx q[0];
rz(-2.0788143) q[0];
sx q[0];
rz(0.99850151) q[0];
x q[1];
rz(0.039723176) q[2];
sx q[2];
rz(-1.2457704) q[2];
sx q[2];
rz(0.15234337) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8833784) q[1];
sx q[1];
rz(-1.6401263) q[1];
sx q[1];
rz(0.99154559) q[1];
x q[2];
rz(2.9028893) q[3];
sx q[3];
rz(-1.2269964) q[3];
sx q[3];
rz(-1.3581004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3945776) q[2];
sx q[2];
rz(-1.0801103) q[2];
sx q[2];
rz(-0.75418312) q[2];
rz(-1.4150103) q[3];
sx q[3];
rz(-0.49608803) q[3];
sx q[3];
rz(-0.32895857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0873782) q[0];
sx q[0];
rz(-2.7890451) q[0];
sx q[0];
rz(1.3341599) q[0];
rz(-1.2913903) q[1];
sx q[1];
rz(-2.1033557) q[1];
sx q[1];
rz(-2.2659567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.401225) q[0];
sx q[0];
rz(-2.2384907) q[0];
sx q[0];
rz(2.4424548) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4562143) q[2];
sx q[2];
rz(-2.6422524) q[2];
sx q[2];
rz(-2.2407766) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8671682) q[1];
sx q[1];
rz(-0.53571415) q[1];
sx q[1];
rz(2.1650326) q[1];
x q[2];
rz(-0.29107901) q[3];
sx q[3];
rz(-1.3716231) q[3];
sx q[3];
rz(-1.2676304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4448173) q[2];
sx q[2];
rz(-0.70187086) q[2];
sx q[2];
rz(0.87265054) q[2];
rz(-2.3526092) q[3];
sx q[3];
rz(-0.9011457) q[3];
sx q[3];
rz(0.22855973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2149684) q[0];
sx q[0];
rz(-1.7890395) q[0];
sx q[0];
rz(-3.1335926) q[0];
rz(-2.3232715) q[1];
sx q[1];
rz(-2.3921831) q[1];
sx q[1];
rz(1.9047033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2674026) q[0];
sx q[0];
rz(-1.2498901) q[0];
sx q[0];
rz(-0.85808922) q[0];
rz(-pi) q[1];
rz(0.34864595) q[2];
sx q[2];
rz(-1.9816035) q[2];
sx q[2];
rz(-2.2553159) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3171009) q[1];
sx q[1];
rz(-2.503509) q[1];
sx q[1];
rz(0.53148766) q[1];
rz(0.45435702) q[3];
sx q[3];
rz(-2.8792692) q[3];
sx q[3];
rz(-1.9087877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41039738) q[2];
sx q[2];
rz(-3.0072913) q[2];
sx q[2];
rz(-2.2006939) q[2];
rz(2.0375552) q[3];
sx q[3];
rz(-1.3528115) q[3];
sx q[3];
rz(-1.997939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10821548) q[0];
sx q[0];
rz(-2.5984851) q[0];
sx q[0];
rz(0.23750842) q[0];
rz(2.4820651) q[1];
sx q[1];
rz(-0.57412761) q[1];
sx q[1];
rz(3.1245756) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40460038) q[0];
sx q[0];
rz(-1.7395818) q[0];
sx q[0];
rz(-3.0927318) q[0];
rz(-pi) q[1];
rz(-0.41827664) q[2];
sx q[2];
rz(-2.661663) q[2];
sx q[2];
rz(-2.0945702) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3238195) q[1];
sx q[1];
rz(-2.5669879) q[1];
sx q[1];
rz(2.9486604) q[1];
x q[2];
rz(2.4260826) q[3];
sx q[3];
rz(-2.1984716) q[3];
sx q[3];
rz(-2.024141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90618769) q[2];
sx q[2];
rz(-1.2013288) q[2];
sx q[2];
rz(2.2271633) q[2];
rz(-0.36214456) q[3];
sx q[3];
rz(-0.80949628) q[3];
sx q[3];
rz(-2.5563498) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2528766) q[0];
sx q[0];
rz(-1.3739561) q[0];
sx q[0];
rz(3.0294898) q[0];
rz(-1.0653227) q[1];
sx q[1];
rz(-2.0707097) q[1];
sx q[1];
rz(-2.0723453) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1493268) q[0];
sx q[0];
rz(-2.7130824) q[0];
sx q[0];
rz(-2.4403768) q[0];
x q[1];
rz(-2.4222071) q[2];
sx q[2];
rz(-0.66443887) q[2];
sx q[2];
rz(1.7094572) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4878975) q[1];
sx q[1];
rz(-1.2733165) q[1];
sx q[1];
rz(1.9031699) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13302044) q[3];
sx q[3];
rz(-0.97571301) q[3];
sx q[3];
rz(1.2979591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4273044) q[2];
sx q[2];
rz(-2.4236743) q[2];
sx q[2];
rz(2.60738) q[2];
rz(-0.30321768) q[3];
sx q[3];
rz(-1.5568308) q[3];
sx q[3];
rz(-1.5155972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41820207) q[0];
sx q[0];
rz(-1.284282) q[0];
sx q[0];
rz(-0.89163017) q[0];
rz(1.3806237) q[1];
sx q[1];
rz(-1.2029519) q[1];
sx q[1];
rz(-2.2183529) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5864201) q[0];
sx q[0];
rz(-1.6229543) q[0];
sx q[0];
rz(-0.098343366) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2963861) q[2];
sx q[2];
rz(-1.1968799) q[2];
sx q[2];
rz(2.3309938) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.90369836) q[1];
sx q[1];
rz(-2.782202) q[1];
sx q[1];
rz(0.87025799) q[1];
rz(-pi) q[2];
rz(-2.4999077) q[3];
sx q[3];
rz(-1.2417792) q[3];
sx q[3];
rz(-1.5050448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.096006958) q[2];
sx q[2];
rz(-0.43232375) q[2];
sx q[2];
rz(1.879479) q[2];
rz(1.1612085) q[3];
sx q[3];
rz(-1.592417) q[3];
sx q[3];
rz(-1.8956634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68814174) q[0];
sx q[0];
rz(-0.83204404) q[0];
sx q[0];
rz(-1.6081109) q[0];
rz(0.43117943) q[1];
sx q[1];
rz(-0.9175514) q[1];
sx q[1];
rz(-1.4088438) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4154177) q[0];
sx q[0];
rz(-1.2841932) q[0];
sx q[0];
rz(2.4917401) q[0];
x q[1];
rz(1.3546014) q[2];
sx q[2];
rz(-0.39297418) q[2];
sx q[2];
rz(0.68902389) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3235493) q[1];
sx q[1];
rz(-2.9406639) q[1];
sx q[1];
rz(2.4646321) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7944144) q[3];
sx q[3];
rz(-2.6287492) q[3];
sx q[3];
rz(-1.7382857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.2048753) q[2];
sx q[2];
rz(-2.016158) q[2];
sx q[2];
rz(1.8219061) q[2];
rz(-0.034505757) q[3];
sx q[3];
rz(-1.5311818) q[3];
sx q[3];
rz(-1.7721133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7337604) q[0];
sx q[0];
rz(-0.5540846) q[0];
sx q[0];
rz(1.9246509) q[0];
rz(-2.2159684) q[1];
sx q[1];
rz(-1.7763014) q[1];
sx q[1];
rz(1.9409174) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2404296) q[0];
sx q[0];
rz(-2.0736742) q[0];
sx q[0];
rz(2.2818185) q[0];
x q[1];
rz(2.578892) q[2];
sx q[2];
rz(-0.68018736) q[2];
sx q[2];
rz(-2.0923751) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7279049) q[1];
sx q[1];
rz(-2.3686045) q[1];
sx q[1];
rz(-0.058900515) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.86494495) q[3];
sx q[3];
rz(-0.53370133) q[3];
sx q[3];
rz(2.7786917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.919148) q[2];
sx q[2];
rz(-1.5347975) q[2];
sx q[2];
rz(-1.0486802) q[2];
rz(2.9564296) q[3];
sx q[3];
rz(-1.984963) q[3];
sx q[3];
rz(-3.0345501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.01934) q[0];
sx q[0];
rz(-0.66146079) q[0];
sx q[0];
rz(-0.81073236) q[0];
rz(-0.73075378) q[1];
sx q[1];
rz(-1.0565051) q[1];
sx q[1];
rz(-2.1810541) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98314032) q[0];
sx q[0];
rz(-3.1380655) q[0];
sx q[0];
rz(-0.28007026) q[0];
x q[1];
rz(0.4308295) q[2];
sx q[2];
rz(-1.1053876) q[2];
sx q[2];
rz(-2.7499466) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62453485) q[1];
sx q[1];
rz(-1.0427999) q[1];
sx q[1];
rz(-2.8370884) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9979565) q[3];
sx q[3];
rz(-1.33516) q[3];
sx q[3];
rz(0.7251265) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9591799) q[2];
sx q[2];
rz(-1.7824087) q[2];
sx q[2];
rz(1.6020927) q[2];
rz(-1.0473853) q[3];
sx q[3];
rz(-1.7644707) q[3];
sx q[3];
rz(-1.2654977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3613116) q[0];
sx q[0];
rz(-0.86532101) q[0];
sx q[0];
rz(-0.58746946) q[0];
rz(2.7777708) q[1];
sx q[1];
rz(-1.9963341) q[1];
sx q[1];
rz(-2.2344373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22940561) q[0];
sx q[0];
rz(-0.36021458) q[0];
sx q[0];
rz(2.9111258) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.67423363) q[2];
sx q[2];
rz(-0.57949726) q[2];
sx q[2];
rz(0.63636875) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.41872596) q[1];
sx q[1];
rz(-1.360824) q[1];
sx q[1];
rz(2.6814815) q[1];
x q[2];
rz(1.7916405) q[3];
sx q[3];
rz(-2.7691602) q[3];
sx q[3];
rz(-0.34162921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.63984799) q[2];
sx q[2];
rz(-3.0526243) q[2];
sx q[2];
rz(-0.42195827) q[2];
rz(-3.1320599) q[3];
sx q[3];
rz(-1.5744753) q[3];
sx q[3];
rz(0.52614051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(-0.52787732) q[0];
sx q[0];
rz(-1.8462702) q[0];
sx q[0];
rz(0.12311002) q[0];
rz(0.70115024) q[1];
sx q[1];
rz(-1.9155365) q[1];
sx q[1];
rz(1.6446) q[1];
rz(0.43177615) q[2];
sx q[2];
rz(-2.3974621) q[2];
sx q[2];
rz(1.5324788) q[2];
rz(-0.99118945) q[3];
sx q[3];
rz(-0.47953987) q[3];
sx q[3];
rz(-2.6399947) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
