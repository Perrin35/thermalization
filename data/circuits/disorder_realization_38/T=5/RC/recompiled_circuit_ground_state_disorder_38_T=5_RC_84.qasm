OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.225086) q[0];
sx q[0];
rz(-0.14296159) q[0];
sx q[0];
rz(11.077865) q[0];
rz(-2.0090964) q[1];
sx q[1];
rz(5.8512591) q[1];
sx q[1];
rz(6.9195256) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637488) q[0];
sx q[0];
rz(-1.9533469) q[0];
sx q[0];
rz(2.2027459) q[0];
x q[1];
rz(1.1875528) q[2];
sx q[2];
rz(-1.9343209) q[2];
sx q[2];
rz(-0.55628796) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.072187034) q[1];
sx q[1];
rz(-1.9652307) q[1];
sx q[1];
rz(3.0071114) q[1];
x q[2];
rz(0.71041469) q[3];
sx q[3];
rz(-2.1542366) q[3];
sx q[3];
rz(-1.5855011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0296313) q[2];
sx q[2];
rz(-1.1250857) q[2];
sx q[2];
rz(-0.65043989) q[2];
rz(1.8247617) q[3];
sx q[3];
rz(-1.3464758) q[3];
sx q[3];
rz(0.10183798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0541075) q[0];
sx q[0];
rz(-0.58687812) q[0];
sx q[0];
rz(2.6869539) q[0];
rz(3.1184323) q[1];
sx q[1];
rz(-1.4440447) q[1];
sx q[1];
rz(-2.815411) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9232193) q[0];
sx q[0];
rz(-1.0207286) q[0];
sx q[0];
rz(1.108842) q[0];
rz(-pi) q[1];
rz(-2.8447215) q[2];
sx q[2];
rz(-1.977747) q[2];
sx q[2];
rz(-2.1708058) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.61039576) q[1];
sx q[1];
rz(-2.463974) q[1];
sx q[1];
rz(0.68790959) q[1];
rz(-1.7214016) q[3];
sx q[3];
rz(-0.66622996) q[3];
sx q[3];
rz(3.1325983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.0252016) q[2];
sx q[2];
rz(-1.2025183) q[2];
sx q[2];
rz(2.1873059) q[2];
rz(-1.8918234) q[3];
sx q[3];
rz(-0.3229177) q[3];
sx q[3];
rz(-0.0013110411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.8702451) q[0];
sx q[0];
rz(-1.988669) q[0];
sx q[0];
rz(-1.9648319) q[0];
rz(-0.36508834) q[1];
sx q[1];
rz(-1.3026214) q[1];
sx q[1];
rz(-1.5286068) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1592401) q[0];
sx q[0];
rz(-1.591658) q[0];
sx q[0];
rz(0.013289159) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0625383) q[2];
sx q[2];
rz(-2.0131265) q[2];
sx q[2];
rz(0.7177663) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6916207) q[1];
sx q[1];
rz(-0.1547389) q[1];
sx q[1];
rz(1.3718894) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69832506) q[3];
sx q[3];
rz(-3.0520682) q[3];
sx q[3];
rz(1.9036628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.24885808) q[2];
sx q[2];
rz(-0.55344075) q[2];
sx q[2];
rz(-1.9888606) q[2];
rz(1.8966127) q[3];
sx q[3];
rz(-0.98074073) q[3];
sx q[3];
rz(-2.8583756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1588441) q[0];
sx q[0];
rz(-11*pi/12) q[0];
sx q[0];
rz(-2.4520279) q[0];
rz(-3.116563) q[1];
sx q[1];
rz(-1.6233147) q[1];
sx q[1];
rz(1.4208581) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.384149) q[0];
sx q[0];
rz(-1.674228) q[0];
sx q[0];
rz(-2.5162656) q[0];
rz(-pi) q[1];
x q[1];
rz(0.17949149) q[2];
sx q[2];
rz(-1.9948261) q[2];
sx q[2];
rz(-1.9551203) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1240439) q[1];
sx q[1];
rz(-2.6376403) q[1];
sx q[1];
rz(0.057226463) q[1];
x q[2];
rz(-0.33735621) q[3];
sx q[3];
rz(-1.3135664) q[3];
sx q[3];
rz(-1.6599865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.20874061) q[2];
sx q[2];
rz(-0.17540652) q[2];
sx q[2];
rz(1.9746732) q[2];
rz(2.6973727) q[3];
sx q[3];
rz(-1.2362213) q[3];
sx q[3];
rz(0.3096295) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2883478) q[0];
sx q[0];
rz(-1.3695559) q[0];
sx q[0];
rz(0.41611588) q[0];
rz(-0.5101282) q[1];
sx q[1];
rz(-2.3704539) q[1];
sx q[1];
rz(-2.3663734) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3907523) q[0];
sx q[0];
rz(-1.3715944) q[0];
sx q[0];
rz(0.61934031) q[0];
rz(-1.3890319) q[2];
sx q[2];
rz(-1.9700288) q[2];
sx q[2];
rz(2.4295074) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0803804) q[1];
sx q[1];
rz(-2.9728372) q[1];
sx q[1];
rz(2.3024998) q[1];
rz(-pi) q[2];
rz(0.59501641) q[3];
sx q[3];
rz(-1.2538246) q[3];
sx q[3];
rz(1.9645312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9267209) q[2];
sx q[2];
rz(-2.1472609) q[2];
sx q[2];
rz(-2.1057687) q[2];
rz(0.18946798) q[3];
sx q[3];
rz(-0.47179705) q[3];
sx q[3];
rz(0.50260472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8871317) q[0];
sx q[0];
rz(-3.0939565) q[0];
sx q[0];
rz(2.7857842) q[0];
rz(-1.6461146) q[1];
sx q[1];
rz(-0.9551841) q[1];
sx q[1];
rz(-1.1411508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8320719) q[0];
sx q[0];
rz(-2.4395925) q[0];
sx q[0];
rz(-1.899748) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5579434) q[2];
sx q[2];
rz(-2.5299978) q[2];
sx q[2];
rz(-2.0919959) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43737632) q[1];
sx q[1];
rz(-2.9659038) q[1];
sx q[1];
rz(0.64196569) q[1];
rz(-pi) q[2];
rz(-1.136888) q[3];
sx q[3];
rz(-1.5818137) q[3];
sx q[3];
rz(0.20996717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5100539) q[2];
sx q[2];
rz(-0.71439356) q[2];
sx q[2];
rz(-1.7047403) q[2];
rz(-2.0884183) q[3];
sx q[3];
rz(-1.5318233) q[3];
sx q[3];
rz(0.50301445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18008867) q[0];
sx q[0];
rz(-0.29619521) q[0];
sx q[0];
rz(-0.12251138) q[0];
rz(2.4854614) q[1];
sx q[1];
rz(-1.8729112) q[1];
sx q[1];
rz(1.8001385) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.487566) q[0];
sx q[0];
rz(-1.6161421) q[0];
sx q[0];
rz(1.5349755) q[0];
rz(-pi) q[1];
rz(2.7846863) q[2];
sx q[2];
rz(-0.30115899) q[2];
sx q[2];
rz(-0.44277175) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16640284) q[1];
sx q[1];
rz(-2.4619048) q[1];
sx q[1];
rz(1.1389334) q[1];
rz(-pi) q[2];
x q[2];
rz(0.44388598) q[3];
sx q[3];
rz(-1.6735014) q[3];
sx q[3];
rz(-3.1004124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.573367) q[2];
sx q[2];
rz(-1.919701) q[2];
sx q[2];
rz(-1.6778256) q[2];
rz(-0.73417869) q[3];
sx q[3];
rz(-2.8575183) q[3];
sx q[3];
rz(-1.3045788) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7952591) q[0];
sx q[0];
rz(-1.7743552) q[0];
sx q[0];
rz(-2.6829868) q[0];
rz(1.0147702) q[1];
sx q[1];
rz(-1.9472803) q[1];
sx q[1];
rz(-1.0626622) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.748007) q[0];
sx q[0];
rz(-2.9655122) q[0];
sx q[0];
rz(-2.8066471) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0222715) q[2];
sx q[2];
rz(-0.8415701) q[2];
sx q[2];
rz(2.4265576) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6148551) q[1];
sx q[1];
rz(-1.4397845) q[1];
sx q[1];
rz(0.06191555) q[1];
rz(-pi) q[2];
rz(1.0414391) q[3];
sx q[3];
rz(-2.5564407) q[3];
sx q[3];
rz(0.063466788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8899272) q[2];
sx q[2];
rz(-1.3883611) q[2];
sx q[2];
rz(-1.7903222) q[2];
rz(1.2190367) q[3];
sx q[3];
rz(-2.7128897) q[3];
sx q[3];
rz(-0.8333227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.079000533) q[0];
sx q[0];
rz(-1.3840249) q[0];
sx q[0];
rz(-0.90743995) q[0];
rz(-2.9145248) q[1];
sx q[1];
rz(-1.0457958) q[1];
sx q[1];
rz(-0.89920941) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2405906) q[0];
sx q[0];
rz(-2.0792897) q[0];
sx q[0];
rz(2.2700538) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3311798) q[2];
sx q[2];
rz(-0.11356662) q[2];
sx q[2];
rz(1.0134987) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.6688706) q[1];
sx q[1];
rz(-2.5963915) q[1];
sx q[1];
rz(-1.9342058) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5158922) q[3];
sx q[3];
rz(-0.74885741) q[3];
sx q[3];
rz(-0.98837534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9565309) q[2];
sx q[2];
rz(-1.9756292) q[2];
sx q[2];
rz(-0.60066191) q[2];
rz(-2.0447842) q[3];
sx q[3];
rz(-0.59022248) q[3];
sx q[3];
rz(-0.87305951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7507062) q[0];
sx q[0];
rz(-2.0426671) q[0];
sx q[0];
rz(0.98439687) q[0];
rz(-0.4920494) q[1];
sx q[1];
rz(-1.4130519) q[1];
sx q[1];
rz(1.1955998) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3541832) q[0];
sx q[0];
rz(-2.4314779) q[0];
sx q[0];
rz(2.5837333) q[0];
rz(-2.1629647) q[2];
sx q[2];
rz(-2.371863) q[2];
sx q[2];
rz(-0.20537381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8859133) q[1];
sx q[1];
rz(-1.1885841) q[1];
sx q[1];
rz(3.0542733) q[1];
rz(-pi) q[2];
rz(-0.77287425) q[3];
sx q[3];
rz(-1.2307271) q[3];
sx q[3];
rz(3.056385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.055858) q[2];
sx q[2];
rz(-1.6207638) q[2];
sx q[2];
rz(-0.98304191) q[2];
rz(2.8899657) q[3];
sx q[3];
rz(-2.7139137) q[3];
sx q[3];
rz(-1.0035286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49878237) q[0];
sx q[0];
rz(-0.97220535) q[0];
sx q[0];
rz(-2.4844949) q[0];
rz(-1.8819173) q[1];
sx q[1];
rz(-1.4114264) q[1];
sx q[1];
rz(-2.2687601) q[1];
rz(-2.9272407) q[2];
sx q[2];
rz(-0.8690693) q[2];
sx q[2];
rz(-2.2590841) q[2];
rz(-2.3771277) q[3];
sx q[3];
rz(-0.44776147) q[3];
sx q[3];
rz(0.54816435) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
