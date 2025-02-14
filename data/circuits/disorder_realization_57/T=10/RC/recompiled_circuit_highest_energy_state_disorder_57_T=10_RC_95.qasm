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
rz(1.1433831) q[0];
sx q[0];
rz(-1.5902061) q[0];
sx q[0];
rz(0.32713148) q[0];
rz(1.7786572) q[1];
sx q[1];
rz(6.8805334) q[1];
sx q[1];
rz(9.7272275) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0898092) q[0];
sx q[0];
rz(-2.8642388) q[0];
sx q[0];
rz(-0.92602317) q[0];
rz(-pi) q[1];
x q[1];
rz(0.53989065) q[2];
sx q[2];
rz(-2.1587025) q[2];
sx q[2];
rz(1.3164275) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.11034549) q[1];
sx q[1];
rz(-2.7953138) q[1];
sx q[1];
rz(-1.9900761) q[1];
rz(-1.6749945) q[3];
sx q[3];
rz(-1.0999354) q[3];
sx q[3];
rz(-2.6764118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.7678307) q[2];
sx q[2];
rz(-3.0475898) q[2];
sx q[2];
rz(-2.8685699) q[2];
rz(0.36429575) q[3];
sx q[3];
rz(-1.7184869) q[3];
sx q[3];
rz(-0.020615904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.453422) q[0];
sx q[0];
rz(-1.0978798) q[0];
sx q[0];
rz(1.4349487) q[0];
rz(-2.9393328) q[1];
sx q[1];
rz(-2.5113998) q[1];
sx q[1];
rz(-2.8369301) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7633491) q[0];
sx q[0];
rz(-1.8190666) q[0];
sx q[0];
rz(2.2909478) q[0];
x q[1];
rz(-0.29918798) q[2];
sx q[2];
rz(-1.9721037) q[2];
sx q[2];
rz(0.68986675) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.92485904) q[1];
sx q[1];
rz(-1.2302515) q[1];
sx q[1];
rz(1.3751283) q[1];
rz(-pi) q[2];
rz(-1.4502787) q[3];
sx q[3];
rz(-1.4637865) q[3];
sx q[3];
rz(1.4329239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8649586) q[2];
sx q[2];
rz(-0.43143299) q[2];
sx q[2];
rz(-1.6032093) q[2];
rz(0.81613427) q[3];
sx q[3];
rz(-0.82681257) q[3];
sx q[3];
rz(2.2468467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64033878) q[0];
sx q[0];
rz(-2.8479939) q[0];
sx q[0];
rz(-1.5119934) q[0];
rz(-1.8005796) q[1];
sx q[1];
rz(-0.87313849) q[1];
sx q[1];
rz(0.27892932) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90401201) q[0];
sx q[0];
rz(-1.4424217) q[0];
sx q[0];
rz(-1.4610941) q[0];
rz(0.52742676) q[2];
sx q[2];
rz(-1.0042388) q[2];
sx q[2];
rz(-1.1184831) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.046849068) q[1];
sx q[1];
rz(-1.2654116) q[1];
sx q[1];
rz(0.8253936) q[1];
rz(-pi) q[2];
rz(1.5749212) q[3];
sx q[3];
rz(-2.4253411) q[3];
sx q[3];
rz(0.28820693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5035847) q[2];
sx q[2];
rz(-1.6562485) q[2];
sx q[2];
rz(-2.8032362) q[2];
rz(1.407297) q[3];
sx q[3];
rz(-2.0140078) q[3];
sx q[3];
rz(-2.1710763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19891837) q[0];
sx q[0];
rz(-0.022138683) q[0];
sx q[0];
rz(2.4778147) q[0];
rz(-0.55555073) q[1];
sx q[1];
rz(-1.5063565) q[1];
sx q[1];
rz(-1.3551855) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595918) q[0];
sx q[0];
rz(-2.2170236) q[0];
sx q[0];
rz(-1.610397) q[0];
rz(-0.023472114) q[2];
sx q[2];
rz(-1.9292574) q[2];
sx q[2];
rz(-1.9752928) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.0993886) q[1];
sx q[1];
rz(-0.70800938) q[1];
sx q[1];
rz(0.96570496) q[1];
rz(1.1045214) q[3];
sx q[3];
rz(-1.5111088) q[3];
sx q[3];
rz(-2.3883143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3163471) q[2];
sx q[2];
rz(-2.1630042) q[2];
sx q[2];
rz(-2.6157731) q[2];
rz(0.60872269) q[3];
sx q[3];
rz(-2.5059097) q[3];
sx q[3];
rz(-2.4354602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(2.4290001) q[0];
sx q[0];
rz(-1.1845931) q[0];
sx q[0];
rz(-0.96479601) q[0];
rz(-0.4110128) q[1];
sx q[1];
rz(-0.95714584) q[1];
sx q[1];
rz(2.029665) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72090688) q[0];
sx q[0];
rz(-1.5941275) q[0];
sx q[0];
rz(-3.1258778) q[0];
rz(-pi) q[1];
rz(0.50639407) q[2];
sx q[2];
rz(-0.7686457) q[2];
sx q[2];
rz(-2.9646693) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7644123) q[1];
sx q[1];
rz(-2.0736338) q[1];
sx q[1];
rz(2.2187114) q[1];
rz(-pi) q[2];
rz(-0.49647496) q[3];
sx q[3];
rz(-0.30204812) q[3];
sx q[3];
rz(0.99629096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3045706) q[2];
sx q[2];
rz(-1.3304173) q[2];
sx q[2];
rz(-2.0988317) q[2];
rz(1.9581095) q[3];
sx q[3];
rz(-0.8684929) q[3];
sx q[3];
rz(2.1709501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17827621) q[0];
sx q[0];
rz(-1.2238418) q[0];
sx q[0];
rz(-2.8455612) q[0];
rz(0.0010679642) q[1];
sx q[1];
rz(-1.8031392) q[1];
sx q[1];
rz(1.0329049) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72275513) q[0];
sx q[0];
rz(-1.5745409) q[0];
sx q[0];
rz(0.0054971183) q[0];
rz(0.11992137) q[2];
sx q[2];
rz(-1.9397221) q[2];
sx q[2];
rz(-0.030583965) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1329676) q[1];
sx q[1];
rz(-1.6264249) q[1];
sx q[1];
rz(0.24764307) q[1];
x q[2];
rz(1.9719129) q[3];
sx q[3];
rz(-0.84904957) q[3];
sx q[3];
rz(-0.37506074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.380015) q[2];
sx q[2];
rz(-1.6753316) q[2];
sx q[2];
rz(1.638691) q[2];
rz(0.6014398) q[3];
sx q[3];
rz(-2.1419339) q[3];
sx q[3];
rz(-2.7433266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10694207) q[0];
sx q[0];
rz(-1.5926188) q[0];
sx q[0];
rz(0.73367992) q[0];
rz(0.94379464) q[1];
sx q[1];
rz(-1.5993886) q[1];
sx q[1];
rz(-2.5234047) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0252903) q[0];
sx q[0];
rz(-2.9553614) q[0];
sx q[0];
rz(1.2982606) q[0];
rz(2.3011977) q[2];
sx q[2];
rz(-1.5777418) q[2];
sx q[2];
rz(2.9027129) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.71431464) q[1];
sx q[1];
rz(-0.66525092) q[1];
sx q[1];
rz(1.6131496) q[1];
x q[2];
rz(2.6691766) q[3];
sx q[3];
rz(-2.4497486) q[3];
sx q[3];
rz(-1.399768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1980629) q[2];
sx q[2];
rz(-1.9518096) q[2];
sx q[2];
rz(2.5243536) q[2];
rz(1.1715568) q[3];
sx q[3];
rz(-2.7660683) q[3];
sx q[3];
rz(2.530781) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.010855762) q[0];
sx q[0];
rz(-0.78748381) q[0];
sx q[0];
rz(0.47295824) q[0];
rz(-1.4594151) q[1];
sx q[1];
rz(-1.727203) q[1];
sx q[1];
rz(0.032729538) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6115295) q[0];
sx q[0];
rz(-0.034688799) q[0];
sx q[0];
rz(0.37337025) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1104537) q[2];
sx q[2];
rz(-0.86058577) q[2];
sx q[2];
rz(-1.6510162) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.53712691) q[1];
sx q[1];
rz(-0.49931881) q[1];
sx q[1];
rz(2.0313655) q[1];
rz(-pi) q[2];
rz(0.030278997) q[3];
sx q[3];
rz(-1.4631728) q[3];
sx q[3];
rz(-0.054747907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.98256436) q[2];
sx q[2];
rz(-0.50243598) q[2];
sx q[2];
rz(-1.095613) q[2];
rz(2.380373) q[3];
sx q[3];
rz(-1.2844362) q[3];
sx q[3];
rz(2.718954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6765167) q[0];
sx q[0];
rz(-0.07971555) q[0];
sx q[0];
rz(-0.7088784) q[0];
rz(2.830016) q[1];
sx q[1];
rz(-2.0496924) q[1];
sx q[1];
rz(0.30837217) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29211048) q[0];
sx q[0];
rz(-1.1774039) q[0];
sx q[0];
rz(1.7500004) q[0];
x q[1];
rz(2.4089085) q[2];
sx q[2];
rz(-0.53712979) q[2];
sx q[2];
rz(-2.2718475) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4235792) q[1];
sx q[1];
rz(-0.96829674) q[1];
sx q[1];
rz(-1.275255) q[1];
rz(-pi) q[2];
rz(0.25777581) q[3];
sx q[3];
rz(-0.80398241) q[3];
sx q[3];
rz(2.3641158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.10820216) q[2];
sx q[2];
rz(-1.6609001) q[2];
sx q[2];
rz(-2.9208276) q[2];
rz(-1.0734142) q[3];
sx q[3];
rz(-2.150841) q[3];
sx q[3];
rz(0.78607917) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6440426) q[0];
sx q[0];
rz(-2.2798517) q[0];
sx q[0];
rz(1.0368689) q[0];
rz(-1.6715624) q[1];
sx q[1];
rz(-2.122888) q[1];
sx q[1];
rz(-1.1933614) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92591015) q[0];
sx q[0];
rz(-1.128845) q[0];
sx q[0];
rz(1.6041167) q[0];
rz(-1.1374129) q[2];
sx q[2];
rz(-2.2221474) q[2];
sx q[2];
rz(-2.1232126) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5037704) q[1];
sx q[1];
rz(-2.42747) q[1];
sx q[1];
rz(-0.28890606) q[1];
rz(-pi) q[2];
rz(-0.12659215) q[3];
sx q[3];
rz(-1.3401573) q[3];
sx q[3];
rz(1.0455719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6679473) q[2];
sx q[2];
rz(-1.8569943) q[2];
sx q[2];
rz(-3.0038707) q[2];
rz(1.600945) q[3];
sx q[3];
rz(-2.9100304) q[3];
sx q[3];
rz(0.93129492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4731428) q[0];
sx q[0];
rz(-1.4746329) q[0];
sx q[0];
rz(2.0869577) q[0];
rz(2.6558381) q[1];
sx q[1];
rz(-1.2243441) q[1];
sx q[1];
rz(-1.4996554) q[1];
rz(-1.3254017) q[2];
sx q[2];
rz(-0.90127887) q[2];
sx q[2];
rz(-0.56294914) q[2];
rz(-2.3790142) q[3];
sx q[3];
rz(-1.1840631) q[3];
sx q[3];
rz(-2.2889591) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
