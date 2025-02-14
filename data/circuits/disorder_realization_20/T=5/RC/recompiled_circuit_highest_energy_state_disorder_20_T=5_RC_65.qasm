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
rz(0.19658495) q[0];
sx q[0];
rz(-2.6994446) q[0];
sx q[0];
rz(-0.86139876) q[0];
rz(1.5098894) q[1];
sx q[1];
rz(-1.8169401) q[1];
sx q[1];
rz(-0.035592508) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7073878) q[0];
sx q[0];
rz(-1.2363096) q[0];
sx q[0];
rz(-1.7777966) q[0];
rz(0.30649779) q[2];
sx q[2];
rz(-1.4671637) q[2];
sx q[2];
rz(2.181884) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5000677) q[1];
sx q[1];
rz(-1.5350122) q[1];
sx q[1];
rz(-0.72662093) q[1];
rz(-pi) q[2];
x q[2];
rz(0.017982131) q[3];
sx q[3];
rz(-2.2507022) q[3];
sx q[3];
rz(0.79538345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5888136) q[2];
sx q[2];
rz(-1.8401044) q[2];
sx q[2];
rz(1.0944875) q[2];
rz(-1.5158481) q[3];
sx q[3];
rz(-0.87259126) q[3];
sx q[3];
rz(0.14357963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.0482408) q[0];
sx q[0];
rz(-0.99026647) q[0];
sx q[0];
rz(2.7675203) q[0];
rz(-0.85809842) q[1];
sx q[1];
rz(-2.2007807) q[1];
sx q[1];
rz(-2.0657952) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8825622) q[0];
sx q[0];
rz(-2.9963065) q[0];
sx q[0];
rz(0.7913211) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7611742) q[2];
sx q[2];
rz(-0.58092344) q[2];
sx q[2];
rz(-0.38512938) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.058671) q[1];
sx q[1];
rz(-1.4545014) q[1];
sx q[1];
rz(-0.40177675) q[1];
x q[2];
rz(-2.4655618) q[3];
sx q[3];
rz(-1.0272678) q[3];
sx q[3];
rz(-1.6401122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51521236) q[2];
sx q[2];
rz(-0.43312803) q[2];
sx q[2];
rz(-2.058378) q[2];
rz(-1.2927239) q[3];
sx q[3];
rz(-1.1336528) q[3];
sx q[3];
rz(2.2129272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.885963) q[0];
sx q[0];
rz(-0.77115458) q[0];
sx q[0];
rz(-0.5249002) q[0];
rz(-1.6820172) q[1];
sx q[1];
rz(-2.4039905) q[1];
sx q[1];
rz(-2.9958013) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6896714) q[0];
sx q[0];
rz(-1.6761685) q[0];
sx q[0];
rz(-1.6622048) q[0];
x q[1];
rz(-0.96876933) q[2];
sx q[2];
rz(-2.0240236) q[2];
sx q[2];
rz(2.9211958) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5151095) q[1];
sx q[1];
rz(-1.4806387) q[1];
sx q[1];
rz(0.43403352) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.57246142) q[3];
sx q[3];
rz(-0.92558555) q[3];
sx q[3];
rz(2.0544586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.24474457) q[2];
sx q[2];
rz(-1.1620099) q[2];
sx q[2];
rz(2.9761918) q[2];
rz(0.8616972) q[3];
sx q[3];
rz(-2.4498037) q[3];
sx q[3];
rz(-2.944788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2794613) q[0];
sx q[0];
rz(-0.52424163) q[0];
sx q[0];
rz(-0.72823802) q[0];
rz(-2.8184452) q[1];
sx q[1];
rz(-2.1073982) q[1];
sx q[1];
rz(2.1023777) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6704039) q[0];
sx q[0];
rz(-1.9946596) q[0];
sx q[0];
rz(2.7039862) q[0];
rz(-2.2065998) q[2];
sx q[2];
rz(-1.3424917) q[2];
sx q[2];
rz(1.993597) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8700712) q[1];
sx q[1];
rz(-1.0876552) q[1];
sx q[1];
rz(-1.6752233) q[1];
rz(2.4432602) q[3];
sx q[3];
rz(-2.192333) q[3];
sx q[3];
rz(2.1752398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.49179658) q[2];
sx q[2];
rz(-1.0504664) q[2];
sx q[2];
rz(-2.6226079) q[2];
rz(1.0745878) q[3];
sx q[3];
rz(-1.2858177) q[3];
sx q[3];
rz(-2.6463267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0786667) q[0];
sx q[0];
rz(-2.2195897) q[0];
sx q[0];
rz(-0.17157383) q[0];
rz(-1.0942787) q[1];
sx q[1];
rz(-1.3010052) q[1];
sx q[1];
rz(-0.23460728) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0594992) q[0];
sx q[0];
rz(-1.6986058) q[0];
sx q[0];
rz(-1.8508085) q[0];
rz(-pi) q[1];
rz(1.6233088) q[2];
sx q[2];
rz(-2.0774475) q[2];
sx q[2];
rz(-1.6135297) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3329016) q[1];
sx q[1];
rz(-2.1655443) q[1];
sx q[1];
rz(-2.4392209) q[1];
rz(2.7127277) q[3];
sx q[3];
rz(-2.3842616) q[3];
sx q[3];
rz(-2.0067818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91308633) q[2];
sx q[2];
rz(-1.191491) q[2];
sx q[2];
rz(1.4027493) q[2];
rz(-2.5168822) q[3];
sx q[3];
rz(-2.1731845) q[3];
sx q[3];
rz(-1.3431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.6828571) q[0];
sx q[0];
rz(-0.0077489297) q[0];
sx q[0];
rz(-2.8140581) q[0];
rz(-3.1134743) q[1];
sx q[1];
rz(-1.997812) q[1];
sx q[1];
rz(0.99172529) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.018143749) q[0];
sx q[0];
rz(-0.94052343) q[0];
sx q[0];
rz(2.2176377) q[0];
rz(1.6794793) q[2];
sx q[2];
rz(-1.3460438) q[2];
sx q[2];
rz(-2.0456184) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2497445) q[1];
sx q[1];
rz(-1.5288903) q[1];
sx q[1];
rz(0.81467198) q[1];
x q[2];
rz(1.5434274) q[3];
sx q[3];
rz(-2.0413766) q[3];
sx q[3];
rz(-1.9675627) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.907054) q[2];
sx q[2];
rz(-1.7639284) q[2];
sx q[2];
rz(-1.482359) q[2];
rz(2.0756857) q[3];
sx q[3];
rz(-1.9612471) q[3];
sx q[3];
rz(2.0445686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7557573) q[0];
sx q[0];
rz(-0.11180728) q[0];
sx q[0];
rz(-2.8461611) q[0];
rz(0.23751986) q[1];
sx q[1];
rz(-2.2078881) q[1];
sx q[1];
rz(-2.3942153) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6791466) q[0];
sx q[0];
rz(-2.0751795) q[0];
sx q[0];
rz(0.64179582) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0355936) q[2];
sx q[2];
rz(-1.0135302) q[2];
sx q[2];
rz(1.1320499) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.805757) q[1];
sx q[1];
rz(-2.0621967) q[1];
sx q[1];
rz(-2.7193428) q[1];
x q[2];
rz(-0.57716863) q[3];
sx q[3];
rz(-2.5169229) q[3];
sx q[3];
rz(-0.28055252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6682917) q[2];
sx q[2];
rz(-1.1857727) q[2];
sx q[2];
rz(0.2549003) q[2];
rz(2.9292246) q[3];
sx q[3];
rz(-2.2771211) q[3];
sx q[3];
rz(-0.00096850639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85106987) q[0];
sx q[0];
rz(-1.9177328) q[0];
sx q[0];
rz(-0.12390027) q[0];
rz(2.2663785) q[1];
sx q[1];
rz(-0.83299914) q[1];
sx q[1];
rz(-2.5828054) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66984719) q[0];
sx q[0];
rz(-2.5515243) q[0];
sx q[0];
rz(2.0315995) q[0];
rz(-pi) q[1];
rz(1.0718078) q[2];
sx q[2];
rz(-1.3200511) q[2];
sx q[2];
rz(-2.8607228) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9001101) q[1];
sx q[1];
rz(-2.0962976) q[1];
sx q[1];
rz(1.8837758) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.41591538) q[3];
sx q[3];
rz(-2.1634566) q[3];
sx q[3];
rz(-0.76502467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.045804068) q[2];
sx q[2];
rz(-2.6498821) q[2];
sx q[2];
rz(2.4208505) q[2];
rz(-2.7785684) q[3];
sx q[3];
rz(-1.9178773) q[3];
sx q[3];
rz(1.629841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95661288) q[0];
sx q[0];
rz(-0.19421254) q[0];
sx q[0];
rz(3.0525364) q[0];
rz(0.36674276) q[1];
sx q[1];
rz(-2.2373503) q[1];
sx q[1];
rz(-1.4595703) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37292591) q[0];
sx q[0];
rz(-0.67250508) q[0];
sx q[0];
rz(-1.9422533) q[0];
rz(-pi) q[1];
rz(-2.1543571) q[2];
sx q[2];
rz(-2.1676807) q[2];
sx q[2];
rz(-1.0839628) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.68946099) q[1];
sx q[1];
rz(-0.88400999) q[1];
sx q[1];
rz(2.6584714) q[1];
rz(1.1330539) q[3];
sx q[3];
rz(-2.1244085) q[3];
sx q[3];
rz(-2.5598524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2064994) q[2];
sx q[2];
rz(-1.7743899) q[2];
sx q[2];
rz(1.8915141) q[2];
rz(1.9153204) q[3];
sx q[3];
rz(-2.5076702) q[3];
sx q[3];
rz(2.5026076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6298237) q[0];
sx q[0];
rz(-1.1331695) q[0];
sx q[0];
rz(1.1400219) q[0];
rz(0.58311588) q[1];
sx q[1];
rz(-1.4864328) q[1];
sx q[1];
rz(-2.4737632) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17530046) q[0];
sx q[0];
rz(-2.326823) q[0];
sx q[0];
rz(0.45481429) q[0];
x q[1];
rz(2.2726353) q[2];
sx q[2];
rz(-1.7545106) q[2];
sx q[2];
rz(0.76317235) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7857504) q[1];
sx q[1];
rz(-1.3763345) q[1];
sx q[1];
rz(2.0114312) q[1];
rz(1.1125426) q[3];
sx q[3];
rz(-2.309121) q[3];
sx q[3];
rz(0.74362459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84393152) q[2];
sx q[2];
rz(-0.49489489) q[2];
sx q[2];
rz(-0.75221357) q[2];
rz(0.080605896) q[3];
sx q[3];
rz(-1.7919431) q[3];
sx q[3];
rz(-0.54752553) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8053631) q[0];
sx q[0];
rz(-1.6148051) q[0];
sx q[0];
rz(-0.057407277) q[0];
rz(0.84857955) q[1];
sx q[1];
rz(-0.72863693) q[1];
sx q[1];
rz(-1.4562664) q[1];
rz(0.66203881) q[2];
sx q[2];
rz(-1.0218191) q[2];
sx q[2];
rz(-1.54984) q[2];
rz(0.46075321) q[3];
sx q[3];
rz(-1.8965773) q[3];
sx q[3];
rz(-2.0544485) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
