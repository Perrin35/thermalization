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
rz(2.6235629) q[0];
sx q[0];
rz(-1.1413483) q[0];
sx q[0];
rz(-3.1014693) q[0];
rz(0.24815458) q[1];
sx q[1];
rz(5.2207898) q[1];
sx q[1];
rz(9.3407486) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3798767) q[0];
sx q[0];
rz(-1.5119702) q[0];
sx q[0];
rz(0.071278871) q[0];
rz(-pi) q[1];
rz(-2.7366287) q[2];
sx q[2];
rz(-2.7128007) q[2];
sx q[2];
rz(2.2909209) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1511615) q[1];
sx q[1];
rz(-0.50499454) q[1];
sx q[1];
rz(-2.7846365) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1575451) q[3];
sx q[3];
rz(-1.0344369) q[3];
sx q[3];
rz(2.9265917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4485126) q[2];
sx q[2];
rz(-1.7757519) q[2];
sx q[2];
rz(2.8903294) q[2];
rz(-1.1723899) q[3];
sx q[3];
rz(-2.0386212) q[3];
sx q[3];
rz(-3.0917061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4548816) q[0];
sx q[0];
rz(-2.659681) q[0];
sx q[0];
rz(1.3739817) q[0];
rz(2.8107367) q[1];
sx q[1];
rz(-1.4127981) q[1];
sx q[1];
rz(-0.27423283) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62456709) q[0];
sx q[0];
rz(-0.32450482) q[0];
sx q[0];
rz(2.9175678) q[0];
rz(-pi) q[1];
x q[1];
rz(3.129446) q[2];
sx q[2];
rz(-0.53330219) q[2];
sx q[2];
rz(0.27517327) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3220904) q[1];
sx q[1];
rz(-1.354447) q[1];
sx q[1];
rz(0.33239969) q[1];
rz(-pi) q[2];
rz(3.0586309) q[3];
sx q[3];
rz(-1.9881691) q[3];
sx q[3];
rz(-0.42835945) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8809044) q[2];
sx q[2];
rz(-1.3284677) q[2];
sx q[2];
rz(-1.0880967) q[2];
rz(-2.0982096) q[3];
sx q[3];
rz(-2.9731396) q[3];
sx q[3];
rz(2.8823631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.008721) q[0];
sx q[0];
rz(-0.50478029) q[0];
sx q[0];
rz(-2.0399427) q[0];
rz(2.3654826) q[1];
sx q[1];
rz(-1.848369) q[1];
sx q[1];
rz(0.64220846) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4896079) q[0];
sx q[0];
rz(-0.63417681) q[0];
sx q[0];
rz(-0.80192566) q[0];
x q[1];
rz(2.8611254) q[2];
sx q[2];
rz(-2.4020122) q[2];
sx q[2];
rz(-2.6698339) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.94564181) q[1];
sx q[1];
rz(-1.6744876) q[1];
sx q[1];
rz(1.1432462) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3989556) q[3];
sx q[3];
rz(-0.81002319) q[3];
sx q[3];
rz(-2.1049316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8316101) q[2];
sx q[2];
rz(-1.9881366) q[2];
sx q[2];
rz(2.2115121) q[2];
rz(0.81405226) q[3];
sx q[3];
rz(-1.1907153) q[3];
sx q[3];
rz(-1.2306151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.880068) q[0];
sx q[0];
rz(-1.9517169) q[0];
sx q[0];
rz(0.50450182) q[0];
rz(-0.91744676) q[1];
sx q[1];
rz(-2.4023299) q[1];
sx q[1];
rz(2.9333072) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84346164) q[0];
sx q[0];
rz(-1.4220474) q[0];
sx q[0];
rz(0.24555969) q[0];
rz(-pi) q[1];
x q[1];
rz(0.31625749) q[2];
sx q[2];
rz(-0.74511792) q[2];
sx q[2];
rz(-0.032501246) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.444297) q[1];
sx q[1];
rz(-2.136886) q[1];
sx q[1];
rz(-2.5725911) q[1];
x q[2];
rz(-0.39842968) q[3];
sx q[3];
rz(-2.5442985) q[3];
sx q[3];
rz(2.3526827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.31671277) q[2];
sx q[2];
rz(-2.0471639) q[2];
sx q[2];
rz(-1.971395) q[2];
rz(-2.1757388) q[3];
sx q[3];
rz(-1.0713157) q[3];
sx q[3];
rz(0.14537183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.589094) q[0];
sx q[0];
rz(-0.67411244) q[0];
sx q[0];
rz(-2.476995) q[0];
rz(2.873114) q[1];
sx q[1];
rz(-1.6842027) q[1];
sx q[1];
rz(0.99884117) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.474668) q[0];
sx q[0];
rz(-2.5380236) q[0];
sx q[0];
rz(-1.3003028) q[0];
rz(-pi) q[1];
rz(1.2810318) q[2];
sx q[2];
rz(-0.43998566) q[2];
sx q[2];
rz(-0.54674613) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4376862) q[1];
sx q[1];
rz(-1.4515299) q[1];
sx q[1];
rz(1.7410623) q[1];
x q[2];
rz(1.5454917) q[3];
sx q[3];
rz(-0.69195834) q[3];
sx q[3];
rz(1.4922752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9773679) q[2];
sx q[2];
rz(-2.2013142) q[2];
sx q[2];
rz(-0.11717907) q[2];
rz(-3.0818648) q[3];
sx q[3];
rz(-1.2195339) q[3];
sx q[3];
rz(-2.3265649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3672459) q[0];
sx q[0];
rz(-2.7278439) q[0];
sx q[0];
rz(1.8814948) q[0];
rz(0.11481181) q[1];
sx q[1];
rz(-1.3453307) q[1];
sx q[1];
rz(0.095349163) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0336825) q[0];
sx q[0];
rz(-2.9386407) q[0];
sx q[0];
rz(-0.072921948) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2749407) q[2];
sx q[2];
rz(-2.618578) q[2];
sx q[2];
rz(-2.366334) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1667249) q[1];
sx q[1];
rz(-2.0151867) q[1];
sx q[1];
rz(-0.4687029) q[1];
rz(0.010706832) q[3];
sx q[3];
rz(-1.043415) q[3];
sx q[3];
rz(-2.7675865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0872588) q[2];
sx q[2];
rz(-2.2237873) q[2];
sx q[2];
rz(-0.66162649) q[2];
rz(0.76964393) q[3];
sx q[3];
rz(-1.6911643) q[3];
sx q[3];
rz(-0.86696398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0041644) q[0];
sx q[0];
rz(-0.34808174) q[0];
sx q[0];
rz(0.75781703) q[0];
rz(3.1163395) q[1];
sx q[1];
rz(-2.7460637) q[1];
sx q[1];
rz(1.1753561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34285082) q[0];
sx q[0];
rz(-1.9254058) q[0];
sx q[0];
rz(-2.6358428) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.25039) q[2];
sx q[2];
rz(-1.9632578) q[2];
sx q[2];
rz(-3.0796555) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0846128) q[1];
sx q[1];
rz(-0.86922042) q[1];
sx q[1];
rz(-2.6210611) q[1];
rz(-pi) q[2];
rz(1.5478915) q[3];
sx q[3];
rz(-1.9996907) q[3];
sx q[3];
rz(1.0064841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5777099) q[2];
sx q[2];
rz(-1.5241728) q[2];
sx q[2];
rz(1.2152524) q[2];
rz(0.75872129) q[3];
sx q[3];
rz(-1.4164378) q[3];
sx q[3];
rz(1.5846213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6587081) q[0];
sx q[0];
rz(-1.8526798) q[0];
sx q[0];
rz(2.7523852) q[0];
rz(-3.0534741) q[1];
sx q[1];
rz(-1.4382818) q[1];
sx q[1];
rz(1.5315936) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0813839) q[0];
sx q[0];
rz(-1.7870149) q[0];
sx q[0];
rz(0.24198089) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5814452) q[2];
sx q[2];
rz(-2.3986926) q[2];
sx q[2];
rz(-2.7882238) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4952505) q[1];
sx q[1];
rz(-1.6227437) q[1];
sx q[1];
rz(2.6990128) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1346899) q[3];
sx q[3];
rz(-0.75911555) q[3];
sx q[3];
rz(-1.5838983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7342928) q[2];
sx q[2];
rz(-0.6137085) q[2];
sx q[2];
rz(-1.3235922) q[2];
rz(-1.343441) q[3];
sx q[3];
rz(-2.3517793) q[3];
sx q[3];
rz(0.33671236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4765428) q[0];
sx q[0];
rz(-0.90710586) q[0];
sx q[0];
rz(0.46515775) q[0];
rz(0.05050412) q[1];
sx q[1];
rz(-0.89021325) q[1];
sx q[1];
rz(-0.49096289) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4024593) q[0];
sx q[0];
rz(-1.8998892) q[0];
sx q[0];
rz(1.2335445) q[0];
x q[1];
rz(-0.30740909) q[2];
sx q[2];
rz(-1.7176065) q[2];
sx q[2];
rz(1.4585782) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.95066324) q[1];
sx q[1];
rz(-2.1219538) q[1];
sx q[1];
rz(-2.3930156) q[1];
rz(-pi) q[2];
rz(-1.4598941) q[3];
sx q[3];
rz(-1.8506128) q[3];
sx q[3];
rz(-0.96074897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.65840536) q[2];
sx q[2];
rz(-2.5453973) q[2];
sx q[2];
rz(2.2779706) q[2];
rz(-2.9577799) q[3];
sx q[3];
rz(-2.176216) q[3];
sx q[3];
rz(-2.1553154) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42613906) q[0];
sx q[0];
rz(-1.621839) q[0];
sx q[0];
rz(-2.5157628) q[0];
rz(-2.4848056) q[1];
sx q[1];
rz(-2.1757809) q[1];
sx q[1];
rz(1.2084557) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55126429) q[0];
sx q[0];
rz(-2.9737824) q[0];
sx q[0];
rz(2.2237334) q[0];
x q[1];
rz(0.65932806) q[2];
sx q[2];
rz(-1.8272597) q[2];
sx q[2];
rz(1.9224482) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.085026) q[1];
sx q[1];
rz(-1.4843656) q[1];
sx q[1];
rz(1.5734476) q[1];
x q[2];
rz(1.3937065) q[3];
sx q[3];
rz(-1.0064126) q[3];
sx q[3];
rz(-2.0616814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5586231) q[2];
sx q[2];
rz(-1.8331336) q[2];
sx q[2];
rz(2.6226131) q[2];
rz(2.6988622) q[3];
sx q[3];
rz(-1.3230007) q[3];
sx q[3];
rz(1.3435266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9291572) q[0];
sx q[0];
rz(-2.1514308) q[0];
sx q[0];
rz(1.9707752) q[0];
rz(-2.9875372) q[1];
sx q[1];
rz(-2.2046721) q[1];
sx q[1];
rz(2.2278723) q[1];
rz(0.51657233) q[2];
sx q[2];
rz(-1.6067098) q[2];
sx q[2];
rz(-1.4756257) q[2];
rz(1.1476573) q[3];
sx q[3];
rz(-0.71865766) q[3];
sx q[3];
rz(-1.3620993) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
