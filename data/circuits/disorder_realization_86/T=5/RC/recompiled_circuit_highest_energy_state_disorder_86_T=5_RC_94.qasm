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
rz(1.6744094) q[0];
sx q[0];
rz(-1.5404584) q[0];
sx q[0];
rz(1.2920446) q[0];
rz(2.0436824) q[1];
sx q[1];
rz(-0.75610375) q[1];
sx q[1];
rz(2.4296711) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24589892) q[0];
sx q[0];
rz(-2.993949) q[0];
sx q[0];
rz(2.9013322) q[0];
rz(2.8991382) q[2];
sx q[2];
rz(-2.3676845) q[2];
sx q[2];
rz(-0.36786181) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.10140534) q[1];
sx q[1];
rz(-2.1899472) q[1];
sx q[1];
rz(-1.2907486) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96613066) q[3];
sx q[3];
rz(-1.1388766) q[3];
sx q[3];
rz(-2.8063065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0277675) q[2];
sx q[2];
rz(-1.0372459) q[2];
sx q[2];
rz(-2.2859331) q[2];
rz(-1.8850108) q[3];
sx q[3];
rz(-2.515007) q[3];
sx q[3];
rz(1.1945061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7282495) q[0];
sx q[0];
rz(-3.0276868) q[0];
sx q[0];
rz(-0.86694992) q[0];
rz(1.4504245) q[1];
sx q[1];
rz(-1.9536641) q[1];
sx q[1];
rz(-0.11122045) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7635083) q[0];
sx q[0];
rz(-1.0279546) q[0];
sx q[0];
rz(1.1397902) q[0];
rz(-pi) q[1];
rz(-0.84663518) q[2];
sx q[2];
rz(-2.4754582) q[2];
sx q[2];
rz(-0.29633477) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4091089) q[1];
sx q[1];
rz(-0.95327158) q[1];
sx q[1];
rz(0.78583048) q[1];
x q[2];
rz(2.9807509) q[3];
sx q[3];
rz(-2.2713158) q[3];
sx q[3];
rz(-0.15047401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1078681) q[2];
sx q[2];
rz(-2.3794231) q[2];
sx q[2];
rz(-0.49125683) q[2];
rz(-1.2563502) q[3];
sx q[3];
rz(-1.7137073) q[3];
sx q[3];
rz(-0.45891941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97378174) q[0];
sx q[0];
rz(-1.351492) q[0];
sx q[0];
rz(-0.0034045086) q[0];
rz(-0.34524521) q[1];
sx q[1];
rz(-0.41989851) q[1];
sx q[1];
rz(-0.87510625) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7686979) q[0];
sx q[0];
rz(-2.5803496) q[0];
sx q[0];
rz(-3.0054557) q[0];
rz(-pi) q[1];
rz(0.11981182) q[2];
sx q[2];
rz(-1.3917599) q[2];
sx q[2];
rz(1.8869635) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.1701974) q[1];
sx q[1];
rz(-1.9985693) q[1];
sx q[1];
rz(2.0712159) q[1];
rz(-3.1140675) q[3];
sx q[3];
rz(-0.77139957) q[3];
sx q[3];
rz(1.9302544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.379999) q[2];
sx q[2];
rz(-0.21445175) q[2];
sx q[2];
rz(2.8144515) q[2];
rz(1.7548148) q[3];
sx q[3];
rz(-1.1969457) q[3];
sx q[3];
rz(-3.0676214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73579329) q[0];
sx q[0];
rz(-1.0013591) q[0];
sx q[0];
rz(-1.5843947) q[0];
rz(0.34218511) q[1];
sx q[1];
rz(-2.1125968) q[1];
sx q[1];
rz(-2.7142966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7988069) q[0];
sx q[0];
rz(-1.5166993) q[0];
sx q[0];
rz(1.6038206) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1085243) q[2];
sx q[2];
rz(-1.8934774) q[2];
sx q[2];
rz(3.0442381) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8646331) q[1];
sx q[1];
rz(-1.0103419) q[1];
sx q[1];
rz(0.87866453) q[1];
rz(-1.8713636) q[3];
sx q[3];
rz(-1.0247314) q[3];
sx q[3];
rz(-0.48496839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7067318) q[2];
sx q[2];
rz(-0.92919436) q[2];
sx q[2];
rz(2.2959607) q[2];
rz(3.1137858) q[3];
sx q[3];
rz(-1.7866725) q[3];
sx q[3];
rz(1.4771247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-2.528275) q[0];
sx q[0];
rz(-2.1149825) q[0];
sx q[0];
rz(-3.0897621) q[0];
rz(-1.9822281) q[1];
sx q[1];
rz(-2.4720188) q[1];
sx q[1];
rz(-0.61028284) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6298415) q[0];
sx q[0];
rz(-1.7423706) q[0];
sx q[0];
rz(0.20769329) q[0];
rz(-pi) q[1];
rz(-1.9699783) q[2];
sx q[2];
rz(-1.295168) q[2];
sx q[2];
rz(2.3720019) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.8598664) q[1];
sx q[1];
rz(-1.7471581) q[1];
sx q[1];
rz(2.6540666) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5334643) q[3];
sx q[3];
rz(-1.593489) q[3];
sx q[3];
rz(-0.80858999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74756885) q[2];
sx q[2];
rz(-1.3538066) q[2];
sx q[2];
rz(0.58689153) q[2];
rz(1.759257) q[3];
sx q[3];
rz(-2.0943677) q[3];
sx q[3];
rz(-2.1658678) q[3];
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
rz(-pi/2) q[0];
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
rz(-1.0368283) q[0];
sx q[0];
rz(-2.6703175) q[0];
sx q[0];
rz(-2.9608744) q[0];
rz(-1.7432927) q[1];
sx q[1];
rz(-1.9075874) q[1];
sx q[1];
rz(-1.3999375) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0796868) q[0];
sx q[0];
rz(-2.3374515) q[0];
sx q[0];
rz(-1.8820394) q[0];
rz(-pi) q[1];
rz(2.7488974) q[2];
sx q[2];
rz(-1.1989856) q[2];
sx q[2];
rz(1.8657045) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.1776947) q[1];
sx q[1];
rz(-2.438058) q[1];
sx q[1];
rz(-1.7046609) q[1];
x q[2];
rz(0.67541452) q[3];
sx q[3];
rz(-2.3590909) q[3];
sx q[3];
rz(-2.7030962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.7515298) q[2];
sx q[2];
rz(-0.84948245) q[2];
sx q[2];
rz(-2.0709822) q[2];
rz(-1.4812329) q[3];
sx q[3];
rz(-2.345572) q[3];
sx q[3];
rz(-1.3720366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91167951) q[0];
sx q[0];
rz(-2.3340618) q[0];
sx q[0];
rz(-2.5813778) q[0];
rz(2.920976) q[1];
sx q[1];
rz(-0.95181528) q[1];
sx q[1];
rz(-0.87699786) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79898873) q[0];
sx q[0];
rz(-1.5681802) q[0];
sx q[0];
rz(1.553234) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7438597) q[2];
sx q[2];
rz(-1.8200201) q[2];
sx q[2];
rz(-0.13532369) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6798579) q[1];
sx q[1];
rz(-1.468424) q[1];
sx q[1];
rz(1.8385248) q[1];
rz(-pi) q[2];
rz(-3.0362066) q[3];
sx q[3];
rz(-2.2544754) q[3];
sx q[3];
rz(-0.87382853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.47473869) q[2];
sx q[2];
rz(-2.2467504) q[2];
sx q[2];
rz(-1.3372927) q[2];
rz(-0.21833359) q[3];
sx q[3];
rz(-2.1209013) q[3];
sx q[3];
rz(1.4214424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5050932) q[0];
sx q[0];
rz(-1.8890843) q[0];
sx q[0];
rz(-3.002758) q[0];
rz(0.88169634) q[1];
sx q[1];
rz(-1.9813184) q[1];
sx q[1];
rz(0.82463157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76060847) q[0];
sx q[0];
rz(-0.09481096) q[0];
sx q[0];
rz(1.3349956) q[0];
rz(3.0443848) q[2];
sx q[2];
rz(-1.6748472) q[2];
sx q[2];
rz(-1.0417633) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.127434) q[1];
sx q[1];
rz(-2.6175272) q[1];
sx q[1];
rz(0.31619831) q[1];
rz(-pi) q[2];
rz(0.91757284) q[3];
sx q[3];
rz(-1.6503449) q[3];
sx q[3];
rz(-3.0419526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3489939) q[2];
sx q[2];
rz(-0.94453347) q[2];
sx q[2];
rz(0.83354956) q[2];
rz(2.870657) q[3];
sx q[3];
rz(-2.0201594) q[3];
sx q[3];
rz(-0.84735316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18868294) q[0];
sx q[0];
rz(-1.2485349) q[0];
sx q[0];
rz(2.6771925) q[0];
rz(2.4480827) q[1];
sx q[1];
rz(-2.214407) q[1];
sx q[1];
rz(-2.4448591) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.369608) q[0];
sx q[0];
rz(-1.4596997) q[0];
sx q[0];
rz(0.055776061) q[0];
x q[1];
rz(1.3905754) q[2];
sx q[2];
rz(-1.5308793) q[2];
sx q[2];
rz(-0.60933622) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3426376) q[1];
sx q[1];
rz(-0.25051719) q[1];
sx q[1];
rz(1.6813075) q[1];
x q[2];
rz(1.5570627) q[3];
sx q[3];
rz(-0.93051118) q[3];
sx q[3];
rz(0.69956799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.36737475) q[2];
sx q[2];
rz(-0.092270277) q[2];
sx q[2];
rz(-0.83887678) q[2];
rz(1.4539666) q[3];
sx q[3];
rz(-1.5980915) q[3];
sx q[3];
rz(1.4723697) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0812747) q[0];
sx q[0];
rz(-1.769279) q[0];
sx q[0];
rz(-1.6171612) q[0];
rz(0.35628191) q[1];
sx q[1];
rz(-1.7226912) q[1];
sx q[1];
rz(1.3391395) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55521783) q[0];
sx q[0];
rz(-1.5065743) q[0];
sx q[0];
rz(-0.97126605) q[0];
x q[1];
rz(-0.3199749) q[2];
sx q[2];
rz(-1.5029482) q[2];
sx q[2];
rz(-2.0463129) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4831743) q[1];
sx q[1];
rz(-0.69026679) q[1];
sx q[1];
rz(-0.19259318) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8175199) q[3];
sx q[3];
rz(-2.9949247) q[3];
sx q[3];
rz(0.038692729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.39849207) q[2];
sx q[2];
rz(-2.3091381) q[2];
sx q[2];
rz(-1.8771646) q[2];
rz(2.8049331) q[3];
sx q[3];
rz(-2.2008937) q[3];
sx q[3];
rz(-1.3338859) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7548512) q[0];
sx q[0];
rz(-1.6254397) q[0];
sx q[0];
rz(2.0425015) q[0];
rz(0.59973888) q[1];
sx q[1];
rz(-0.44282423) q[1];
sx q[1];
rz(2.5437358) q[1];
rz(-1.1588396) q[2];
sx q[2];
rz(-1.2745274) q[2];
sx q[2];
rz(1.6714255) q[2];
rz(2.609008) q[3];
sx q[3];
rz(-1.5556794) q[3];
sx q[3];
rz(-1.5928768) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
