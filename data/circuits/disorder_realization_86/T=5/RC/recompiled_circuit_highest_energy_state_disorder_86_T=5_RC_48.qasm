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
rz(-1.4671833) q[0];
sx q[0];
rz(-1.6011342) q[0];
sx q[0];
rz(1.849548) q[0];
rz(-1.0979103) q[1];
sx q[1];
rz(-2.3854889) q[1];
sx q[1];
rz(0.71192157) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24589892) q[0];
sx q[0];
rz(-0.14764365) q[0];
sx q[0];
rz(-2.9013322) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3825157) q[2];
sx q[2];
rz(-1.7393987) q[2];
sx q[2];
rz(2.1137266) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56150964) q[1];
sx q[1];
rz(-0.67189184) q[1];
sx q[1];
rz(-0.36999934) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.175462) q[3];
sx q[3];
rz(-1.1388766) q[3];
sx q[3];
rz(2.8063065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0277675) q[2];
sx q[2];
rz(-1.0372459) q[2];
sx q[2];
rz(-2.2859331) q[2];
rz(1.8850108) q[3];
sx q[3];
rz(-2.515007) q[3];
sx q[3];
rz(-1.1945061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.41334316) q[0];
sx q[0];
rz(-0.11390587) q[0];
sx q[0];
rz(2.2746427) q[0];
rz(-1.4504245) q[1];
sx q[1];
rz(-1.1879286) q[1];
sx q[1];
rz(-0.11122045) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3780843) q[0];
sx q[0];
rz(-2.1136381) q[0];
sx q[0];
rz(-1.1397902) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2949575) q[2];
sx q[2];
rz(-2.4754582) q[2];
sx q[2];
rz(2.8452579) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.31383816) q[1];
sx q[1];
rz(-0.95664531) q[1];
sx q[1];
rz(-0.78734963) q[1];
rz(-pi) q[2];
rz(2.9807509) q[3];
sx q[3];
rz(-2.2713158) q[3];
sx q[3];
rz(2.9911186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0337246) q[2];
sx q[2];
rz(-0.76216951) q[2];
sx q[2];
rz(2.6503358) q[2];
rz(1.8852425) q[3];
sx q[3];
rz(-1.4278853) q[3];
sx q[3];
rz(-2.6826732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1678109) q[0];
sx q[0];
rz(-1.351492) q[0];
sx q[0];
rz(3.1381881) q[0];
rz(-2.7963474) q[1];
sx q[1];
rz(-2.7216941) q[1];
sx q[1];
rz(2.2664864) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3728947) q[0];
sx q[0];
rz(-2.5803496) q[0];
sx q[0];
rz(0.13613693) q[0];
rz(-pi) q[1];
rz(-2.1544863) q[2];
sx q[2];
rz(-0.21506992) q[2];
sx q[2];
rz(-1.849086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.1701974) q[1];
sx q[1];
rz(-1.1430234) q[1];
sx q[1];
rz(-2.0712159) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77121021) q[3];
sx q[3];
rz(-1.589984) q[3];
sx q[3];
rz(2.8018708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.7615937) q[2];
sx q[2];
rz(-2.9271409) q[2];
sx q[2];
rz(0.3271412) q[2];
rz(-1.7548148) q[3];
sx q[3];
rz(-1.1969457) q[3];
sx q[3];
rz(3.0676214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73579329) q[0];
sx q[0];
rz(-2.1402335) q[0];
sx q[0];
rz(1.5843947) q[0];
rz(2.7994075) q[1];
sx q[1];
rz(-2.1125968) q[1];
sx q[1];
rz(2.7142966) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9153684) q[0];
sx q[0];
rz(-1.6037723) q[0];
sx q[0];
rz(-3.0874662) q[0];
rz(-2.1491545) q[2];
sx q[2];
rz(-0.61884862) q[2];
sx q[2];
rz(-2.1568795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4179743) q[1];
sx q[1];
rz(-2.2811541) q[1];
sx q[1];
rz(-0.79377004) q[1];
rz(1.8713636) q[3];
sx q[3];
rz(-2.1168613) q[3];
sx q[3];
rz(-0.48496839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4348609) q[2];
sx q[2];
rz(-2.2123983) q[2];
sx q[2];
rz(-2.2959607) q[2];
rz(3.1137858) q[3];
sx q[3];
rz(-1.3549201) q[3];
sx q[3];
rz(1.6644679) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.528275) q[0];
sx q[0];
rz(-1.0266101) q[0];
sx q[0];
rz(0.051830526) q[0];
rz(-1.1593646) q[1];
sx q[1];
rz(-2.4720188) q[1];
sx q[1];
rz(0.61028284) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1185111) q[0];
sx q[0];
rz(-1.7753965) q[0];
sx q[0];
rz(-1.3955297) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29783037) q[2];
sx q[2];
rz(-1.9541085) q[2];
sx q[2];
rz(-0.9155067) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1111022) q[1];
sx q[1];
rz(-0.51601948) q[1];
sx q[1];
rz(-0.36353021) q[1];
rz(1.5984437) q[3];
sx q[3];
rz(-2.1787454) q[3];
sx q[3];
rz(-0.77800084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.74756885) q[2];
sx q[2];
rz(-1.3538066) q[2];
sx q[2];
rz(-0.58689153) q[2];
rz(1.3823357) q[3];
sx q[3];
rz(-2.0943677) q[3];
sx q[3];
rz(-0.97572485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0368283) q[0];
sx q[0];
rz(-2.6703175) q[0];
sx q[0];
rz(2.9608744) q[0];
rz(-1.7432927) q[1];
sx q[1];
rz(-1.9075874) q[1];
sx q[1];
rz(-1.3999375) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0619059) q[0];
sx q[0];
rz(-0.8041412) q[0];
sx q[0];
rz(-1.8820394) q[0];
rz(2.7488974) q[2];
sx q[2];
rz(-1.1989856) q[2];
sx q[2];
rz(1.8657045) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8508262) q[1];
sx q[1];
rz(-1.6572448) q[1];
sx q[1];
rz(-2.269901) q[1];
x q[2];
rz(-2.1269538) q[3];
sx q[3];
rz(-0.98812166) q[3];
sx q[3];
rz(-1.2847054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7515298) q[2];
sx q[2];
rz(-2.2921102) q[2];
sx q[2];
rz(-2.0709822) q[2];
rz(-1.6603598) q[3];
sx q[3];
rz(-2.345572) q[3];
sx q[3];
rz(-1.7695561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91167951) q[0];
sx q[0];
rz(-2.3340618) q[0];
sx q[0];
rz(0.56021488) q[0];
rz(-0.22061661) q[1];
sx q[1];
rz(-0.95181528) q[1];
sx q[1];
rz(-0.87699786) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3426039) q[0];
sx q[0];
rz(-1.5681802) q[0];
sx q[0];
rz(1.553234) q[0];
rz(2.7438597) q[2];
sx q[2];
rz(-1.3215725) q[2];
sx q[2];
rz(-3.006269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6798579) q[1];
sx q[1];
rz(-1.468424) q[1];
sx q[1];
rz(-1.3030679) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6991922) q[3];
sx q[3];
rz(-0.69046016) q[3];
sx q[3];
rz(-0.70790509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.666854) q[2];
sx q[2];
rz(-2.2467504) q[2];
sx q[2];
rz(1.8043) q[2];
rz(-0.21833359) q[3];
sx q[3];
rz(-1.0206914) q[3];
sx q[3];
rz(-1.4214424) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6364994) q[0];
sx q[0];
rz(-1.2525083) q[0];
sx q[0];
rz(-3.002758) q[0];
rz(-2.2598963) q[1];
sx q[1];
rz(-1.9813184) q[1];
sx q[1];
rz(-2.3169611) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0449681) q[0];
sx q[0];
rz(-1.5929149) q[0];
sx q[0];
rz(1.6629987) q[0];
x q[1];
rz(0.82197676) q[2];
sx q[2];
rz(-0.14227371) q[2];
sx q[2];
rz(1.345944) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.6528553) q[1];
sx q[1];
rz(-1.0751546) q[1];
sx q[1];
rz(1.3929709) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-1.3489939) q[2];
sx q[2];
rz(-0.94453347) q[2];
sx q[2];
rz(-0.83354956) q[2];
rz(2.870657) q[3];
sx q[3];
rz(-1.1214333) q[3];
sx q[3];
rz(0.84735316) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9529097) q[0];
sx q[0];
rz(-1.2485349) q[0];
sx q[0];
rz(2.6771925) q[0];
rz(-0.69350997) q[1];
sx q[1];
rz(-2.214407) q[1];
sx q[1];
rz(0.69673353) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90307921) q[0];
sx q[0];
rz(-3.0173324) q[0];
sx q[0];
rz(2.0342229) q[0];
rz(1.3905754) q[2];
sx q[2];
rz(-1.6107133) q[2];
sx q[2];
rz(-2.5322564) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9129968) q[1];
sx q[1];
rz(-1.8197528) q[1];
sx q[1];
rz(0.028214022) q[1];
rz(2.5012623) q[3];
sx q[3];
rz(-1.5818095) q[3];
sx q[3];
rz(2.2785694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.7742179) q[2];
sx q[2];
rz(-3.0493224) q[2];
sx q[2];
rz(-2.3027159) q[2];
rz(-1.687626) q[3];
sx q[3];
rz(-1.5980915) q[3];
sx q[3];
rz(-1.669223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0812747) q[0];
sx q[0];
rz(-1.769279) q[0];
sx q[0];
rz(1.5244315) q[0];
rz(-2.7853107) q[1];
sx q[1];
rz(-1.4189015) q[1];
sx q[1];
rz(1.8024532) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55521783) q[0];
sx q[0];
rz(-1.5065743) q[0];
sx q[0];
rz(-2.1703266) q[0];
rz(1.4993323) q[2];
sx q[2];
rz(-1.8900089) q[2];
sx q[2];
rz(-0.45305529) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.061627541) q[1];
sx q[1];
rz(-1.6929757) q[1];
sx q[1];
rz(-2.4604718) q[1];
rz(-pi) q[2];
rz(-1.7130835) q[3];
sx q[3];
rz(-1.5350966) q[3];
sx q[3];
rz(-1.3653099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7431006) q[2];
sx q[2];
rz(-2.3091381) q[2];
sx q[2];
rz(-1.8771646) q[2];
rz(0.33665952) q[3];
sx q[3];
rz(-2.2008937) q[3];
sx q[3];
rz(1.3338859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7548512) q[0];
sx q[0];
rz(-1.6254397) q[0];
sx q[0];
rz(2.0425015) q[0];
rz(-0.59973888) q[1];
sx q[1];
rz(-2.6987684) q[1];
sx q[1];
rz(-0.59785688) q[1];
rz(-2.8200321) q[2];
sx q[2];
rz(-1.9637932) q[2];
sx q[2];
rz(3.1153352) q[2];
rz(1.5883432) q[3];
sx q[3];
rz(-1.038279) q[3];
sx q[3];
rz(-0.013171218) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
