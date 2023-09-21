OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34713137) q[0];
sx q[0];
rz(5.2678582) q[0];
sx q[0];
rz(9.8922748) q[0];
rz(-0.52019083) q[1];
sx q[1];
rz(-1.3462892) q[1];
sx q[1];
rz(-0.68038124) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7079948) q[0];
sx q[0];
rz(-2.2888219) q[0];
sx q[0];
rz(-0.17311592) q[0];
rz(-0.19619588) q[2];
sx q[2];
rz(-1.2805403) q[2];
sx q[2];
rz(0.44858518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8531224) q[1];
sx q[1];
rz(-0.56325699) q[1];
sx q[1];
rz(-2.017615) q[1];
x q[2];
rz(-2.8075571) q[3];
sx q[3];
rz(-1.4035657) q[3];
sx q[3];
rz(-0.86080307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9900069) q[2];
sx q[2];
rz(-0.9881343) q[2];
sx q[2];
rz(-0.087466784) q[2];
rz(0.72922373) q[3];
sx q[3];
rz(-0.37344033) q[3];
sx q[3];
rz(0.27174404) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4366348) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(1.0013642) q[0];
rz(0.17240605) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66602089) q[0];
sx q[0];
rz(-2.4040939) q[0];
sx q[0];
rz(-1.1230418) q[0];
x q[1];
rz(2.0662202) q[2];
sx q[2];
rz(-2.3790857) q[2];
sx q[2];
rz(2.5410595) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9973035) q[1];
sx q[1];
rz(-1.3376437) q[1];
sx q[1];
rz(-2.7551329) q[1];
rz(-pi) q[2];
rz(1.4161795) q[3];
sx q[3];
rz(-1.039045) q[3];
sx q[3];
rz(-2.1173409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-0.45980644) q[2];
sx q[2];
rz(1.8015507) q[2];
rz(-0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-0.036858233) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3443417) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(-0.62477338) q[0];
rz(2.1773188) q[1];
sx q[1];
rz(-0.48502973) q[1];
sx q[1];
rz(-2.952081) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0441168) q[0];
sx q[0];
rz(-2.2783845) q[0];
sx q[0];
rz(-1.9301027) q[0];
rz(-pi) q[1];
rz(-2.8787829) q[2];
sx q[2];
rz(-2.7542979) q[2];
sx q[2];
rz(0.58236052) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.5513969) q[1];
sx q[1];
rz(-1.8101748) q[1];
sx q[1];
rz(-1.6281284) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8241006) q[3];
sx q[3];
rz(-2.5251303) q[3];
sx q[3];
rz(-0.14641078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0345962) q[2];
sx q[2];
rz(-2.2980289) q[2];
sx q[2];
rz(1.3872046) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(-2.0534024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65961924) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(1.5455998) q[0];
rz(2.003147) q[1];
sx q[1];
rz(-0.7754511) q[1];
sx q[1];
rz(-2.0514354) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2290303) q[0];
sx q[0];
rz(-1.9711442) q[0];
sx q[0];
rz(1.1926665) q[0];
x q[1];
rz(-2.3218669) q[2];
sx q[2];
rz(-1.1851386) q[2];
sx q[2];
rz(-2.0476598) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.7294457) q[1];
sx q[1];
rz(-2.4033961) q[1];
sx q[1];
rz(0.63936887) q[1];
rz(1.3642717) q[3];
sx q[3];
rz(-1.1518475) q[3];
sx q[3];
rz(-1.2057613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.99889341) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(2.2857655) q[2];
rz(1.9479729) q[3];
sx q[3];
rz(-1.6222298) q[3];
sx q[3];
rz(-0.91317552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49044931) q[0];
sx q[0];
rz(-0.97390807) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(-0.99114746) q[1];
sx q[1];
rz(-1.2065572) q[1];
sx q[1];
rz(1.1700464) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67384185) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(-0.90676102) q[0];
rz(-1.9636376) q[2];
sx q[2];
rz(-2.422214) q[2];
sx q[2];
rz(0.17578416) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7157117) q[1];
sx q[1];
rz(-2.5320146) q[1];
sx q[1];
rz(0.73519911) q[1];
rz(-pi) q[2];
rz(-2.823579) q[3];
sx q[3];
rz(-2.5807096) q[3];
sx q[3];
rz(-2.6624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2698007) q[2];
sx q[2];
rz(-1.5503927) q[2];
sx q[2];
rz(3.0043547) q[2];
rz(-1.8042701) q[3];
sx q[3];
rz(-0.58115712) q[3];
sx q[3];
rz(-1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18734922) q[0];
sx q[0];
rz(-1.3784778) q[0];
sx q[0];
rz(-1.3504008) q[0];
rz(-2.2987135) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(2.4687016) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.580299) q[0];
sx q[0];
rz(-1.467388) q[0];
sx q[0];
rz(2.1456477) q[0];
rz(-pi) q[1];
rz(-1.5422836) q[2];
sx q[2];
rz(-2.4270504) q[2];
sx q[2];
rz(1.6376094) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.37398794) q[1];
sx q[1];
rz(-1.7786221) q[1];
sx q[1];
rz(-0.61088224) q[1];
x q[2];
rz(-1.4937431) q[3];
sx q[3];
rz(-1.9814081) q[3];
sx q[3];
rz(-1.3239469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.3488397) q[2];
sx q[2];
rz(-1.2092084) q[2];
sx q[2];
rz(-2.7065281) q[2];
rz(1.3600291) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(-2.8939261) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.172794) q[0];
sx q[0];
rz(-2.2491169) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(-1.1116213) q[1];
sx q[1];
rz(-1.9042791) q[1];
sx q[1];
rz(-1.7395082) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60915011) q[0];
sx q[0];
rz(-0.46450588) q[0];
sx q[0];
rz(0.64446489) q[0];
rz(-pi) q[1];
rz(2.3767396) q[2];
sx q[2];
rz(-1.7424889) q[2];
sx q[2];
rz(-2.9219251) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40560383) q[1];
sx q[1];
rz(-2.7466008) q[1];
sx q[1];
rz(-2.7337381) q[1];
x q[2];
rz(-0.98251179) q[3];
sx q[3];
rz(-1.4559828) q[3];
sx q[3];
rz(-1.3548152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7193675) q[2];
sx q[2];
rz(-2.5740467) q[2];
sx q[2];
rz(-2.4613703) q[2];
rz(-2.7133572) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(-1.359882) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.632804) q[0];
sx q[0];
rz(-1.4161685) q[0];
sx q[0];
rz(1.592214) q[0];
rz(0.24958615) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(0.54135281) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68529785) q[0];
sx q[0];
rz(-2.1349499) q[0];
sx q[0];
rz(-1.5840522) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63033732) q[2];
sx q[2];
rz(-0.54140831) q[2];
sx q[2];
rz(-0.88592096) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7975446) q[1];
sx q[1];
rz(-2.27564) q[1];
sx q[1];
rz(2.123358) q[1];
rz(-2.7612711) q[3];
sx q[3];
rz(-1.5725279) q[3];
sx q[3];
rz(-2.1960432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.044518746) q[2];
sx q[2];
rz(-2.7842583) q[2];
sx q[2];
rz(2.3642335) q[2];
rz(0.85123953) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(1.8204934) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7580496) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(-3.1316277) q[0];
rz(2.126157) q[1];
sx q[1];
rz(-0.76428691) q[1];
sx q[1];
rz(-1.4452971) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1709258) q[0];
sx q[0];
rz(-2.2609841) q[0];
sx q[0];
rz(0.70113457) q[0];
rz(-pi) q[1];
rz(-0.5814914) q[2];
sx q[2];
rz(-0.98052374) q[2];
sx q[2];
rz(1.7828538) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.5053133) q[1];
sx q[1];
rz(-1.5552102) q[1];
sx q[1];
rz(-0.63025766) q[1];
rz(-pi) q[2];
x q[2];
rz(0.61049283) q[3];
sx q[3];
rz(-0.99184147) q[3];
sx q[3];
rz(2.0905153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4426667) q[2];
sx q[2];
rz(-0.93124229) q[2];
sx q[2];
rz(0.60738579) q[2];
rz(-1.4390885) q[3];
sx q[3];
rz(-1.3834229) q[3];
sx q[3];
rz(2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8059175) q[0];
sx q[0];
rz(-1.6573925) q[0];
sx q[0];
rz(2.5277396) q[0];
rz(1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(-0.39224958) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0961571) q[0];
sx q[0];
rz(-2.3491612) q[0];
sx q[0];
rz(2.5655454) q[0];
x q[1];
rz(-0.73108436) q[2];
sx q[2];
rz(-0.77027551) q[2];
sx q[2];
rz(-2.0172271) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0543921) q[1];
sx q[1];
rz(-1.1322081) q[1];
sx q[1];
rz(2.2589576) q[1];
x q[2];
rz(1.5997821) q[3];
sx q[3];
rz(-1.3658804) q[3];
sx q[3];
rz(0.25587413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.3839046) q[2];
sx q[2];
rz(2.5411141) q[2];
rz(1.0673808) q[3];
sx q[3];
rz(-1.821527) q[3];
sx q[3];
rz(-1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7077211) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(2.1451163) q[1];
sx q[1];
rz(-1.6468208) q[1];
sx q[1];
rz(-1.5368808) q[1];
rz(-0.36841064) q[2];
sx q[2];
rz(-2.340292) q[2];
sx q[2];
rz(3.1030263) q[2];
rz(2.2157833) q[3];
sx q[3];
rz(-2.3498597) q[3];
sx q[3];
rz(1.3510977) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
