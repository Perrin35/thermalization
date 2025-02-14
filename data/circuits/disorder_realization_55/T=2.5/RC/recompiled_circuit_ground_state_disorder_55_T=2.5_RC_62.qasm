OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.80630535) q[0];
sx q[0];
rz(-2.1016313) q[0];
sx q[0];
rz(0.34997532) q[0];
rz(1.5179874) q[1];
sx q[1];
rz(-2.1972158) q[1];
sx q[1];
rz(1.1854393) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3458334) q[0];
sx q[0];
rz(-1.1544219) q[0];
sx q[0];
rz(-1.9365947) q[0];
x q[1];
rz(1.8030422) q[2];
sx q[2];
rz(-1.9793538) q[2];
sx q[2];
rz(0.76257715) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6373316) q[1];
sx q[1];
rz(-2.0954014) q[1];
sx q[1];
rz(-0.27224147) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8270742) q[3];
sx q[3];
rz(-2.5102237) q[3];
sx q[3];
rz(2.6169962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8650633) q[2];
sx q[2];
rz(-0.41215602) q[2];
sx q[2];
rz(-2.7318447) q[2];
rz(-0.77541882) q[3];
sx q[3];
rz(-1.5284458) q[3];
sx q[3];
rz(0.30737901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2839117) q[0];
sx q[0];
rz(-2.6429521) q[0];
sx q[0];
rz(-0.45251244) q[0];
rz(2.9127938) q[1];
sx q[1];
rz(-0.50140536) q[1];
sx q[1];
rz(-2.9948044) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89167021) q[0];
sx q[0];
rz(-1.9324713) q[0];
sx q[0];
rz(-0.2043496) q[0];
x q[1];
rz(0.37756672) q[2];
sx q[2];
rz(-1.7930328) q[2];
sx q[2];
rz(1.3106032) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1498857) q[1];
sx q[1];
rz(-3.001431) q[1];
sx q[1];
rz(2.8542942) q[1];
rz(-pi) q[2];
x q[2];
rz(0.66007075) q[3];
sx q[3];
rz(-1.9418896) q[3];
sx q[3];
rz(1.1255655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.16810922) q[2];
sx q[2];
rz(-0.3349458) q[2];
sx q[2];
rz(-2.2701263) q[2];
rz(-2.787309) q[3];
sx q[3];
rz(-0.93577093) q[3];
sx q[3];
rz(-1.833545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.5062434) q[0];
sx q[0];
rz(-0.62017089) q[0];
sx q[0];
rz(2.035602) q[0];
rz(-0.94331074) q[1];
sx q[1];
rz(-0.7990852) q[1];
sx q[1];
rz(-1.6518637) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7794253) q[0];
sx q[0];
rz(-0.12813103) q[0];
sx q[0];
rz(1.4904899) q[0];
x q[1];
rz(1.4749799) q[2];
sx q[2];
rz(-1.7783217) q[2];
sx q[2];
rz(2.735161) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.572444) q[1];
sx q[1];
rz(-0.97578543) q[1];
sx q[1];
rz(1.5250564) q[1];
rz(-pi) q[2];
rz(1.6451646) q[3];
sx q[3];
rz(-1.8083124) q[3];
sx q[3];
rz(0.46922153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.6279383) q[2];
sx q[2];
rz(-1.9460121) q[2];
sx q[2];
rz(-2.1453843) q[2];
rz(0.57322383) q[3];
sx q[3];
rz(-0.51155353) q[3];
sx q[3];
rz(2.5073124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117677) q[0];
sx q[0];
rz(-1.2720164) q[0];
sx q[0];
rz(1.5559394) q[0];
rz(0.32444435) q[1];
sx q[1];
rz(-0.8322081) q[1];
sx q[1];
rz(-1.9437887) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65069136) q[0];
sx q[0];
rz(-1.8244716) q[0];
sx q[0];
rz(-0.59532292) q[0];
x q[1];
rz(2.156331) q[2];
sx q[2];
rz(-1.6949495) q[2];
sx q[2];
rz(0.92701605) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50655925) q[1];
sx q[1];
rz(-2.6406619) q[1];
sx q[1];
rz(-1.8440203) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3088552) q[3];
sx q[3];
rz(-1.4125694) q[3];
sx q[3];
rz(-0.55486521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0962301) q[2];
sx q[2];
rz(-2.5924293) q[2];
sx q[2];
rz(-1.9722923) q[2];
rz(0.62396389) q[3];
sx q[3];
rz(-1.4493161) q[3];
sx q[3];
rz(-0.72871488) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8600334) q[0];
sx q[0];
rz(-2.7427854) q[0];
sx q[0];
rz(1.1997724) q[0];
rz(-1.3072321) q[1];
sx q[1];
rz(-0.93991005) q[1];
sx q[1];
rz(-1.7050381) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2685674) q[0];
sx q[0];
rz(-1.6313261) q[0];
sx q[0];
rz(1.1824754) q[0];
x q[1];
rz(-1.8017112) q[2];
sx q[2];
rz(-1.4428939) q[2];
sx q[2];
rz(-0.85258871) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.58545326) q[1];
sx q[1];
rz(-1.2188497) q[1];
sx q[1];
rz(-1.8896249) q[1];
rz(-pi) q[2];
rz(-0.41884274) q[3];
sx q[3];
rz(-0.9658893) q[3];
sx q[3];
rz(-2.1640282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.78937078) q[2];
sx q[2];
rz(-1.7033966) q[2];
sx q[2];
rz(2.5686725) q[2];
rz(2.2574183) q[3];
sx q[3];
rz(-2.1638162) q[3];
sx q[3];
rz(2.6463215) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5837412) q[0];
sx q[0];
rz(-1.8889677) q[0];
sx q[0];
rz(2.0649233) q[0];
rz(0.50648266) q[1];
sx q[1];
rz(-0.61512893) q[1];
sx q[1];
rz(1.8411676) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.444435) q[0];
sx q[0];
rz(-1.5845254) q[0];
sx q[0];
rz(-1.5221217) q[0];
rz(2.9175148) q[2];
sx q[2];
rz(-2.5505318) q[2];
sx q[2];
rz(-3.0853809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4757931) q[1];
sx q[1];
rz(-2.8912326) q[1];
sx q[1];
rz(-2.1264892) q[1];
rz(-pi) q[2];
rz(-1.0664942) q[3];
sx q[3];
rz(-0.39736727) q[3];
sx q[3];
rz(0.25542828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41636813) q[2];
sx q[2];
rz(-2.29795) q[2];
sx q[2];
rz(1.9492524) q[2];
rz(-0.013817712) q[3];
sx q[3];
rz(-2.4255224) q[3];
sx q[3];
rz(2.1684087) q[3];
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
rz(-pi) q[0];
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
rz(-1.248812) q[0];
sx q[0];
rz(-2.2154494) q[0];
sx q[0];
rz(-2.7698351) q[0];
rz(-0.88428503) q[1];
sx q[1];
rz(-2.6339032) q[1];
sx q[1];
rz(0.54381347) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9659198) q[0];
sx q[0];
rz(-1.9750711) q[0];
sx q[0];
rz(-0.61575295) q[0];
rz(-pi) q[1];
rz(1.061166) q[2];
sx q[2];
rz(-2.3544534) q[2];
sx q[2];
rz(0.9796046) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2993967) q[1];
sx q[1];
rz(-0.37494117) q[1];
sx q[1];
rz(3.123466) q[1];
rz(-pi) q[2];
rz(0.94856222) q[3];
sx q[3];
rz(-0.18634054) q[3];
sx q[3];
rz(-2.2295913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.030293327) q[2];
sx q[2];
rz(-2.0798422) q[2];
sx q[2];
rz(0.57478762) q[2];
rz(-2.287367) q[3];
sx q[3];
rz(-1.4762907) q[3];
sx q[3];
rz(-2.0265719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92693555) q[0];
sx q[0];
rz(-0.68100005) q[0];
sx q[0];
rz(0.32447234) q[0];
rz(2.532161) q[1];
sx q[1];
rz(-0.84171265) q[1];
sx q[1];
rz(-2.6268974) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99583944) q[0];
sx q[0];
rz(-1.1356562) q[0];
sx q[0];
rz(0.010220411) q[0];
x q[1];
rz(1.1616522) q[2];
sx q[2];
rz(-2.3858641) q[2];
sx q[2];
rz(-0.31015304) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.62186723) q[1];
sx q[1];
rz(-0.61841209) q[1];
sx q[1];
rz(-2.7296978) q[1];
rz(-pi) q[2];
rz(-0.10917337) q[3];
sx q[3];
rz(-1.2965373) q[3];
sx q[3];
rz(1.7133939) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1631761) q[2];
sx q[2];
rz(-1.5055483) q[2];
sx q[2];
rz(0.27463883) q[2];
rz(-0.93242532) q[3];
sx q[3];
rz(-2.7115188) q[3];
sx q[3];
rz(2.8453804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68359971) q[0];
sx q[0];
rz(-1.6807115) q[0];
sx q[0];
rz(3.006111) q[0];
rz(0.97700351) q[1];
sx q[1];
rz(-0.75825399) q[1];
sx q[1];
rz(-0.0010842222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8196306) q[0];
sx q[0];
rz(-2.062383) q[0];
sx q[0];
rz(1.7474773) q[0];
rz(2.094467) q[2];
sx q[2];
rz(-1.6738679) q[2];
sx q[2];
rz(1.4541986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7338176) q[1];
sx q[1];
rz(-2.1548591) q[1];
sx q[1];
rz(-2.255846) q[1];
rz(-0.67768456) q[3];
sx q[3];
rz(-2.2608058) q[3];
sx q[3];
rz(-0.023493903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.5929426) q[2];
sx q[2];
rz(-1.0385916) q[2];
sx q[2];
rz(-0.20498094) q[2];
rz(-0.18260469) q[3];
sx q[3];
rz(-0.20668106) q[3];
sx q[3];
rz(2.4038147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36271998) q[0];
sx q[0];
rz(-0.39411476) q[0];
sx q[0];
rz(0.47738281) q[0];
rz(0.1167156) q[1];
sx q[1];
rz(-1.621403) q[1];
sx q[1];
rz(0.83888549) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.981659) q[0];
sx q[0];
rz(-0.75272876) q[0];
sx q[0];
rz(-1.9636867) q[0];
x q[1];
rz(2.7493189) q[2];
sx q[2];
rz(-2.3773411) q[2];
sx q[2];
rz(1.6583673) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.5236256) q[1];
sx q[1];
rz(-2.2516421) q[1];
sx q[1];
rz(1.3300072) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.1118591) q[3];
sx q[3];
rz(-1.8095867) q[3];
sx q[3];
rz(2.5532171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.78187752) q[2];
sx q[2];
rz(-2.8350267) q[2];
sx q[2];
rz(0.26620418) q[2];
rz(-0.5075469) q[3];
sx q[3];
rz(-0.72269732) q[3];
sx q[3];
rz(-0.44107309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97892852) q[0];
sx q[0];
rz(-1.6829818) q[0];
sx q[0];
rz(-0.80243954) q[0];
rz(-0.91213999) q[1];
sx q[1];
rz(-1.1911387) q[1];
sx q[1];
rz(-2.3008507) q[1];
rz(2.5723614) q[2];
sx q[2];
rz(-2.8494159) q[2];
sx q[2];
rz(-0.055622774) q[2];
rz(0.75822266) q[3];
sx q[3];
rz(-0.52715404) q[3];
sx q[3];
rz(-0.84670443) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
