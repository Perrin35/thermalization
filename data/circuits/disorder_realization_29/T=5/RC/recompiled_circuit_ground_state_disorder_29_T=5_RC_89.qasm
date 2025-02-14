OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(3.0512515) q[0];
sx q[0];
rz(-2.0988965) q[0];
sx q[0];
rz(-0.27867499) q[0];
rz(-1.9947808) q[1];
sx q[1];
rz(-2.6169701) q[1];
sx q[1];
rz(-1.8705179) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88408121) q[0];
sx q[0];
rz(-1.8313098) q[0];
sx q[0];
rz(-1.5896912) q[0];
rz(-0.15058168) q[2];
sx q[2];
rz(-1.4408334) q[2];
sx q[2];
rz(-0.58565631) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8184389) q[1];
sx q[1];
rz(-0.85994321) q[1];
sx q[1];
rz(3.0368785) q[1];
rz(0.60802976) q[3];
sx q[3];
rz(-1.078106) q[3];
sx q[3];
rz(-1.0580491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.110454) q[2];
sx q[2];
rz(-1.9240856) q[2];
sx q[2];
rz(-2.8793867) q[2];
rz(-2.0860705) q[3];
sx q[3];
rz(-1.2473829) q[3];
sx q[3];
rz(3.054255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6156886) q[0];
sx q[0];
rz(-2.6380802) q[0];
sx q[0];
rz(0.18158922) q[0];
rz(-1.4008105) q[1];
sx q[1];
rz(-1.9488275) q[1];
sx q[1];
rz(2.3006732) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5530295) q[0];
sx q[0];
rz(-2.8210777) q[0];
sx q[0];
rz(3.0856087) q[0];
rz(-pi) q[1];
rz(-0.96581755) q[2];
sx q[2];
rz(-0.73558319) q[2];
sx q[2];
rz(-3.0515665) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.052472981) q[1];
sx q[1];
rz(-1.5346255) q[1];
sx q[1];
rz(-1.7391886) q[1];
rz(2.0729468) q[3];
sx q[3];
rz(-0.6752033) q[3];
sx q[3];
rz(-1.6197551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.759364) q[2];
sx q[2];
rz(-0.21491924) q[2];
sx q[2];
rz(-1.0648897) q[2];
rz(-0.061428849) q[3];
sx q[3];
rz(-1.3353142) q[3];
sx q[3];
rz(-1.6399062) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29838022) q[0];
sx q[0];
rz(-1.9623663) q[0];
sx q[0];
rz(2.9425353) q[0];
rz(2.3152323) q[1];
sx q[1];
rz(-2.2233456) q[1];
sx q[1];
rz(0.77702776) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3129757) q[0];
sx q[0];
rz(-1.6369633) q[0];
sx q[0];
rz(1.7278759) q[0];
rz(-pi) q[1];
rz(0.92240693) q[2];
sx q[2];
rz(-2.4925579) q[2];
sx q[2];
rz(0.051608406) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9183082) q[1];
sx q[1];
rz(-1.8294083) q[1];
sx q[1];
rz(2.1373939) q[1];
x q[2];
rz(-0.97717173) q[3];
sx q[3];
rz(-1.3389753) q[3];
sx q[3];
rz(-2.1013732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0394773) q[2];
sx q[2];
rz(-1.9934883) q[2];
sx q[2];
rz(-2.6285505) q[2];
rz(0.78220621) q[3];
sx q[3];
rz(-1.2057883) q[3];
sx q[3];
rz(-1.3768844) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2143329) q[0];
sx q[0];
rz(-1.5434649) q[0];
sx q[0];
rz(-1.4904892) q[0];
rz(-0.53228846) q[1];
sx q[1];
rz(-2.5106301) q[1];
sx q[1];
rz(1.5375563) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2435115) q[0];
sx q[0];
rz(-1.4644721) q[0];
sx q[0];
rz(2.9778019) q[0];
rz(3.0481993) q[2];
sx q[2];
rz(-2.1150781) q[2];
sx q[2];
rz(1.5285847) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2080431) q[1];
sx q[1];
rz(-1.6416104) q[1];
sx q[1];
rz(2.2068702) q[1];
rz(-pi) q[2];
x q[2];
rz(0.38942899) q[3];
sx q[3];
rz(-0.5618605) q[3];
sx q[3];
rz(1.495468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.535546) q[2];
sx q[2];
rz(-0.54680768) q[2];
sx q[2];
rz(-0.01288506) q[2];
rz(-1.9384725) q[3];
sx q[3];
rz(-1.4667526) q[3];
sx q[3];
rz(-1.0217246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0032229) q[0];
sx q[0];
rz(-2.5165181) q[0];
sx q[0];
rz(-1.3997929) q[0];
rz(-2.7895582) q[1];
sx q[1];
rz(-1.0540009) q[1];
sx q[1];
rz(-0.83784109) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.261917) q[0];
sx q[0];
rz(-2.3220372) q[0];
sx q[0];
rz(2.9330334) q[0];
rz(-pi) q[1];
rz(2.3532392) q[2];
sx q[2];
rz(-2.3559795) q[2];
sx q[2];
rz(2.2503302) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.074118) q[1];
sx q[1];
rz(-2.5606541) q[1];
sx q[1];
rz(0.23260637) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9772163) q[3];
sx q[3];
rz(-1.7225186) q[3];
sx q[3];
rz(-2.3325932) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0248523) q[2];
sx q[2];
rz(-0.96472538) q[2];
sx q[2];
rz(-1.9579197) q[2];
rz(-0.27070326) q[3];
sx q[3];
rz(-0.94047061) q[3];
sx q[3];
rz(1.5326327) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1149513) q[0];
sx q[0];
rz(-2.5178435) q[0];
sx q[0];
rz(0.57914105) q[0];
rz(-1.9193513) q[1];
sx q[1];
rz(-1.3074343) q[1];
sx q[1];
rz(2.0319669) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5752788) q[0];
sx q[0];
rz(-1.405861) q[0];
sx q[0];
rz(1.724302) q[0];
rz(-pi) q[1];
rz(-2.9640084) q[2];
sx q[2];
rz(-2.8320667) q[2];
sx q[2];
rz(-0.9847275) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0640969) q[1];
sx q[1];
rz(-1.9019097) q[1];
sx q[1];
rz(1.3084133) q[1];
rz(-pi) q[2];
rz(-0.51646502) q[3];
sx q[3];
rz(-1.0871917) q[3];
sx q[3];
rz(-0.91396871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1960725) q[2];
sx q[2];
rz(-1.8438945) q[2];
sx q[2];
rz(-0.35745364) q[2];
rz(-1.9516021) q[3];
sx q[3];
rz(-1.9414732) q[3];
sx q[3];
rz(-1.801871) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59915197) q[0];
sx q[0];
rz(-1.0170794) q[0];
sx q[0];
rz(-0.79208148) q[0];
rz(0.7041086) q[1];
sx q[1];
rz(-0.45583615) q[1];
sx q[1];
rz(1.5584996) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9607199) q[0];
sx q[0];
rz(-1.5737434) q[0];
sx q[0];
rz(2.6774027) q[0];
rz(0.59410166) q[2];
sx q[2];
rz(-2.1081446) q[2];
sx q[2];
rz(1.0460451) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.0593554) q[1];
sx q[1];
rz(-1.6052639) q[1];
sx q[1];
rz(-0.42621526) q[1];
rz(-pi) q[2];
x q[2];
rz(0.75140636) q[3];
sx q[3];
rz(-1.4045241) q[3];
sx q[3];
rz(2.5010482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.030297) q[2];
sx q[2];
rz(-2.3534677) q[2];
sx q[2];
rz(-3.1374068) q[2];
rz(-2.8748416) q[3];
sx q[3];
rz(-1.5240074) q[3];
sx q[3];
rz(-0.72223103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58077145) q[0];
sx q[0];
rz(-0.72951356) q[0];
sx q[0];
rz(-0.15604493) q[0];
rz(-1.5849628) q[1];
sx q[1];
rz(-0.86235756) q[1];
sx q[1];
rz(0.59476605) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9879919) q[0];
sx q[0];
rz(-0.56454851) q[0];
sx q[0];
rz(-0.49086824) q[0];
rz(-pi) q[1];
rz(-2.4492599) q[2];
sx q[2];
rz(-1.344948) q[2];
sx q[2];
rz(-0.69162265) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.76164272) q[1];
sx q[1];
rz(-2.159286) q[1];
sx q[1];
rz(2.2698927) q[1];
rz(-pi) q[2];
rz(0.72799369) q[3];
sx q[3];
rz(-2.1001951) q[3];
sx q[3];
rz(1.4024783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5742089) q[2];
sx q[2];
rz(-1.5625861) q[2];
sx q[2];
rz(-2.1136005) q[2];
rz(1.7993641) q[3];
sx q[3];
rz(-2.4864311) q[3];
sx q[3];
rz(-1.3348234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8780355) q[0];
sx q[0];
rz(-1.4565383) q[0];
sx q[0];
rz(0.79826075) q[0];
rz(2.3410666) q[1];
sx q[1];
rz(-2.6526178) q[1];
sx q[1];
rz(2.5953603) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6123264) q[0];
sx q[0];
rz(-2.5298425) q[0];
sx q[0];
rz(-1.3806369) q[0];
rz(-pi) q[1];
rz(0.58878657) q[2];
sx q[2];
rz(-1.9336039) q[2];
sx q[2];
rz(-2.6984359) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45193397) q[1];
sx q[1];
rz(-2.0469189) q[1];
sx q[1];
rz(1.0187638) q[1];
rz(-pi) q[2];
rz(1.6437552) q[3];
sx q[3];
rz(-1.7966509) q[3];
sx q[3];
rz(2.7166004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74786782) q[2];
sx q[2];
rz(-1.0900898) q[2];
sx q[2];
rz(0.063035034) q[2];
rz(-2.4437599) q[3];
sx q[3];
rz(-0.2581667) q[3];
sx q[3];
rz(2.2285018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49552712) q[0];
sx q[0];
rz(-0.87764469) q[0];
sx q[0];
rz(-0.416042) q[0];
rz(-1.7795732) q[1];
sx q[1];
rz(-1.1241309) q[1];
sx q[1];
rz(-1.8488041) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.775499) q[0];
sx q[0];
rz(-1.3028569) q[0];
sx q[0];
rz(1.3331205) q[0];
x q[1];
rz(-2.6854158) q[2];
sx q[2];
rz(-1.8229228) q[2];
sx q[2];
rz(1.1751564) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4626235) q[1];
sx q[1];
rz(-2.0324283) q[1];
sx q[1];
rz(-2.1831475) q[1];
rz(-pi) q[2];
x q[2];
rz(1.148801) q[3];
sx q[3];
rz(-1.3316324) q[3];
sx q[3];
rz(2.0668427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4585939) q[2];
sx q[2];
rz(-1.8051882) q[2];
sx q[2];
rz(2.2340753) q[2];
rz(1.8806184) q[3];
sx q[3];
rz(-2.5670299) q[3];
sx q[3];
rz(-1.449409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64611971) q[0];
sx q[0];
rz(-1.4541805) q[0];
sx q[0];
rz(-2.4021586) q[0];
rz(-0.44838913) q[1];
sx q[1];
rz(-1.3403475) q[1];
sx q[1];
rz(1.1833804) q[1];
rz(2.7896055) q[2];
sx q[2];
rz(-1.864973) q[2];
sx q[2];
rz(-0.61054338) q[2];
rz(-0.84239324) q[3];
sx q[3];
rz(-1.471023) q[3];
sx q[3];
rz(-1.4820549) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
