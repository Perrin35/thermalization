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
rz(-2.6851299) q[0];
sx q[0];
rz(-0.87378341) q[0];
sx q[0];
rz(-1.1377347) q[0];
rz(-3.0955834) q[1];
sx q[1];
rz(0.53541056) q[1];
sx q[1];
rz(14.169197) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41661501) q[0];
sx q[0];
rz(-3.0862245) q[0];
sx q[0];
rz(2.2105818) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0191604) q[2];
sx q[2];
rz(-2.2847567) q[2];
sx q[2];
rz(-3.0453504) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5933758) q[1];
sx q[1];
rz(-1.1187828) q[1];
sx q[1];
rz(-1.7991245) q[1];
rz(-pi) q[2];
rz(1.7102431) q[3];
sx q[3];
rz(-1.6260423) q[3];
sx q[3];
rz(-1.6193438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5753182) q[2];
sx q[2];
rz(-2.7988269) q[2];
sx q[2];
rz(-2.9052367) q[2];
rz(-2.6929839) q[3];
sx q[3];
rz(-1.4834504) q[3];
sx q[3];
rz(-2.1382704) q[3];
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
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3818632) q[0];
sx q[0];
rz(-2.6549082) q[0];
sx q[0];
rz(2.3554262) q[0];
rz(0.78474125) q[1];
sx q[1];
rz(-1.6678383) q[1];
sx q[1];
rz(-0.10202185) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1624296) q[0];
sx q[0];
rz(-1.3693045) q[0];
sx q[0];
rz(0.54690806) q[0];
x q[1];
rz(-1.2743852) q[2];
sx q[2];
rz(-2.3326932) q[2];
sx q[2];
rz(-2.1975225) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0077304) q[1];
sx q[1];
rz(-1.8419767) q[1];
sx q[1];
rz(0.31409632) q[1];
rz(2.3524093) q[3];
sx q[3];
rz(-1.7459646) q[3];
sx q[3];
rz(-3.0449245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.027779) q[2];
sx q[2];
rz(-1.6543829) q[2];
sx q[2];
rz(0.2629183) q[2];
rz(1.8391838) q[3];
sx q[3];
rz(-2.0263367) q[3];
sx q[3];
rz(-2.3579679) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1394434) q[0];
sx q[0];
rz(-3.1378916) q[0];
sx q[0];
rz(-1.8078467) q[0];
rz(-1.3847146) q[1];
sx q[1];
rz(-2.2127071) q[1];
sx q[1];
rz(-0.76470107) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0686808) q[0];
sx q[0];
rz(-0.31038767) q[0];
sx q[0];
rz(-0.68279065) q[0];
rz(-pi) q[1];
rz(-0.5506634) q[2];
sx q[2];
rz(-2.705859) q[2];
sx q[2];
rz(-2.7274362) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4746318) q[1];
sx q[1];
rz(-1.3159915) q[1];
sx q[1];
rz(-0.17837015) q[1];
rz(-pi) q[2];
rz(1.9165048) q[3];
sx q[3];
rz(-2.3937493) q[3];
sx q[3];
rz(-1.2338828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1769522) q[2];
sx q[2];
rz(-1.2531589) q[2];
sx q[2];
rz(-2.9580252) q[2];
rz(-0.14487264) q[3];
sx q[3];
rz(-1.8795857) q[3];
sx q[3];
rz(0.22411331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1862828) q[0];
sx q[0];
rz(-1.1306385) q[0];
sx q[0];
rz(-2.8593707) q[0];
rz(1.9918293) q[1];
sx q[1];
rz(-1.7825922) q[1];
sx q[1];
rz(1.5001635) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.041752664) q[0];
sx q[0];
rz(-0.40099335) q[0];
sx q[0];
rz(-2.6869511) q[0];
x q[1];
rz(2.8871782) q[2];
sx q[2];
rz(-2.1037648) q[2];
sx q[2];
rz(-2.0526469) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.73622656) q[1];
sx q[1];
rz(-0.42441503) q[1];
sx q[1];
rz(-2.2176803) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.30434609) q[3];
sx q[3];
rz(-1.658421) q[3];
sx q[3];
rz(0.1803785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.942261) q[2];
sx q[2];
rz(-1.4456238) q[2];
sx q[2];
rz(1.65421) q[2];
rz(0.88996327) q[3];
sx q[3];
rz(-1.7421937) q[3];
sx q[3];
rz(2.9922805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99059659) q[0];
sx q[0];
rz(-2.6166333) q[0];
sx q[0];
rz(-1.6816444) q[0];
rz(1.8563942) q[1];
sx q[1];
rz(-1.3672914) q[1];
sx q[1];
rz(-0.45305124) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.734402) q[0];
sx q[0];
rz(-1.2670867) q[0];
sx q[0];
rz(0.056945368) q[0];
rz(2.1531851) q[2];
sx q[2];
rz(-1.526579) q[2];
sx q[2];
rz(-2.1427936) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.510524) q[1];
sx q[1];
rz(-1.9434526) q[1];
sx q[1];
rz(1.7700397) q[1];
rz(-pi) q[2];
rz(-2.3052081) q[3];
sx q[3];
rz(-0.57698133) q[3];
sx q[3];
rz(1.0671237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.65752658) q[2];
sx q[2];
rz(-0.47407293) q[2];
sx q[2];
rz(1.5849812) q[2];
rz(-1.5786242) q[3];
sx q[3];
rz(-1.2789187) q[3];
sx q[3];
rz(-1.0274308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22353657) q[0];
sx q[0];
rz(-0.99128381) q[0];
sx q[0];
rz(-1.1728485) q[0];
rz(0.46074834) q[1];
sx q[1];
rz(-1.2604424) q[1];
sx q[1];
rz(-1.6430829) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7448085) q[0];
sx q[0];
rz(-2.5288594) q[0];
sx q[0];
rz(-1.5721129) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3042816) q[2];
sx q[2];
rz(-1.892748) q[2];
sx q[2];
rz(2.6139174) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9869796) q[1];
sx q[1];
rz(-1.8334532) q[1];
sx q[1];
rz(3.0820552) q[1];
rz(2.2727358) q[3];
sx q[3];
rz(-1.464499) q[3];
sx q[3];
rz(0.93206638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0074761) q[2];
sx q[2];
rz(-1.4223998) q[2];
sx q[2];
rz(2.9134992) q[2];
rz(2.5988233) q[3];
sx q[3];
rz(-2.1398862) q[3];
sx q[3];
rz(-0.91226474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0246564) q[0];
sx q[0];
rz(-1.0448562) q[0];
sx q[0];
rz(0.60633099) q[0];
rz(1.8527276) q[1];
sx q[1];
rz(-1.5325129) q[1];
sx q[1];
rz(-0.85561633) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1042306) q[0];
sx q[0];
rz(-2.3920569) q[0];
sx q[0];
rz(-0.61893344) q[0];
rz(-0.087302999) q[2];
sx q[2];
rz(-2.3663967) q[2];
sx q[2];
rz(-2.6784507) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2860171) q[1];
sx q[1];
rz(-0.97980503) q[1];
sx q[1];
rz(-1.819064) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0511914) q[3];
sx q[3];
rz(-1.3452936) q[3];
sx q[3];
rz(2.4236039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3768846) q[2];
sx q[2];
rz(-0.43262425) q[2];
sx q[2];
rz(-3.063859) q[2];
rz(3.0149095) q[3];
sx q[3];
rz(-1.0418016) q[3];
sx q[3];
rz(-2.3104987) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7938101) q[0];
sx q[0];
rz(-0.3599444) q[0];
sx q[0];
rz(2.1600294) q[0];
rz(1.3579824) q[1];
sx q[1];
rz(-0.49109083) q[1];
sx q[1];
rz(-0.29409274) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61035778) q[0];
sx q[0];
rz(-1.2623275) q[0];
sx q[0];
rz(1.6334794) q[0];
x q[1];
rz(-2.7693439) q[2];
sx q[2];
rz(-1.8342092) q[2];
sx q[2];
rz(1.5836704) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8998125) q[1];
sx q[1];
rz(-0.93437785) q[1];
sx q[1];
rz(-2.7905725) q[1];
rz(-pi) q[2];
rz(0.96947396) q[3];
sx q[3];
rz(-3.0532928) q[3];
sx q[3];
rz(1.2754319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.10245094) q[2];
sx q[2];
rz(-0.65639085) q[2];
sx q[2];
rz(1.0605109) q[2];
rz(-0.55824009) q[3];
sx q[3];
rz(-0.57359901) q[3];
sx q[3];
rz(3.1225045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7790262) q[0];
sx q[0];
rz(-2.9988852) q[0];
sx q[0];
rz(-2.3574164) q[0];
rz(1.342429) q[1];
sx q[1];
rz(-1.8524086) q[1];
sx q[1];
rz(-2.4450891) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7292228) q[0];
sx q[0];
rz(-1.7984838) q[0];
sx q[0];
rz(-0.21842893) q[0];
x q[1];
rz(2.0663459) q[2];
sx q[2];
rz(-2.181567) q[2];
sx q[2];
rz(1.8004168) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.9827658) q[1];
sx q[1];
rz(-0.67081314) q[1];
sx q[1];
rz(2.4300086) q[1];
rz(2.302565) q[3];
sx q[3];
rz(-2.9425042) q[3];
sx q[3];
rz(-1.8323048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.8037618) q[2];
sx q[2];
rz(-2.3904114) q[2];
sx q[2];
rz(1.8776228) q[2];
rz(1.3761282) q[3];
sx q[3];
rz(-2.2270146) q[3];
sx q[3];
rz(-1.5862563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(-0.92397583) q[0];
sx q[0];
rz(-0.042984977) q[0];
sx q[0];
rz(-1.6643583) q[0];
rz(-2.2337275) q[1];
sx q[1];
rz(-1.6198747) q[1];
sx q[1];
rz(-0.97000617) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15031397) q[0];
sx q[0];
rz(-1.3320001) q[0];
sx q[0];
rz(-2.1697609) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97899784) q[2];
sx q[2];
rz(-1.3221957) q[2];
sx q[2];
rz(-0.73038855) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1906311) q[1];
sx q[1];
rz(-0.64563692) q[1];
sx q[1];
rz(-0.26107045) q[1];
rz(2.4360551) q[3];
sx q[3];
rz(-0.85800401) q[3];
sx q[3];
rz(2.4487011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9475391) q[2];
sx q[2];
rz(-0.87475646) q[2];
sx q[2];
rz(0.27210316) q[2];
rz(-0.84723204) q[3];
sx q[3];
rz(-2.6341485) q[3];
sx q[3];
rz(2.0525172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67615164) q[0];
sx q[0];
rz(-1.1347102) q[0];
sx q[0];
rz(1.4665428) q[0];
rz(-0.33879406) q[1];
sx q[1];
rz(-1.6962961) q[1];
sx q[1];
rz(1.2512339) q[1];
rz(-0.39739824) q[2];
sx q[2];
rz(-0.76042475) q[2];
sx q[2];
rz(-2.3971192) q[2];
rz(0.083214464) q[3];
sx q[3];
rz(-1.1490001) q[3];
sx q[3];
rz(0.84925539) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
