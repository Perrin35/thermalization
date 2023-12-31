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
rz(2.4612114) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9677833) q[0];
sx q[0];
rz(-0.73497226) q[0];
sx q[0];
rz(1.3761139) q[0];
rz(0.19619588) q[2];
sx q[2];
rz(-1.8610524) q[2];
sx q[2];
rz(0.44858518) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.89741325) q[1];
sx q[1];
rz(-1.8036098) q[1];
sx q[1];
rz(2.0884872) q[1];
rz(-2.8075571) q[3];
sx q[3];
rz(-1.4035657) q[3];
sx q[3];
rz(-0.86080307) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9900069) q[2];
sx q[2];
rz(-2.1534584) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7049578) q[0];
sx q[0];
rz(-2.3925662) q[0];
sx q[0];
rz(-1.0013642) q[0];
rz(-2.9691866) q[1];
sx q[1];
rz(-2.0253851) q[1];
sx q[1];
rz(0.52406812) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8952626) q[0];
sx q[0];
rz(-1.2753914) q[0];
sx q[0];
rz(2.2569879) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0662202) q[2];
sx q[2];
rz(-2.3790857) q[2];
sx q[2];
rz(0.60053315) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3327649) q[1];
sx q[1];
rz(-1.9462703) q[1];
sx q[1];
rz(1.3198225) q[1];
rz(-1.4161795) q[3];
sx q[3];
rz(-1.039045) q[3];
sx q[3];
rz(-1.0242517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9841763) q[2];
sx q[2];
rz(-2.6817862) q[2];
sx q[2];
rz(-1.8015507) q[2];
rz(0.79483461) q[3];
sx q[3];
rz(-2.0017616) q[3];
sx q[3];
rz(-3.1047344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79725093) q[0];
sx q[0];
rz(-2.4448555) q[0];
sx q[0];
rz(2.5168193) q[0];
rz(0.96427381) q[1];
sx q[1];
rz(-2.6565629) q[1];
sx q[1];
rz(-2.952081) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7127842) q[0];
sx q[0];
rz(-1.841294) q[0];
sx q[0];
rz(0.74032797) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37522845) q[2];
sx q[2];
rz(-1.669075) q[2];
sx q[2];
rz(-2.3972942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1358007) q[1];
sx q[1];
rz(-1.6264919) q[1];
sx q[1];
rz(-2.9018351) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.965851) q[3];
sx q[3];
rz(-0.97676859) q[3];
sx q[3];
rz(-2.9807846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1069964) q[2];
sx q[2];
rz(-0.84356374) q[2];
sx q[2];
rz(1.3872046) q[2];
rz(-2.710279) q[3];
sx q[3];
rz(-1.8547736) q[3];
sx q[3];
rz(1.0881902) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4819734) q[0];
sx q[0];
rz(-2.1932333) q[0];
sx q[0];
rz(1.5455998) q[0];
rz(-2.003147) q[1];
sx q[1];
rz(-2.3661416) q[1];
sx q[1];
rz(1.0901573) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91256234) q[0];
sx q[0];
rz(-1.1704485) q[0];
sx q[0];
rz(1.9489261) q[0];
x q[1];
rz(2.6345989) q[2];
sx q[2];
rz(-2.255313) q[2];
sx q[2];
rz(-0.13912858) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.480099) q[1];
sx q[1];
rz(-1.9839994) q[1];
sx q[1];
rz(0.6306298) q[1];
x q[2];
rz(0.4269883) q[3];
sx q[3];
rz(-1.3823576) q[3];
sx q[3];
rz(-2.6915336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.99889341) q[2];
sx q[2];
rz(-2.6901851) q[2];
sx q[2];
rz(-2.2857655) q[2];
rz(-1.1936197) q[3];
sx q[3];
rz(-1.5193628) q[3];
sx q[3];
rz(-2.2284171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49044931) q[0];
sx q[0];
rz(-2.1676846) q[0];
sx q[0];
rz(-0.14973101) q[0];
rz(-0.99114746) q[1];
sx q[1];
rz(-1.9350355) q[1];
sx q[1];
rz(1.9715462) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67384185) q[0];
sx q[0];
rz(-1.8918599) q[0];
sx q[0];
rz(0.90676102) q[0];
rz(-pi) q[1];
rz(-0.32355002) q[2];
sx q[2];
rz(-2.2252482) q[2];
sx q[2];
rz(0.3277342) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7331446) q[1];
sx q[1];
rz(-2.0093579) q[1];
sx q[1];
rz(-2.0088197) q[1];
rz(2.6036161) q[3];
sx q[3];
rz(-1.7378983) q[3];
sx q[3];
rz(0.81975421) q[3];
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
rz(-1.3373226) q[3];
sx q[3];
rz(-2.5604355) q[3];
sx q[3];
rz(-1.8491245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9542434) q[0];
sx q[0];
rz(-1.7631148) q[0];
sx q[0];
rz(1.3504008) q[0];
rz(-0.84287914) q[1];
sx q[1];
rz(-0.73892361) q[1];
sx q[1];
rz(0.67289105) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0653172) q[0];
sx q[0];
rz(-0.99940171) q[0];
sx q[0];
rz(3.0185643) q[0];
x q[1];
rz(1.599309) q[2];
sx q[2];
rz(-2.4270504) q[2];
sx q[2];
rz(-1.5039832) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.37398794) q[1];
sx q[1];
rz(-1.3629706) q[1];
sx q[1];
rz(0.61088224) q[1];
rz(-pi) q[2];
rz(-2.7298922) q[3];
sx q[3];
rz(-1.6414335) q[3];
sx q[3];
rz(-0.27765805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.3488397) q[2];
sx q[2];
rz(-1.9323843) q[2];
sx q[2];
rz(2.7065281) q[2];
rz(1.7815636) q[3];
sx q[3];
rz(-2.3924148) q[3];
sx q[3];
rz(2.8939261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9687987) q[0];
sx q[0];
rz(-0.89247576) q[0];
sx q[0];
rz(-2.0794179) q[0];
rz(-2.0299714) q[1];
sx q[1];
rz(-1.2373135) q[1];
sx q[1];
rz(1.4020845) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5324425) q[0];
sx q[0];
rz(-0.46450588) q[0];
sx q[0];
rz(-2.4971278) q[0];
rz(-1.8066605) q[2];
sx q[2];
rz(-0.81996041) q[2];
sx q[2];
rz(1.6279398) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.84343108) q[1];
sx q[1];
rz(-1.2097675) q[1];
sx q[1];
rz(-1.4069188) q[1];
rz(-pi) q[2];
rz(1.3659031) q[3];
sx q[3];
rz(-2.5435102) q[3];
sx q[3];
rz(0.3860592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.7193675) q[2];
sx q[2];
rz(-0.56754595) q[2];
sx q[2];
rz(-2.4613703) q[2];
rz(0.42823544) q[3];
sx q[3];
rz(-1.8869583) q[3];
sx q[3];
rz(1.7817106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.5087886) q[0];
sx q[0];
rz(-1.7254242) q[0];
sx q[0];
rz(1.592214) q[0];
rz(-2.8920065) q[1];
sx q[1];
rz(-1.9892178) q[1];
sx q[1];
rz(0.54135281) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.66051018) q[0];
sx q[0];
rz(-2.5773002) q[0];
sx q[0];
rz(-3.1206467) q[0];
rz(-pi) q[1];
rz(-2.5112553) q[2];
sx q[2];
rz(-2.6001843) q[2];
sx q[2];
rz(-0.88592096) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3440481) q[1];
sx q[1];
rz(-0.86595264) q[1];
sx q[1];
rz(-1.0182347) q[1];
x q[2];
rz(3.1369282) q[3];
sx q[3];
rz(-0.38032535) q[3];
sx q[3];
rz(2.5120146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.044518746) q[2];
sx q[2];
rz(-0.35733435) q[2];
sx q[2];
rz(-0.77735916) q[2];
rz(-2.2903531) q[3];
sx q[3];
rz(-1.0612396) q[3];
sx q[3];
rz(-1.3210993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38354307) q[0];
sx q[0];
rz(-1.8109011) q[0];
sx q[0];
rz(-3.1316277) q[0];
rz(2.126157) q[1];
sx q[1];
rz(-2.3773057) q[1];
sx q[1];
rz(-1.6962956) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9706668) q[0];
sx q[0];
rz(-2.2609841) q[0];
sx q[0];
rz(2.4404581) q[0];
rz(0.5814914) q[2];
sx q[2];
rz(-2.1610689) q[2];
sx q[2];
rz(-1.3587388) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.0547486) q[1];
sx q[1];
rz(-0.63042414) q[1];
sx q[1];
rz(0.026442095) q[1];
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
sx q[1];
rz(pi/2) q[1];
rz(-0.69892591) q[2];
sx q[2];
rz(-2.2103504) q[2];
sx q[2];
rz(-0.60738579) q[2];
rz(1.7025042) q[3];
sx q[3];
rz(-1.7581698) q[3];
sx q[3];
rz(-2.3506892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8059175) q[0];
sx q[0];
rz(-1.4842002) q[0];
sx q[0];
rz(-0.6138531) q[0];
rz(1.0461668) q[1];
sx q[1];
rz(-0.26509735) q[1];
sx q[1];
rz(-0.39224958) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.046767226) q[0];
sx q[0];
rz(-1.1724768) q[0];
sx q[0];
rz(-2.4368068) q[0];
rz(-pi) q[1];
rz(-2.1456111) q[2];
sx q[2];
rz(-1.0258342) q[2];
sx q[2];
rz(-2.0202707) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0543921) q[1];
sx q[1];
rz(-2.0093845) q[1];
sx q[1];
rz(0.88263504) q[1];
rz(-pi) q[2];
rz(1.5418105) q[3];
sx q[3];
rz(-1.7757123) q[3];
sx q[3];
rz(0.25587413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0739416) q[2];
sx q[2];
rz(-1.757688) q[2];
sx q[2];
rz(-2.5411141) q[2];
rz(-2.0742119) q[3];
sx q[3];
rz(-1.3200656) q[3];
sx q[3];
rz(1.0935121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4338715) q[0];
sx q[0];
rz(-0.43294551) q[0];
sx q[0];
rz(-1.659163) q[0];
rz(-2.1451163) q[1];
sx q[1];
rz(-1.4947718) q[1];
sx q[1];
rz(1.6047118) q[1];
rz(-1.2148576) q[2];
sx q[2];
rz(-2.3050953) q[2];
sx q[2];
rz(-0.54511025) q[2];
rz(-0.89041238) q[3];
sx q[3];
rz(-2.012841) q[3];
sx q[3];
rz(-2.8750318) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
