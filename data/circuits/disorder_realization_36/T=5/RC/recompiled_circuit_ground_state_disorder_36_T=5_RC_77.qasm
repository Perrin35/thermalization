OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1898243) q[0];
sx q[0];
rz(-0.70280743) q[0];
sx q[0];
rz(-1.9195317) q[0];
rz(3.089978) q[1];
sx q[1];
rz(-1.6366704) q[1];
sx q[1];
rz(1.6610891) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5334789) q[0];
sx q[0];
rz(-2.2256269) q[0];
sx q[0];
rz(-2.9346908) q[0];
rz(-1.0009345) q[2];
sx q[2];
rz(-2.9409932) q[2];
sx q[2];
rz(-2.4127485) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6340629) q[1];
sx q[1];
rz(-0.46137091) q[1];
sx q[1];
rz(2.7410236) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0546255) q[3];
sx q[3];
rz(-1.9619521) q[3];
sx q[3];
rz(0.018875518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9606955) q[2];
sx q[2];
rz(-0.81520671) q[2];
sx q[2];
rz(1.0308456) q[2];
rz(2.893462) q[3];
sx q[3];
rz(-1.7957567) q[3];
sx q[3];
rz(0.26721755) q[3];
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
rz(0.88028518) q[0];
sx q[0];
rz(-1.1957276) q[0];
sx q[0];
rz(-2.0174568) q[0];
rz(2.0696438) q[1];
sx q[1];
rz(-1.5644466) q[1];
sx q[1];
rz(-2.7148278) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2248707) q[0];
sx q[0];
rz(-1.568075) q[0];
sx q[0];
rz(0.82758521) q[0];
rz(0.34554225) q[2];
sx q[2];
rz(-1.8176518) q[2];
sx q[2];
rz(-3.0902362) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76600915) q[1];
sx q[1];
rz(-1.462862) q[1];
sx q[1];
rz(-0.69891478) q[1];
rz(-pi) q[2];
rz(1.9598448) q[3];
sx q[3];
rz(-1.3836244) q[3];
sx q[3];
rz(-2.1120918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.65853226) q[2];
sx q[2];
rz(-2.6250562) q[2];
sx q[2];
rz(-0.10022441) q[2];
rz(1.077486) q[3];
sx q[3];
rz(-1.5651549) q[3];
sx q[3];
rz(1.235435) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.539262) q[0];
sx q[0];
rz(-2.1644008) q[0];
sx q[0];
rz(-1.7072898) q[0];
rz(2.3151248) q[1];
sx q[1];
rz(-1.4974599) q[1];
sx q[1];
rz(-0.43063393) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9481107) q[0];
sx q[0];
rz(-1.8405582) q[0];
sx q[0];
rz(1.4007691) q[0];
x q[1];
rz(-1.1330092) q[2];
sx q[2];
rz(-2.8418645) q[2];
sx q[2];
rz(-2.9309811) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.34738008) q[1];
sx q[1];
rz(-2.2092704) q[1];
sx q[1];
rz(0.49974167) q[1];
x q[2];
rz(0.39726302) q[3];
sx q[3];
rz(-1.0618342) q[3];
sx q[3];
rz(-2.404117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.87381252) q[2];
sx q[2];
rz(-1.9766821) q[2];
sx q[2];
rz(0.57360345) q[2];
rz(-2.7576647) q[3];
sx q[3];
rz(-1.826518) q[3];
sx q[3];
rz(-2.4887776) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7225994) q[0];
sx q[0];
rz(-2.3498131) q[0];
sx q[0];
rz(2.0844039) q[0];
rz(0.81622299) q[1];
sx q[1];
rz(-2.1482601) q[1];
sx q[1];
rz(-0.23424558) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1907983) q[0];
sx q[0];
rz(-1.4613348) q[0];
sx q[0];
rz(1.9881161) q[0];
rz(0.90936898) q[2];
sx q[2];
rz(-0.63440914) q[2];
sx q[2];
rz(2.6480791) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.034742753) q[1];
sx q[1];
rz(-2.2487246) q[1];
sx q[1];
rz(-2.2449916) q[1];
x q[2];
rz(-2.0755505) q[3];
sx q[3];
rz(-1.1467993) q[3];
sx q[3];
rz(-3.0847065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8546042) q[2];
sx q[2];
rz(-1.8082666) q[2];
sx q[2];
rz(-2.2387779) q[2];
rz(1.3826987) q[3];
sx q[3];
rz(-2.8185676) q[3];
sx q[3];
rz(-0.84669101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6803902) q[0];
sx q[0];
rz(-1.8620551) q[0];
sx q[0];
rz(2.4317106) q[0];
rz(-2.6704125) q[1];
sx q[1];
rz(-1.9147562) q[1];
sx q[1];
rz(3.0772298) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.956139) q[0];
sx q[0];
rz(-0.50489932) q[0];
sx q[0];
rz(1.288815) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9856312) q[2];
sx q[2];
rz(-2.6081787) q[2];
sx q[2];
rz(-0.87863011) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2379265) q[1];
sx q[1];
rz(-1.1519379) q[1];
sx q[1];
rz(-2.5861857) q[1];
x q[2];
rz(0.086802089) q[3];
sx q[3];
rz(-1.5613114) q[3];
sx q[3];
rz(-2.529151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.71089661) q[2];
sx q[2];
rz(-0.5527834) q[2];
sx q[2];
rz(0.95332471) q[2];
rz(0.55771762) q[3];
sx q[3];
rz(-1.6744813) q[3];
sx q[3];
rz(2.6127156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8271269) q[0];
sx q[0];
rz(-2.0966661) q[0];
sx q[0];
rz(-0.99660981) q[0];
rz(-3.0820471) q[1];
sx q[1];
rz(-0.98809067) q[1];
sx q[1];
rz(-1.7893808) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4137581) q[0];
sx q[0];
rz(-1.0737077) q[0];
sx q[0];
rz(0.97151206) q[0];
rz(0.77162051) q[2];
sx q[2];
rz(-1.811027) q[2];
sx q[2];
rz(0.26773237) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.775573) q[1];
sx q[1];
rz(-1.8692353) q[1];
sx q[1];
rz(1.2828903) q[1];
rz(-2.7963403) q[3];
sx q[3];
rz(-0.94178978) q[3];
sx q[3];
rz(-2.0382413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.754564) q[2];
sx q[2];
rz(-1.8661541) q[2];
sx q[2];
rz(-2.582029) q[2];
rz(-0.83545056) q[3];
sx q[3];
rz(-0.42568922) q[3];
sx q[3];
rz(-1.3028418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68592042) q[0];
sx q[0];
rz(-2.4113825) q[0];
sx q[0];
rz(-2.7749104) q[0];
rz(0.8815676) q[1];
sx q[1];
rz(-2.5036948) q[1];
sx q[1];
rz(1.4311904) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0019497) q[0];
sx q[0];
rz(-0.77660034) q[0];
sx q[0];
rz(0.5712255) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7273799) q[2];
sx q[2];
rz(-0.72986947) q[2];
sx q[2];
rz(-0.038686631) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8578355) q[1];
sx q[1];
rz(-1.5997442) q[1];
sx q[1];
rz(3.0250103) q[1];
rz(0.65815417) q[3];
sx q[3];
rz(-1.7580716) q[3];
sx q[3];
rz(-1.3211471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7876579) q[2];
sx q[2];
rz(-0.5160318) q[2];
sx q[2];
rz(0.77989522) q[2];
rz(0.53906131) q[3];
sx q[3];
rz(-1.5122248) q[3];
sx q[3];
rz(-0.61162925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1270545) q[0];
sx q[0];
rz(-0.67033613) q[0];
sx q[0];
rz(-0.77914733) q[0];
rz(0.51891333) q[1];
sx q[1];
rz(-1.8592368) q[1];
sx q[1];
rz(-2.7386477) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76501827) q[0];
sx q[0];
rz(-2.2857762) q[0];
sx q[0];
rz(-0.9299703) q[0];
x q[1];
rz(1.0874426) q[2];
sx q[2];
rz(-1.8515808) q[2];
sx q[2];
rz(-1.7710242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.0190474) q[1];
sx q[1];
rz(-1.9662939) q[1];
sx q[1];
rz(0.67087733) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9698079) q[3];
sx q[3];
rz(-0.61099377) q[3];
sx q[3];
rz(-2.4965167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.11401033) q[2];
sx q[2];
rz(-1.1677914) q[2];
sx q[2];
rz(-2.0243417) q[2];
rz(-2.4452325) q[3];
sx q[3];
rz(-2.0317234) q[3];
sx q[3];
rz(1.2874359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47525147) q[0];
sx q[0];
rz(-2.9069558) q[0];
sx q[0];
rz(-2.7269205) q[0];
rz(-2.0798202) q[1];
sx q[1];
rz(-0.88645187) q[1];
sx q[1];
rz(0.83962238) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1621203) q[0];
sx q[0];
rz(-1.1651255) q[0];
sx q[0];
rz(1.0825054) q[0];
x q[1];
rz(2.1346852) q[2];
sx q[2];
rz(-2.208935) q[2];
sx q[2];
rz(-1.2096805) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.4794107) q[1];
sx q[1];
rz(-2.4644797) q[1];
sx q[1];
rz(0.22775316) q[1];
rz(-pi) q[2];
rz(-0.35653161) q[3];
sx q[3];
rz(-1.339669) q[3];
sx q[3];
rz(0.69510776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7825369) q[2];
sx q[2];
rz(-1.3691207) q[2];
sx q[2];
rz(-0.071652023) q[2];
rz(3.0960633) q[3];
sx q[3];
rz(-1.9436911) q[3];
sx q[3];
rz(-0.45904747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5422106) q[0];
sx q[0];
rz(-0.57294232) q[0];
sx q[0];
rz(1.0546767) q[0];
rz(-0.032546267) q[1];
sx q[1];
rz(-0.97144214) q[1];
sx q[1];
rz(2.6194825) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.78799) q[0];
sx q[0];
rz(-2.4756469) q[0];
sx q[0];
rz(0.5086201) q[0];
rz(-pi) q[1];
rz(-0.63154659) q[2];
sx q[2];
rz(-0.8435404) q[2];
sx q[2];
rz(2.2377559) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5473816) q[1];
sx q[1];
rz(-1.949233) q[1];
sx q[1];
rz(-0.40706472) q[1];
rz(-pi) q[2];
rz(2.9154112) q[3];
sx q[3];
rz(-1.0345248) q[3];
sx q[3];
rz(2.8482343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9944428) q[2];
sx q[2];
rz(-2.9394737) q[2];
sx q[2];
rz(-1.408255) q[2];
rz(-1.5857006) q[3];
sx q[3];
rz(-0.2318016) q[3];
sx q[3];
rz(-2.2304992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0178575) q[0];
sx q[0];
rz(-1.3029079) q[0];
sx q[0];
rz(-2.2444176) q[0];
rz(-0.97902117) q[1];
sx q[1];
rz(-1.9985825) q[1];
sx q[1];
rz(-0.40252007) q[1];
rz(-1.903309) q[2];
sx q[2];
rz(-1.3855743) q[2];
sx q[2];
rz(-0.34045548) q[2];
rz(-0.66313638) q[3];
sx q[3];
rz(-1.4226154) q[3];
sx q[3];
rz(-2.0285574) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
