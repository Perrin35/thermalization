OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8146347) q[0];
sx q[0];
rz(-0.0184514) q[0];
sx q[0];
rz(0.17395328) q[0];
rz(-1.3375018) q[1];
sx q[1];
rz(-0.87556535) q[1];
sx q[1];
rz(-2.8101885) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50179304) q[0];
sx q[0];
rz(-1.26609) q[0];
sx q[0];
rz(1.7250389) q[0];
rz(-1.061756) q[2];
sx q[2];
rz(-0.94359829) q[2];
sx q[2];
rz(-0.043832253) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.2690588) q[1];
sx q[1];
rz(-0.89947723) q[1];
sx q[1];
rz(-1.3289245) q[1];
x q[2];
rz(-2.9865206) q[3];
sx q[3];
rz(-1.6033246) q[3];
sx q[3];
rz(-0.28964466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5212253) q[2];
sx q[2];
rz(-1.2646893) q[2];
sx q[2];
rz(-0.59086665) q[2];
rz(-1.8477731) q[3];
sx q[3];
rz(-1.7655244) q[3];
sx q[3];
rz(0.53904343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(2.6227601) q[0];
sx q[0];
rz(-2.9157214) q[0];
sx q[0];
rz(-0.87509218) q[0];
rz(-3.0488293) q[1];
sx q[1];
rz(-1.0567254) q[1];
sx q[1];
rz(-2.1327532) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0419249) q[0];
sx q[0];
rz(-2.0001051) q[0];
sx q[0];
rz(-0.03089182) q[0];
rz(-2.5549803) q[2];
sx q[2];
rz(-1.5924708) q[2];
sx q[2];
rz(2.7119111) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31263217) q[1];
sx q[1];
rz(-0.66753879) q[1];
sx q[1];
rz(1.8132736) q[1];
x q[2];
rz(1.6289234) q[3];
sx q[3];
rz(-1.584313) q[3];
sx q[3];
rz(-2.0292387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.86866093) q[2];
sx q[2];
rz(-1.2596143) q[2];
sx q[2];
rz(-1.3966365) q[2];
rz(-2.7331288) q[3];
sx q[3];
rz(-0.68578458) q[3];
sx q[3];
rz(0.43533152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.246966) q[0];
sx q[0];
rz(-2.8530687) q[0];
sx q[0];
rz(0.95046473) q[0];
rz(1.2735927) q[1];
sx q[1];
rz(-1.7968977) q[1];
sx q[1];
rz(-1.3709995) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34370658) q[0];
sx q[0];
rz(-2.4737353) q[0];
sx q[0];
rz(-2.9178828) q[0];
rz(2.7080405) q[2];
sx q[2];
rz(-1.6759652) q[2];
sx q[2];
rz(-2.1702332) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.85670722) q[1];
sx q[1];
rz(-2.6213985) q[1];
sx q[1];
rz(1.3150645) q[1];
rz(-pi) q[2];
rz(0.14601918) q[3];
sx q[3];
rz(-0.64490025) q[3];
sx q[3];
rz(-2.9420946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3889968) q[2];
sx q[2];
rz(-1.7277371) q[2];
sx q[2];
rz(1.5163806) q[2];
rz(-0.73836941) q[3];
sx q[3];
rz(-1.5117896) q[3];
sx q[3];
rz(0.71877688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3589288) q[0];
sx q[0];
rz(-1.5072701) q[0];
sx q[0];
rz(1.2231109) q[0];
rz(-0.75543985) q[1];
sx q[1];
rz(-2.0016045) q[1];
sx q[1];
rz(-3.0208407) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2837388) q[0];
sx q[0];
rz(-1.4072937) q[0];
sx q[0];
rz(0.16450926) q[0];
rz(1.3322387) q[2];
sx q[2];
rz(-1.0285796) q[2];
sx q[2];
rz(-1.2921367) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.11343458) q[1];
sx q[1];
rz(-2.096281) q[1];
sx q[1];
rz(2.3281086) q[1];
rz(1.6183245) q[3];
sx q[3];
rz(-1.589425) q[3];
sx q[3];
rz(0.19028529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.52189031) q[2];
sx q[2];
rz(-1.4554687) q[2];
sx q[2];
rz(-1.5187368) q[2];
rz(-1.3269199) q[3];
sx q[3];
rz(-0.36553317) q[3];
sx q[3];
rz(-1.1091308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13490747) q[0];
sx q[0];
rz(-1.6970072) q[0];
sx q[0];
rz(-2.4476442) q[0];
rz(0.42849439) q[1];
sx q[1];
rz(-2.4887812) q[1];
sx q[1];
rz(1.8983715) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.080414198) q[0];
sx q[0];
rz(-1.5251056) q[0];
sx q[0];
rz(-1.3752328) q[0];
rz(0.17036338) q[2];
sx q[2];
rz(-1.3479832) q[2];
sx q[2];
rz(-0.52056584) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32520884) q[1];
sx q[1];
rz(-1.3482136) q[1];
sx q[1];
rz(-1.9680262) q[1];
rz(-pi) q[2];
rz(-1.6749773) q[3];
sx q[3];
rz(-1.0910463) q[3];
sx q[3];
rz(2.097985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.99416939) q[2];
sx q[2];
rz(-1.8403515) q[2];
sx q[2];
rz(-0.48747882) q[2];
rz(-3.0726037) q[3];
sx q[3];
rz(-2.4853849) q[3];
sx q[3];
rz(-2.5330353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2212806) q[0];
sx q[0];
rz(-0.99928904) q[0];
sx q[0];
rz(-2.6699303) q[0];
rz(-2.3645875) q[1];
sx q[1];
rz(-2.0174556) q[1];
sx q[1];
rz(-1.5832925) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64028231) q[0];
sx q[0];
rz(-0.030661432) q[0];
sx q[0];
rz(2.8778381) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.110811) q[2];
sx q[2];
rz(-1.938463) q[2];
sx q[2];
rz(0.45258507) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1780381) q[1];
sx q[1];
rz(-2.6500672) q[1];
sx q[1];
rz(2.6705533) q[1];
rz(1.8040001) q[3];
sx q[3];
rz(-2.2836855) q[3];
sx q[3];
rz(1.9645129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.16288699) q[2];
sx q[2];
rz(-1.4894314) q[2];
sx q[2];
rz(1.9784652) q[2];
rz(-0.85824054) q[3];
sx q[3];
rz(-1.6909928) q[3];
sx q[3];
rz(0.8391909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6100886) q[0];
sx q[0];
rz(-1.6781582) q[0];
sx q[0];
rz(-2.8186744) q[0];
rz(-1.8220176) q[1];
sx q[1];
rz(-1.1406621) q[1];
sx q[1];
rz(-1.3986157) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8136637) q[0];
sx q[0];
rz(-1.3998056) q[0];
sx q[0];
rz(1.405262) q[0];
rz(-pi) q[1];
rz(2.1115913) q[2];
sx q[2];
rz(-0.68533603) q[2];
sx q[2];
rz(0.50619307) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4109282) q[1];
sx q[1];
rz(-0.81984869) q[1];
sx q[1];
rz(-1.5265205) q[1];
rz(-pi) q[2];
rz(0.62039726) q[3];
sx q[3];
rz(-2.982938) q[3];
sx q[3];
rz(-0.59172025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4601595) q[2];
sx q[2];
rz(-0.56542772) q[2];
sx q[2];
rz(0.73271712) q[2];
rz(-3.1198464) q[3];
sx q[3];
rz(-1.6201868) q[3];
sx q[3];
rz(2.963613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9916423) q[0];
sx q[0];
rz(-3.0198779) q[0];
sx q[0];
rz(2.6974086) q[0];
rz(-1.3581879) q[1];
sx q[1];
rz(-1.5636179) q[1];
sx q[1];
rz(-1.1202728) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2884707) q[0];
sx q[0];
rz(-0.74672304) q[0];
sx q[0];
rz(1.7225368) q[0];
rz(-pi) q[1];
rz(0.6702126) q[2];
sx q[2];
rz(-2.5465917) q[2];
sx q[2];
rz(-0.52192823) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.6669338) q[1];
sx q[1];
rz(-0.28945112) q[1];
sx q[1];
rz(0.54929026) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1577663) q[3];
sx q[3];
rz(-2.0011229) q[3];
sx q[3];
rz(2.4347507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.84919763) q[2];
sx q[2];
rz(-2.2870543) q[2];
sx q[2];
rz(1.2348068) q[2];
rz(0.97709996) q[3];
sx q[3];
rz(-1.2950725) q[3];
sx q[3];
rz(0.60404122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4567756) q[0];
sx q[0];
rz(-3.1148532) q[0];
sx q[0];
rz(-2.7944414) q[0];
rz(0.30255643) q[1];
sx q[1];
rz(-1.1056933) q[1];
sx q[1];
rz(-0.0060630719) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13841471) q[0];
sx q[0];
rz(-0.18269953) q[0];
sx q[0];
rz(1.4783786) q[0];
rz(-pi) q[1];
rz(2.8441695) q[2];
sx q[2];
rz(-2.2280577) q[2];
sx q[2];
rz(-1.707945) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2711941) q[1];
sx q[1];
rz(-2.4358303) q[1];
sx q[1];
rz(1.0949542) q[1];
rz(-2.361287) q[3];
sx q[3];
rz(-1.6041451) q[3];
sx q[3];
rz(1.0409955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.2431295) q[2];
sx q[2];
rz(-0.96776882) q[2];
sx q[2];
rz(-1.751162) q[2];
rz(0.27680382) q[3];
sx q[3];
rz(-1.0172458) q[3];
sx q[3];
rz(-1.0677342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7368363) q[0];
sx q[0];
rz(-2.3911349) q[0];
sx q[0];
rz(-0.87338895) q[0];
rz(2.5923173) q[1];
sx q[1];
rz(-2.4875689) q[1];
sx q[1];
rz(2.3347847) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72708541) q[0];
sx q[0];
rz(-0.9555409) q[0];
sx q[0];
rz(-1.8919957) q[0];
x q[1];
rz(-2.4960356) q[2];
sx q[2];
rz(-0.19536138) q[2];
sx q[2];
rz(-0.33249172) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.52376952) q[1];
sx q[1];
rz(-1.6865936) q[1];
sx q[1];
rz(1.8207184) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8397614) q[3];
sx q[3];
rz(-2.1595397) q[3];
sx q[3];
rz(-1.4422688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4775382) q[2];
sx q[2];
rz(-1.3128023) q[2];
sx q[2];
rz(-0.99267268) q[2];
rz(-0.70332518) q[3];
sx q[3];
rz(-2.0086918) q[3];
sx q[3];
rz(-2.4542184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.138388) q[0];
sx q[0];
rz(-2.3532372) q[0];
sx q[0];
rz(-1.1113356) q[0];
rz(-2.5913024) q[1];
sx q[1];
rz(-2.7688409) q[1];
sx q[1];
rz(0.54584835) q[1];
rz(-1.7334123) q[2];
sx q[2];
rz(-2.2615216) q[2];
sx q[2];
rz(0.79590454) q[2];
rz(1.1543858) q[3];
sx q[3];
rz(-1.9999123) q[3];
sx q[3];
rz(-2.6725389) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
