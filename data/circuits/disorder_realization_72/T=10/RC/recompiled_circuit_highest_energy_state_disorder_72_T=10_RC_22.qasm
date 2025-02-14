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
rz(-1.8347972) q[0];
sx q[0];
rz(-2.2936294) q[0];
sx q[0];
rz(-3.1041978) q[0];
rz(2.1283863) q[1];
sx q[1];
rz(-2.6599045) q[1];
sx q[1];
rz(-1.4451292) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.405482) q[0];
sx q[0];
rz(-2.6713712) q[0];
sx q[0];
rz(2.7313563) q[0];
x q[1];
rz(3.0710241) q[2];
sx q[2];
rz(-1.6733992) q[2];
sx q[2];
rz(3.0377394) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2593258) q[1];
sx q[1];
rz(-2.1082467) q[1];
sx q[1];
rz(-2.5718556) q[1];
x q[2];
rz(-2.1358498) q[3];
sx q[3];
rz(-1.6673207) q[3];
sx q[3];
rz(2.9526476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41363132) q[2];
sx q[2];
rz(-3.0484338) q[2];
sx q[2];
rz(-1.5067345) q[2];
rz(2.9480751) q[3];
sx q[3];
rz(-2.3573124) q[3];
sx q[3];
rz(2.1957652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0993318) q[0];
sx q[0];
rz(-1.3657382) q[0];
sx q[0];
rz(-0.17901626) q[0];
rz(-1.5366813) q[1];
sx q[1];
rz(-0.34114006) q[1];
sx q[1];
rz(-1.590033) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1354438) q[0];
sx q[0];
rz(-1.28147) q[0];
sx q[0];
rz(-0.10162123) q[0];
x q[1];
rz(-2.8423772) q[2];
sx q[2];
rz(-2.4096074) q[2];
sx q[2];
rz(-2.3927286) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.034160567) q[1];
sx q[1];
rz(-1.9692268) q[1];
sx q[1];
rz(0.28657367) q[1];
rz(2.1158982) q[3];
sx q[3];
rz(-1.1602931) q[3];
sx q[3];
rz(0.59343597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4701074) q[2];
sx q[2];
rz(-1.6652197) q[2];
sx q[2];
rz(3.1191077) q[2];
rz(0.89618987) q[3];
sx q[3];
rz(-2.360207) q[3];
sx q[3];
rz(1.5145068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80980587) q[0];
sx q[0];
rz(-2.9349194) q[0];
sx q[0];
rz(-1.608954) q[0];
rz(-1.1398075) q[1];
sx q[1];
rz(-2.8228357) q[1];
sx q[1];
rz(-0.87361139) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5667359) q[0];
sx q[0];
rz(-2.2169161) q[0];
sx q[0];
rz(0.91888756) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.457946) q[2];
sx q[2];
rz(-0.77221671) q[2];
sx q[2];
rz(1.5600086) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1908299) q[1];
sx q[1];
rz(-2.0520981) q[1];
sx q[1];
rz(1.7442622) q[1];
rz(1.4954044) q[3];
sx q[3];
rz(-1.8429379) q[3];
sx q[3];
rz(-0.69313245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5736299) q[2];
sx q[2];
rz(-2.6831388) q[2];
sx q[2];
rz(-0.47631329) q[2];
rz(-2.1161946) q[3];
sx q[3];
rz(-1.3566596) q[3];
sx q[3];
rz(0.29388139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77182257) q[0];
sx q[0];
rz(-0.85404587) q[0];
sx q[0];
rz(-0.59637946) q[0];
rz(-2.3884933) q[1];
sx q[1];
rz(-1.3558847) q[1];
sx q[1];
rz(2.7900043) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4575949) q[0];
sx q[0];
rz(-0.54349923) q[0];
sx q[0];
rz(0.25590055) q[0];
rz(-0.48248191) q[2];
sx q[2];
rz(-2.1957779) q[2];
sx q[2];
rz(-0.04908726) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5758613) q[1];
sx q[1];
rz(-1.8003776) q[1];
sx q[1];
rz(-1.9750392) q[1];
rz(-pi) q[2];
rz(-0.63372374) q[3];
sx q[3];
rz(-1.0175287) q[3];
sx q[3];
rz(0.80851698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4190462) q[2];
sx q[2];
rz(-1.6435577) q[2];
sx q[2];
rz(-2.2124115) q[2];
rz(1.2292713) q[3];
sx q[3];
rz(-0.036673948) q[3];
sx q[3];
rz(1.2693955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96255985) q[0];
sx q[0];
rz(-2.5293009) q[0];
sx q[0];
rz(-2.6137733) q[0];
rz(-0.57012308) q[1];
sx q[1];
rz(-2.5716883) q[1];
sx q[1];
rz(-2.747587) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.67585) q[0];
sx q[0];
rz(-2.619474) q[0];
sx q[0];
rz(2.7444849) q[0];
rz(-pi) q[1];
x q[1];
rz(0.5929596) q[2];
sx q[2];
rz(-1.2622132) q[2];
sx q[2];
rz(1.9401996) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.24007639) q[1];
sx q[1];
rz(-1.4528767) q[1];
sx q[1];
rz(1.7941956) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3218003) q[3];
sx q[3];
rz(-0.91320437) q[3];
sx q[3];
rz(-0.70943225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3923308) q[2];
sx q[2];
rz(-1.2886084) q[2];
sx q[2];
rz(-1.1345081) q[2];
rz(0.025731651) q[3];
sx q[3];
rz(-2.0870356) q[3];
sx q[3];
rz(-2.499685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5662017) q[0];
sx q[0];
rz(-1.8113149) q[0];
sx q[0];
rz(2.6824685) q[0];
rz(-2.3205914) q[1];
sx q[1];
rz(-2.2586925) q[1];
sx q[1];
rz(2.6301036) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2393414) q[0];
sx q[0];
rz(-0.56846877) q[0];
sx q[0];
rz(1.0064378) q[0];
rz(-0.66368565) q[2];
sx q[2];
rz(-1.7567217) q[2];
sx q[2];
rz(2.3472361) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.59907179) q[1];
sx q[1];
rz(-2.5498646) q[1];
sx q[1];
rz(0.27246957) q[1];
rz(-1.7396163) q[3];
sx q[3];
rz(-0.69069117) q[3];
sx q[3];
rz(-0.21571479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5795827) q[2];
sx q[2];
rz(-2.0333813) q[2];
sx q[2];
rz(-2.3902334) q[2];
rz(0.56685081) q[3];
sx q[3];
rz(-1.6040498) q[3];
sx q[3];
rz(-1.3273299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.421627) q[0];
sx q[0];
rz(-0.1816853) q[0];
sx q[0];
rz(3.1340461) q[0];
rz(-2.201572) q[1];
sx q[1];
rz(-0.66497856) q[1];
sx q[1];
rz(0.13253458) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.433411) q[0];
sx q[0];
rz(-1.1601747) q[0];
sx q[0];
rz(0.93958409) q[0];
rz(-pi) q[1];
rz(1.437625) q[2];
sx q[2];
rz(-1.2107163) q[2];
sx q[2];
rz(-3.0293426) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.30250994) q[1];
sx q[1];
rz(-1.9614451) q[1];
sx q[1];
rz(-1.7724228) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2279195) q[3];
sx q[3];
rz(-1.4063121) q[3];
sx q[3];
rz(-0.13492385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.14564766) q[2];
sx q[2];
rz(-1.4952679) q[2];
sx q[2];
rz(-3.0403467) q[2];
rz(-2.4920987) q[3];
sx q[3];
rz(-2.5443304) q[3];
sx q[3];
rz(0.35972843) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7262909) q[0];
sx q[0];
rz(-1.2969718) q[0];
sx q[0];
rz(-1.2671965) q[0];
rz(1.8750809) q[1];
sx q[1];
rz(-2.9837065) q[1];
sx q[1];
rz(-0.74323851) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4871536) q[0];
sx q[0];
rz(-1.3689353) q[0];
sx q[0];
rz(2.5623723) q[0];
rz(-pi) q[1];
rz(0.27979346) q[2];
sx q[2];
rz(-0.81088561) q[2];
sx q[2];
rz(1.5992129) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.083588138) q[1];
sx q[1];
rz(-0.28540719) q[1];
sx q[1];
rz(-2.1607375) q[1];
rz(-pi) q[2];
rz(0.8193797) q[3];
sx q[3];
rz(-1.0257105) q[3];
sx q[3];
rz(-1.2631288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0172952) q[2];
sx q[2];
rz(-2.9399019) q[2];
sx q[2];
rz(-1.4264433) q[2];
rz(2.1258449) q[3];
sx q[3];
rz(-1.7044715) q[3];
sx q[3];
rz(2.3118741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.5468686) q[0];
sx q[0];
rz(-2.4877553) q[0];
sx q[0];
rz(0.015901707) q[0];
rz(-1.6390027) q[1];
sx q[1];
rz(-2.465261) q[1];
sx q[1];
rz(-2.1471088) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54399207) q[0];
sx q[0];
rz(-0.74403896) q[0];
sx q[0];
rz(-0.97980325) q[0];
rz(-pi) q[1];
rz(-2.9839462) q[2];
sx q[2];
rz(-0.38262832) q[2];
sx q[2];
rz(-1.1747109) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3164036) q[1];
sx q[1];
rz(-1.087552) q[1];
sx q[1];
rz(-0.8512599) q[1];
rz(-pi) q[2];
rz(-0.1511622) q[3];
sx q[3];
rz(-1.0619319) q[3];
sx q[3];
rz(0.14948949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9859621) q[2];
sx q[2];
rz(-1.7113643) q[2];
sx q[2];
rz(0.17781167) q[2];
rz(-1.6832247) q[3];
sx q[3];
rz(-0.42176133) q[3];
sx q[3];
rz(-1.873707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.6630702) q[0];
sx q[0];
rz(-1.8926184) q[0];
sx q[0];
rz(-1.2836237) q[0];
rz(-2.3578857) q[1];
sx q[1];
rz(-0.77373928) q[1];
sx q[1];
rz(-1.8342038) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4625774) q[0];
sx q[0];
rz(-2.0082012) q[0];
sx q[0];
rz(-1.4694197) q[0];
rz(-pi) q[1];
rz(-1.1828912) q[2];
sx q[2];
rz(-1.5426621) q[2];
sx q[2];
rz(1.6468587) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.1124005) q[1];
sx q[1];
rz(-1.629273) q[1];
sx q[1];
rz(1.7450208) q[1];
rz(-pi) q[2];
rz(2.1204167) q[3];
sx q[3];
rz(-1.9299639) q[3];
sx q[3];
rz(0.041426126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.54138294) q[2];
sx q[2];
rz(-1.7317438) q[2];
sx q[2];
rz(2.6516338) q[2];
rz(2.0269035) q[3];
sx q[3];
rz(-1.9652818) q[3];
sx q[3];
rz(0.32976845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.020323044) q[0];
sx q[0];
rz(-1.9135495) q[0];
sx q[0];
rz(1.3187153) q[0];
rz(-1.2976788) q[1];
sx q[1];
rz(-0.73328016) q[1];
sx q[1];
rz(2.7136623) q[1];
rz(-1.3135459) q[2];
sx q[2];
rz(-1.3197109) q[2];
sx q[2];
rz(-0.33725658) q[2];
rz(-0.19965996) q[3];
sx q[3];
rz(-0.87423751) q[3];
sx q[3];
rz(-2.8578491) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
