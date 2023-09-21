OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.3553319) q[0];
sx q[0];
rz(-3.0769899) q[0];
sx q[0];
rz(-0.021615418) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3992213) q[0];
sx q[0];
rz(-2.135015) q[0];
sx q[0];
rz(-2.7447634) q[0];
rz(-1.1515491) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(-2.575945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6380438) q[1];
sx q[1];
rz(-0.81144864) q[1];
sx q[1];
rz(-1.391295) q[1];
x q[2];
rz(1.1159775) q[3];
sx q[3];
rz(-2.4664719) q[3];
sx q[3];
rz(-1.7636553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.33064476) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(2.5374106) q[2];
rz(-1.0243246) q[3];
sx q[3];
rz(-1.1573135) q[3];
sx q[3];
rz(-2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8169096) q[0];
sx q[0];
rz(-3.1284101) q[0];
sx q[0];
rz(1.0634134) q[0];
rz(-0.88513199) q[1];
sx q[1];
rz(-1.5848031) q[1];
sx q[1];
rz(0.0016454776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7749274) q[0];
sx q[0];
rz(-1.5257376) q[0];
sx q[0];
rz(-0.010618322) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4992141) q[2];
sx q[2];
rz(-1.3628236) q[2];
sx q[2];
rz(-0.013052879) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.26047036) q[1];
sx q[1];
rz(-0.80658856) q[1];
sx q[1];
rz(1.2299728) q[1];
rz(0.9886338) q[3];
sx q[3];
rz(-2.8391264) q[3];
sx q[3];
rz(-0.31103381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.985618) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(1.0937141) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(-2.1480907) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.920632) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(2.3143342) q[0];
rz(0.0050841252) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(2.0522096) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8436962) q[0];
sx q[0];
rz(-1.9934405) q[0];
sx q[0];
rz(0.98590132) q[0];
rz(-pi) q[1];
rz(0.58615746) q[2];
sx q[2];
rz(-2.456089) q[2];
sx q[2];
rz(3.0229178) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5685801) q[1];
sx q[1];
rz(-1.9553361) q[1];
sx q[1];
rz(0.09210715) q[1];
x q[2];
rz(2.5855519) q[3];
sx q[3];
rz(-2.2696113) q[3];
sx q[3];
rz(-1.8204821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0777145) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(0.88469488) q[2];
rz(1.1832773) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(-1.6515091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(2.4147721) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(2.7817536) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.115288) q[0];
sx q[0];
rz(-1.3867154) q[0];
sx q[0];
rz(-2.227965) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4904576) q[2];
sx q[2];
rz(-1.480181) q[2];
sx q[2];
rz(3.0504984) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6105146) q[1];
sx q[1];
rz(-0.58137608) q[1];
sx q[1];
rz(0.73552144) q[1];
rz(-pi) q[2];
rz(0.8538586) q[3];
sx q[3];
rz(-2.6918908) q[3];
sx q[3];
rz(0.85292294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5741253) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(0.7652258) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-0.59195834) q[3];
sx q[3];
rz(2.9005907) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9969479) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(-1.416052) q[0];
rz(-0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(-1.0353154) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.028776289) q[0];
sx q[0];
rz(-1.1338286) q[0];
sx q[0];
rz(-0.90941888) q[0];
rz(-pi) q[1];
rz(1.5315227) q[2];
sx q[2];
rz(-1.0796667) q[2];
sx q[2];
rz(2.479535) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.1977735) q[1];
sx q[1];
rz(-0.10995956) q[1];
sx q[1];
rz(2.6206559) q[1];
rz(-pi) q[2];
rz(-1.3607929) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(-0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7845903) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(2.3727097) q[2];
rz(0.33603493) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.4679573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96200213) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(-1.996421) q[0];
rz(1.1046474) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(-0.2072269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80400318) q[0];
sx q[0];
rz(-2.1875256) q[0];
sx q[0];
rz(-1.3998652) q[0];
rz(-0.60116641) q[2];
sx q[2];
rz(-1.1784369) q[2];
sx q[2];
rz(-0.62077921) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2726278) q[1];
sx q[1];
rz(-0.3513063) q[1];
sx q[1];
rz(2.3428194) q[1];
x q[2];
rz(-1.7436142) q[3];
sx q[3];
rz(-2.380905) q[3];
sx q[3];
rz(-0.93483227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-2.2613566) q[2];
sx q[2];
rz(-0.72171372) q[2];
rz(1.2747814) q[3];
sx q[3];
rz(-1.5721679) q[3];
sx q[3];
rz(-0.26003626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-2.8969144) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(2.4095643) q[0];
rz(0.0094982068) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(-0.2917372) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4371944) q[0];
sx q[0];
rz(-1.6219553) q[0];
sx q[0];
rz(0.41914661) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8740011) q[2];
sx q[2];
rz(-1.469194) q[2];
sx q[2];
rz(0.50819699) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.26475036) q[1];
sx q[1];
rz(-1.2683588) q[1];
sx q[1];
rz(-1.7251863) q[1];
rz(2.6187906) q[3];
sx q[3];
rz(-1.8409981) q[3];
sx q[3];
rz(-0.064388007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8774524) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(2.3042802) q[2];
rz(-1.1710179) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(-0.062019197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432805) q[0];
sx q[0];
rz(-2.4536112) q[0];
sx q[0];
rz(0.6638546) q[0];
rz(-0.10617667) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(-0.95867872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230969) q[0];
sx q[0];
rz(-2.3765058) q[0];
sx q[0];
rz(-0.78406783) q[0];
x q[1];
rz(-2.0197228) q[2];
sx q[2];
rz(-2.556483) q[2];
sx q[2];
rz(1.4657071) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9740323) q[1];
sx q[1];
rz(-1.6775727) q[1];
sx q[1];
rz(2.0978931) q[1];
rz(-pi) q[2];
rz(-0.58416768) q[3];
sx q[3];
rz(-1.6730047) q[3];
sx q[3];
rz(2.971873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.133193) q[2];
sx q[2];
rz(-1.9063176) q[2];
sx q[2];
rz(0.91910249) q[2];
rz(1.7637926) q[3];
sx q[3];
rz(-1.2104687) q[3];
sx q[3];
rz(1.1834043) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(-2.7046955) q[0];
rz(-2.4412952) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(2.6760496) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31637329) q[0];
sx q[0];
rz(-1.4857978) q[0];
sx q[0];
rz(-2.4337016) q[0];
rz(0.12840694) q[2];
sx q[2];
rz(-1.1616716) q[2];
sx q[2];
rz(-2.1985334) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.1411966) q[1];
sx q[1];
rz(-1.9138412) q[1];
sx q[1];
rz(0.47980206) q[1];
rz(-2.1699379) q[3];
sx q[3];
rz(-0.42704901) q[3];
sx q[3];
rz(0.077582434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.643606) q[2];
rz(-2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(2.8695316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4964504) q[0];
sx q[0];
rz(-0.56149879) q[0];
sx q[0];
rz(1.1219332) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.5263999) q[1];
sx q[1];
rz(-2.5591992) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96888992) q[0];
sx q[0];
rz(-2.4318998) q[0];
sx q[0];
rz(1.2560647) q[0];
x q[1];
rz(1.2920612) q[2];
sx q[2];
rz(-2.6346452) q[2];
sx q[2];
rz(-2.7915733) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.460234) q[1];
sx q[1];
rz(-1.6885307) q[1];
sx q[1];
rz(-1.387499) q[1];
rz(-1.598761) q[3];
sx q[3];
rz(-1.3536605) q[3];
sx q[3];
rz(2.12487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1311243) q[2];
sx q[2];
rz(-0.98103395) q[2];
sx q[2];
rz(2.3790512) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-2.2236731) q[3];
sx q[3];
rz(1.5677174) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0653771) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(1.3394042) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(-1.273524) q[2];
sx q[2];
rz(-2.3625629) q[2];
sx q[2];
rz(0.071803781) q[2];
rz(-1.0675666) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
