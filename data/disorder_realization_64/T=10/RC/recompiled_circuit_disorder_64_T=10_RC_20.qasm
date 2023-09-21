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
rz(-1.3316766) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0488779) q[0];
sx q[0];
rz(-1.2381011) q[0];
sx q[0];
rz(-0.96941745) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9900436) q[2];
sx q[2];
rz(-0.24818072) q[2];
sx q[2];
rz(0.56564769) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6380438) q[1];
sx q[1];
rz(-2.330144) q[1];
sx q[1];
rz(1.7502977) q[1];
rz(-pi) q[2];
x q[2];
rz(0.33820037) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(2.3232834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8109479) q[2];
sx q[2];
rz(-1.4928879) q[2];
sx q[2];
rz(-2.5374106) q[2];
rz(-2.1172681) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(-2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8169096) q[0];
sx q[0];
rz(-0.01318251) q[0];
sx q[0];
rz(2.0781793) q[0];
rz(-2.2564607) q[1];
sx q[1];
rz(-1.5567895) q[1];
sx q[1];
rz(0.0016454776) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3666653) q[0];
sx q[0];
rz(-1.6158551) q[0];
sx q[0];
rz(3.1309743) q[0];
x q[1];
rz(-2.4992141) q[2];
sx q[2];
rz(-1.7787691) q[2];
sx q[2];
rz(-0.013052879) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.26047036) q[1];
sx q[1];
rz(-2.3350041) q[1];
sx q[1];
rz(-1.2299728) q[1];
x q[2];
rz(-2.9716773) q[3];
sx q[3];
rz(-1.3193469) q[3];
sx q[3];
rz(-2.2268695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.985618) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(2.0478785) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(-0.99350199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.920632) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(0.0050841252) q[1];
sx q[1];
rz(-1.2132443) q[1];
sx q[1];
rz(2.0522096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1339061) q[0];
sx q[0];
rz(-2.0984762) q[0];
sx q[0];
rz(0.49467996) q[0];
x q[1];
rz(2.543534) q[2];
sx q[2];
rz(-1.9285678) q[2];
sx q[2];
rz(-1.9269112) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8100064) q[1];
sx q[1];
rz(-0.3948822) q[1];
sx q[1];
rz(1.7942795) q[1];
x q[2];
rz(0.55604071) q[3];
sx q[3];
rz(-0.87198139) q[3];
sx q[3];
rz(1.3211105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0638782) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(2.2568978) q[2];
rz(1.1832773) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5692212) q[0];
sx q[0];
rz(-0.5643934) q[0];
sx q[0];
rz(-2.4147721) q[0];
rz(-0.80365333) q[1];
sx q[1];
rz(-2.0539961) q[1];
sx q[1];
rz(0.35983905) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3641748) q[0];
sx q[0];
rz(-0.67876498) q[0];
sx q[0];
rz(-1.8666408) q[0];
x q[1];
rz(-1.4570518) q[2];
sx q[2];
rz(-2.2188088) q[2];
sx q[2];
rz(1.5930454) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.53107809) q[1];
sx q[1];
rz(-2.5602166) q[1];
sx q[1];
rz(2.4060712) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8538586) q[3];
sx q[3];
rz(-0.44970185) q[3];
sx q[3];
rz(0.85292294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.56746733) q[2];
sx q[2];
rz(-1.987477) q[2];
sx q[2];
rz(-0.7652258) q[2];
rz(-2.3848173) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(-2.9005907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1446447) q[0];
sx q[0];
rz(-0.47676555) q[0];
sx q[0];
rz(1.416052) q[0];
rz(-0.36711806) q[1];
sx q[1];
rz(-1.3857625) q[1];
sx q[1];
rz(2.1062772) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1015548) q[0];
sx q[0];
rz(-2.3674175) q[0];
sx q[0];
rz(2.220962) q[0];
x q[1];
rz(3.0683124) q[2];
sx q[2];
rz(-2.6490232) q[2];
sx q[2];
rz(2.5626593) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.94381911) q[1];
sx q[1];
rz(-3.0316331) q[1];
sx q[1];
rz(2.6206559) q[1];
rz(-pi) q[2];
rz(1.3607929) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.7845903) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(0.76888293) q[2];
rz(2.8055577) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(1.6736354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1795905) q[0];
sx q[0];
rz(-2.5080894) q[0];
sx q[0];
rz(-1.996421) q[0];
rz(-2.0369453) q[1];
sx q[1];
rz(-1.9112174) q[1];
sx q[1];
rz(2.9343658) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.094039) q[0];
sx q[0];
rz(-0.63699603) q[0];
sx q[0];
rz(-2.9061222) q[0];
rz(-0.63164288) q[2];
sx q[2];
rz(-2.4372059) q[2];
sx q[2];
rz(2.7001675) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6991899) q[1];
sx q[1];
rz(-1.3283722) q[1];
sx q[1];
rz(-1.8276023) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3239922) q[3];
sx q[3];
rz(-1.451965) q[3];
sx q[3];
rz(0.51018754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3383639) q[2];
sx q[2];
rz(-0.88023606) q[2];
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
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24467829) q[0];
sx q[0];
rz(-1.5748064) q[0];
sx q[0];
rz(-2.4095643) q[0];
rz(-3.1320944) q[1];
sx q[1];
rz(-2.5932725) q[1];
sx q[1];
rz(-0.2917372) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7043982) q[0];
sx q[0];
rz(-1.5196374) q[0];
sx q[0];
rz(-2.722446) q[0];
rz(-1.8998434) q[2];
sx q[2];
rz(-2.8223158) q[2];
sx q[2];
rz(-0.74908756) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8818672) q[1];
sx q[1];
rz(-1.7181267) q[1];
sx q[1];
rz(-2.8357361) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6187906) q[3];
sx q[3];
rz(-1.8409981) q[3];
sx q[3];
rz(0.064388007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8774524) q[2];
sx q[2];
rz(-1.9984657) q[2];
sx q[2];
rz(-2.3042802) q[2];
rz(-1.9705747) q[3];
sx q[3];
rz(-1.5191017) q[3];
sx q[3];
rz(-3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0983122) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(-2.4777381) q[0];
rz(3.035416) q[1];
sx q[1];
rz(-0.60634923) q[1];
sx q[1];
rz(2.1829139) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52866919) q[0];
sx q[0];
rz(-2.0818424) q[0];
sx q[0];
rz(-0.59707609) q[0];
rz(1.1218698) q[2];
sx q[2];
rz(-0.58510963) q[2];
sx q[2];
rz(1.6758855) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.465185) q[1];
sx q[1];
rz(-2.0945815) q[1];
sx q[1];
rz(-0.1233867) q[1];
x q[2];
rz(-0.18387353) q[3];
sx q[3];
rz(-0.5920147) q[3];
sx q[3];
rz(1.5873991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0083996) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(-2.2224902) q[2];
rz(1.3778) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(-1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.33525) q[0];
sx q[0];
rz(-2.2882473) q[0];
sx q[0];
rz(-0.43689716) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.7179787) q[1];
sx q[1];
rz(-2.6760496) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1818905) q[0];
sx q[0];
rz(-2.2756016) q[0];
sx q[0];
rz(-1.4591135) q[0];
x q[1];
rz(-1.1586458) q[2];
sx q[2];
rz(-1.6885542) q[2];
sx q[2];
rz(2.462537) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.25632875) q[1];
sx q[1];
rz(-1.1210821) q[1];
sx q[1];
rz(1.9535669) q[1];
x q[2];
rz(-0.97165473) q[3];
sx q[3];
rz(-2.7145436) q[3];
sx q[3];
rz(0.077582434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(1.4979866) q[2];
rz(-2.9368029) q[3];
sx q[3];
rz(-1.0054761) q[3];
sx q[3];
rz(-0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4964504) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(-2.0196594) q[0];
rz(2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(-0.5823935) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1727027) q[0];
sx q[0];
rz(-2.4318998) q[0];
sx q[0];
rz(-1.2560647) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15162823) q[2];
sx q[2];
rz(-2.0564338) q[2];
sx q[2];
rz(0.033657311) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0092587) q[1];
sx q[1];
rz(-1.7528105) q[1];
sx q[1];
rz(0.11972129) q[1];
x q[2];
rz(0.12606975) q[3];
sx q[3];
rz(-0.21890103) q[3];
sx q[3];
rz(2.253988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1311243) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(-0.76254145) q[2];
rz(1.7307581) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0762155) q[0];
sx q[0];
rz(-2.1552754) q[0];
sx q[0];
rz(-1.7393204) q[0];
rz(1.3394042) q[1];
sx q[1];
rz(-1.4504455) q[1];
sx q[1];
rz(1.6557678) q[1];
rz(0.81417685) q[2];
sx q[2];
rz(-1.3635175) q[2];
sx q[2];
rz(-1.2843532) q[2];
rz(2.0740261) q[3];
sx q[3];
rz(-1.9964841) q[3];
sx q[3];
rz(-1.869429) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];