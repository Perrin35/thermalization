OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7862608) q[0];
sx q[0];
rz(-0.064602764) q[0];
sx q[0];
rz(0.021615418) q[0];
rz(-0.99524438) q[1];
sx q[1];
rz(-1.3270451) q[1];
sx q[1];
rz(1.8099161) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0639122) q[0];
sx q[0];
rz(-2.4644116) q[0];
sx q[0];
rz(-1.022524) q[0];
x q[1];
rz(1.3433427) q[2];
sx q[2];
rz(-1.6709575) q[2];
sx q[2];
rz(-1.4129461) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3803346) q[1];
sx q[1];
rz(-0.77612703) q[1];
sx q[1];
rz(0.18591979) q[1];
rz(-2.8033923) q[3];
sx q[3];
rz(-0.97465289) q[3];
sx q[3];
rz(-0.81830922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8109479) q[2];
sx q[2];
rz(-1.6487048) q[2];
sx q[2];
rz(-2.5374106) q[2];
rz(1.0243246) q[3];
sx q[3];
rz(-1.9842792) q[3];
sx q[3];
rz(-2.0143051) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
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
rz(-1.5848031) q[1];
sx q[1];
rz(-0.0016454776) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9379399) q[0];
sx q[0];
rz(-1.5601888) q[0];
sx q[0];
rz(-1.525735) q[0];
rz(-2.4992141) q[2];
sx q[2];
rz(-1.7787691) q[2];
sx q[2];
rz(-0.013052879) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.590608) q[1];
sx q[1];
rz(-1.3270757) q[1];
sx q[1];
rz(-2.3477712) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1529589) q[3];
sx q[3];
rz(-2.8391264) q[3];
sx q[3];
rz(0.31103381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1559747) q[2];
sx q[2];
rz(-1.5094455) q[2];
sx q[2];
rz(-2.4334811) q[2];
rz(2.0478785) q[3];
sx q[3];
rz(-1.3392859) q[3];
sx q[3];
rz(2.1480907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.920632) q[0];
sx q[0];
rz(-1.7376124) q[0];
sx q[0];
rz(-0.8272585) q[0];
rz(-3.1365085) q[1];
sx q[1];
rz(-1.9283483) q[1];
sx q[1];
rz(1.089383) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8436962) q[0];
sx q[0];
rz(-1.9934405) q[0];
sx q[0];
rz(-0.98590132) q[0];
rz(-pi) q[1];
rz(1.9956279) q[2];
sx q[2];
rz(-2.1263188) q[2];
sx q[2];
rz(-0.59031634) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8100064) q[1];
sx q[1];
rz(-0.3948822) q[1];
sx q[1];
rz(-1.3473131) q[1];
rz(-pi) q[2];
x q[2];
rz(0.55604071) q[3];
sx q[3];
rz(-0.87198139) q[3];
sx q[3];
rz(-1.8204821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0638782) q[2];
sx q[2];
rz(-0.87746799) q[2];
sx q[2];
rz(-0.88469488) q[2];
rz(-1.9583154) q[3];
sx q[3];
rz(-1.8235455) q[3];
sx q[3];
rz(1.4900835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5692212) q[0];
sx q[0];
rz(-2.5771993) q[0];
sx q[0];
rz(-2.4147721) q[0];
rz(2.3379393) q[1];
sx q[1];
rz(-1.0875965) q[1];
sx q[1];
rz(-0.35983905) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7774178) q[0];
sx q[0];
rz(-2.4628277) q[0];
sx q[0];
rz(-1.8666408) q[0];
rz(-pi) q[1];
x q[1];
rz(0.14881046) q[2];
sx q[2];
rz(-2.4850922) q[2];
sx q[2];
rz(1.3615001) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.53107809) q[1];
sx q[1];
rz(-2.5602166) q[1];
sx q[1];
rz(2.4060712) q[1];
rz(-pi) q[2];
rz(1.2218277) q[3];
sx q[3];
rz(-1.8604606) q[3];
sx q[3];
rz(1.7581913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.56746733) q[2];
sx q[2];
rz(-1.1541157) q[2];
sx q[2];
rz(-2.3763669) q[2];
rz(0.75677538) q[3];
sx q[3];
rz(-2.5496343) q[3];
sx q[3];
rz(0.24100196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9969479) q[0];
sx q[0];
rz(-2.6648271) q[0];
sx q[0];
rz(1.416052) q[0];
rz(2.7744746) q[1];
sx q[1];
rz(-1.7558302) q[1];
sx q[1];
rz(1.0353154) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1015548) q[0];
sx q[0];
rz(-2.3674175) q[0];
sx q[0];
rz(0.9206307) q[0];
rz(2.6501422) q[2];
sx q[2];
rz(-1.6054258) q[2];
sx q[2];
rz(-2.2143242) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.6742179) q[1];
sx q[1];
rz(-1.6661223) q[1];
sx q[1];
rz(1.6256888) q[1];
x q[2];
rz(-1.3607929) q[3];
sx q[3];
rz(-1.8404507) q[3];
sx q[3];
rz(-0.3013914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7845903) q[2];
sx q[2];
rz(-1.720153) q[2];
sx q[2];
rz(0.76888293) q[2];
rz(-2.8055577) q[3];
sx q[3];
rz(-0.78032812) q[3];
sx q[3];
rz(-1.6736354) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96200213) q[0];
sx q[0];
rz(-0.63350326) q[0];
sx q[0];
rz(1.1451716) q[0];
rz(-1.1046474) q[1];
sx q[1];
rz(-1.2303753) q[1];
sx q[1];
rz(2.9343658) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4743054) q[0];
sx q[0];
rz(-1.7100088) q[0];
sx q[0];
rz(2.5179203) q[0];
x q[1];
rz(1.1057165) q[2];
sx q[2];
rz(-1.0208703) q[2];
sx q[2];
rz(1.9351026) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2726278) q[1];
sx q[1];
rz(-0.3513063) q[1];
sx q[1];
rz(2.3428194) q[1];
rz(-pi) q[2];
rz(-2.3239922) q[3];
sx q[3];
rz(-1.6896276) q[3];
sx q[3];
rz(0.51018754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3383639) q[2];
sx q[2];
rz(-0.88023606) q[2];
sx q[2];
rz(-2.4198789) q[2];
rz(-1.8668113) q[3];
sx q[3];
rz(-1.5694247) q[3];
sx q[3];
rz(-2.8815564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8969144) q[0];
sx q[0];
rz(-1.5667863) q[0];
sx q[0];
rz(2.4095643) q[0];
rz(3.1320944) q[1];
sx q[1];
rz(-0.54832012) q[1];
sx q[1];
rz(2.8498555) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1222629) q[0];
sx q[0];
rz(-0.42207345) q[0];
sx q[0];
rz(3.0164369) q[0];
rz(-pi) q[1];
rz(1.2675915) q[2];
sx q[2];
rz(-1.469194) q[2];
sx q[2];
rz(-0.50819699) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.74624324) q[1];
sx q[1];
rz(-0.33848539) q[1];
sx q[1];
rz(-2.6836718) q[1];
rz(-pi) q[2];
rz(-0.50645701) q[3];
sx q[3];
rz(-0.582687) q[3];
sx q[3];
rz(-1.9401693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26414028) q[2];
sx q[2];
rz(-1.143127) q[2];
sx q[2];
rz(2.3042802) q[2];
rz(-1.9705747) q[3];
sx q[3];
rz(-1.6224909) q[3];
sx q[3];
rz(3.0795735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0432805) q[0];
sx q[0];
rz(-0.68798143) q[0];
sx q[0];
rz(-0.6638546) q[0];
rz(3.035416) q[1];
sx q[1];
rz(-2.5352434) q[1];
sx q[1];
rz(0.95867872) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230969) q[0];
sx q[0];
rz(-2.3765058) q[0];
sx q[0];
rz(2.3575248) q[0];
rz(-1.1218698) q[2];
sx q[2];
rz(-2.556483) q[2];
sx q[2];
rz(-1.4657071) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.6764076) q[1];
sx q[1];
rz(-1.0470111) q[1];
sx q[1];
rz(-0.1233867) q[1];
rz(-pi) q[2];
rz(-2.9577191) q[3];
sx q[3];
rz(-2.549578) q[3];
sx q[3];
rz(1.5873991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.133193) q[2];
sx q[2];
rz(-1.2352751) q[2];
sx q[2];
rz(2.2224902) q[2];
rz(-1.7637926) q[3];
sx q[3];
rz(-1.931124) q[3];
sx q[3];
rz(-1.9581883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8063426) q[0];
sx q[0];
rz(-0.85334539) q[0];
sx q[0];
rz(-2.7046955) q[0];
rz(0.70029744) q[1];
sx q[1];
rz(-1.4236139) q[1];
sx q[1];
rz(2.6760496) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8252194) q[0];
sx q[0];
rz(-1.4857978) q[0];
sx q[0];
rz(0.70789106) q[0];
x q[1];
rz(1.1586458) q[2];
sx q[2];
rz(-1.4530384) q[2];
sx q[2];
rz(2.462537) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.8852639) q[1];
sx q[1];
rz(-2.0205106) q[1];
sx q[1];
rz(1.9535669) q[1];
x q[2];
rz(-2.1699379) q[3];
sx q[3];
rz(-0.42704901) q[3];
sx q[3];
rz(-3.0640102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.41032252) q[2];
sx q[2];
rz(-2.6022544) q[2];
sx q[2];
rz(-1.643606) q[2];
rz(0.20478976) q[3];
sx q[3];
rz(-2.1361165) q[3];
sx q[3];
rz(0.27206102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64514226) q[0];
sx q[0];
rz(-2.5800939) q[0];
sx q[0];
rz(-2.0196594) q[0];
rz(-2.3902068) q[1];
sx q[1];
rz(-1.6151927) q[1];
sx q[1];
rz(-2.5591992) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1727027) q[0];
sx q[0];
rz(-0.70969289) q[0];
sx q[0];
rz(-1.2560647) q[0];
rz(1.0803797) q[2];
sx q[2];
rz(-1.4368125) q[2];
sx q[2];
rz(1.4659363) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6813587) q[1];
sx q[1];
rz(-1.6885307) q[1];
sx q[1];
rz(1.7540936) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.12606975) q[3];
sx q[3];
rz(-0.21890103) q[3];
sx q[3];
rz(-2.253988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0104684) q[2];
sx q[2];
rz(-2.1605587) q[2];
sx q[2];
rz(2.3790512) q[2];
rz(1.4108346) q[3];
sx q[3];
rz(-0.91791955) q[3];
sx q[3];
rz(-1.5738752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0762155) q[0];
sx q[0];
rz(-0.98631728) q[0];
sx q[0];
rz(1.4022723) q[0];
rz(-1.3394042) q[1];
sx q[1];
rz(-1.6911472) q[1];
sx q[1];
rz(-1.4858248) q[1];
rz(-1.8680686) q[2];
sx q[2];
rz(-0.77902972) q[2];
sx q[2];
rz(-3.0697889) q[2];
rz(2.6639832) q[3];
sx q[3];
rz(-2.0255247) q[3];
sx q[3];
rz(-0.075102641) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
