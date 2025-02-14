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
rz(-1.4827363) q[0];
sx q[0];
rz(4.1229376) q[0];
sx q[0];
rz(11.468588) q[0];
rz(-0.78805796) q[1];
sx q[1];
rz(-2.0907953) q[1];
sx q[1];
rz(3.1346336) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6774782) q[0];
sx q[0];
rz(-0.65271806) q[0];
sx q[0];
rz(1.663289) q[0];
rz(1.8780776) q[2];
sx q[2];
rz(-0.6699962) q[2];
sx q[2];
rz(0.43596632) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3462249) q[1];
sx q[1];
rz(-2.9034894) q[1];
sx q[1];
rz(1.7134929) q[1];
x q[2];
rz(2.8473749) q[3];
sx q[3];
rz(-0.3539043) q[3];
sx q[3];
rz(-0.82415165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9030582) q[2];
sx q[2];
rz(-1.7944585) q[2];
sx q[2];
rz(1.4702338) q[2];
rz(-3.0155731) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(-2.457705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11494342) q[0];
sx q[0];
rz(-2.6925955) q[0];
sx q[0];
rz(2.5057416) q[0];
rz(-2.3682829) q[1];
sx q[1];
rz(-0.67716235) q[1];
sx q[1];
rz(1.4453452) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4659368) q[0];
sx q[0];
rz(-1.1268864) q[0];
sx q[0];
rz(2.2919284) q[0];
x q[1];
rz(-2.6670806) q[2];
sx q[2];
rz(-2.2998428) q[2];
sx q[2];
rz(3.0345033) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4559926) q[1];
sx q[1];
rz(-1.7280792) q[1];
sx q[1];
rz(-2.434751) q[1];
rz(1.0926756) q[3];
sx q[3];
rz(-2.993686) q[3];
sx q[3];
rz(-2.2151102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.14629743) q[2];
sx q[2];
rz(-1.3714906) q[2];
sx q[2];
rz(-3.0038317) q[2];
rz(-0.45608258) q[3];
sx q[3];
rz(-0.59499732) q[3];
sx q[3];
rz(-0.47541398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.544203) q[0];
sx q[0];
rz(-1.5578288) q[0];
sx q[0];
rz(2.5624045) q[0];
rz(3.0384565) q[1];
sx q[1];
rz(-1.3214279) q[1];
sx q[1];
rz(-0.78027049) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0294757) q[0];
sx q[0];
rz(-1.5519841) q[0];
sx q[0];
rz(3.0464059) q[0];
rz(-pi) q[1];
rz(0.19635503) q[2];
sx q[2];
rz(-2.4508424) q[2];
sx q[2];
rz(0.029411246) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5835598) q[1];
sx q[1];
rz(-1.3509024) q[1];
sx q[1];
rz(1.1660485) q[1];
rz(-pi) q[2];
rz(0.22091945) q[3];
sx q[3];
rz(-0.88730371) q[3];
sx q[3];
rz(-1.4190159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4812193) q[2];
sx q[2];
rz(-1.4619091) q[2];
sx q[2];
rz(-0.090864651) q[2];
rz(-1.6615435) q[3];
sx q[3];
rz(-1.816498) q[3];
sx q[3];
rz(-0.81502325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12260967) q[0];
sx q[0];
rz(-0.18773395) q[0];
sx q[0];
rz(0.51938272) q[0];
rz(1.9620365) q[1];
sx q[1];
rz(-0.68125454) q[1];
sx q[1];
rz(3.093241) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3607442) q[0];
sx q[0];
rz(-2.0791302) q[0];
sx q[0];
rz(-0.8698277) q[0];
rz(-pi) q[1];
rz(0.3860571) q[2];
sx q[2];
rz(-1.1163082) q[2];
sx q[2];
rz(0.41815652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.0569339) q[1];
sx q[1];
rz(-2.6531583) q[1];
sx q[1];
rz(0.059367511) q[1];
rz(1.8353668) q[3];
sx q[3];
rz(-2.3481784) q[3];
sx q[3];
rz(-0.95246623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4590596) q[2];
sx q[2];
rz(-2.8352663) q[2];
sx q[2];
rz(1.4502067) q[2];
rz(1.1059443) q[3];
sx q[3];
rz(-1.9927497) q[3];
sx q[3];
rz(2.2753184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3806216) q[0];
sx q[0];
rz(-2.3411317) q[0];
sx q[0];
rz(-2.3642484) q[0];
rz(0.49579534) q[1];
sx q[1];
rz(-2.7340041) q[1];
sx q[1];
rz(-0.19283238) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0986493) q[0];
sx q[0];
rz(-1.0355562) q[0];
sx q[0];
rz(-1.213206) q[0];
rz(-pi) q[1];
rz(1.8106789) q[2];
sx q[2];
rz(-0.81776202) q[2];
sx q[2];
rz(-2.8155934) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1950394) q[1];
sx q[1];
rz(-1.6583539) q[1];
sx q[1];
rz(1.273543) q[1];
rz(-pi) q[2];
rz(-1.939658) q[3];
sx q[3];
rz(-0.4274803) q[3];
sx q[3];
rz(2.1457246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.58482802) q[2];
sx q[2];
rz(-2.3636221) q[2];
sx q[2];
rz(-0.83621109) q[2];
rz(2.1861475) q[3];
sx q[3];
rz(-1.7183869) q[3];
sx q[3];
rz(-2.6098765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8504976) q[0];
sx q[0];
rz(-1.9170772) q[0];
sx q[0];
rz(3.1078872) q[0];
rz(-0.91824245) q[1];
sx q[1];
rz(-2.4372209) q[1];
sx q[1];
rz(-2.4423626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3501773) q[0];
sx q[0];
rz(-0.79424131) q[0];
sx q[0];
rz(0.82483952) q[0];
rz(-pi) q[1];
rz(-2.9543058) q[2];
sx q[2];
rz(-1.5178871) q[2];
sx q[2];
rz(-1.0165745) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4650326) q[1];
sx q[1];
rz(-2.0102215) q[1];
sx q[1];
rz(1.7656754) q[1];
rz(-pi) q[2];
rz(-2.8630303) q[3];
sx q[3];
rz(-1.2435561) q[3];
sx q[3];
rz(2.183941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.08192) q[2];
sx q[2];
rz(-2.4884188) q[2];
sx q[2];
rz(3.0461404) q[2];
rz(-2.0898315) q[3];
sx q[3];
rz(-1.2442518) q[3];
sx q[3];
rz(-0.89422798) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8555701) q[0];
sx q[0];
rz(-0.48208553) q[0];
sx q[0];
rz(1.9739738) q[0];
rz(-2.7883912) q[1];
sx q[1];
rz(-2.628852) q[1];
sx q[1];
rz(2.8599427) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62171157) q[0];
sx q[0];
rz(-0.48292749) q[0];
sx q[0];
rz(1.3819225) q[0];
x q[1];
rz(3.0007576) q[2];
sx q[2];
rz(-0.65208921) q[2];
sx q[2];
rz(-1.2229133) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.901501) q[1];
sx q[1];
rz(-2.9411843) q[1];
sx q[1];
rz(2.3095381) q[1];
rz(-1.0135039) q[3];
sx q[3];
rz(-1.9907328) q[3];
sx q[3];
rz(0.98972631) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99792751) q[2];
sx q[2];
rz(-1.9196271) q[2];
sx q[2];
rz(-1.3227051) q[2];
rz(1.5023242) q[3];
sx q[3];
rz(-2.0570677) q[3];
sx q[3];
rz(0.098202078) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1439576) q[0];
sx q[0];
rz(-1.2459545) q[0];
sx q[0];
rz(2.8676046) q[0];
rz(-0.59448376) q[1];
sx q[1];
rz(-1.4130054) q[1];
sx q[1];
rz(-0.53057539) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.119736) q[0];
sx q[0];
rz(-1.5314191) q[0];
sx q[0];
rz(1.574613) q[0];
rz(3.0184348) q[2];
sx q[2];
rz(-2.5854857) q[2];
sx q[2];
rz(1.7983939) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.952632) q[1];
sx q[1];
rz(-2.7285693) q[1];
sx q[1];
rz(-0.4255336) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85803558) q[3];
sx q[3];
rz(-1.4467561) q[3];
sx q[3];
rz(-0.8197561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24255594) q[2];
sx q[2];
rz(-1.9344923) q[2];
sx q[2];
rz(-0.25699082) q[2];
rz(-1.1993923) q[3];
sx q[3];
rz(-1.2446087) q[3];
sx q[3];
rz(-1.6283584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5523819) q[0];
sx q[0];
rz(-0.95443812) q[0];
sx q[0];
rz(-0.5300262) q[0];
rz(-1.5026622) q[1];
sx q[1];
rz(-1.9866147) q[1];
sx q[1];
rz(-2.2030305) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28806557) q[0];
sx q[0];
rz(-0.99366436) q[0];
sx q[0];
rz(-1.0017299) q[0];
rz(-1.6907127) q[2];
sx q[2];
rz(-1.7133822) q[2];
sx q[2];
rz(0.8220807) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5512498) q[1];
sx q[1];
rz(-1.9971202) q[1];
sx q[1];
rz(0.23386441) q[1];
x q[2];
rz(0.13646941) q[3];
sx q[3];
rz(-2.1519063) q[3];
sx q[3];
rz(0.438353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3878801) q[2];
sx q[2];
rz(-2.2215999) q[2];
sx q[2];
rz(-1.0803427) q[2];
rz(-2.3267817) q[3];
sx q[3];
rz(-1.1956513) q[3];
sx q[3];
rz(-1.4174392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1431047) q[0];
sx q[0];
rz(-1.657722) q[0];
sx q[0];
rz(-2.0527573) q[0];
rz(-3.0740652) q[1];
sx q[1];
rz(-1.8636401) q[1];
sx q[1];
rz(0.75928226) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1812173) q[0];
sx q[0];
rz(-1.3160254) q[0];
sx q[0];
rz(0.1228469) q[0];
x q[1];
rz(1.4467054) q[2];
sx q[2];
rz(-1.795009) q[2];
sx q[2];
rz(-1.2808764) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2226505) q[1];
sx q[1];
rz(-2.3098051) q[1];
sx q[1];
rz(-1.2779425) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7966657) q[3];
sx q[3];
rz(-1.7874679) q[3];
sx q[3];
rz(-0.39122736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33409432) q[2];
sx q[2];
rz(-1.0958025) q[2];
sx q[2];
rz(2.5661772) q[2];
rz(0.56257644) q[3];
sx q[3];
rz(-2.176599) q[3];
sx q[3];
rz(-1.3728728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0811049) q[0];
sx q[0];
rz(-1.2860379) q[0];
sx q[0];
rz(2.2338569) q[0];
rz(-1.0870712) q[1];
sx q[1];
rz(-2.0094951) q[1];
sx q[1];
rz(3.1241945) q[1];
rz(-2.9870177) q[2];
sx q[2];
rz(-1.7107202) q[2];
sx q[2];
rz(-1.163027) q[2];
rz(1.1997052) q[3];
sx q[3];
rz(-2.1368847) q[3];
sx q[3];
rz(0.18659244) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
