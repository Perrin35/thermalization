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
rz(-2.5950522) q[0];
sx q[0];
rz(-0.15931436) q[0];
sx q[0];
rz(-1.8848609) q[0];
rz(0.15428421) q[1];
sx q[1];
rz(4.2255314) q[1];
sx q[1];
rz(10.780938) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6124466) q[0];
sx q[0];
rz(-2.6137798) q[0];
sx q[0];
rz(1.1124658) q[0];
rz(-pi) q[1];
rz(1.5370813) q[2];
sx q[2];
rz(-1.4178992) q[2];
sx q[2];
rz(-0.70399415) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4843085) q[1];
sx q[1];
rz(-2.2806232) q[1];
sx q[1];
rz(1.6508938) q[1];
x q[2];
rz(-2.431193) q[3];
sx q[3];
rz(-0.59962873) q[3];
sx q[3];
rz(0.3223294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95947444) q[2];
sx q[2];
rz(-2.0017767) q[2];
sx q[2];
rz(-3.0457022) q[2];
rz(2.828756) q[3];
sx q[3];
rz(-1.0597119) q[3];
sx q[3];
rz(-2.634341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7392015) q[0];
sx q[0];
rz(-1.5132138) q[0];
sx q[0];
rz(1.2865404) q[0];
rz(1.8402428) q[1];
sx q[1];
rz(-1.8845314) q[1];
sx q[1];
rz(-1.4305065) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5848918) q[0];
sx q[0];
rz(-1.624361) q[0];
sx q[0];
rz(2.0840097) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.78708055) q[2];
sx q[2];
rz(-2.5445017) q[2];
sx q[2];
rz(-2.7379089) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.41005653) q[1];
sx q[1];
rz(-1.6997196) q[1];
sx q[1];
rz(1.6881264) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8589571) q[3];
sx q[3];
rz(-0.728033) q[3];
sx q[3];
rz(-1.9835492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3133189) q[2];
sx q[2];
rz(-2.2602849) q[2];
sx q[2];
rz(2.4271097) q[2];
rz(-2.1099527) q[3];
sx q[3];
rz(-2.1561626) q[3];
sx q[3];
rz(-2.6825582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76630509) q[0];
sx q[0];
rz(-2.3630688) q[0];
sx q[0];
rz(-2.9189723) q[0];
rz(0.34046945) q[1];
sx q[1];
rz(-1.6727996) q[1];
sx q[1];
rz(-0.31390831) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8720377) q[0];
sx q[0];
rz(-2.0282708) q[0];
sx q[0];
rz(-0.78367718) q[0];
x q[1];
rz(0.84157518) q[2];
sx q[2];
rz(-1.944659) q[2];
sx q[2];
rz(2.8453927) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7063278) q[1];
sx q[1];
rz(-1.8774596) q[1];
sx q[1];
rz(-2.4198662) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.39775689) q[3];
sx q[3];
rz(-0.53103775) q[3];
sx q[3];
rz(0.49285313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.24806222) q[2];
sx q[2];
rz(-1.5241104) q[2];
sx q[2];
rz(2.1171872) q[2];
rz(1.8223358) q[3];
sx q[3];
rz(-2.8596467) q[3];
sx q[3];
rz(-2.8446741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9803186) q[0];
sx q[0];
rz(-1.4264822) q[0];
sx q[0];
rz(2.1153765) q[0];
rz(2.9769492) q[1];
sx q[1];
rz(-0.34521979) q[1];
sx q[1];
rz(-2.3688597) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29870009) q[0];
sx q[0];
rz(-0.94988454) q[0];
sx q[0];
rz(1.1197394) q[0];
x q[1];
rz(2.7838559) q[2];
sx q[2];
rz(-1.5669723) q[2];
sx q[2];
rz(1.8568904) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3983359) q[1];
sx q[1];
rz(-0.90088049) q[1];
sx q[1];
rz(-0.27781133) q[1];
rz(-3.1373259) q[3];
sx q[3];
rz(-1.6223309) q[3];
sx q[3];
rz(-0.70642614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2356977) q[2];
sx q[2];
rz(-2.5696281) q[2];
sx q[2];
rz(0.31615654) q[2];
rz(2.9302127) q[3];
sx q[3];
rz(-1.6930765) q[3];
sx q[3];
rz(1.0630509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.137977) q[0];
sx q[0];
rz(-1.6162385) q[0];
sx q[0];
rz(0.77950087) q[0];
rz(-1.7100854) q[1];
sx q[1];
rz(-1.6110438) q[1];
sx q[1];
rz(-0.12758189) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0128827) q[0];
sx q[0];
rz(-1.6054376) q[0];
sx q[0];
rz(1.4564464) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1912135) q[2];
sx q[2];
rz(-2.1135672) q[2];
sx q[2];
rz(-1.448878) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0605088) q[1];
sx q[1];
rz(-1.2049872) q[1];
sx q[1];
rz(2.8660956) q[1];
rz(-pi) q[2];
rz(-2.1688858) q[3];
sx q[3];
rz(-1.0216301) q[3];
sx q[3];
rz(1.8795667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.40464193) q[2];
sx q[2];
rz(-1.5751795) q[2];
sx q[2];
rz(2.289782) q[2];
rz(0.24438721) q[3];
sx q[3];
rz(-2.4308725) q[3];
sx q[3];
rz(-1.3062668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7459625) q[0];
sx q[0];
rz(-1.4816477) q[0];
sx q[0];
rz(-0.10598824) q[0];
rz(1.7164187) q[1];
sx q[1];
rz(-1.1499848) q[1];
sx q[1];
rz(2.0521767) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9734151) q[0];
sx q[0];
rz(-1.3333048) q[0];
sx q[0];
rz(3.0701145) q[0];
x q[1];
rz(2.6174215) q[2];
sx q[2];
rz(-0.75850111) q[2];
sx q[2];
rz(0.82265178) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9346843) q[1];
sx q[1];
rz(-1.6977773) q[1];
sx q[1];
rz(1.8909353) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1852538) q[3];
sx q[3];
rz(-1.4439665) q[3];
sx q[3];
rz(-2.9143132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.085422903) q[2];
sx q[2];
rz(-1.9499754) q[2];
sx q[2];
rz(-0.58758152) q[2];
rz(-1.9116631) q[3];
sx q[3];
rz(-0.39837024) q[3];
sx q[3];
rz(2.5689094) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040445) q[0];
sx q[0];
rz(-1.2296822) q[0];
sx q[0];
rz(0.44573927) q[0];
rz(-0.84272376) q[1];
sx q[1];
rz(-2.3363967) q[1];
sx q[1];
rz(-0.79949784) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.580502) q[0];
sx q[0];
rz(-1.1691302) q[0];
sx q[0];
rz(-2.7343844) q[0];
x q[1];
rz(-0.72321524) q[2];
sx q[2];
rz(-2.2951381) q[2];
sx q[2];
rz(0.81956132) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.362613) q[1];
sx q[1];
rz(-1.1396761) q[1];
sx q[1];
rz(1.1289454) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1506924) q[3];
sx q[3];
rz(-1.5154723) q[3];
sx q[3];
rz(-0.25159454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.94109702) q[2];
sx q[2];
rz(-2.7636187) q[2];
sx q[2];
rz(-0.7824347) q[2];
rz(0.93067074) q[3];
sx q[3];
rz(-1.2717671) q[3];
sx q[3];
rz(2.8455287) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74872577) q[0];
sx q[0];
rz(-2.1935232) q[0];
sx q[0];
rz(-1.123708) q[0];
rz(0.8404845) q[1];
sx q[1];
rz(-1.1437462) q[1];
sx q[1];
rz(-1.3927654) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18735841) q[0];
sx q[0];
rz(-1.4692882) q[0];
sx q[0];
rz(-0.049075394) q[0];
rz(-pi) q[1];
rz(-2.5729449) q[2];
sx q[2];
rz(-0.37867448) q[2];
sx q[2];
rz(-0.17998047) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.44237524) q[1];
sx q[1];
rz(-1.3997119) q[1];
sx q[1];
rz(1.2198001) q[1];
rz(-pi) q[2];
x q[2];
rz(1.487562) q[3];
sx q[3];
rz(-0.95531598) q[3];
sx q[3];
rz(2.7546101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8380518) q[2];
sx q[2];
rz(-2.3736062) q[2];
sx q[2];
rz(-1.5194019) q[2];
rz(2.0766025) q[3];
sx q[3];
rz(-1.0737373) q[3];
sx q[3];
rz(2.2860693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72724718) q[0];
sx q[0];
rz(-2.1942744) q[0];
sx q[0];
rz(1.1668246) q[0];
rz(-0.30544063) q[1];
sx q[1];
rz(-2.0716397) q[1];
sx q[1];
rz(-0.50183141) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.017664488) q[0];
sx q[0];
rz(-0.2166324) q[0];
sx q[0];
rz(0.44938918) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73955543) q[2];
sx q[2];
rz(-2.2408492) q[2];
sx q[2];
rz(0.53436992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.045446) q[1];
sx q[1];
rz(-1.6442361) q[1];
sx q[1];
rz(1.9364843) q[1];
rz(-pi) q[2];
rz(2.3324729) q[3];
sx q[3];
rz(-0.71965766) q[3];
sx q[3];
rz(-1.1728731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66494232) q[2];
sx q[2];
rz(-2.3822337) q[2];
sx q[2];
rz(-0.13258983) q[2];
rz(0.29738033) q[3];
sx q[3];
rz(-1.6011651) q[3];
sx q[3];
rz(-2.4471349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7731758) q[0];
sx q[0];
rz(-1.6650124) q[0];
sx q[0];
rz(-1.52894) q[0];
rz(-1.348314) q[1];
sx q[1];
rz(-1.7650083) q[1];
sx q[1];
rz(-2.7682159) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38988977) q[0];
sx q[0];
rz(-1.9284964) q[0];
sx q[0];
rz(2.974984) q[0];
x q[1];
rz(2.6033295) q[2];
sx q[2];
rz(-0.54780947) q[2];
sx q[2];
rz(0.73073587) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58907813) q[1];
sx q[1];
rz(-2.2158872) q[1];
sx q[1];
rz(2.6721775) q[1];
rz(-1.3530497) q[3];
sx q[3];
rz(-1.3065803) q[3];
sx q[3];
rz(1.7154201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1030964) q[2];
sx q[2];
rz(-3.0615443) q[2];
sx q[2];
rz(-2.6888964) q[2];
rz(-2.2321841) q[3];
sx q[3];
rz(-2.1287287) q[3];
sx q[3];
rz(0.67210853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4368923) q[0];
sx q[0];
rz(-1.9743275) q[0];
sx q[0];
rz(-2.0335017) q[0];
rz(-2.3469901) q[1];
sx q[1];
rz(-1.4789076) q[1];
sx q[1];
rz(2.2750003) q[1];
rz(-3.1187155) q[2];
sx q[2];
rz(-2.2616017) q[2];
sx q[2];
rz(0.2837413) q[2];
rz(-2.0183715) q[3];
sx q[3];
rz(-0.84240693) q[3];
sx q[3];
rz(0.014159023) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
