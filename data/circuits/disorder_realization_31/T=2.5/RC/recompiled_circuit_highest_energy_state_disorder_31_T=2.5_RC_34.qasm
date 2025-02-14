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
rz(-2.9873084) q[1];
sx q[1];
rz(-1.0839387) q[1];
sx q[1];
rz(-1.3561603) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0480373) q[0];
sx q[0];
rz(-1.1021656) q[0];
sx q[0];
rz(0.2524391) q[0];
x q[1];
rz(2.9886099) q[2];
sx q[2];
rz(-1.6041179) q[2];
sx q[2];
rz(-2.2696537) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6068597) q[1];
sx q[1];
rz(-0.71355017) q[1];
sx q[1];
rz(0.092852863) q[1];
rz(1.1514444) q[3];
sx q[3];
rz(-1.1287125) q[3];
sx q[3];
rz(2.6579554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1821182) q[2];
sx q[2];
rz(-2.0017767) q[2];
sx q[2];
rz(3.0457022) q[2];
rz(2.828756) q[3];
sx q[3];
rz(-1.0597119) q[3];
sx q[3];
rz(0.50725168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4023912) q[0];
sx q[0];
rz(-1.5132138) q[0];
sx q[0];
rz(-1.2865404) q[0];
rz(-1.3013499) q[1];
sx q[1];
rz(-1.2570612) q[1];
sx q[1];
rz(1.4305065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1576618) q[0];
sx q[0];
rz(-2.0832015) q[0];
sx q[0];
rz(-0.061467193) q[0];
x q[1];
rz(-2.6941259) q[2];
sx q[2];
rz(-1.1612084) q[2];
sx q[2];
rz(1.2818532) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.41005653) q[1];
sx q[1];
rz(-1.441873) q[1];
sx q[1];
rz(-1.4534662) q[1];
x q[2];
rz(-2.4336171) q[3];
sx q[3];
rz(-1.7574508) q[3];
sx q[3];
rz(2.942323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8282738) q[2];
sx q[2];
rz(-0.88130772) q[2];
sx q[2];
rz(0.71448294) q[2];
rz(-2.1099527) q[3];
sx q[3];
rz(-2.1561626) q[3];
sx q[3];
rz(0.45903444) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3752876) q[0];
sx q[0];
rz(-2.3630688) q[0];
sx q[0];
rz(0.2226204) q[0];
rz(0.34046945) q[1];
sx q[1];
rz(-1.4687931) q[1];
sx q[1];
rz(0.31390831) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1134046) q[0];
sx q[0];
rz(-2.2566099) q[0];
sx q[0];
rz(0.96341204) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.48431335) q[2];
sx q[2];
rz(-2.2398758) q[2];
sx q[2];
rz(-0.95916699) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7063278) q[1];
sx q[1];
rz(-1.8774596) q[1];
sx q[1];
rz(0.72172647) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7438358) q[3];
sx q[3];
rz(-0.53103775) q[3];
sx q[3];
rz(-0.49285313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8935304) q[2];
sx q[2];
rz(-1.5241104) q[2];
sx q[2];
rz(2.1171872) q[2];
rz(1.3192568) q[3];
sx q[3];
rz(-2.8596467) q[3];
sx q[3];
rz(-0.2969186) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.1612741) q[0];
sx q[0];
rz(-1.4264822) q[0];
sx q[0];
rz(1.0262161) q[0];
rz(0.16464344) q[1];
sx q[1];
rz(-2.7963729) q[1];
sx q[1];
rz(-2.3688597) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29870009) q[0];
sx q[0];
rz(-2.1917081) q[0];
sx q[0];
rz(-2.0218532) q[0];
rz(-pi) q[1];
x q[1];
rz(0.010920694) q[2];
sx q[2];
rz(-2.7838363) q[2];
sx q[2];
rz(-2.8657279) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.0028006) q[1];
sx q[1];
rz(-1.7874663) q[1];
sx q[1];
rz(-2.2598221) q[1];
x q[2];
rz(3.1373259) q[3];
sx q[3];
rz(-1.6223309) q[3];
sx q[3];
rz(0.70642614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2356977) q[2];
sx q[2];
rz(-2.5696281) q[2];
sx q[2];
rz(-2.8254361) q[2];
rz(0.21137992) q[3];
sx q[3];
rz(-1.4485161) q[3];
sx q[3];
rz(-2.0785418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0036156) q[0];
sx q[0];
rz(-1.6162385) q[0];
sx q[0];
rz(-2.3620918) q[0];
rz(1.7100854) q[1];
sx q[1];
rz(-1.5305488) q[1];
sx q[1];
rz(3.0140108) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12870991) q[0];
sx q[0];
rz(-1.6054376) q[0];
sx q[0];
rz(1.4564464) q[0];
rz(-0.63795264) q[2];
sx q[2];
rz(-1.049713) q[2];
sx q[2];
rz(2.6661154) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7292916) q[1];
sx q[1];
rz(-2.6874098) q[1];
sx q[1];
rz(0.95328625) q[1];
rz(-pi) q[2];
rz(0.63739325) q[3];
sx q[3];
rz(-2.0717607) q[3];
sx q[3];
rz(2.4911309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7369507) q[2];
sx q[2];
rz(-1.5751795) q[2];
sx q[2];
rz(2.289782) q[2];
rz(2.8972054) q[3];
sx q[3];
rz(-0.71072018) q[3];
sx q[3];
rz(-1.3062668) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39563018) q[0];
sx q[0];
rz(-1.6599449) q[0];
sx q[0];
rz(-3.0356044) q[0];
rz(1.7164187) q[1];
sx q[1];
rz(-1.1499848) q[1];
sx q[1];
rz(-1.0894159) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46360923) q[0];
sx q[0];
rz(-0.24782032) q[0];
sx q[0];
rz(-1.8576966) q[0];
rz(-pi) q[1];
rz(-2.4545499) q[2];
sx q[2];
rz(-1.9222447) q[2];
sx q[2];
rz(2.7907651) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8196667) q[1];
sx q[1];
rz(-1.2533256) q[1];
sx q[1];
rz(-0.13369932) q[1];
rz(-pi) q[2];
rz(1.7884729) q[3];
sx q[3];
rz(-0.62575118) q[3];
sx q[3];
rz(1.9754639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0561698) q[2];
sx q[2];
rz(-1.9499754) q[2];
sx q[2];
rz(0.58758152) q[2];
rz(1.2299296) q[3];
sx q[3];
rz(-0.39837024) q[3];
sx q[3];
rz(2.5689094) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2040445) q[0];
sx q[0];
rz(-1.9119104) q[0];
sx q[0];
rz(-0.44573927) q[0];
rz(2.2988689) q[1];
sx q[1];
rz(-0.80519599) q[1];
sx q[1];
rz(-2.3420948) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84265316) q[0];
sx q[0];
rz(-1.1976722) q[0];
sx q[0];
rz(1.1375269) q[0];
rz(-2.4183774) q[2];
sx q[2];
rz(-0.84645459) q[2];
sx q[2];
rz(0.81956132) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7789797) q[1];
sx q[1];
rz(-1.1396761) q[1];
sx q[1];
rz(-2.0126473) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0754856) q[3];
sx q[3];
rz(-0.99190208) q[3];
sx q[3];
rz(-1.2829977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2004956) q[2];
sx q[2];
rz(-2.7636187) q[2];
sx q[2];
rz(-2.359158) q[2];
rz(-0.93067074) q[3];
sx q[3];
rz(-1.8698255) q[3];
sx q[3];
rz(2.8455287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3928669) q[0];
sx q[0];
rz(-0.94806945) q[0];
sx q[0];
rz(-2.0178846) q[0];
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
rz(-0.26396093) q[0];
sx q[0];
rz(-3.0288806) q[0];
sx q[0];
rz(2.0196223) q[0];
rz(-2.8181067) q[2];
sx q[2];
rz(-1.7712103) q[2];
sx q[2];
rz(2.2866164) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.9509134) q[1];
sx q[1];
rz(-1.2251405) q[1];
sx q[1];
rz(2.9596364) q[1];
rz(-pi) q[2];
rz(1.6540307) q[3];
sx q[3];
rz(-2.1862767) q[3];
sx q[3];
rz(-0.38698254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.8380518) q[2];
sx q[2];
rz(-2.3736062) q[2];
sx q[2];
rz(-1.5194019) q[2];
rz(-1.0649902) q[3];
sx q[3];
rz(-2.0678554) q[3];
sx q[3];
rz(-2.2860693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4143455) q[0];
sx q[0];
rz(-0.94731826) q[0];
sx q[0];
rz(-1.974768) q[0];
rz(-2.836152) q[1];
sx q[1];
rz(-2.0716397) q[1];
sx q[1];
rz(0.50183141) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1129285) q[0];
sx q[0];
rz(-1.477286) q[0];
sx q[0];
rz(2.9458955) q[0];
rz(-pi) q[1];
rz(-2.2756349) q[2];
sx q[2];
rz(-0.9532477) q[2];
sx q[2];
rz(2.7028529) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.045446) q[1];
sx q[1];
rz(-1.4973565) q[1];
sx q[1];
rz(1.2051084) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54399854) q[3];
sx q[3];
rz(-2.0680313) q[3];
sx q[3];
rz(-1.0656644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.66494232) q[2];
sx q[2];
rz(-2.3822337) q[2];
sx q[2];
rz(3.0090028) q[2];
rz(2.8442123) q[3];
sx q[3];
rz(-1.6011651) q[3];
sx q[3];
rz(-0.69445777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7731758) q[0];
sx q[0];
rz(-1.4765803) q[0];
sx q[0];
rz(-1.52894) q[0];
rz(-1.7932786) q[1];
sx q[1];
rz(-1.3765843) q[1];
sx q[1];
rz(0.37337676) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7517029) q[0];
sx q[0];
rz(-1.9284964) q[0];
sx q[0];
rz(-0.16660868) q[0];
rz(-pi) q[1];
rz(-1.2676722) q[2];
sx q[2];
rz(-1.1071919) q[2];
sx q[2];
rz(1.3410717) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5525145) q[1];
sx q[1];
rz(-2.2158872) q[1];
sx q[1];
rz(0.46941514) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7885429) q[3];
sx q[3];
rz(-1.3065803) q[3];
sx q[3];
rz(1.7154201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1030964) q[2];
sx q[2];
rz(-0.080048397) q[2];
sx q[2];
rz(-2.6888964) q[2];
rz(0.90940851) q[3];
sx q[3];
rz(-1.012864) q[3];
sx q[3];
rz(-0.67210853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(1.7047003) q[0];
sx q[0];
rz(-1.1672651) q[0];
sx q[0];
rz(1.1080909) q[0];
rz(-0.7946026) q[1];
sx q[1];
rz(-1.662685) q[1];
sx q[1];
rz(-0.86659238) q[1];
rz(-2.2617302) q[2];
sx q[2];
rz(-1.5884279) q[2];
sx q[2];
rz(1.8691155) q[2];
rz(2.0183715) q[3];
sx q[3];
rz(-2.2991857) q[3];
sx q[3];
rz(-3.1274336) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
