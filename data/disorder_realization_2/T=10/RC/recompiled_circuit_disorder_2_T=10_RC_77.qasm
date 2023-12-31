OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4484654) q[0];
sx q[0];
rz(-2.6187596) q[0];
sx q[0];
rz(-2.5180106) q[0];
rz(0.29016718) q[1];
sx q[1];
rz(-2.4224412) q[1];
sx q[1];
rz(2.6410988) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61636415) q[0];
sx q[0];
rz(-0.59983569) q[0];
sx q[0];
rz(-3.1347549) q[0];
rz(-pi) q[1];
rz(3.0481553) q[2];
sx q[2];
rz(-1.4215901) q[2];
sx q[2];
rz(-2.0051533) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3574672) q[1];
sx q[1];
rz(-2.3298414) q[1];
sx q[1];
rz(-2.0264506) q[1];
x q[2];
rz(1.2455363) q[3];
sx q[3];
rz(-0.45913011) q[3];
sx q[3];
rz(-1.1577275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3502675) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(1.9187437) q[2];
rz(1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-2.1712415) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.37110776) q[0];
sx q[0];
rz(-1.4828232) q[0];
sx q[0];
rz(2.0626542) q[0];
rz(1.7547912) q[1];
sx q[1];
rz(-2.3290122) q[1];
sx q[1];
rz(-2.4761377) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6404214) q[0];
sx q[0];
rz(-1.2873532) q[0];
sx q[0];
rz(0.47407504) q[0];
rz(-pi) q[1];
x q[1];
rz(0.32174343) q[2];
sx q[2];
rz(-1.7330568) q[2];
sx q[2];
rz(-2.4010047) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7484819) q[1];
sx q[1];
rz(-2.0512274) q[1];
sx q[1];
rz(1.837681) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5311702) q[3];
sx q[3];
rz(-2.431776) q[3];
sx q[3];
rz(0.94330793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60454303) q[2];
sx q[2];
rz(-1.7384572) q[2];
sx q[2];
rz(0.075142168) q[2];
rz(1.6710619) q[3];
sx q[3];
rz(-2.0139147) q[3];
sx q[3];
rz(-2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55494088) q[0];
sx q[0];
rz(-1.2723158) q[0];
sx q[0];
rz(-2.3828322) q[0];
rz(1.8485908) q[1];
sx q[1];
rz(-1.8483775) q[1];
sx q[1];
rz(-1.4000777) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6436359) q[0];
sx q[0];
rz(-0.47753497) q[0];
sx q[0];
rz(2.5316373) q[0];
rz(-2.5163469) q[2];
sx q[2];
rz(-0.22908224) q[2];
sx q[2];
rz(0.89876995) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.61664596) q[1];
sx q[1];
rz(-2.4541306) q[1];
sx q[1];
rz(-2.0928659) q[1];
rz(-pi) q[2];
rz(-2.9754144) q[3];
sx q[3];
rz(-1.8718534) q[3];
sx q[3];
rz(2.500246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.489958) q[2];
sx q[2];
rz(-2.659446) q[2];
sx q[2];
rz(2.4839694) q[2];
rz(1.1714606) q[3];
sx q[3];
rz(-1.6572584) q[3];
sx q[3];
rz(1.4191779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79384971) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(-1.4452274) q[0];
rz(1.6943278) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(0.34805527) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9075273) q[0];
sx q[0];
rz(-2.4287927) q[0];
sx q[0];
rz(-0.45760052) q[0];
x q[1];
rz(-0.15050998) q[2];
sx q[2];
rz(-1.9410656) q[2];
sx q[2];
rz(-0.18266695) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.5665633) q[1];
sx q[1];
rz(-1.6600779) q[1];
sx q[1];
rz(-0.95972285) q[1];
x q[2];
rz(0.69339852) q[3];
sx q[3];
rz(-2.9326673) q[3];
sx q[3];
rz(2.8914176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.7248914) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(-0.42281881) q[2];
rz(-0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87930644) q[0];
sx q[0];
rz(-1.3129741) q[0];
sx q[0];
rz(-0.53043956) q[0];
rz(0.92492217) q[1];
sx q[1];
rz(-1.568012) q[1];
sx q[1];
rz(1.2984498) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0153774) q[0];
sx q[0];
rz(-1.7324289) q[0];
sx q[0];
rz(0.14590185) q[0];
rz(0.41579397) q[2];
sx q[2];
rz(-2.4928164) q[2];
sx q[2];
rz(-2.6598038) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7309155) q[1];
sx q[1];
rz(-1.2123322) q[1];
sx q[1];
rz(-1.4645542) q[1];
rz(2.2523746) q[3];
sx q[3];
rz(-2.4272356) q[3];
sx q[3];
rz(-2.0444972) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.986074) q[2];
sx q[2];
rz(0.47362622) q[2];
rz(3.04223) q[3];
sx q[3];
rz(-1.8615581) q[3];
sx q[3];
rz(2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6376003) q[0];
sx q[0];
rz(-1.3562599) q[0];
sx q[0];
rz(2.1437058) q[0];
rz(2.2672794) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(-0.46674892) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1186819) q[0];
sx q[0];
rz(-0.70436275) q[0];
sx q[0];
rz(0.33606152) q[0];
rz(-1.7919962) q[2];
sx q[2];
rz(-1.9890519) q[2];
sx q[2];
rz(2.7115371) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.1280061) q[1];
sx q[1];
rz(-0.15833536) q[1];
sx q[1];
rz(2.8360785) q[1];
rz(-pi) q[2];
rz(1.0420226) q[3];
sx q[3];
rz(-0.18897945) q[3];
sx q[3];
rz(2.6721862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2816887) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(-1.5768645) q[2];
rz(2.5148897) q[3];
sx q[3];
rz(-1.7765216) q[3];
sx q[3];
rz(-2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3307813) q[0];
sx q[0];
rz(-2.2450876) q[0];
sx q[0];
rz(-0.41982857) q[0];
rz(0.22142521) q[1];
sx q[1];
rz(-0.47859335) q[1];
sx q[1];
rz(1.5484757) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67328582) q[0];
sx q[0];
rz(-1.5791897) q[0];
sx q[0];
rz(-3.054137) q[0];
rz(-pi) q[1];
rz(0.19303796) q[2];
sx q[2];
rz(-1.8503354) q[2];
sx q[2];
rz(2.2045731) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.78649) q[1];
sx q[1];
rz(-0.10663248) q[1];
sx q[1];
rz(2.998888) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1057304) q[3];
sx q[3];
rz(-1.8814058) q[3];
sx q[3];
rz(0.066699337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.55591136) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(1.7377724) q[2];
rz(-0.79706556) q[3];
sx q[3];
rz(-2.6795487) q[3];
sx q[3];
rz(1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4306915) q[0];
sx q[0];
rz(-0.71390188) q[0];
sx q[0];
rz(-0.28924334) q[0];
rz(-2.5166683) q[1];
sx q[1];
rz(-1.9754675) q[1];
sx q[1];
rz(1.3141059) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4411366) q[0];
sx q[0];
rz(-1.4608129) q[0];
sx q[0];
rz(-2.3814047) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7881175) q[2];
sx q[2];
rz(-0.77958737) q[2];
sx q[2];
rz(-1.0970864) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.96771679) q[1];
sx q[1];
rz(-1.6248676) q[1];
sx q[1];
rz(1.6915583) q[1];
rz(-0.81359158) q[3];
sx q[3];
rz(-1.7980669) q[3];
sx q[3];
rz(3.1384625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4079995) q[2];
sx q[2];
rz(-1.5780129) q[2];
sx q[2];
rz(-1.4286263) q[2];
rz(2.1777878) q[3];
sx q[3];
rz(-2.0644085) q[3];
sx q[3];
rz(1.0296286) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6190417) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(-1.8956986) q[0];
rz(-0.11101162) q[1];
sx q[1];
rz(-1.9440034) q[1];
sx q[1];
rz(0.54661173) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3205991) q[0];
sx q[0];
rz(-2.5490767) q[0];
sx q[0];
rz(-2.2792363) q[0];
rz(1.1999646) q[2];
sx q[2];
rz(-1.2002581) q[2];
sx q[2];
rz(1.3818936) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4709028) q[1];
sx q[1];
rz(-2.7109475) q[1];
sx q[1];
rz(1.9914658) q[1];
rz(-0.9162174) q[3];
sx q[3];
rz(-1.1715874) q[3];
sx q[3];
rz(2.6938714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1116011) q[2];
sx q[2];
rz(-2.2932055) q[2];
sx q[2];
rz(-1.6938422) q[2];
rz(-2.0041806) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(2.494273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85957134) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(-0.71722537) q[0];
rz(-1.2099129) q[1];
sx q[1];
rz(-0.33214339) q[1];
sx q[1];
rz(2.4338914) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.087698547) q[0];
sx q[0];
rz(-2.2367034) q[0];
sx q[0];
rz(0.44184394) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3466481) q[2];
sx q[2];
rz(-1.4321616) q[2];
sx q[2];
rz(2.5773406) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3860491) q[1];
sx q[1];
rz(-1.886133) q[1];
sx q[1];
rz(-2.3029033) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1901312) q[3];
sx q[3];
rz(-1.8639495) q[3];
sx q[3];
rz(-2.4952863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.5132961) q[2];
sx q[2];
rz(-2.4847023) q[2];
sx q[2];
rz(0.90325242) q[2];
rz(-1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(-2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.306504) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(2.3256336) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(2.6782398) q[2];
sx q[2];
rz(-1.1506766) q[2];
sx q[2];
rz(1.2406032) q[2];
rz(0.054374183) q[3];
sx q[3];
rz(-1.6118703) q[3];
sx q[3];
rz(0.76606228) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
