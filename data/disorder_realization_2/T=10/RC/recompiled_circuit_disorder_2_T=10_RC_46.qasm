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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5335124) q[0];
sx q[0];
rz(-0.97097662) q[0];
sx q[0];
rz(1.5661201) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0481553) q[2];
sx q[2];
rz(-1.4215901) q[2];
sx q[2];
rz(-1.1364394) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3574672) q[1];
sx q[1];
rz(-2.3298414) q[1];
sx q[1];
rz(1.1151421) q[1];
x q[2];
rz(-1.1327098) q[3];
sx q[3];
rz(-1.7128908) q[3];
sx q[3];
rz(-2.4349468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7913251) q[2];
sx q[2];
rz(-1.9359549) q[2];
sx q[2];
rz(-1.9187437) q[2];
rz(1.4482927) q[3];
sx q[3];
rz(-2.1494614) q[3];
sx q[3];
rz(-2.1712415) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37110776) q[0];
sx q[0];
rz(-1.6587695) q[0];
sx q[0];
rz(-2.0626542) q[0];
rz(-1.3868015) q[1];
sx q[1];
rz(-0.81258041) q[1];
sx q[1];
rz(-0.66545495) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21214813) q[0];
sx q[0];
rz(-2.0245027) q[0];
sx q[0];
rz(1.8871904) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47768728) q[2];
sx q[2];
rz(-0.3590695) q[2];
sx q[2];
rz(1.2815086) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39311073) q[1];
sx q[1];
rz(-2.0512274) q[1];
sx q[1];
rz(-1.3039116) q[1];
rz(-pi) q[2];
rz(-0.86136787) q[3];
sx q[3];
rz(-1.5449761) q[3];
sx q[3];
rz(2.4840419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5370496) q[2];
sx q[2];
rz(-1.4031354) q[2];
sx q[2];
rz(-0.075142168) q[2];
rz(-1.4705307) q[3];
sx q[3];
rz(-1.1276779) q[3];
sx q[3];
rz(2.4501734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55494088) q[0];
sx q[0];
rz(-1.8692769) q[0];
sx q[0];
rz(0.75876045) q[0];
rz(-1.8485908) q[1];
sx q[1];
rz(-1.8483775) q[1];
sx q[1];
rz(-1.741515) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4979567) q[0];
sx q[0];
rz(-0.47753497) q[0];
sx q[0];
rz(-2.5316373) q[0];
rz(-pi) q[1];
rz(-1.4351575) q[2];
sx q[2];
rz(-1.3856158) q[2];
sx q[2];
rz(2.8806825) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6058265) q[1];
sx q[1];
rz(-1.2488135) q[1];
sx q[1];
rz(2.1893326) q[1];
rz(-pi) q[2];
rz(-1.2657884) q[3];
sx q[3];
rz(-1.4121571) q[3];
sx q[3];
rz(0.97914417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.489958) q[2];
sx q[2];
rz(-0.48214665) q[2];
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
rz(pi/2) q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3477429) q[0];
sx q[0];
rz(-2.1934953) q[0];
sx q[0];
rz(-1.4452274) q[0];
rz(-1.6943278) q[1];
sx q[1];
rz(-1.4935962) q[1];
sx q[1];
rz(2.7935374) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8349583) q[0];
sx q[0];
rz(-1.277703) q[0];
sx q[0];
rz(2.4819863) q[0];
rz(-pi) q[1];
rz(-2.9910827) q[2];
sx q[2];
rz(-1.9410656) q[2];
sx q[2];
rz(-2.9589257) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0192249) q[1];
sx q[1];
rz(-2.5248563) q[1];
sx q[1];
rz(1.4160181) q[1];
rz(-pi) q[2];
rz(2.4481941) q[3];
sx q[3];
rz(-2.9326673) q[3];
sx q[3];
rz(-2.8914176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.41670123) q[2];
sx q[2];
rz(-1.8323703) q[2];
sx q[2];
rz(0.42281881) q[2];
rz(0.73741284) q[3];
sx q[3];
rz(-2.335572) q[3];
sx q[3];
rz(-0.084658682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87930644) q[0];
sx q[0];
rz(-1.8286185) q[0];
sx q[0];
rz(-2.6111531) q[0];
rz(-2.2166705) q[1];
sx q[1];
rz(-1.568012) q[1];
sx q[1];
rz(1.2984498) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0153774) q[0];
sx q[0];
rz(-1.7324289) q[0];
sx q[0];
rz(2.9956908) q[0];
rz(-pi) q[1];
rz(2.7257987) q[2];
sx q[2];
rz(-0.64877629) q[2];
sx q[2];
rz(-2.6598038) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.4106771) q[1];
sx q[1];
rz(-1.9292604) q[1];
sx q[1];
rz(-1.6770384) q[1];
rz(2.1634444) q[3];
sx q[3];
rz(-1.9962629) q[3];
sx q[3];
rz(-3.0654207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2003145) q[2];
sx q[2];
rz(-1.1555187) q[2];
sx q[2];
rz(2.6679664) q[2];
rz(-0.099362699) q[3];
sx q[3];
rz(-1.2800346) q[3];
sx q[3];
rz(-2.30106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6376003) q[0];
sx q[0];
rz(-1.7853328) q[0];
sx q[0];
rz(-0.9978869) q[0];
rz(2.2672794) q[1];
sx q[1];
rz(-2.120178) q[1];
sx q[1];
rz(2.6748437) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.022910718) q[0];
sx q[0];
rz(-2.4372299) q[0];
sx q[0];
rz(2.8055311) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42748638) q[2];
sx q[2];
rz(-1.7726521) q[2];
sx q[2];
rz(2.0919378) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0007243) q[1];
sx q[1];
rz(-1.5233526) q[1];
sx q[1];
rz(-0.15111698) q[1];
rz(2.09957) q[3];
sx q[3];
rz(-2.9526132) q[3];
sx q[3];
rz(2.6721862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2816887) q[2];
sx q[2];
rz(-2.6624661) q[2];
sx q[2];
rz(1.5768645) q[2];
rz(-0.62670296) q[3];
sx q[3];
rz(-1.3650711) q[3];
sx q[3];
rz(2.4600162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3307813) q[0];
sx q[0];
rz(-0.89650506) q[0];
sx q[0];
rz(-0.41982857) q[0];
rz(-0.22142521) q[1];
sx q[1];
rz(-2.6629993) q[1];
sx q[1];
rz(-1.5931169) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3395183) q[0];
sx q[0];
rz(-0.087856494) q[0];
sx q[0];
rz(-0.095803424) q[0];
x q[1];
rz(1.2862455) q[2];
sx q[2];
rz(-1.3853405) q[2];
sx q[2];
rz(-0.68765771) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.49861136) q[1];
sx q[1];
rz(-1.4652518) q[1];
sx q[1];
rz(-1.5555744) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7967954) q[3];
sx q[3];
rz(-1.1296009) q[3];
sx q[3];
rz(-1.7896717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.5856813) q[2];
sx q[2];
rz(-2.6101117) q[2];
sx q[2];
rz(1.7377724) q[2];
rz(-2.3445271) q[3];
sx q[3];
rz(-0.46204391) q[3];
sx q[3];
rz(1.2169303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4306915) q[0];
sx q[0];
rz(-2.4276908) q[0];
sx q[0];
rz(0.28924334) q[0];
rz(-2.5166683) q[1];
sx q[1];
rz(-1.1661252) q[1];
sx q[1];
rz(-1.3141059) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24459141) q[0];
sx q[0];
rz(-0.76652157) q[0];
sx q[0];
rz(-2.9826829) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7881175) q[2];
sx q[2];
rz(-0.77958737) q[2];
sx q[2];
rz(1.0970864) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.96771679) q[1];
sx q[1];
rz(-1.6248676) q[1];
sx q[1];
rz(1.6915583) q[1];
rz(1.8955599) q[3];
sx q[3];
rz(-0.78403463) q[3];
sx q[3];
rz(1.8079545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4079995) q[2];
sx q[2];
rz(-1.5635798) q[2];
sx q[2];
rz(1.4286263) q[2];
rz(-2.1777878) q[3];
sx q[3];
rz(-1.0771841) q[3];
sx q[3];
rz(1.0296286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.52255094) q[0];
sx q[0];
rz(-1.4551117) q[0];
sx q[0];
rz(1.8956986) q[0];
rz(3.030581) q[1];
sx q[1];
rz(-1.1975892) q[1];
sx q[1];
rz(2.5949809) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1221736) q[0];
sx q[0];
rz(-2.008736) q[0];
sx q[0];
rz(0.41284783) q[0];
rz(-pi) q[1];
x q[1];
rz(1.941628) q[2];
sx q[2];
rz(-1.9413345) q[2];
sx q[2];
rz(-1.759699) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.48601549) q[1];
sx q[1];
rz(-1.3994819) q[1];
sx q[1];
rz(1.1737215) q[1];
rz(-pi) q[2];
rz(2.6528477) q[3];
sx q[3];
rz(-2.1663323) q[3];
sx q[3];
rz(2.3084156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.1116011) q[2];
sx q[2];
rz(-0.84838715) q[2];
sx q[2];
rz(1.4477504) q[2];
rz(2.0041806) q[3];
sx q[3];
rz(-1.3237938) q[3];
sx q[3];
rz(0.64731961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2820213) q[0];
sx q[0];
rz(-0.30680007) q[0];
sx q[0];
rz(0.71722537) q[0];
rz(1.2099129) q[1];
sx q[1];
rz(-2.8094493) q[1];
sx q[1];
rz(-0.70770121) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0538941) q[0];
sx q[0];
rz(-0.90488926) q[0];
sx q[0];
rz(0.44184394) q[0];
x q[1];
rz(1.3466481) q[2];
sx q[2];
rz(-1.7094311) q[2];
sx q[2];
rz(-0.56425205) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.3860491) q[1];
sx q[1];
rz(-1.886133) q[1];
sx q[1];
rz(-2.3029033) q[1];
rz(-pi) q[2];
rz(-0.95146146) q[3];
sx q[3];
rz(-1.8639495) q[3];
sx q[3];
rz(2.4952863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6282965) q[2];
sx q[2];
rz(-0.6568903) q[2];
sx q[2];
rz(-0.90325242) q[2];
rz(1.603027) q[3];
sx q[3];
rz(-2.2730946) q[3];
sx q[3];
rz(2.2911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83508867) q[0];
sx q[0];
rz(-0.36515129) q[0];
sx q[0];
rz(-0.93602244) q[0];
rz(2.3256336) q[1];
sx q[1];
rz(-0.42146704) q[1];
sx q[1];
rz(-2.0889919) q[1];
rz(2.03394) q[2];
sx q[2];
rz(-1.1504428) q[2];
sx q[2];
rz(3.0124315) q[2];
rz(-1.5296616) q[3];
sx q[3];
rz(-1.516468) q[3];
sx q[3];
rz(-0.80249912) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
