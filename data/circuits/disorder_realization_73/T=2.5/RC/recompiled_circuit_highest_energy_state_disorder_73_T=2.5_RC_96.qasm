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
rz(1.3923378) q[0];
sx q[0];
rz(-2.7633986) q[0];
sx q[0];
rz(2.7910772) q[0];
rz(-2.2502083) q[1];
sx q[1];
rz(-0.79217029) q[1];
sx q[1];
rz(-1.1729191) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.354686) q[0];
sx q[0];
rz(-1.3808492) q[0];
sx q[0];
rz(-2.9482916) q[0];
x q[1];
rz(0.15526659) q[2];
sx q[2];
rz(-1.7874996) q[2];
sx q[2];
rz(-2.7086176) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.71095548) q[1];
sx q[1];
rz(-1.4078377) q[1];
sx q[1];
rz(-0.32117543) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8137742) q[3];
sx q[3];
rz(-2.647305) q[3];
sx q[3];
rz(1.0456628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17060908) q[2];
sx q[2];
rz(-1.8656518) q[2];
sx q[2];
rz(-0.21273908) q[2];
rz(0.085266026) q[3];
sx q[3];
rz(-1.250896) q[3];
sx q[3];
rz(-1.2712449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3474715) q[0];
sx q[0];
rz(-0.81231064) q[0];
sx q[0];
rz(-1.043327) q[0];
rz(3.0087545) q[1];
sx q[1];
rz(-2.9265407) q[1];
sx q[1];
rz(-1.1588233) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83860676) q[0];
sx q[0];
rz(-1.9667454) q[0];
sx q[0];
rz(1.4141757) q[0];
x q[1];
rz(0.01387502) q[2];
sx q[2];
rz(-1.4094947) q[2];
sx q[2];
rz(-0.32729766) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.68090428) q[1];
sx q[1];
rz(-1.0863606) q[1];
sx q[1];
rz(-1.2699732) q[1];
rz(-pi) q[2];
rz(2.0015249) q[3];
sx q[3];
rz(-1.8716806) q[3];
sx q[3];
rz(-1.6611851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.1966689) q[2];
sx q[2];
rz(-2.479574) q[2];
sx q[2];
rz(0.72466737) q[2];
rz(-2.479018) q[3];
sx q[3];
rz(-0.96810883) q[3];
sx q[3];
rz(-0.29964963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.9579983) q[0];
sx q[0];
rz(-1.8738926) q[0];
sx q[0];
rz(-0.9147574) q[0];
rz(2.5547408) q[1];
sx q[1];
rz(-1.4205168) q[1];
sx q[1];
rz(-0.14952001) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0256074) q[0];
sx q[0];
rz(-1.43029) q[0];
sx q[0];
rz(1.0464267) q[0];
rz(-pi) q[1];
rz(-2.3517866) q[2];
sx q[2];
rz(-0.49587223) q[2];
sx q[2];
rz(-0.42988955) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.65626493) q[1];
sx q[1];
rz(-1.0228923) q[1];
sx q[1];
rz(-2.5119971) q[1];
x q[2];
rz(-2.3749659) q[3];
sx q[3];
rz(-2.2671923) q[3];
sx q[3];
rz(1.6625787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.84825039) q[2];
sx q[2];
rz(-1.9363656) q[2];
sx q[2];
rz(0.80043522) q[2];
rz(0.46659255) q[3];
sx q[3];
rz(-1.1060017) q[3];
sx q[3];
rz(2.485937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34224299) q[0];
sx q[0];
rz(-1.290134) q[0];
sx q[0];
rz(-0.85860646) q[0];
rz(2.730992) q[1];
sx q[1];
rz(-0.20315367) q[1];
sx q[1];
rz(0.98791775) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3666375) q[0];
sx q[0];
rz(-1.8995842) q[0];
sx q[0];
rz(-2.9820697) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0201077) q[2];
sx q[2];
rz(-2.0628235) q[2];
sx q[2];
rz(2.55388) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.021182755) q[1];
sx q[1];
rz(-1.2721918) q[1];
sx q[1];
rz(2.8575767) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32295042) q[3];
sx q[3];
rz(-1.3846372) q[3];
sx q[3];
rz(1.8017839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9504488) q[2];
sx q[2];
rz(-1.3942275) q[2];
sx q[2];
rz(0.44551715) q[2];
rz(2.4281003) q[3];
sx q[3];
rz(-2.4092509) q[3];
sx q[3];
rz(-2.2215686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54818654) q[0];
sx q[0];
rz(-1.0819409) q[0];
sx q[0];
rz(1.5997546) q[0];
rz(-1.8707188) q[1];
sx q[1];
rz(-1.4022695) q[1];
sx q[1];
rz(0.52938968) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75225509) q[0];
sx q[0];
rz(-1.9175954) q[0];
sx q[0];
rz(2.8777468) q[0];
rz(-pi) q[1];
rz(-0.28044706) q[2];
sx q[2];
rz(-0.63647645) q[2];
sx q[2];
rz(2.150879) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3965949) q[1];
sx q[1];
rz(-1.0675808) q[1];
sx q[1];
rz(-1.2835591) q[1];
rz(-0.1481109) q[3];
sx q[3];
rz(-1.221773) q[3];
sx q[3];
rz(2.7452041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.14220898) q[2];
sx q[2];
rz(-1.9500407) q[2];
sx q[2];
rz(-1.3231529) q[2];
rz(0.38149825) q[3];
sx q[3];
rz(-2.0254717) q[3];
sx q[3];
rz(2.748446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8489654) q[0];
sx q[0];
rz(-2.3830074) q[0];
sx q[0];
rz(1.0211771) q[0];
rz(-1.609833) q[1];
sx q[1];
rz(-2.1949218) q[1];
sx q[1];
rz(2.3381332) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6670584) q[0];
sx q[0];
rz(-0.34862384) q[0];
sx q[0];
rz(-2.6571996) q[0];
rz(-pi) q[1];
x q[1];
rz(0.19101269) q[2];
sx q[2];
rz(-0.41897853) q[2];
sx q[2];
rz(1.524817) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.72587126) q[1];
sx q[1];
rz(-2.5588003) q[1];
sx q[1];
rz(-2.9225213) q[1];
rz(2.5147076) q[3];
sx q[3];
rz(-1.6520368) q[3];
sx q[3];
rz(2.8089942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6959186) q[2];
sx q[2];
rz(-0.58738223) q[2];
sx q[2];
rz(-0.43061259) q[2];
rz(-0.63498354) q[3];
sx q[3];
rz(-3.1228784) q[3];
sx q[3];
rz(-1.2682605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(1.1465313) q[0];
sx q[0];
rz(-0.54556161) q[0];
sx q[0];
rz(-1.5978285) q[0];
rz(0.68430463) q[1];
sx q[1];
rz(-1.6938208) q[1];
sx q[1];
rz(0.79944557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7921126) q[0];
sx q[0];
rz(-1.4941915) q[0];
sx q[0];
rz(1.576206) q[0];
rz(-0.075887738) q[2];
sx q[2];
rz(-2.2846756) q[2];
sx q[2];
rz(-2.3347428) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.909976) q[1];
sx q[1];
rz(-1.5607395) q[1];
sx q[1];
rz(0.46044989) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0897899) q[3];
sx q[3];
rz(-2.3481927) q[3];
sx q[3];
rz(-0.97360669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.53576175) q[2];
sx q[2];
rz(-0.83767086) q[2];
sx q[2];
rz(-2.3568995) q[2];
rz(0.11624087) q[3];
sx q[3];
rz(-1.616547) q[3];
sx q[3];
rz(-1.8893265) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5410974) q[0];
sx q[0];
rz(-1.9116115) q[0];
sx q[0];
rz(-1.0294718) q[0];
rz(0.12044278) q[1];
sx q[1];
rz(-1.30013) q[1];
sx q[1];
rz(-2.5708503) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.127205) q[0];
sx q[0];
rz(-0.61349166) q[0];
sx q[0];
rz(-0.86742371) q[0];
rz(-pi) q[1];
rz(1.458857) q[2];
sx q[2];
rz(-2.2862391) q[2];
sx q[2];
rz(-1.2067522) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2201669) q[1];
sx q[1];
rz(-1.1296185) q[1];
sx q[1];
rz(0.36018546) q[1];
rz(1.086497) q[3];
sx q[3];
rz(-0.37621337) q[3];
sx q[3];
rz(-2.0813326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8735147) q[2];
sx q[2];
rz(-2.2669078) q[2];
sx q[2];
rz(2.6178005) q[2];
rz(-1.9889471) q[3];
sx q[3];
rz(-0.58061424) q[3];
sx q[3];
rz(-1.7136278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6607587) q[0];
sx q[0];
rz(-1.9402215) q[0];
sx q[0];
rz(-2.7177366) q[0];
rz(-2.2747874) q[1];
sx q[1];
rz(-2.1173756) q[1];
sx q[1];
rz(-2.9383235) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05860672) q[0];
sx q[0];
rz(-1.5241677) q[0];
sx q[0];
rz(-2.4871422) q[0];
rz(-pi) q[1];
rz(0.31000455) q[2];
sx q[2];
rz(-1.5192724) q[2];
sx q[2];
rz(2.132694) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8228559) q[1];
sx q[1];
rz(-1.7137495) q[1];
sx q[1];
rz(-2.8566384) q[1];
rz(1.8762693) q[3];
sx q[3];
rz(-1.8169699) q[3];
sx q[3];
rz(2.9873029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.086143494) q[2];
sx q[2];
rz(-1.2148427) q[2];
sx q[2];
rz(-0.27935585) q[2];
rz(-1.2934359) q[3];
sx q[3];
rz(-1.3118298) q[3];
sx q[3];
rz(-1.9259341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9781037) q[0];
sx q[0];
rz(-0.80253974) q[0];
sx q[0];
rz(-1.7247024) q[0];
rz(0.92719999) q[1];
sx q[1];
rz(-1.5798774) q[1];
sx q[1];
rz(1.7318116) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5997084) q[0];
sx q[0];
rz(-2.4711631) q[0];
sx q[0];
rz(1.7267358) q[0];
rz(-pi) q[1];
rz(2.1883606) q[2];
sx q[2];
rz(-2.5698863) q[2];
sx q[2];
rz(2.1505411) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.77892196) q[1];
sx q[1];
rz(-0.92423981) q[1];
sx q[1];
rz(-2.126466) q[1];
rz(1.4019484) q[3];
sx q[3];
rz(-2.8956684) q[3];
sx q[3];
rz(-2.3000474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.37567821) q[2];
sx q[2];
rz(-1.3216852) q[2];
sx q[2];
rz(-0.970617) q[2];
rz(2.4961903) q[3];
sx q[3];
rz(-0.97964764) q[3];
sx q[3];
rz(-2.1360883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8476625) q[0];
sx q[0];
rz(-1.4959338) q[0];
sx q[0];
rz(-1.3069859) q[0];
rz(-0.38181276) q[1];
sx q[1];
rz(-1.3419071) q[1];
sx q[1];
rz(0.76795427) q[1];
rz(2.1099595) q[2];
sx q[2];
rz(-2.3066386) q[2];
sx q[2];
rz(2.5021449) q[2];
rz(-2.0313203) q[3];
sx q[3];
rz(-1.9721748) q[3];
sx q[3];
rz(1.0655793) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
