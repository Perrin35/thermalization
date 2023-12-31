OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.62087286) q[0];
sx q[0];
rz(-1.3735266) q[0];
sx q[0];
rz(1.5078478) q[0];
rz(-3.0942492) q[1];
sx q[1];
rz(-0.77818692) q[1];
sx q[1];
rz(-0.49931061) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0190174) q[0];
sx q[0];
rz(-2.8589773) q[0];
sx q[0];
rz(-2.9119888) q[0];
x q[1];
rz(-0.17272858) q[2];
sx q[2];
rz(-0.85731259) q[2];
sx q[2];
rz(-1.2781065) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5159113) q[1];
sx q[1];
rz(-1.7209007) q[1];
sx q[1];
rz(2.4019269) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3189089) q[3];
sx q[3];
rz(-3.0348274) q[3];
sx q[3];
rz(-1.3085384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9177861) q[2];
sx q[2];
rz(-0.97057682) q[2];
sx q[2];
rz(1.0144368) q[2];
rz(0.23400083) q[3];
sx q[3];
rz(-0.52105415) q[3];
sx q[3];
rz(-2.8570989) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6681799) q[0];
sx q[0];
rz(-1.4665335) q[0];
sx q[0];
rz(2.1372674) q[0];
rz(-1.6218119) q[1];
sx q[1];
rz(-0.92679778) q[1];
sx q[1];
rz(1.0027592) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28264499) q[0];
sx q[0];
rz(-0.98416057) q[0];
sx q[0];
rz(-2.100201) q[0];
rz(1.3833212) q[2];
sx q[2];
rz(-1.7406929) q[2];
sx q[2];
rz(-1.0152917) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7111295) q[1];
sx q[1];
rz(-1.7257479) q[1];
sx q[1];
rz(-2.1610545) q[1];
rz(-1.2172132) q[3];
sx q[3];
rz(-2.1401322) q[3];
sx q[3];
rz(0.47488892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26943794) q[2];
sx q[2];
rz(-2.1472011) q[2];
sx q[2];
rz(0.15094748) q[2];
rz(-2.7271467) q[3];
sx q[3];
rz(-2.5413385) q[3];
sx q[3];
rz(3.0533561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8931483) q[0];
sx q[0];
rz(-1.3131161) q[0];
sx q[0];
rz(0.77899581) q[0];
rz(2.7422854) q[1];
sx q[1];
rz(-1.8931959) q[1];
sx q[1];
rz(-0.88358203) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6813864) q[0];
sx q[0];
rz(-1.9708512) q[0];
sx q[0];
rz(-1.7494739) q[0];
x q[1];
rz(2.1951139) q[2];
sx q[2];
rz(-2.6636071) q[2];
sx q[2];
rz(1.1849272) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2391501) q[1];
sx q[1];
rz(-0.28581866) q[1];
sx q[1];
rz(1.6509389) q[1];
rz(-0.15150841) q[3];
sx q[3];
rz(-2.476859) q[3];
sx q[3];
rz(2.1745149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.3016004) q[2];
sx q[2];
rz(-1.2568544) q[2];
sx q[2];
rz(0.42923129) q[2];
rz(0.99003506) q[3];
sx q[3];
rz(-1.0777377) q[3];
sx q[3];
rz(-0.86301962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(0.4317959) q[0];
sx q[0];
rz(-1.3732095) q[0];
sx q[0];
rz(0.61169949) q[0];
rz(1.1071831) q[1];
sx q[1];
rz(-0.8586084) q[1];
sx q[1];
rz(-0.5501737) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.600425) q[0];
sx q[0];
rz(-0.74099243) q[0];
sx q[0];
rz(-3.0279972) q[0];
rz(-0.29454622) q[2];
sx q[2];
rz(-1.6606332) q[2];
sx q[2];
rz(-1.7812658) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.53766996) q[1];
sx q[1];
rz(-2.1426755) q[1];
sx q[1];
rz(3.0043383) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9228966) q[3];
sx q[3];
rz(-2.9069206) q[3];
sx q[3];
rz(0.72363879) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.52175534) q[2];
sx q[2];
rz(-0.48626128) q[2];
sx q[2];
rz(0.27553976) q[2];
rz(-3.0299305) q[3];
sx q[3];
rz(-1.2005946) q[3];
sx q[3];
rz(-2.6707941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4438641) q[0];
sx q[0];
rz(-2.0621018) q[0];
sx q[0];
rz(-2.2648947) q[0];
rz(-0.69119167) q[1];
sx q[1];
rz(-0.87379876) q[1];
sx q[1];
rz(-2.2263288) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0278496) q[0];
sx q[0];
rz(-0.62999524) q[0];
sx q[0];
rz(-1.6508474) q[0];
rz(-pi) q[1];
x q[1];
rz(0.077276262) q[2];
sx q[2];
rz(-2.6151267) q[2];
sx q[2];
rz(-0.72409814) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.76449672) q[1];
sx q[1];
rz(-0.15863523) q[1];
sx q[1];
rz(1.3678958) q[1];
x q[2];
rz(0.81766537) q[3];
sx q[3];
rz(-2.2125803) q[3];
sx q[3];
rz(-1.4828009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9203732) q[2];
sx q[2];
rz(-1.0125786) q[2];
sx q[2];
rz(0.4635703) q[2];
rz(2.5772337) q[3];
sx q[3];
rz(-0.99269358) q[3];
sx q[3];
rz(-0.92818964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
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
rz(-2.1148465) q[0];
sx q[0];
rz(-0.62018728) q[0];
sx q[0];
rz(-2.0625431) q[0];
rz(-2.5462529) q[1];
sx q[1];
rz(-0.7535615) q[1];
sx q[1];
rz(2.7450096) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.840344) q[0];
sx q[0];
rz(-0.52919555) q[0];
sx q[0];
rz(1.0521786) q[0];
rz(-pi) q[1];
rz(1.2601389) q[2];
sx q[2];
rz(-1.1144181) q[2];
sx q[2];
rz(3.0820456) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.10627667) q[1];
sx q[1];
rz(-2.4748487) q[1];
sx q[1];
rz(1.4186526) q[1];
x q[2];
rz(2.5572705) q[3];
sx q[3];
rz(-1.7752247) q[3];
sx q[3];
rz(2.1401329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.1534999) q[2];
sx q[2];
rz(-2.2160539) q[2];
sx q[2];
rz(-1.5552103) q[2];
rz(1.68613) q[3];
sx q[3];
rz(-0.60417914) q[3];
sx q[3];
rz(-0.82715183) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39847386) q[0];
sx q[0];
rz(-1.159659) q[0];
sx q[0];
rz(0.78480762) q[0];
rz(1.2706884) q[1];
sx q[1];
rz(-1.7653468) q[1];
sx q[1];
rz(1.126359) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2369909) q[0];
sx q[0];
rz(-2.8303574) q[0];
sx q[0];
rz(0.29445946) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9992505) q[2];
sx q[2];
rz(-1.0668313) q[2];
sx q[2];
rz(-3.1359283) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69648751) q[1];
sx q[1];
rz(-1.7587874) q[1];
sx q[1];
rz(-0.45447116) q[1];
rz(1.2231636) q[3];
sx q[3];
rz(-1.4813444) q[3];
sx q[3];
rz(-1.4060494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9324947) q[2];
sx q[2];
rz(-1.1065437) q[2];
sx q[2];
rz(-2.9879925) q[2];
rz(-0.30512729) q[3];
sx q[3];
rz(-2.025445) q[3];
sx q[3];
rz(1.7677527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7711733) q[0];
sx q[0];
rz(-0.53806794) q[0];
sx q[0];
rz(-3.0287108) q[0];
rz(1.0007292) q[1];
sx q[1];
rz(-0.74644867) q[1];
sx q[1];
rz(-0.032756068) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93452867) q[0];
sx q[0];
rz(-1.3158187) q[0];
sx q[0];
rz(1.9051139) q[0];
rz(-2.6997386) q[2];
sx q[2];
rz(-2.0896857) q[2];
sx q[2];
rz(-1.9372802) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0748464) q[1];
sx q[1];
rz(-2.7201338) q[1];
sx q[1];
rz(1.7154109) q[1];
rz(-pi) q[2];
rz(1.350358) q[3];
sx q[3];
rz(-0.38808295) q[3];
sx q[3];
rz(-2.2409852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.96413606) q[2];
sx q[2];
rz(-0.63825858) q[2];
sx q[2];
rz(-2.5194871) q[2];
rz(-1.9744251) q[3];
sx q[3];
rz(-1.111235) q[3];
sx q[3];
rz(0.39045236) q[3];
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
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099982925) q[0];
sx q[0];
rz(-0.5031302) q[0];
sx q[0];
rz(1.5266248) q[0];
rz(-2.408662) q[1];
sx q[1];
rz(-0.65299487) q[1];
sx q[1];
rz(-0.48535767) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5068685) q[0];
sx q[0];
rz(-1.6244495) q[0];
sx q[0];
rz(0.13454484) q[0];
rz(2.543407) q[2];
sx q[2];
rz(-1.4387812) q[2];
sx q[2];
rz(-2.8796632) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7945054) q[1];
sx q[1];
rz(-2.3350888) q[1];
sx q[1];
rz(-1.1361213) q[1];
x q[2];
rz(2.5076809) q[3];
sx q[3];
rz(-0.68448193) q[3];
sx q[3];
rz(-0.67542911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1099403) q[2];
sx q[2];
rz(-1.8376708) q[2];
sx q[2];
rz(0.1594485) q[2];
rz(1.738328) q[3];
sx q[3];
rz(-0.38968971) q[3];
sx q[3];
rz(3.1260417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6195246) q[0];
sx q[0];
rz(-2.5043026) q[0];
sx q[0];
rz(1.7074701) q[0];
rz(1.2592978) q[1];
sx q[1];
rz(-0.95016304) q[1];
sx q[1];
rz(-2.5433345) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3529417) q[0];
sx q[0];
rz(-1.7876248) q[0];
sx q[0];
rz(-1.3593332) q[0];
x q[1];
rz(-1.7190476) q[2];
sx q[2];
rz(-1.7808) q[2];
sx q[2];
rz(1.460618) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4926589) q[1];
sx q[1];
rz(-1.4760541) q[1];
sx q[1];
rz(-1.9365063) q[1];
rz(-1.7353021) q[3];
sx q[3];
rz(-1.4146155) q[3];
sx q[3];
rz(0.7848878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0593807) q[2];
sx q[2];
rz(-1.2450612) q[2];
sx q[2];
rz(-0.65199488) q[2];
rz(-0.59371289) q[3];
sx q[3];
rz(-1.9605325) q[3];
sx q[3];
rz(1.1317071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.839529) q[0];
sx q[0];
rz(-2.4840214) q[0];
sx q[0];
rz(-2.1485463) q[0];
rz(1.2364173) q[1];
sx q[1];
rz(-2.1724783) q[1];
sx q[1];
rz(1.9289) q[1];
rz(-2.9700206) q[2];
sx q[2];
rz(-2.5149859) q[2];
sx q[2];
rz(1.7563216) q[2];
rz(-2.0486352) q[3];
sx q[3];
rz(-2.408705) q[3];
sx q[3];
rz(-0.68791289) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
