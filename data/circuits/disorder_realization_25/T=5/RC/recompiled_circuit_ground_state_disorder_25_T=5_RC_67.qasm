OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8018262) q[0];
sx q[0];
rz(-2.5751994) q[0];
sx q[0];
rz(-1.1993778) q[0];
rz(-0.83130032) q[1];
sx q[1];
rz(-1.5273233) q[1];
sx q[1];
rz(-0.49973127) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2466149) q[0];
sx q[0];
rz(-0.76789373) q[0];
sx q[0];
rz(1.270986) q[0];
x q[1];
rz(-1.0775282) q[2];
sx q[2];
rz(-2.7367711) q[2];
sx q[2];
rz(-2.2641393) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.079165212) q[1];
sx q[1];
rz(-1.8621596) q[1];
sx q[1];
rz(2.3645762) q[1];
x q[2];
rz(-1.6263032) q[3];
sx q[3];
rz(-0.71904564) q[3];
sx q[3];
rz(2.2628502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.1685593) q[2];
sx q[2];
rz(-1.3848105) q[2];
sx q[2];
rz(-2.2200269) q[2];
rz(-1.5714931) q[3];
sx q[3];
rz(-2.470033) q[3];
sx q[3];
rz(2.4422372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2498995) q[0];
sx q[0];
rz(-2.1136916) q[0];
sx q[0];
rz(1.4693042) q[0];
rz(-1.4768614) q[1];
sx q[1];
rz(-1.7988484) q[1];
sx q[1];
rz(2.6254168) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75137532) q[0];
sx q[0];
rz(-1.1773407) q[0];
sx q[0];
rz(0.1050001) q[0];
rz(-pi) q[1];
rz(0.24476972) q[2];
sx q[2];
rz(-1.8932071) q[2];
sx q[2];
rz(-0.66157326) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6440638) q[1];
sx q[1];
rz(-1.6402771) q[1];
sx q[1];
rz(-2.890088) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2616181) q[3];
sx q[3];
rz(-1.954292) q[3];
sx q[3];
rz(-0.10552204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76253647) q[2];
sx q[2];
rz(-1.9409337) q[2];
sx q[2];
rz(-2.8442247) q[2];
rz(0.51658806) q[3];
sx q[3];
rz(-1.9818431) q[3];
sx q[3];
rz(-0.62492257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.49279889) q[0];
sx q[0];
rz(-0.53557098) q[0];
sx q[0];
rz(0.18185644) q[0];
rz(2.6975373) q[1];
sx q[1];
rz(-2.2477138) q[1];
sx q[1];
rz(0.081238834) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96222197) q[0];
sx q[0];
rz(-1.5419382) q[0];
sx q[0];
rz(-2.7462237) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2368579) q[2];
sx q[2];
rz(-0.71839625) q[2];
sx q[2];
rz(-0.12018724) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7162345) q[1];
sx q[1];
rz(-0.50941804) q[1];
sx q[1];
rz(-1.1822027) q[1];
rz(-pi) q[2];
rz(1.0100097) q[3];
sx q[3];
rz(-1.4727691) q[3];
sx q[3];
rz(-1.3658226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2150725) q[2];
sx q[2];
rz(-2.7762065) q[2];
sx q[2];
rz(0.60058769) q[2];
rz(1.8961204) q[3];
sx q[3];
rz(-0.54491091) q[3];
sx q[3];
rz(3.1020402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1323701) q[0];
sx q[0];
rz(-0.70766574) q[0];
sx q[0];
rz(3.1170377) q[0];
rz(2.764616) q[1];
sx q[1];
rz(-1.2908649) q[1];
sx q[1];
rz(2.7797508) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013600969) q[0];
sx q[0];
rz(-0.34843081) q[0];
sx q[0];
rz(-2.6763787) q[0];
rz(1.76775) q[2];
sx q[2];
rz(-2.7598689) q[2];
sx q[2];
rz(-2.2854855) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3256073) q[1];
sx q[1];
rz(-2.4802842) q[1];
sx q[1];
rz(-2.8490752) q[1];
x q[2];
rz(-2.0482158) q[3];
sx q[3];
rz(-1.9481812) q[3];
sx q[3];
rz(1.6536825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.61818608) q[2];
sx q[2];
rz(-1.8768825) q[2];
sx q[2];
rz(2.3183909) q[2];
rz(2.9472255) q[3];
sx q[3];
rz(-0.40209642) q[3];
sx q[3];
rz(2.226734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.119656) q[0];
sx q[0];
rz(-1.5384262) q[0];
sx q[0];
rz(-0.62171474) q[0];
rz(0.74553982) q[1];
sx q[1];
rz(-2.3108683) q[1];
sx q[1];
rz(3.0527557) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5977166) q[0];
sx q[0];
rz(-1.5041699) q[0];
sx q[0];
rz(-3.0140829) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.88314573) q[2];
sx q[2];
rz(-2.7252203) q[2];
sx q[2];
rz(1.0036381) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.93394132) q[1];
sx q[1];
rz(-0.7867032) q[1];
sx q[1];
rz(0.3982597) q[1];
rz(-pi) q[2];
rz(-2.0114548) q[3];
sx q[3];
rz(-2.7719636) q[3];
sx q[3];
rz(-2.2274048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2290153) q[2];
sx q[2];
rz(-1.503016) q[2];
sx q[2];
rz(-1.5360606) q[2];
rz(-0.80794263) q[3];
sx q[3];
rz(-0.84351051) q[3];
sx q[3];
rz(-1.1725175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.262893) q[0];
sx q[0];
rz(-1.7182588) q[0];
sx q[0];
rz(-2.189157) q[0];
rz(-1.7772504) q[1];
sx q[1];
rz(-1.9352501) q[1];
sx q[1];
rz(-2.2946024) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95147371) q[0];
sx q[0];
rz(-1.4757809) q[0];
sx q[0];
rz(2.998178) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79361992) q[2];
sx q[2];
rz(-1.334903) q[2];
sx q[2];
rz(1.2737887) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27375633) q[1];
sx q[1];
rz(-2.057192) q[1];
sx q[1];
rz(2.0517906) q[1];
rz(1.5853508) q[3];
sx q[3];
rz(-1.4676009) q[3];
sx q[3];
rz(1.1905014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.18818894) q[2];
sx q[2];
rz(-2.0764515) q[2];
sx q[2];
rz(-2.7597001) q[2];
rz(2.0415107) q[3];
sx q[3];
rz(-0.81172687) q[3];
sx q[3];
rz(-0.41128099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5683658) q[0];
sx q[0];
rz(-1.8170284) q[0];
sx q[0];
rz(0.49760094) q[0];
rz(-0.35500232) q[1];
sx q[1];
rz(-1.4811131) q[1];
sx q[1];
rz(-0.76638806) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3610182) q[0];
sx q[0];
rz(-1.0583335) q[0];
sx q[0];
rz(1.6105525) q[0];
rz(-2.1649579) q[2];
sx q[2];
rz(-2.4554002) q[2];
sx q[2];
rz(2.2399069) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16535266) q[1];
sx q[1];
rz(-1.8596949) q[1];
sx q[1];
rz(-1.5084748) q[1];
rz(-pi) q[2];
rz(-1.3192407) q[3];
sx q[3];
rz(-0.90853229) q[3];
sx q[3];
rz(-1.4317625) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4088722) q[2];
sx q[2];
rz(-1.6072105) q[2];
sx q[2];
rz(1.4917779) q[2];
rz(1.07897) q[3];
sx q[3];
rz(-1.8337199) q[3];
sx q[3];
rz(2.5274247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33477467) q[0];
sx q[0];
rz(-0.65009999) q[0];
sx q[0];
rz(-0.41710576) q[0];
rz(-0.053622309) q[1];
sx q[1];
rz(-1.6157179) q[1];
sx q[1];
rz(-2.0448304) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2963144) q[0];
sx q[0];
rz(-0.9818378) q[0];
sx q[0];
rz(0.14473578) q[0];
rz(-1.1100329) q[2];
sx q[2];
rz(-1.5189324) q[2];
sx q[2];
rz(3.101609) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5442526) q[1];
sx q[1];
rz(-1.5778004) q[1];
sx q[1];
rz(-0.16503741) q[1];
rz(-pi) q[2];
rz(-1.4194059) q[3];
sx q[3];
rz(-2.243379) q[3];
sx q[3];
rz(-1.1248433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9216448) q[2];
sx q[2];
rz(-1.767445) q[2];
sx q[2];
rz(2.5448223) q[2];
rz(0.88045949) q[3];
sx q[3];
rz(-1.4245278) q[3];
sx q[3];
rz(2.5318291) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6841458) q[0];
sx q[0];
rz(-2.5992751) q[0];
sx q[0];
rz(0.31998262) q[0];
rz(0.34842247) q[1];
sx q[1];
rz(-2.6967144) q[1];
sx q[1];
rz(0.50331032) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1419174) q[0];
sx q[0];
rz(-1.4223232) q[0];
sx q[0];
rz(0.093080892) q[0];
x q[1];
rz(-1.032802) q[2];
sx q[2];
rz(-2.3176486) q[2];
sx q[2];
rz(1.8593553) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1784672) q[1];
sx q[1];
rz(-1.1390634) q[1];
sx q[1];
rz(0.086612062) q[1];
rz(-1.6942523) q[3];
sx q[3];
rz(-2.0887825) q[3];
sx q[3];
rz(-0.65604612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.2549501) q[2];
sx q[2];
rz(-0.4526259) q[2];
sx q[2];
rz(0.63549271) q[2];
rz(1.5358216) q[3];
sx q[3];
rz(-1.7255892) q[3];
sx q[3];
rz(0.087285727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6814293) q[0];
sx q[0];
rz(-0.59578139) q[0];
sx q[0];
rz(1.3687362) q[0];
rz(-2.6112556) q[1];
sx q[1];
rz(-1.6531205) q[1];
sx q[1];
rz(0.77290767) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0470974) q[0];
sx q[0];
rz(-0.90029923) q[0];
sx q[0];
rz(1.0675724) q[0];
rz(-pi) q[1];
rz(-0.79276104) q[2];
sx q[2];
rz(-1.498641) q[2];
sx q[2];
rz(-2.4506086) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0046282) q[1];
sx q[1];
rz(-1.1234979) q[1];
sx q[1];
rz(-2.9699259) q[1];
rz(1.0152994) q[3];
sx q[3];
rz(-0.51768983) q[3];
sx q[3];
rz(0.836984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3717926) q[2];
sx q[2];
rz(-2.1830406) q[2];
sx q[2];
rz(2.8690763) q[2];
rz(-2.6767327) q[3];
sx q[3];
rz(-1.9404989) q[3];
sx q[3];
rz(2.4242937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7335325) q[0];
sx q[0];
rz(-1.536674) q[0];
sx q[0];
rz(1.6996171) q[0];
rz(1.6170665) q[1];
sx q[1];
rz(-2.1433612) q[1];
sx q[1];
rz(-2.8467766) q[1];
rz(-2.4649302) q[2];
sx q[2];
rz(-0.39702111) q[2];
sx q[2];
rz(0.71297356) q[2];
rz(-2.0130264) q[3];
sx q[3];
rz(-1.7719405) q[3];
sx q[3];
rz(2.8476261) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
