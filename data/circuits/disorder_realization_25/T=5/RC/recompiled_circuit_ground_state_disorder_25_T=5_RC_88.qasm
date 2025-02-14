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
rz(1.9422148) q[0];
rz(-0.83130032) q[1];
sx q[1];
rz(4.755862) q[1];
sx q[1];
rz(8.9250467) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.840755) q[0];
sx q[0];
rz(-0.84512701) q[0];
sx q[0];
rz(0.27780224) q[0];
rz(-pi) q[1];
rz(-1.9316767) q[2];
sx q[2];
rz(-1.7583876) q[2];
sx q[2];
rz(-2.9071992) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.079165212) q[1];
sx q[1];
rz(-1.2794331) q[1];
sx q[1];
rz(-2.3645762) q[1];
rz(-pi) q[2];
rz(3.0930661) q[3];
sx q[3];
rz(-2.2884946) q[3];
sx q[3];
rz(-0.80503073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1685593) q[2];
sx q[2];
rz(-1.3848105) q[2];
sx q[2];
rz(0.92156571) q[2];
rz(1.5700995) q[3];
sx q[3];
rz(-2.470033) q[3];
sx q[3];
rz(-0.69935548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2498995) q[0];
sx q[0];
rz(-2.1136916) q[0];
sx q[0];
rz(1.4693042) q[0];
rz(1.4768614) q[1];
sx q[1];
rz(-1.3427443) q[1];
sx q[1];
rz(-0.5161759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0196386) q[0];
sx q[0];
rz(-2.7350744) q[0];
sx q[0];
rz(1.818114) q[0];
x q[1];
rz(-2.198368) q[2];
sx q[2];
rz(-0.4021968) q[2];
sx q[2];
rz(3.135596) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.055430977) q[1];
sx q[1];
rz(-1.3199116) q[1];
sx q[1];
rz(-1.4990663) q[1];
x q[2];
rz(1.0063051) q[3];
sx q[3];
rz(-2.3670475) q[3];
sx q[3];
rz(-1.2513127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.76253647) q[2];
sx q[2];
rz(-1.2006589) q[2];
sx q[2];
rz(-2.8442247) q[2];
rz(-2.6250046) q[3];
sx q[3];
rz(-1.9818431) q[3];
sx q[3];
rz(-0.62492257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49279889) q[0];
sx q[0];
rz(-2.6060217) q[0];
sx q[0];
rz(2.9597362) q[0];
rz(-0.44405538) q[1];
sx q[1];
rz(-2.2477138) q[1];
sx q[1];
rz(0.081238834) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67760175) q[0];
sx q[0];
rz(-0.39636546) q[0];
sx q[0];
rz(-0.074808077) q[0];
x q[1];
rz(0.88043682) q[2];
sx q[2];
rz(-1.7882344) q[2];
sx q[2];
rz(-1.7060929) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42535815) q[1];
sx q[1];
rz(-0.50941804) q[1];
sx q[1];
rz(1.1822027) q[1];
rz(2.131583) q[3];
sx q[3];
rz(-1.6688235) q[3];
sx q[3];
rz(1.7757701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2150725) q[2];
sx q[2];
rz(-2.7762065) q[2];
sx q[2];
rz(-0.60058769) q[2];
rz(1.8961204) q[3];
sx q[3];
rz(-2.5966817) q[3];
sx q[3];
rz(0.039552461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1323701) q[0];
sx q[0];
rz(-0.70766574) q[0];
sx q[0];
rz(-3.1170377) q[0];
rz(2.764616) q[1];
sx q[1];
rz(-1.8507277) q[1];
sx q[1];
rz(-2.7797508) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.637476) q[0];
sx q[0];
rz(-1.8808805) q[0];
sx q[0];
rz(-1.4092567) q[0];
x q[1];
rz(0.078388647) q[2];
sx q[2];
rz(-1.9447717) q[2];
sx q[2];
rz(2.0736935) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1805318) q[1];
sx q[1];
rz(-2.1994563) q[1];
sx q[1];
rz(1.791545) q[1];
x q[2];
rz(-1.0933769) q[3];
sx q[3];
rz(-1.1934115) q[3];
sx q[3];
rz(1.6536825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.61818608) q[2];
sx q[2];
rz(-1.8768825) q[2];
sx q[2];
rz(-2.3183909) q[2];
rz(-0.19436714) q[3];
sx q[3];
rz(-2.7394962) q[3];
sx q[3];
rz(-2.226734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.119656) q[0];
sx q[0];
rz(-1.6031665) q[0];
sx q[0];
rz(0.62171474) q[0];
rz(0.74553982) q[1];
sx q[1];
rz(-2.3108683) q[1];
sx q[1];
rz(3.0527557) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54387605) q[0];
sx q[0];
rz(-1.5041699) q[0];
sx q[0];
rz(-0.1275098) q[0];
x q[1];
rz(-2.2584469) q[2];
sx q[2];
rz(-2.7252203) q[2];
sx q[2];
rz(2.1379545) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93394132) q[1];
sx q[1];
rz(-0.7867032) q[1];
sx q[1];
rz(-0.3982597) q[1];
x q[2];
rz(-0.16377512) q[3];
sx q[3];
rz(-1.9036674) q[3];
sx q[3];
rz(2.6956357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2290153) q[2];
sx q[2];
rz(-1.6385767) q[2];
sx q[2];
rz(-1.5360606) q[2];
rz(-2.33365) q[3];
sx q[3];
rz(-0.84351051) q[3];
sx q[3];
rz(-1.9690751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.262893) q[0];
sx q[0];
rz(-1.4233339) q[0];
sx q[0];
rz(2.189157) q[0];
rz(1.3643422) q[1];
sx q[1];
rz(-1.9352501) q[1];
sx q[1];
rz(0.84699026) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1035396) q[0];
sx q[0];
rz(-2.969739) q[0];
sx q[0];
rz(2.553493) q[0];
rz(1.2405841) q[2];
sx q[2];
rz(-2.3366513) q[2];
sx q[2];
rz(3.0778468) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8678363) q[1];
sx q[1];
rz(-2.057192) q[1];
sx q[1];
rz(-2.0517906) q[1];
rz(-pi) q[2];
rz(-3.0019747) q[3];
sx q[3];
rz(-3.0373796) q[3];
sx q[3];
rz(-2.0914608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.18818894) q[2];
sx q[2];
rz(-2.0764515) q[2];
sx q[2];
rz(0.38189253) q[2];
rz(-1.1000819) q[3];
sx q[3];
rz(-0.81172687) q[3];
sx q[3];
rz(-0.41128099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
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
rz(-1.6604796) q[1];
sx q[1];
rz(-2.3752046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9123133) q[0];
sx q[0];
rz(-1.5361495) q[0];
sx q[0];
rz(-2.628792) q[0];
x q[1];
rz(0.42986912) q[2];
sx q[2];
rz(-2.1235222) q[2];
sx q[2];
rz(-1.6195219) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38098225) q[1];
sx q[1];
rz(-2.8462324) q[1];
sx q[1];
rz(2.9350314) q[1];
x q[2];
rz(2.4637632) q[3];
sx q[3];
rz(-1.7683709) q[3];
sx q[3];
rz(0.017700087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7327205) q[2];
sx q[2];
rz(-1.5343821) q[2];
sx q[2];
rz(-1.4917779) q[2];
rz(-1.07897) q[3];
sx q[3];
rz(-1.8337199) q[3];
sx q[3];
rz(0.61416793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.33477467) q[0];
sx q[0];
rz(-2.4914927) q[0];
sx q[0];
rz(-0.41710576) q[0];
rz(-0.053622309) q[1];
sx q[1];
rz(-1.5258748) q[1];
sx q[1];
rz(2.0448304) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3552719) q[0];
sx q[0];
rz(-1.4505761) q[0];
sx q[0];
rz(2.1646196) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.057889537) q[2];
sx q[2];
rz(-1.1107003) q[2];
sx q[2];
rz(-1.636508) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.027710304) q[1];
sx q[1];
rz(-1.7358297) q[1];
sx q[1];
rz(-1.5636958) q[1];
rz(0.67819579) q[3];
sx q[3];
rz(-1.6890397) q[3];
sx q[3];
rz(-0.54071301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9216448) q[2];
sx q[2];
rz(-1.3741477) q[2];
sx q[2];
rz(2.5448223) q[2];
rz(-0.88045949) q[3];
sx q[3];
rz(-1.4245278) q[3];
sx q[3];
rz(-2.5318291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45744687) q[0];
sx q[0];
rz(-2.5992751) q[0];
sx q[0];
rz(2.82161) q[0];
rz(0.34842247) q[1];
sx q[1];
rz(-2.6967144) q[1];
sx q[1];
rz(-2.6382823) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5849294) q[0];
sx q[0];
rz(-1.4787424) q[0];
sx q[0];
rz(1.7199055) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.50554363) q[2];
sx q[2];
rz(-2.2525666) q[2];
sx q[2];
rz(-0.56150061) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7585246) q[1];
sx q[1];
rz(-0.43979859) q[1];
sx q[1];
rz(-1.7563933) q[1];
rz(-pi) q[2];
rz(1.4473404) q[3];
sx q[3];
rz(-1.0528101) q[3];
sx q[3];
rz(0.65604612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8866426) q[2];
sx q[2];
rz(-2.6889668) q[2];
sx q[2];
rz(-2.5060999) q[2];
rz(-1.5358216) q[3];
sx q[3];
rz(-1.7255892) q[3];
sx q[3];
rz(-0.087285727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6814293) q[0];
sx q[0];
rz(-0.59578139) q[0];
sx q[0];
rz(1.3687362) q[0];
rz(-0.5303371) q[1];
sx q[1];
rz(-1.6531205) q[1];
sx q[1];
rz(2.368685) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0470974) q[0];
sx q[0];
rz(-2.2412934) q[0];
sx q[0];
rz(-2.0740202) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3488316) q[2];
sx q[2];
rz(-1.6429516) q[2];
sx q[2];
rz(2.4506086) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7557395) q[1];
sx q[1];
rz(-2.6645711) q[1];
sx q[1];
rz(1.2287089) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1201376) q[3];
sx q[3];
rz(-1.3067596) q[3];
sx q[3];
rz(2.9024189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3717926) q[2];
sx q[2];
rz(-0.95855203) q[2];
sx q[2];
rz(0.27251631) q[2];
rz(0.46485999) q[3];
sx q[3];
rz(-1.9404989) q[3];
sx q[3];
rz(-0.71729898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40806017) q[0];
sx q[0];
rz(-1.6049186) q[0];
sx q[0];
rz(-1.4419755) q[0];
rz(-1.6170665) q[1];
sx q[1];
rz(-0.99823141) q[1];
sx q[1];
rz(0.29481606) q[1];
rz(0.67666247) q[2];
sx q[2];
rz(-0.39702111) q[2];
sx q[2];
rz(0.71297356) q[2];
rz(-2.0154304) q[3];
sx q[3];
rz(-2.6585326) q[3];
sx q[3];
rz(-1.4654893) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
