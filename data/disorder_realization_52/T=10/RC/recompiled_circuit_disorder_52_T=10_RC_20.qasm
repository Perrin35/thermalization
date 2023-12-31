OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.62587005) q[0];
sx q[0];
rz(6.8318879) q[0];
sx q[0];
rz(5.3988342) q[0];
rz(1.4305152) q[1];
sx q[1];
rz(-2.1880452) q[1];
sx q[1];
rz(1.5024827) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9031154) q[0];
sx q[0];
rz(-1.3912541) q[0];
sx q[0];
rz(1.205501) q[0];
x q[1];
rz(-1.5760954) q[2];
sx q[2];
rz(-1.0056579) q[2];
sx q[2];
rz(-2.0084755) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.2734387) q[1];
sx q[1];
rz(-1.8167129) q[1];
sx q[1];
rz(-1.6383365) q[1];
rz(2.2803335) q[3];
sx q[3];
rz(-1.2951295) q[3];
sx q[3];
rz(1.1112801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3551336) q[2];
sx q[2];
rz(-0.81374514) q[2];
sx q[2];
rz(0.65594977) q[2];
rz(1.2077228) q[3];
sx q[3];
rz(-1.9658807) q[3];
sx q[3];
rz(2.1470127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8939963) q[0];
sx q[0];
rz(-2.7212454) q[0];
sx q[0];
rz(-2.7080652) q[0];
rz(2.9128089) q[1];
sx q[1];
rz(-0.42963916) q[1];
sx q[1];
rz(3.1343592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318856) q[0];
sx q[0];
rz(-1.9623555) q[0];
sx q[0];
rz(-0.10086985) q[0];
rz(1.8859293) q[2];
sx q[2];
rz(-1.2909781) q[2];
sx q[2];
rz(-1.8880106) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4870287) q[1];
sx q[1];
rz(-1.9332758) q[1];
sx q[1];
rz(1.3351721) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.238027) q[3];
sx q[3];
rz(-1.1840608) q[3];
sx q[3];
rz(2.8500593) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5269512) q[2];
sx q[2];
rz(-2.3336637) q[2];
sx q[2];
rz(-2.4439404) q[2];
rz(-3.0200322) q[3];
sx q[3];
rz(-1.2391042) q[3];
sx q[3];
rz(-2.8377623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4678629) q[0];
sx q[0];
rz(-0.91600743) q[0];
sx q[0];
rz(-1.3695705) q[0];
rz(1.2415775) q[1];
sx q[1];
rz(-1.4135655) q[1];
sx q[1];
rz(-0.26161584) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8825892) q[0];
sx q[0];
rz(-3.1212174) q[0];
sx q[0];
rz(1.0783608) q[0];
x q[1];
rz(3.0171266) q[2];
sx q[2];
rz(-2.3111768) q[2];
sx q[2];
rz(2.894573) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.97754543) q[1];
sx q[1];
rz(-2.3512212) q[1];
sx q[1];
rz(0.17069823) q[1];
x q[2];
rz(0.51298489) q[3];
sx q[3];
rz(-1.5981711) q[3];
sx q[3];
rz(0.56752906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.6802784) q[2];
sx q[2];
rz(-1.6363982) q[2];
sx q[2];
rz(-0.20351163) q[2];
rz(2.2198548) q[3];
sx q[3];
rz(-1.8739871) q[3];
sx q[3];
rz(-2.8620499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97776425) q[0];
sx q[0];
rz(-1.5777359) q[0];
sx q[0];
rz(-1.571636) q[0];
rz(1.0034026) q[1];
sx q[1];
rz(-1.3137716) q[1];
sx q[1];
rz(-1.2483695) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49444775) q[0];
sx q[0];
rz(-1.6361423) q[0];
sx q[0];
rz(1.4739743) q[0];
rz(-0.7123956) q[2];
sx q[2];
rz(-0.27370307) q[2];
sx q[2];
rz(1.7143539) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.215938) q[1];
sx q[1];
rz(-0.92080322) q[1];
sx q[1];
rz(-1.710612) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.29420935) q[3];
sx q[3];
rz(-0.74511408) q[3];
sx q[3];
rz(1.7780768) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0288329) q[2];
sx q[2];
rz(-1.8303822) q[2];
sx q[2];
rz(0.83703414) q[2];
rz(1.933243) q[3];
sx q[3];
rz(-1.874606) q[3];
sx q[3];
rz(0.78554955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1355302) q[0];
sx q[0];
rz(-1.0536138) q[0];
sx q[0];
rz(-0.77520448) q[0];
rz(-0.40183055) q[1];
sx q[1];
rz(-0.95087516) q[1];
sx q[1];
rz(-0.90243375) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2941023) q[0];
sx q[0];
rz(-0.93203629) q[0];
sx q[0];
rz(-1.627648) q[0];
x q[1];
rz(2.6685153) q[2];
sx q[2];
rz(-2.6396751) q[2];
sx q[2];
rz(2.9830473) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8946998) q[1];
sx q[1];
rz(-2.1293853) q[1];
sx q[1];
rz(0.65319368) q[1];
rz(-pi) q[2];
rz(-1.8875214) q[3];
sx q[3];
rz(-1.5734908) q[3];
sx q[3];
rz(-2.5731784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19501413) q[2];
sx q[2];
rz(-2.6533551) q[2];
sx q[2];
rz(-1.1966594) q[2];
rz(-1.6992016) q[3];
sx q[3];
rz(-1.2714352) q[3];
sx q[3];
rz(-1.6825914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3301795) q[0];
sx q[0];
rz(-2.9646962) q[0];
sx q[0];
rz(0.50317558) q[0];
rz(1.4563837) q[1];
sx q[1];
rz(-2.0676985) q[1];
sx q[1];
rz(0.20176372) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0148894) q[0];
sx q[0];
rz(-0.12173437) q[0];
sx q[0];
rz(-0.99862167) q[0];
x q[1];
rz(-2.3927275) q[2];
sx q[2];
rz(-0.61486926) q[2];
sx q[2];
rz(0.90235898) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4469874) q[1];
sx q[1];
rz(-1.1519377) q[1];
sx q[1];
rz(-0.5254196) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4089036) q[3];
sx q[3];
rz(-1.2142039) q[3];
sx q[3];
rz(2.4458812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.21489828) q[2];
sx q[2];
rz(-1.639067) q[2];
sx q[2];
rz(-2.8743437) q[2];
rz(2.3184508) q[3];
sx q[3];
rz(-0.079113364) q[3];
sx q[3];
rz(-0.87583035) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.293752) q[0];
sx q[0];
rz(-0.1593312) q[0];
sx q[0];
rz(3.0840432) q[0];
rz(1.6607704) q[1];
sx q[1];
rz(-1.7819504) q[1];
sx q[1];
rz(2.1988791) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7563815) q[0];
sx q[0];
rz(-2.4009973) q[0];
sx q[0];
rz(2.5368607) q[0];
rz(1.7985293) q[2];
sx q[2];
rz(-1.7274389) q[2];
sx q[2];
rz(-1.1802955) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6616933) q[1];
sx q[1];
rz(-0.53240314) q[1];
sx q[1];
rz(-0.9049306) q[1];
x q[2];
rz(-1.058217) q[3];
sx q[3];
rz(-1.3430542) q[3];
sx q[3];
rz(0.77373576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0223579) q[2];
sx q[2];
rz(-2.0998349) q[2];
sx q[2];
rz(0.38267246) q[2];
rz(2.102397) q[3];
sx q[3];
rz(-1.305205) q[3];
sx q[3];
rz(0.97810811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4246178) q[0];
sx q[0];
rz(-3.0452403) q[0];
sx q[0];
rz(2.8714645) q[0];
rz(0.62942901) q[1];
sx q[1];
rz(-0.71297485) q[1];
sx q[1];
rz(0.28392917) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1536381) q[0];
sx q[0];
rz(-2.0001912) q[0];
sx q[0];
rz(-0.9149787) q[0];
x q[1];
rz(-2.9291199) q[2];
sx q[2];
rz(-2.5310235) q[2];
sx q[2];
rz(0.072710466) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.200951) q[1];
sx q[1];
rz(-1.1108228) q[1];
sx q[1];
rz(-0.56401395) q[1];
rz(-pi) q[2];
rz(3.1030032) q[3];
sx q[3];
rz(-0.66017294) q[3];
sx q[3];
rz(0.71036464) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5876864) q[2];
sx q[2];
rz(-0.17051414) q[2];
sx q[2];
rz(1.2109057) q[2];
rz(-2.8816913) q[3];
sx q[3];
rz(-0.62429684) q[3];
sx q[3];
rz(-0.62817812) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0652086) q[0];
sx q[0];
rz(-2.5807091) q[0];
sx q[0];
rz(0.22924766) q[0];
rz(0.30300888) q[1];
sx q[1];
rz(-1.3906994) q[1];
sx q[1];
rz(-1.680826) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54366771) q[0];
sx q[0];
rz(-2.857508) q[0];
sx q[0];
rz(0.062505917) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42713366) q[2];
sx q[2];
rz(-2.7286227) q[2];
sx q[2];
rz(2.4954748) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8742121) q[1];
sx q[1];
rz(-2.0955718) q[1];
sx q[1];
rz(-1.5449779) q[1];
rz(-pi) q[2];
rz(0.38735729) q[3];
sx q[3];
rz(-0.82327561) q[3];
sx q[3];
rz(2.2203317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.71172697) q[2];
sx q[2];
rz(-1.2167598) q[2];
sx q[2];
rz(2.4460068) q[2];
rz(-2.7097278) q[3];
sx q[3];
rz(-2.6769107) q[3];
sx q[3];
rz(2.4263884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.573695) q[0];
sx q[0];
rz(-1.9298113) q[0];
sx q[0];
rz(-0.33690548) q[0];
rz(2.9341872) q[1];
sx q[1];
rz(-1.0284871) q[1];
sx q[1];
rz(2.7609603) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5779553) q[0];
sx q[0];
rz(-1.6407688) q[0];
sx q[0];
rz(-1.3689343) q[0];
x q[1];
rz(0.62217525) q[2];
sx q[2];
rz(-2.2969349) q[2];
sx q[2];
rz(-0.13571339) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.8457348) q[1];
sx q[1];
rz(-1.0746135) q[1];
sx q[1];
rz(2.9198398) q[1];
rz(0.57772824) q[3];
sx q[3];
rz(-1.2168222) q[3];
sx q[3];
rz(0.83565328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1404861) q[2];
sx q[2];
rz(-1.9753186) q[2];
sx q[2];
rz(1.1364737) q[2];
rz(-3.100637) q[3];
sx q[3];
rz(-2.332873) q[3];
sx q[3];
rz(-1.3142746) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3363591) q[0];
sx q[0];
rz(-2.2362066) q[0];
sx q[0];
rz(2.8295828) q[0];
rz(-2.1144755) q[1];
sx q[1];
rz(-1.849091) q[1];
sx q[1];
rz(-1.0277933) q[1];
rz(-2.3958191) q[2];
sx q[2];
rz(-0.33200982) q[2];
sx q[2];
rz(0.40340323) q[2];
rz(-2.0102262) q[3];
sx q[3];
rz(-1.7527179) q[3];
sx q[3];
rz(-2.6408165) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
