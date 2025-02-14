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
rz(-2.1498635) q[0];
sx q[0];
rz(-2.5340134) q[0];
sx q[0];
rz(2.1460331) q[0];
rz(0.52752703) q[1];
sx q[1];
rz(4.1726569) q[1];
sx q[1];
rz(9.9095681) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3178892) q[0];
sx q[0];
rz(-1.7021588) q[0];
sx q[0];
rz(-1.0571558) q[0];
x q[1];
rz(-2.0071335) q[2];
sx q[2];
rz(-0.84056652) q[2];
sx q[2];
rz(0.98524994) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5875549) q[1];
sx q[1];
rz(-1.5120737) q[1];
sx q[1];
rz(-0.96129932) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23956128) q[3];
sx q[3];
rz(-0.29920855) q[3];
sx q[3];
rz(-2.9394958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.506044) q[2];
sx q[2];
rz(-2.1940993) q[2];
sx q[2];
rz(0.49723899) q[2];
rz(0.55073589) q[3];
sx q[3];
rz(-0.68998718) q[3];
sx q[3];
rz(-2.0042888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-2.9841442) q[0];
sx q[0];
rz(-1.7406311) q[0];
sx q[0];
rz(2.3781811) q[0];
rz(-0.58452559) q[1];
sx q[1];
rz(-1.7598563) q[1];
sx q[1];
rz(1.0113299) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7147192) q[0];
sx q[0];
rz(-1.4397286) q[0];
sx q[0];
rz(3.0253973) q[0];
rz(-0.23638921) q[2];
sx q[2];
rz(-1.4589785) q[2];
sx q[2];
rz(-0.37416247) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3015131) q[1];
sx q[1];
rz(-2.7373145) q[1];
sx q[1];
rz(-2.1188645) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2143651) q[3];
sx q[3];
rz(-2.5303741) q[3];
sx q[3];
rz(-0.54655308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4804907) q[2];
sx q[2];
rz(-1.9384408) q[2];
sx q[2];
rz(1.0002381) q[2];
rz(-0.58146042) q[3];
sx q[3];
rz(-2.4088819) q[3];
sx q[3];
rz(-1.2015517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8345399) q[0];
sx q[0];
rz(-2.5426799) q[0];
sx q[0];
rz(-0.55759984) q[0];
rz(2.8343976) q[1];
sx q[1];
rz(-0.59978825) q[1];
sx q[1];
rz(-2.6934326) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4059999) q[0];
sx q[0];
rz(-1.4383398) q[0];
sx q[0];
rz(0.51901061) q[0];
rz(-1.9314172) q[2];
sx q[2];
rz(-0.54979392) q[2];
sx q[2];
rz(-2.2711693) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3002172) q[1];
sx q[1];
rz(-0.63134495) q[1];
sx q[1];
rz(1.6168827) q[1];
x q[2];
rz(-2.7072315) q[3];
sx q[3];
rz(-1.5944527) q[3];
sx q[3];
rz(-2.0659741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3192516) q[2];
sx q[2];
rz(-2.4288869) q[2];
sx q[2];
rz(2.5437497) q[2];
rz(0.44761014) q[3];
sx q[3];
rz(-1.2667344) q[3];
sx q[3];
rz(-3.083057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7287801) q[0];
sx q[0];
rz(-0.030981177) q[0];
sx q[0];
rz(2.431751) q[0];
rz(2.1624883) q[1];
sx q[1];
rz(-0.85936463) q[1];
sx q[1];
rz(-1.8202579) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13793547) q[0];
sx q[0];
rz(-2.4383579) q[0];
sx q[0];
rz(3.0556995) q[0];
rz(-0.15726451) q[2];
sx q[2];
rz(-0.7866106) q[2];
sx q[2];
rz(-0.32341126) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7497341) q[1];
sx q[1];
rz(-2.1987111) q[1];
sx q[1];
rz(2.8086752) q[1];
rz(-1.0274326) q[3];
sx q[3];
rz(-1.8413944) q[3];
sx q[3];
rz(2.4222514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9909624) q[2];
sx q[2];
rz(-1.4272828) q[2];
sx q[2];
rz(-0.70685753) q[2];
rz(1.3458716) q[3];
sx q[3];
rz(-1.462505) q[3];
sx q[3];
rz(-0.70485392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3020346) q[0];
sx q[0];
rz(-0.78080451) q[0];
sx q[0];
rz(-1.6395521) q[0];
rz(1.5976277) q[1];
sx q[1];
rz(-2.2815956) q[1];
sx q[1];
rz(1.647321) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65022138) q[0];
sx q[0];
rz(-1.5724535) q[0];
sx q[0];
rz(2.1412388) q[0];
rz(1.5635033) q[2];
sx q[2];
rz(-2.1661758) q[2];
sx q[2];
rz(-3.1356168) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0922904) q[1];
sx q[1];
rz(-0.43152025) q[1];
sx q[1];
rz(2.3445361) q[1];
rz(-pi) q[2];
rz(-0.022541209) q[3];
sx q[3];
rz(-2.1858741) q[3];
sx q[3];
rz(-0.32451567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.26022252) q[2];
sx q[2];
rz(-1.2726731) q[2];
sx q[2];
rz(-0.55848813) q[2];
rz(-0.50885606) q[3];
sx q[3];
rz(-1.5285834) q[3];
sx q[3];
rz(1.2671027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59638554) q[0];
sx q[0];
rz(-1.8983497) q[0];
sx q[0];
rz(-0.16045706) q[0];
rz(-2.4682553) q[1];
sx q[1];
rz(-1.9247749) q[1];
sx q[1];
rz(-2.014726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25803963) q[0];
sx q[0];
rz(-0.48012381) q[0];
sx q[0];
rz(-0.86999805) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8062988) q[2];
sx q[2];
rz(-0.79293434) q[2];
sx q[2];
rz(-1.060134) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7221795) q[1];
sx q[1];
rz(-0.99687885) q[1];
sx q[1];
rz(-0.24146059) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0288208) q[3];
sx q[3];
rz(-1.752697) q[3];
sx q[3];
rz(2.101909) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1274073) q[2];
sx q[2];
rz(-1.5421474) q[2];
sx q[2];
rz(0.7737774) q[2];
rz(-2.1524147) q[3];
sx q[3];
rz(-0.4072322) q[3];
sx q[3];
rz(2.0686843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2268552) q[0];
sx q[0];
rz(-1.4829153) q[0];
sx q[0];
rz(-2.3684655) q[0];
rz(-1.5559366) q[1];
sx q[1];
rz(-1.2890041) q[1];
sx q[1];
rz(1.1770491) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53922916) q[0];
sx q[0];
rz(-1.6022816) q[0];
sx q[0];
rz(-0.8336153) q[0];
x q[1];
rz(0.99780166) q[2];
sx q[2];
rz(-1.8353476) q[2];
sx q[2];
rz(-0.71392347) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.29876626) q[1];
sx q[1];
rz(-0.93888679) q[1];
sx q[1];
rz(0.34159391) q[1];
x q[2];
rz(0.1775976) q[3];
sx q[3];
rz(-1.7081722) q[3];
sx q[3];
rz(0.053973764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3561463) q[2];
sx q[2];
rz(-2.7768504) q[2];
sx q[2];
rz(2.6981603) q[2];
rz(-2.7555079) q[3];
sx q[3];
rz(-1.723879) q[3];
sx q[3];
rz(2.940787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3767553) q[0];
sx q[0];
rz(-1.5489464) q[0];
sx q[0];
rz(-3.1078597) q[0];
rz(-2.4434166) q[1];
sx q[1];
rz(-0.59711421) q[1];
sx q[1];
rz(2.4765292) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8830983) q[0];
sx q[0];
rz(-1.8332714) q[0];
sx q[0];
rz(1.1061944) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2970398) q[2];
sx q[2];
rz(-1.204927) q[2];
sx q[2];
rz(-1.6576115) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6463975) q[1];
sx q[1];
rz(-0.84719275) q[1];
sx q[1];
rz(-0.88314914) q[1];
x q[2];
rz(-3.0785776) q[3];
sx q[3];
rz(-2.821453) q[3];
sx q[3];
rz(0.053016114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.72208059) q[2];
sx q[2];
rz(-0.91700143) q[2];
sx q[2];
rz(1.7783995) q[2];
rz(0.52115399) q[3];
sx q[3];
rz(-1.3241973) q[3];
sx q[3];
rz(0.49191973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9687965) q[0];
sx q[0];
rz(-1.3461312) q[0];
sx q[0];
rz(-2.3222493) q[0];
rz(-0.68483812) q[1];
sx q[1];
rz(-1.2872773) q[1];
sx q[1];
rz(-0.12231621) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8259895) q[0];
sx q[0];
rz(-1.1513183) q[0];
sx q[0];
rz(-2.2537838) q[0];
rz(2.7610131) q[2];
sx q[2];
rz(-1.0099482) q[2];
sx q[2];
rz(2.3083936) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.12625369) q[1];
sx q[1];
rz(-2.1623134) q[1];
sx q[1];
rz(1.9638318) q[1];
rz(-pi) q[2];
rz(0.58845406) q[3];
sx q[3];
rz(-1.1752442) q[3];
sx q[3];
rz(-1.5992529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.704432) q[2];
sx q[2];
rz(-1.8081343) q[2];
sx q[2];
rz(-2.2740347) q[2];
rz(-1.7433172) q[3];
sx q[3];
rz(-0.38845348) q[3];
sx q[3];
rz(-0.5790264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45359465) q[0];
sx q[0];
rz(-1.9168251) q[0];
sx q[0];
rz(1.0892185) q[0];
rz(0.047529686) q[1];
sx q[1];
rz(-2.0952974) q[1];
sx q[1];
rz(-0.68738031) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4682161) q[0];
sx q[0];
rz(-0.49096732) q[0];
sx q[0];
rz(-1.6223905) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0105404) q[2];
sx q[2];
rz(-2.6446335) q[2];
sx q[2];
rz(0.69253773) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2694305) q[1];
sx q[1];
rz(-0.80842847) q[1];
sx q[1];
rz(2.4926659) q[1];
x q[2];
rz(-3.0484588) q[3];
sx q[3];
rz(-2.2054005) q[3];
sx q[3];
rz(-2.840691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.99531949) q[2];
sx q[2];
rz(-2.5226888) q[2];
sx q[2];
rz(-1.6590365) q[2];
rz(-2.4847374) q[3];
sx q[3];
rz(-1.992179) q[3];
sx q[3];
rz(-0.38203865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.252608) q[0];
sx q[0];
rz(-1.1532619) q[0];
sx q[0];
rz(1.1070195) q[0];
rz(0.57869115) q[1];
sx q[1];
rz(-2.3042669) q[1];
sx q[1];
rz(-0.28490983) q[1];
rz(-2.9215521) q[2];
sx q[2];
rz(-2.7773428) q[2];
sx q[2];
rz(-0.11014948) q[2];
rz(2.9978776) q[3];
sx q[3];
rz(-0.52142177) q[3];
sx q[3];
rz(1.8038255) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
