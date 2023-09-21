OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8945583) q[0];
sx q[0];
rz(-1.2556827) q[0];
sx q[0];
rz(-0.32796252) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(-3.0501563) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5782195) q[0];
sx q[0];
rz(-1.1735998) q[0];
sx q[0];
rz(0.63011516) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.435906) q[2];
sx q[2];
rz(-1.3139259) q[2];
sx q[2];
rz(0.97193064) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5603197) q[1];
sx q[1];
rz(-0.25174192) q[1];
sx q[1];
rz(1.9074349) q[1];
rz(-pi) q[2];
rz(-2.1655733) q[3];
sx q[3];
rz(-1.3888437) q[3];
sx q[3];
rz(0.44427696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.4771007) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(-1.1260024) q[2];
rz(-0.27515718) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-0.90308213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0098410957) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(2.9512761) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8961401) q[0];
sx q[0];
rz(-1.7691233) q[0];
sx q[0];
rz(-1.9818927) q[0];
x q[1];
rz(1.9637738) q[2];
sx q[2];
rz(-2.071638) q[2];
sx q[2];
rz(-0.21265342) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.43860501) q[1];
sx q[1];
rz(-0.72241306) q[1];
sx q[1];
rz(1.6045251) q[1];
rz(-2.042949) q[3];
sx q[3];
rz(-2.241579) q[3];
sx q[3];
rz(-2.9428279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(1.2197536) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(0.90863168) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.398657) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(-2.3213342) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(-1.2480199) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6662212) q[0];
sx q[0];
rz(-0.9733805) q[0];
sx q[0];
rz(1.682838) q[0];
x q[1];
rz(-1.6349995) q[2];
sx q[2];
rz(-1.777613) q[2];
sx q[2];
rz(-1.1484255) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65590811) q[1];
sx q[1];
rz(-0.67042065) q[1];
sx q[1];
rz(0.1078492) q[1];
rz(2.8415801) q[3];
sx q[3];
rz(-2.1999151) q[3];
sx q[3];
rz(-2.574207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(1.9281663) q[2];
rz(0.16472566) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320025) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(-0.18606342) q[0];
rz(0.20446725) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.8444555) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.064917795) q[0];
sx q[0];
rz(-3.0041822) q[0];
sx q[0];
rz(1.2491559) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2816216) q[2];
sx q[2];
rz(-1.4099858) q[2];
sx q[2];
rz(-2.5922054) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.90721005) q[1];
sx q[1];
rz(-1.0115336) q[1];
sx q[1];
rz(-2.2940966) q[1];
rz(-pi) q[2];
x q[2];
rz(0.013838776) q[3];
sx q[3];
rz(-1.7503563) q[3];
sx q[3];
rz(-1.1650712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.90594784) q[2];
sx q[2];
rz(-2.3489958) q[2];
sx q[2];
rz(2.2223991) q[2];
rz(2.8202608) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(-1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(0.17386757) q[0];
sx q[0];
rz(-2.6350832) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(1.8163266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(1.7153046) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25709897) q[0];
sx q[0];
rz(-1.6528659) q[0];
sx q[0];
rz(0.9135855) q[0];
x q[1];
rz(2.9039188) q[2];
sx q[2];
rz(-1.875068) q[2];
sx q[2];
rz(-1.8397457) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6225699) q[1];
sx q[1];
rz(-1.4624603) q[1];
sx q[1];
rz(0.83801724) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8910847) q[3];
sx q[3];
rz(-2.3158145) q[3];
sx q[3];
rz(-0.68009963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2720126) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(3.0997979) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-0.085952856) q[0];
sx q[0];
rz(-1.0304931) q[0];
rz(0.73973918) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(-2.5700263) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58493462) q[0];
sx q[0];
rz(-0.99781636) q[0];
sx q[0];
rz(2.1893188) q[0];
x q[1];
rz(-1.4095441) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(1.6598998) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.9532277) q[1];
sx q[1];
rz(-2.0247012) q[1];
sx q[1];
rz(0.26785775) q[1];
rz(-2.3267641) q[3];
sx q[3];
rz(-1.0199254) q[3];
sx q[3];
rz(-1.8967472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.968154) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(3.0155904) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512222) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(-0.5258711) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(0.66551048) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3411361) q[0];
sx q[0];
rz(-1.3285713) q[0];
sx q[0];
rz(-0.0010629396) q[0];
x q[1];
rz(1.7475142) q[2];
sx q[2];
rz(-0.70509796) q[2];
sx q[2];
rz(1.9213898) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29384089) q[1];
sx q[1];
rz(-1.3517225) q[1];
sx q[1];
rz(-1.1024544) q[1];
rz(-2.9969278) q[3];
sx q[3];
rz(-1.33107) q[3];
sx q[3];
rz(0.02558115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9843288) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(-1.7187913) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(-0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0916864) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-0.19083047) q[0];
rz(-0.62675369) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(0.33871067) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3074293) q[0];
sx q[0];
rz(-1.9142262) q[0];
sx q[0];
rz(0.15983454) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3391101) q[2];
sx q[2];
rz(-1.1294239) q[2];
sx q[2];
rz(1.0362792) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6845477) q[1];
sx q[1];
rz(-1.0984048) q[1];
sx q[1];
rz(-0.12462516) q[1];
rz(-0.2770822) q[3];
sx q[3];
rz(-1.8859366) q[3];
sx q[3];
rz(1.806123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-2.4411566) q[2];
rz(-2.2436079) q[3];
sx q[3];
rz(-1.286641) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(3.0019965) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1698648) q[0];
sx q[0];
rz(-1.5507878) q[0];
sx q[0];
rz(1.5097029) q[0];
rz(-pi) q[1];
rz(3.082824) q[2];
sx q[2];
rz(-2.3261056) q[2];
sx q[2];
rz(-1.2440484) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.86233625) q[1];
sx q[1];
rz(-1.5039872) q[1];
sx q[1];
rz(-0.52211296) q[1];
rz(-1.6406052) q[3];
sx q[3];
rz(-2.4959292) q[3];
sx q[3];
rz(-2.0696236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5247941) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(-0.33111462) q[2];
rz(-2.3838499) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(-0.087879114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8452633) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(0.18558003) q[0];
rz(-1.0962076) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.4846444) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0508704) q[0];
sx q[0];
rz(-1.8403247) q[0];
sx q[0];
rz(-1.1140633) q[0];
rz(-pi) q[1];
rz(2.2514258) q[2];
sx q[2];
rz(-2.7919127) q[2];
sx q[2];
rz(-0.65469757) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.1440925) q[1];
sx q[1];
rz(-1.1966685) q[1];
sx q[1];
rz(-3.0348026) q[1];
rz(-0.86587571) q[3];
sx q[3];
rz(-1.5745592) q[3];
sx q[3];
rz(-0.83434425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.93402702) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(-0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26375297) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(0.18763018) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(0.4734584) q[2];
sx q[2];
rz(-0.31242328) q[2];
sx q[2];
rz(-1.7996126) q[2];
rz(-0.7836624) q[3];
sx q[3];
rz(-0.85859921) q[3];
sx q[3];
rz(2.8006299) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
