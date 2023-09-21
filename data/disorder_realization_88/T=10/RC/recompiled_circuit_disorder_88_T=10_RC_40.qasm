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
rz(2.8136301) q[0];
rz(-0.2344996) q[1];
sx q[1];
rz(-2.9405138) q[1];
sx q[1];
rz(-0.091436401) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8592005) q[0];
sx q[0];
rz(-0.99635591) q[0];
sx q[0];
rz(2.0496856) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7056866) q[2];
sx q[2];
rz(-1.3139259) q[2];
sx q[2];
rz(2.169662) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.5603197) q[1];
sx q[1];
rz(-0.25174192) q[1];
sx q[1];
rz(-1.2341577) q[1];
x q[2];
rz(1.8880635) q[3];
sx q[3];
rz(-0.61875611) q[3];
sx q[3];
rz(1.3878824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4771007) q[2];
sx q[2];
rz(-2.1680919) q[2];
sx q[2];
rz(-2.0155902) q[2];
rz(0.27515718) q[3];
sx q[3];
rz(-2.5313009) q[3];
sx q[3];
rz(-0.90308213) q[3];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(2.4480208) q[0];
rz(2.0454848) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(0.19031659) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41102558) q[0];
sx q[0];
rz(-1.9733631) q[0];
sx q[0];
rz(-0.21582027) q[0];
rz(-1.1778189) q[2];
sx q[2];
rz(-1.0699546) q[2];
sx q[2];
rz(-2.9289392) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.39365921) q[1];
sx q[1];
rz(-0.84888443) q[1];
sx q[1];
rz(0.02971239) q[1];
rz(-pi) q[2];
rz(-0.72782794) q[3];
sx q[3];
rz(-1.2065294) q[3];
sx q[3];
rz(1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6341614) q[2];
sx q[2];
rz(-2.4298411) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(0.1427342) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(0.90863168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74293566) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(2.3213342) q[0];
rz(-0.29207686) q[1];
sx q[1];
rz(-2.0675817) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0322198) q[0];
sx q[0];
rz(-1.4782227) q[0];
sx q[0];
rz(2.5412482) q[0];
rz(2.8448361) q[2];
sx q[2];
rz(-2.9251758) q[2];
sx q[2];
rz(1.4518472) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6229912) q[1];
sx q[1];
rz(-0.90497436) q[1];
sx q[1];
rz(-1.4856505) q[1];
x q[2];
rz(-0.91979023) q[3];
sx q[3];
rz(-1.3295104) q[3];
sx q[3];
rz(1.1834708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-1.0504861) q[2];
sx q[2];
rz(1.2134264) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(-0.081610672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7320025) q[0];
sx q[0];
rz(-1.859917) q[0];
sx q[0];
rz(2.9555292) q[0];
rz(-0.20446725) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0766749) q[0];
sx q[0];
rz(-3.0041822) q[0];
sx q[0];
rz(-1.8924367) q[0];
x q[1];
rz(-0.16764955) q[2];
sx q[2];
rz(-1.8561346) q[2];
sx q[2];
rz(-0.97380762) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.90721005) q[1];
sx q[1];
rz(-2.130059) q[1];
sx q[1];
rz(-0.84749605) q[1];
x q[2];
rz(-1.3912195) q[3];
sx q[3];
rz(-1.5571801) q[3];
sx q[3];
rz(0.40819685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2356448) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17386757) q[0];
sx q[0];
rz(-2.6350832) q[0];
sx q[0];
rz(2.2633973) q[0];
rz(-1.8163266) q[1];
sx q[1];
rz(-1.7293431) q[1];
sx q[1];
rz(-1.426288) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844937) q[0];
sx q[0];
rz(-1.6528659) q[0];
sx q[0];
rz(2.2280072) q[0];
rz(-pi) q[1];
x q[1];
rz(0.23767383) q[2];
sx q[2];
rz(-1.875068) q[2];
sx q[2];
rz(1.8397457) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1713472) q[1];
sx q[1];
rz(-2.4023224) q[1];
sx q[1];
rz(-1.7319748) q[1];
rz(-0.77107314) q[3];
sx q[3];
rz(-1.3372476) q[3];
sx q[3];
rz(2.0296823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2720126) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(0.041794725) q[2];
rz(3.0801008) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(2.7048236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7866661) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(1.0304931) q[0];
rz(-0.73973918) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(-0.57156634) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8067236) q[0];
sx q[0];
rz(-0.81672251) q[0];
sx q[0];
rz(-0.73210324) q[0];
x q[1];
rz(1.7320485) q[2];
sx q[2];
rz(-2.469392) q[2];
sx q[2];
rz(1.4816928) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.4042582) q[1];
sx q[1];
rz(-1.3306276) q[1];
sx q[1];
rz(2.0391697) q[1];
rz(-pi) q[2];
x q[2];
rz(0.70116331) q[3];
sx q[3];
rz(-0.94651604) q[3];
sx q[3];
rz(-3.0091156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.968154) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(2.5773933) q[2];
rz(-0.12600222) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(2.8469767) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0512222) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(-2.4293161) q[0];
rz(-2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-0.66551048) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3411361) q[0];
sx q[0];
rz(-1.8130214) q[0];
sx q[0];
rz(3.1405297) q[0];
rz(-pi) q[1];
x q[1];
rz(2.268157) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(0.21542491) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.682714) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(-1.1125803) q[1];
rz(-pi) q[2];
rz(-0.14466488) q[3];
sx q[3];
rz(-1.33107) q[3];
sx q[3];
rz(-0.02558115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.9843288) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(-3.026399) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(-2.9490024) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0916864) q[0];
sx q[0];
rz(-2.005907) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(-0.62675369) q[1];
sx q[1];
rz(-2.125506) q[1];
sx q[1];
rz(2.802882) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3074293) q[0];
sx q[0];
rz(-1.2273664) q[0];
sx q[0];
rz(-2.9817581) q[0];
rz(-pi) q[1];
x q[1];
rz(2.68967) q[2];
sx q[2];
rz(-1.3616614) q[2];
sx q[2];
rz(2.506633) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6845477) q[1];
sx q[1];
rz(-1.0984048) q[1];
sx q[1];
rz(-3.0169675) q[1];
rz(-pi) q[2];
rz(-1.2440153) q[3];
sx q[3];
rz(-1.3076926) q[3];
sx q[3];
rz(-0.14740482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-0.70043606) q[2];
rz(-0.8979848) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.609628) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(-2.530653) q[1];
sx q[1];
rz(-1.2555723) q[1];
sx q[1];
rz(3.0019965) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0567719) q[0];
sx q[0];
rz(-0.064282566) q[0];
sx q[0];
rz(-1.8875185) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81462461) q[2];
sx q[2];
rz(-1.6135718) q[2];
sx q[2];
rz(2.8551561) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4715251) q[1];
sx q[1];
rz(-1.0499665) q[1];
sx q[1];
rz(1.6478369) q[1];
x q[2];
rz(0.052501909) q[3];
sx q[3];
rz(-0.92696654) q[3];
sx q[3];
rz(2.1569463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5247941) q[2];
sx q[2];
rz(-2.1308265) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(2.3838499) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(3.0537135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8452633) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(0.18558003) q[0];
rz(1.0962076) q[1];
sx q[1];
rz(-2.9269693) q[1];
sx q[1];
rz(1.6569482) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34996841) q[0];
sx q[0];
rz(-2.0098643) q[0];
sx q[0];
rz(-2.8429948) q[0];
rz(-0.89016685) q[2];
sx q[2];
rz(-0.3496799) q[2];
sx q[2];
rz(0.65469757) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2828196) q[1];
sx q[1];
rz(-0.38837896) q[1];
sx q[1];
rz(1.3057083) q[1];
x q[2];
rz(-2.2757169) q[3];
sx q[3];
rz(-1.5745592) q[3];
sx q[3];
rz(0.83434425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.2075656) q[2];
sx q[2];
rz(-0.88946122) q[2];
sx q[2];
rz(-0.55220848) q[2];
rz(2.3637555) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(0.51789969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8778397) q[0];
sx q[0];
rz(-1.629566) q[0];
sx q[0];
rz(-1.4019479) q[0];
rz(-2.9539625) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(-1.4245695) q[2];
sx q[2];
rz(-1.8478827) q[2];
sx q[2];
rz(-2.293496) q[2];
rz(0.68708146) q[3];
sx q[3];
rz(-2.1344746) q[3];
sx q[3];
rz(1.8070756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];