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
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(0.091436401) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592005) q[0];
sx q[0];
rz(-2.1452367) q[0];
sx q[0];
rz(-1.0919071) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47317998) q[2];
sx q[2];
rz(-0.28943974) q[2];
sx q[2];
rz(2.6602886) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.58127296) q[1];
sx q[1];
rz(-0.25174192) q[1];
sx q[1];
rz(1.2341577) q[1];
x q[2];
rz(-0.9760194) q[3];
sx q[3];
rz(-1.3888437) q[3];
sx q[3];
rz(2.6973157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4771007) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(2.0155902) q[2];
rz(-2.8664355) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1317516) q[0];
sx q[0];
rz(-2.4480902) q[0];
sx q[0];
rz(0.69357187) q[0];
rz(-2.0454848) q[1];
sx q[1];
rz(-0.98384905) q[1];
sx q[1];
rz(-2.9512761) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0424461) q[0];
sx q[0];
rz(-2.6876246) q[0];
sx q[0];
rz(1.1048261) q[0];
rz(2.6066577) q[2];
sx q[2];
rz(-1.2282279) q[2];
sx q[2];
rz(-1.554622) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.1574993) q[1];
sx q[1];
rz(-1.5484973) q[1];
sx q[1];
rz(0.84866546) q[1];
rz(-pi) q[2];
rz(-1.0986436) q[3];
sx q[3];
rz(-0.90001366) q[3];
sx q[3];
rz(-2.9428279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.50743121) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(1.921839) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-2.1295363) q[3];
sx q[3];
rz(2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74293566) q[0];
sx q[0];
rz(-1.0899028) q[0];
sx q[0];
rz(2.3213342) q[0];
rz(0.29207686) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67278636) q[0];
sx q[0];
rz(-0.60657036) q[0];
sx q[0];
rz(-0.16288217) q[0];
rz(0.2072316) q[2];
sx q[2];
rz(-1.5079632) q[2];
sx q[2];
rz(-0.40916967) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6229912) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(-1.6559421) q[1];
rz(-pi) q[2];
rz(2.8415801) q[3];
sx q[3];
rz(-2.1999151) q[3];
sx q[3];
rz(-2.574207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(-1.2134264) q[2];
rz(-2.976867) q[3];
sx q[3];
rz(-0.89403331) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.40959013) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(-2.9555292) q[0];
rz(-0.20446725) q[1];
sx q[1];
rz(-0.47195131) q[1];
sx q[1];
rz(-1.8444555) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0766749) q[0];
sx q[0];
rz(-0.13741048) q[0];
sx q[0];
rz(1.2491559) q[0];
rz(-pi) q[1];
rz(1.2816216) q[2];
sx q[2];
rz(-1.7316069) q[2];
sx q[2];
rz(-0.5493872) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.12249494) q[1];
sx q[1];
rz(-2.2593453) q[1];
sx q[1];
rz(-0.81329878) q[1];
x q[2];
rz(0.013838776) q[3];
sx q[3];
rz(-1.3912364) q[3];
sx q[3];
rz(-1.9765215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2356448) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-2.2223991) q[2];
rz(-2.8202608) q[3];
sx q[3];
rz(-2.0733757) q[3];
sx q[3];
rz(1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
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
rz(-1.7293431) q[1];
sx q[1];
rz(1.426288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.250524) q[0];
sx q[0];
rz(-2.2254125) q[0];
sx q[0];
rz(3.0380681) q[0];
x q[1];
rz(-1.2582785) q[2];
sx q[2];
rz(-1.3442355) q[2];
sx q[2];
rz(0.3414008) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1713472) q[1];
sx q[1];
rz(-0.73927021) q[1];
sx q[1];
rz(1.4096178) q[1];
x q[2];
rz(-2.812643) q[3];
sx q[3];
rz(-2.3429686) q[3];
sx q[3];
rz(0.22508276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8695801) q[2];
sx q[2];
rz(-2*pi/13) q[2];
sx q[2];
rz(3.0997979) q[2];
rz(3.0801008) q[3];
sx q[3];
rz(-1.0890591) q[3];
sx q[3];
rz(0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3549266) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(-1.0304931) q[0];
rz(0.73973918) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(0.57156634) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33486903) q[0];
sx q[0];
rz(-0.81672251) q[0];
sx q[0];
rz(-0.73210324) q[0];
rz(-pi) q[1];
rz(-0.90494855) q[2];
sx q[2];
rz(-1.6709423) q[2];
sx q[2];
rz(0.03749321) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7476269) q[1];
sx q[1];
rz(-0.52226258) q[1];
sx q[1];
rz(2.0678492) q[1];
rz(-pi) q[2];
rz(2.3267641) q[3];
sx q[3];
rz(-1.0199254) q[3];
sx q[3];
rz(-1.2448454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.17343865) q[2];
sx q[2];
rz(-2.3996694) q[2];
sx q[2];
rz(0.56419939) q[2];
rz(0.12600222) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(-2.8469767) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0903704) q[0];
sx q[0];
rz(-0.19989656) q[0];
sx q[0];
rz(2.4293161) q[0];
rz(2.6157216) q[1];
sx q[1];
rz(-2.7253175) q[1];
sx q[1];
rz(-2.4760822) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3411361) q[0];
sx q[0];
rz(-1.8130214) q[0];
sx q[0];
rz(-0.0010629396) q[0];
rz(-0.14851103) q[2];
sx q[2];
rz(-2.2627137) q[2];
sx q[2];
rz(1.4505381) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1674541) q[1];
sx q[1];
rz(-1.1145089) q[1];
sx q[1];
rz(-0.2445226) q[1];
x q[2];
rz(0.14466488) q[3];
sx q[3];
rz(-1.33107) q[3];
sx q[3];
rz(-3.1160115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.15726382) q[2];
sx q[2];
rz(-2.4815346) q[2];
sx q[2];
rz(-1.7187913) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-0.51518232) q[3];
sx q[3];
rz(2.9490024) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.049906235) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(-2.9507622) q[0];
rz(-0.62675369) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(0.33871067) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8341634) q[0];
sx q[0];
rz(-1.9142262) q[0];
sx q[0];
rz(-0.15983454) q[0];
rz(-1.8024826) q[2];
sx q[2];
rz(-2.0121687) q[2];
sx q[2];
rz(1.0362792) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.72570669) q[1];
sx q[1];
rz(-0.48735122) q[1];
sx q[1];
rz(-1.8094256) q[1];
rz(-0.2770822) q[3];
sx q[3];
rz(-1.8859366) q[3];
sx q[3];
rz(-1.3354697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(0.70043606) q[2];
rz(-2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(-2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53196466) q[0];
sx q[0];
rz(-0.48159972) q[0];
sx q[0];
rz(-0.65834808) q[0];
rz(-0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-0.13959612) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4021554) q[0];
sx q[0];
rz(-1.5097152) q[0];
sx q[0];
rz(0.020045965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81462461) q[2];
sx q[2];
rz(-1.5280208) q[2];
sx q[2];
rz(2.8551561) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3176206) q[1];
sx q[1];
rz(-2.6156153) q[1];
sx q[1];
rz(0.13336639) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6406052) q[3];
sx q[3];
rz(-2.4959292) q[3];
sx q[3];
rz(1.071969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(-2.3838499) q[3];
sx q[3];
rz(-0.38882935) q[3];
sx q[3];
rz(3.0537135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963294) q[0];
sx q[0];
rz(-1.0110649) q[0];
sx q[0];
rz(0.18558003) q[0];
rz(2.045385) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(1.6569482) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.164924) q[0];
sx q[0];
rz(-2.6161368) q[0];
sx q[0];
rz(-2.1303961) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89016685) q[2];
sx q[2];
rz(-2.7919127) q[2];
sx q[2];
rz(-2.4868951) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.1440925) q[1];
sx q[1];
rz(-1.9449241) q[1];
sx q[1];
rz(0.10679006) q[1];
rz(2.2757169) q[3];
sx q[3];
rz(-1.5745592) q[3];
sx q[3];
rz(-0.83434425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.93402702) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(-2.3637555) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(1.4245695) q[2];
sx q[2];
rz(-1.29371) q[2];
sx q[2];
rz(0.84809662) q[2];
rz(0.88541661) q[3];
sx q[3];
rz(-2.1366742) q[3];
sx q[3];
rz(-2.4921806) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
