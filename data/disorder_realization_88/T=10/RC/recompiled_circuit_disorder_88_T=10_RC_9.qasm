OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.24703439) q[0];
sx q[0];
rz(4.3972754) q[0];
sx q[0];
rz(9.7527405) q[0];
rz(2.9070931) q[1];
sx q[1];
rz(-0.20107888) q[1];
sx q[1];
rz(-3.0501563) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2823921) q[0];
sx q[0];
rz(-2.1452367) q[0];
sx q[0];
rz(1.0919071) q[0];
rz(-pi) q[1];
rz(1.435906) q[2];
sx q[2];
rz(-1.3139259) q[2];
sx q[2];
rz(-0.97193064) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.2135703) q[1];
sx q[1];
rz(-1.3334647) q[1];
sx q[1];
rz(-0.084753239) q[1];
rz(-pi) q[2];
rz(0.9760194) q[3];
sx q[3];
rz(-1.752749) q[3];
sx q[3];
rz(-0.44427696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.66449195) q[2];
sx q[2];
rz(-0.97350073) q[2];
sx q[2];
rz(-2.0155902) q[2];
rz(-2.8664355) q[3];
sx q[3];
rz(-0.61029172) q[3];
sx q[3];
rz(-2.2385105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0098410957) q[0];
sx q[0];
rz(-0.69350243) q[0];
sx q[0];
rz(-0.69357187) q[0];
rz(-1.0961078) q[1];
sx q[1];
rz(-2.1577436) q[1];
sx q[1];
rz(-2.9512761) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099146518) q[0];
sx q[0];
rz(-2.6876246) q[0];
sx q[0];
rz(2.0367665) q[0];
rz(-0.53493494) q[2];
sx q[2];
rz(-1.2282279) q[2];
sx q[2];
rz(1.5869706) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9840934) q[1];
sx q[1];
rz(-1.5484973) q[1];
sx q[1];
rz(0.84866546) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72782794) q[3];
sx q[3];
rz(-1.9350633) q[3];
sx q[3];
rz(1.4621853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6341614) q[2];
sx q[2];
rz(-0.71175152) q[2];
sx q[2];
rz(-1.921839) q[2];
rz(-2.9988585) q[3];
sx q[3];
rz(-1.0120564) q[3];
sx q[3];
rz(-2.232961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-0.74293566) q[0];
sx q[0];
rz(-2.0516899) q[0];
sx q[0];
rz(-0.8202585) q[0];
rz(-2.8495158) q[1];
sx q[1];
rz(-1.074011) q[1];
sx q[1];
rz(-1.8935727) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4753715) q[0];
sx q[0];
rz(-2.1682122) q[0];
sx q[0];
rz(-1.682838) q[0];
rz(-pi) q[1];
x q[1];
rz(0.29675656) q[2];
sx q[2];
rz(-2.9251758) q[2];
sx q[2];
rz(-1.4518472) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.51860147) q[1];
sx q[1];
rz(-2.2366183) q[1];
sx q[1];
rz(1.4856505) q[1];
rz(-2.2218024) q[3];
sx q[3];
rz(-1.8120822) q[3];
sx q[3];
rz(1.1834708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.8212006) q[2];
sx q[2];
rz(-2.0911066) q[2];
sx q[2];
rz(1.9281663) q[2];
rz(2.976867) q[3];
sx q[3];
rz(-2.2475593) q[3];
sx q[3];
rz(3.059982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40959013) q[0];
sx q[0];
rz(-1.2816757) q[0];
sx q[0];
rz(0.18606342) q[0];
rz(-2.9371254) q[1];
sx q[1];
rz(-2.6696413) q[1];
sx q[1];
rz(-1.8444555) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0766749) q[0];
sx q[0];
rz(-3.0041822) q[0];
sx q[0];
rz(-1.2491559) q[0];
rz(-1.859971) q[2];
sx q[2];
rz(-1.4099858) q[2];
sx q[2];
rz(-2.5922054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0399196) q[1];
sx q[1];
rz(-2.1664201) q[1];
sx q[1];
rz(-2.4458829) q[1];
rz(-pi) q[2];
rz(0.013838776) q[3];
sx q[3];
rz(-1.7503563) q[3];
sx q[3];
rz(1.9765215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2356448) q[2];
sx q[2];
rz(-0.79259688) q[2];
sx q[2];
rz(-0.91919351) q[2];
rz(-0.32133189) q[3];
sx q[3];
rz(-1.0682169) q[3];
sx q[3];
rz(1.2478158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9677251) q[0];
sx q[0];
rz(-0.50650948) q[0];
sx q[0];
rz(0.87819535) q[0];
rz(-1.325266) q[1];
sx q[1];
rz(-1.4122496) q[1];
sx q[1];
rz(1.7153046) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25709897) q[0];
sx q[0];
rz(-1.4887267) q[0];
sx q[0];
rz(-0.9135855) q[0];
x q[1];
rz(0.23767383) q[2];
sx q[2];
rz(-1.2665247) q[2];
sx q[2];
rz(-1.8397457) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1713472) q[1];
sx q[1];
rz(-0.73927021) q[1];
sx q[1];
rz(1.4096178) q[1];
x q[2];
rz(-0.32894965) q[3];
sx q[3];
rz(-2.3429686) q[3];
sx q[3];
rz(-0.22508276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2720126) q[2];
sx q[2];
rz(-11*pi/13) q[2];
sx q[2];
rz(-3.0997979) q[2];
rz(0.061491866) q[3];
sx q[3];
rz(-2.0525335) q[3];
sx q[3];
rz(0.4367691) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7866661) q[0];
sx q[0];
rz(-3.0556398) q[0];
sx q[0];
rz(2.1110995) q[0];
rz(0.73973918) q[1];
sx q[1];
rz(-1.61295) q[1];
sx q[1];
rz(0.57156634) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3540928) q[0];
sx q[0];
rz(-1.0618853) q[0];
sx q[0];
rz(-0.66977588) q[0];
rz(-pi) q[1];
rz(-1.7320485) q[2];
sx q[2];
rz(-0.67220062) q[2];
sx q[2];
rz(-1.6598998) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.188365) q[1];
sx q[1];
rz(-1.1168915) q[1];
sx q[1];
rz(-0.26785775) q[1];
rz(2.3011175) q[3];
sx q[3];
rz(-2.2395036) q[3];
sx q[3];
rz(-0.83276444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.968154) q[2];
sx q[2];
rz(-0.74192321) q[2];
sx q[2];
rz(-2.5773933) q[2];
rz(-0.12600222) q[3];
sx q[3];
rz(-1.4583476) q[3];
sx q[3];
rz(2.8469767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0512222) q[0];
sx q[0];
rz(-2.9416961) q[0];
sx q[0];
rz(2.4293161) q[0];
rz(-0.5258711) q[1];
sx q[1];
rz(-0.41627517) q[1];
sx q[1];
rz(2.4760822) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3367046) q[0];
sx q[0];
rz(-0.24222736) q[0];
sx q[0];
rz(1.5750984) q[0];
x q[1];
rz(-2.268157) q[2];
sx q[2];
rz(-1.4566112) q[2];
sx q[2];
rz(-0.21542491) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4588786) q[1];
sx q[1];
rz(-2.6280118) q[1];
sx q[1];
rz(1.1125803) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9969278) q[3];
sx q[3];
rz(-1.8105227) q[3];
sx q[3];
rz(0.02558115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.15726382) q[2];
sx q[2];
rz(-0.66005808) q[2];
sx q[2];
rz(-1.4228014) q[2];
rz(-0.11519365) q[3];
sx q[3];
rz(-2.6264103) q[3];
sx q[3];
rz(0.19259024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0916864) q[0];
sx q[0];
rz(-1.1356857) q[0];
sx q[0];
rz(2.9507622) q[0];
rz(-0.62675369) q[1];
sx q[1];
rz(-1.0160867) q[1];
sx q[1];
rz(-2.802882) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8341634) q[0];
sx q[0];
rz(-1.9142262) q[0];
sx q[0];
rz(0.15983454) q[0];
rz(-pi) q[1];
rz(2.6892002) q[2];
sx q[2];
rz(-0.49491844) q[2];
sx q[2];
rz(0.53168833) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0847816) q[1];
sx q[1];
rz(-1.4598795) q[1];
sx q[1];
rz(-2.0463498) q[1];
rz(-2.268928) q[3];
sx q[3];
rz(-2.7250395) q[3];
sx q[3];
rz(2.0779028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4954341) q[2];
sx q[2];
rz(-0.63445264) q[2];
sx q[2];
rz(-0.70043606) q[2];
rz(2.2436079) q[3];
sx q[3];
rz(-1.8549517) q[3];
sx q[3];
rz(2.9680796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.609628) q[0];
sx q[0];
rz(-2.6599929) q[0];
sx q[0];
rz(2.4832446) q[0];
rz(0.61093962) q[1];
sx q[1];
rz(-1.8860203) q[1];
sx q[1];
rz(-3.0019965) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1698648) q[0];
sx q[0];
rz(-1.5908049) q[0];
sx q[0];
rz(1.5097029) q[0];
x q[1];
rz(2.326968) q[2];
sx q[2];
rz(-1.5280208) q[2];
sx q[2];
rz(-0.2864366) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.3176206) q[1];
sx q[1];
rz(-2.6156153) q[1];
sx q[1];
rz(-0.13336639) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5009874) q[3];
sx q[3];
rz(-0.6456635) q[3];
sx q[3];
rz(-2.0696236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5247941) q[2];
sx q[2];
rz(-1.0107661) q[2];
sx q[2];
rz(-2.810478) q[2];
rz(0.75774276) q[3];
sx q[3];
rz(-2.7527633) q[3];
sx q[3];
rz(0.087879114) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2963294) q[0];
sx q[0];
rz(-2.1305278) q[0];
sx q[0];
rz(2.9560126) q[0];
rz(2.045385) q[1];
sx q[1];
rz(-0.21462333) q[1];
sx q[1];
rz(-1.4846444) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0907222) q[0];
sx q[0];
rz(-1.301268) q[0];
sx q[0];
rz(1.1140633) q[0];
rz(-pi) q[1];
rz(-1.2946285) q[2];
sx q[2];
rz(-1.3535) q[2];
sx q[2];
rz(1.5751788) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.6757322) q[1];
sx q[1];
rz(-1.670174) q[1];
sx q[1];
rz(-1.194721) q[1];
rz(-1.5766034) q[3];
sx q[3];
rz(-0.70492893) q[3];
sx q[3];
rz(0.73202902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2075656) q[2];
sx q[2];
rz(-2.2521314) q[2];
sx q[2];
rz(-2.5893842) q[2];
rz(0.77783716) q[3];
sx q[3];
rz(-0.85770291) q[3];
sx q[3];
rz(2.623693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778397) q[0];
sx q[0];
rz(-1.5120266) q[0];
sx q[0];
rz(1.7396447) q[0];
rz(0.18763018) q[1];
sx q[1];
rz(-1.7788806) q[1];
sx q[1];
rz(2.3685041) q[1];
rz(-1.7170231) q[2];
sx q[2];
rz(-1.29371) q[2];
sx q[2];
rz(0.84809662) q[2];
rz(-2.4545112) q[3];
sx q[3];
rz(-2.1344746) q[3];
sx q[3];
rz(1.8070756) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
