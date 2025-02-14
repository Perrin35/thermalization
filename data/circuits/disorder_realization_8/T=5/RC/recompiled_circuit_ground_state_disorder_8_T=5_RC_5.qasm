OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.79787624) q[0];
sx q[0];
rz(-0.98348445) q[0];
sx q[0];
rz(2.9498192) q[0];
rz(-2.6311488) q[1];
sx q[1];
rz(-1.3671083) q[1];
sx q[1];
rz(-1.7459315) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5001427) q[0];
sx q[0];
rz(-3.1063188) q[0];
sx q[0];
rz(-2.7462237) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7566827) q[2];
sx q[2];
rz(-1.7966401) q[2];
sx q[2];
rz(-2.5511044) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6916445) q[1];
sx q[1];
rz(-0.80300036) q[1];
sx q[1];
rz(-1.8464441) q[1];
rz(-1.9463948) q[3];
sx q[3];
rz(-1.7197242) q[3];
sx q[3];
rz(-2.4747839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50769606) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(-0.92602777) q[2];
rz(-1.5121459) q[3];
sx q[3];
rz(-0.67073268) q[3];
sx q[3];
rz(-1.9984455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40874261) q[0];
sx q[0];
rz(-0.73246211) q[0];
sx q[0];
rz(2.9090885) q[0];
rz(-1.1393503) q[1];
sx q[1];
rz(-1.573223) q[1];
sx q[1];
rz(0.92612902) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0588603) q[0];
sx q[0];
rz(-1.7960567) q[0];
sx q[0];
rz(0.20998574) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.861169) q[2];
sx q[2];
rz(-0.33992514) q[2];
sx q[2];
rz(-2.6715793) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.35891102) q[1];
sx q[1];
rz(-0.57702438) q[1];
sx q[1];
rz(0.38459528) q[1];
rz(-pi) q[2];
rz(-0.6517209) q[3];
sx q[3];
rz(-2.4656396) q[3];
sx q[3];
rz(-0.60871688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.53655475) q[2];
sx q[2];
rz(-1.522541) q[2];
sx q[2];
rz(2.3489595) q[2];
rz(0.41147453) q[3];
sx q[3];
rz(-0.94235197) q[3];
sx q[3];
rz(2.3365432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8234392) q[0];
sx q[0];
rz(-2.1397488) q[0];
sx q[0];
rz(-0.26082984) q[0];
rz(2.8584495) q[1];
sx q[1];
rz(-1.6915551) q[1];
sx q[1];
rz(2.0859437) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973141) q[0];
sx q[0];
rz(-1.3998654) q[0];
sx q[0];
rz(-2.2173082) q[0];
x q[1];
rz(-1.5670384) q[2];
sx q[2];
rz(-2.2039287) q[2];
sx q[2];
rz(0.95852493) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.187589) q[1];
sx q[1];
rz(-0.84792811) q[1];
sx q[1];
rz(1.4299117) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.061402873) q[3];
sx q[3];
rz(-2.0656584) q[3];
sx q[3];
rz(1.7375378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5259214) q[2];
sx q[2];
rz(-2.0522223) q[2];
sx q[2];
rz(-1.8833359) q[2];
rz(1.0387756) q[3];
sx q[3];
rz(-0.98826161) q[3];
sx q[3];
rz(1.725215) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1309758) q[0];
sx q[0];
rz(-1.6058141) q[0];
sx q[0];
rz(-0.11949874) q[0];
rz(-0.15826982) q[1];
sx q[1];
rz(-0.4522849) q[1];
sx q[1];
rz(1.4403042) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8058469) q[0];
sx q[0];
rz(-2.3290538) q[0];
sx q[0];
rz(2.0772019) q[0];
rz(-3.0624448) q[2];
sx q[2];
rz(-1.9047577) q[2];
sx q[2];
rz(-2.3862267) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0021683) q[1];
sx q[1];
rz(-1.0989215) q[1];
sx q[1];
rz(1.2743239) q[1];
rz(-2.1705008) q[3];
sx q[3];
rz(-2.5673742) q[3];
sx q[3];
rz(2.3410839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4629472) q[2];
sx q[2];
rz(-3.0061649) q[2];
sx q[2];
rz(0.43369183) q[2];
rz(-0.67484754) q[3];
sx q[3];
rz(-2.0516472) q[3];
sx q[3];
rz(1.1564144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-2.7334412) q[0];
sx q[0];
rz(-1.0557405) q[0];
sx q[0];
rz(-0.59930402) q[0];
rz(2.377548) q[1];
sx q[1];
rz(-1.0300449) q[1];
sx q[1];
rz(-2.5291671) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5500887) q[0];
sx q[0];
rz(-1.7180802) q[0];
sx q[0];
rz(0.96571897) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1407498) q[2];
sx q[2];
rz(-0.71686059) q[2];
sx q[2];
rz(1.5886024) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.45578411) q[1];
sx q[1];
rz(-2.1791885) q[1];
sx q[1];
rz(1.366282) q[1];
x q[2];
rz(-2.3633943) q[3];
sx q[3];
rz(-2.457387) q[3];
sx q[3];
rz(-2.9915031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6969488) q[2];
sx q[2];
rz(-2.3418043) q[2];
sx q[2];
rz(-2.6714228) q[2];
rz(-2.7759975) q[3];
sx q[3];
rz(-1.1919034) q[3];
sx q[3];
rz(-0.61075926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.325901) q[0];
sx q[0];
rz(-2.3453562) q[0];
sx q[0];
rz(-0.66194397) q[0];
rz(-3.0603307) q[1];
sx q[1];
rz(-2.6632402) q[1];
sx q[1];
rz(-3.001396) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.79665444) q[0];
sx q[0];
rz(-1.7895763) q[0];
sx q[0];
rz(-0.47056912) q[0];
rz(1.1344994) q[2];
sx q[2];
rz(-0.66777523) q[2];
sx q[2];
rz(-0.89950022) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6684456) q[1];
sx q[1];
rz(-0.47159099) q[1];
sx q[1];
rz(-1.5525329) q[1];
rz(-pi) q[2];
rz(0.34937691) q[3];
sx q[3];
rz(-1.3991068) q[3];
sx q[3];
rz(0.44595815) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7774272) q[2];
sx q[2];
rz(-1.3157996) q[2];
sx q[2];
rz(-2.8167456) q[2];
rz(0.15803629) q[3];
sx q[3];
rz(-2.7285748) q[3];
sx q[3];
rz(2.1260927) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(2.3196816) q[0];
sx q[0];
rz(-0.076229036) q[0];
sx q[0];
rz(0.88777375) q[0];
rz(2.7582788) q[1];
sx q[1];
rz(-2.2593468) q[1];
sx q[1];
rz(2.9820138) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6066642) q[0];
sx q[0];
rz(-1.3114358) q[0];
sx q[0];
rz(1.242437) q[0];
rz(2.331686) q[2];
sx q[2];
rz(-0.72945539) q[2];
sx q[2];
rz(-2.4021961) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3899052) q[1];
sx q[1];
rz(-1.4426822) q[1];
sx q[1];
rz(-2.7634578) q[1];
rz(-pi) q[2];
rz(1.4042312) q[3];
sx q[3];
rz(-1.8565531) q[3];
sx q[3];
rz(0.33657246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5923656) q[2];
sx q[2];
rz(-1.9840019) q[2];
sx q[2];
rz(1.7249736) q[2];
rz(-1.8727411) q[3];
sx q[3];
rz(-0.78059355) q[3];
sx q[3];
rz(-2.1835073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45605993) q[0];
sx q[0];
rz(-2.0262418) q[0];
sx q[0];
rz(-2.2400895) q[0];
rz(-0.69425663) q[1];
sx q[1];
rz(-1.6777104) q[1];
sx q[1];
rz(1.136397) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0944871) q[0];
sx q[0];
rz(-1.4102954) q[0];
sx q[0];
rz(-0.032139924) q[0];
x q[1];
rz(-1.198223) q[2];
sx q[2];
rz(-1.0508453) q[2];
sx q[2];
rz(-2.6432248) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7982805) q[1];
sx q[1];
rz(-1.5250051) q[1];
sx q[1];
rz(-2.2664323) q[1];
x q[2];
rz(-2.5707022) q[3];
sx q[3];
rz(-0.81557214) q[3];
sx q[3];
rz(-0.81787986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.14094341) q[2];
sx q[2];
rz(-1.8724172) q[2];
sx q[2];
rz(0.81306523) q[2];
rz(-0.55417577) q[3];
sx q[3];
rz(-1.7289836) q[3];
sx q[3];
rz(-0.36177844) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5799705) q[0];
sx q[0];
rz(-1.1142092) q[0];
sx q[0];
rz(-2.9680874) q[0];
rz(0.42501998) q[1];
sx q[1];
rz(-1.605502) q[1];
sx q[1];
rz(0.82957155) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7873142) q[0];
sx q[0];
rz(-1.715938) q[0];
sx q[0];
rz(0.73696846) q[0];
rz(-pi) q[1];
rz(0.7220975) q[2];
sx q[2];
rz(-0.061826454) q[2];
sx q[2];
rz(2.823148) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6704055) q[1];
sx q[1];
rz(-2.1758432) q[1];
sx q[1];
rz(1.3208645) q[1];
x q[2];
rz(0.33631992) q[3];
sx q[3];
rz(-1.631853) q[3];
sx q[3];
rz(-1.8674191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.35036626) q[2];
sx q[2];
rz(-1.841505) q[2];
sx q[2];
rz(-1.680797) q[2];
rz(2.1469877) q[3];
sx q[3];
rz(-2.1456783) q[3];
sx q[3];
rz(2.2554956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17722546) q[0];
sx q[0];
rz(-2.565964) q[0];
sx q[0];
rz(-0.10203578) q[0];
rz(2.521934) q[1];
sx q[1];
rz(-1.6435868) q[1];
sx q[1];
rz(2.5443351) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5383947) q[0];
sx q[0];
rz(-1.704633) q[0];
sx q[0];
rz(-0.60026818) q[0];
x q[1];
rz(1.4240392) q[2];
sx q[2];
rz(-1.2191311) q[2];
sx q[2];
rz(0.41837012) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.44587943) q[1];
sx q[1];
rz(-1.6076325) q[1];
sx q[1];
rz(2.1365154) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.064205211) q[3];
sx q[3];
rz(-0.47229015) q[3];
sx q[3];
rz(-0.20087584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9336046) q[2];
sx q[2];
rz(-2.7808069) q[2];
sx q[2];
rz(-2.2130845) q[2];
rz(1.1601296) q[3];
sx q[3];
rz(-1.7103633) q[3];
sx q[3];
rz(2.3742356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3930897) q[0];
sx q[0];
rz(-0.7001644) q[0];
sx q[0];
rz(-2.9099332) q[0];
rz(1.1275935) q[1];
sx q[1];
rz(-2.6078106) q[1];
sx q[1];
rz(-0.003905205) q[1];
rz(-1.4814346) q[2];
sx q[2];
rz(-2.4632774) q[2];
sx q[2];
rz(2.6954271) q[2];
rz(0.87455187) q[3];
sx q[3];
rz(-1.8378432) q[3];
sx q[3];
rz(1.7624553) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
