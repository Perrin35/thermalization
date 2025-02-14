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
rz(-0.62491971) q[0];
sx q[0];
rz(-1.3345557) q[0];
sx q[0];
rz(-0.053243756) q[0];
rz(0.48625311) q[1];
sx q[1];
rz(-0.12737218) q[1];
sx q[1];
rz(1.706634) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93293982) q[0];
sx q[0];
rz(-0.31626748) q[0];
sx q[0];
rz(0.44751944) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7289624) q[2];
sx q[2];
rz(-2.8511471) q[2];
sx q[2];
rz(-1.7163139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.26027825) q[1];
sx q[1];
rz(-1.6242172) q[1];
sx q[1];
rz(-0.39398663) q[1];
rz(-0.94287431) q[3];
sx q[3];
rz(-0.74340313) q[3];
sx q[3];
rz(-1.8441895) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3367553) q[2];
sx q[2];
rz(-1.7982322) q[2];
sx q[2];
rz(-0.27377823) q[2];
rz(-1.9233507) q[3];
sx q[3];
rz(-0.67353407) q[3];
sx q[3];
rz(-0.28606733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.022920595) q[0];
sx q[0];
rz(-0.76347041) q[0];
sx q[0];
rz(-2.3619695) q[0];
rz(-0.66462213) q[1];
sx q[1];
rz(-1.9326991) q[1];
sx q[1];
rz(1.8962616) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1561403) q[0];
sx q[0];
rz(-0.13053556) q[0];
sx q[0];
rz(-2.0640316) q[0];
x q[1];
rz(-1.4636366) q[2];
sx q[2];
rz(-2.4273211) q[2];
sx q[2];
rz(0.34253866) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.14549637) q[1];
sx q[1];
rz(-0.22095535) q[1];
sx q[1];
rz(-1.2902106) q[1];
x q[2];
rz(-0.20540463) q[3];
sx q[3];
rz(-1.0213791) q[3];
sx q[3];
rz(2.7618264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.43210426) q[2];
sx q[2];
rz(-2.4435142) q[2];
sx q[2];
rz(3.091605) q[2];
rz(1.6677808) q[3];
sx q[3];
rz(-1.4155017) q[3];
sx q[3];
rz(2.2059691) q[3];
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
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98899984) q[0];
sx q[0];
rz(-2.279156) q[0];
sx q[0];
rz(-1.1784026) q[0];
rz(-1.9940469) q[1];
sx q[1];
rz(-0.18745628) q[1];
sx q[1];
rz(-2.5386834) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41571799) q[0];
sx q[0];
rz(-2.9511669) q[0];
sx q[0];
rz(1.003153) q[0];
rz(-pi) q[1];
rz(0.96430594) q[2];
sx q[2];
rz(-0.78769257) q[2];
sx q[2];
rz(2.380013) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6711839) q[1];
sx q[1];
rz(-2.3215331) q[1];
sx q[1];
rz(2.5633308) q[1];
rz(-pi) q[2];
rz(2.1985377) q[3];
sx q[3];
rz(-2.3585417) q[3];
sx q[3];
rz(0.32367009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9511562) q[2];
sx q[2];
rz(-0.78593212) q[2];
sx q[2];
rz(2.7583165) q[2];
rz(1.8303309) q[3];
sx q[3];
rz(-1.0878599) q[3];
sx q[3];
rz(-3.0200628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7168147) q[0];
sx q[0];
rz(-0.99262339) q[0];
sx q[0];
rz(-0.92856652) q[0];
rz(0.83944744) q[1];
sx q[1];
rz(-1.8239559) q[1];
sx q[1];
rz(-2.3323257) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.033484785) q[0];
sx q[0];
rz(-1.2039469) q[0];
sx q[0];
rz(-0.97292329) q[0];
rz(-pi) q[1];
rz(-0.42994606) q[2];
sx q[2];
rz(-1.4166178) q[2];
sx q[2];
rz(-1.4232612) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4397341) q[1];
sx q[1];
rz(-1.8161521) q[1];
sx q[1];
rz(-0.071003242) q[1];
rz(2.2146642) q[3];
sx q[3];
rz(-1.9044821) q[3];
sx q[3];
rz(0.22367553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.75260085) q[2];
sx q[2];
rz(-1.5316803) q[2];
sx q[2];
rz(-1.4468225) q[2];
rz(-2.0145448) q[3];
sx q[3];
rz(-0.73486745) q[3];
sx q[3];
rz(0.62895044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9160354) q[0];
sx q[0];
rz(-1.1981523) q[0];
sx q[0];
rz(1.4300562) q[0];
rz(0.23722181) q[1];
sx q[1];
rz(-0.96313852) q[1];
sx q[1];
rz(-2.846834) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4859516) q[0];
sx q[0];
rz(-1.427934) q[0];
sx q[0];
rz(0.44012196) q[0];
x q[1];
rz(1.1072269) q[2];
sx q[2];
rz(-0.71706248) q[2];
sx q[2];
rz(-1.1087369) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.7511616) q[1];
sx q[1];
rz(-2.0195578) q[1];
sx q[1];
rz(0.46350355) q[1];
x q[2];
rz(-2.9137524) q[3];
sx q[3];
rz(-0.58594847) q[3];
sx q[3];
rz(2.4333352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1416867) q[2];
sx q[2];
rz(-0.20076951) q[2];
sx q[2];
rz(0.13299175) q[2];
rz(-0.76987949) q[3];
sx q[3];
rz(-1.2917638) q[3];
sx q[3];
rz(-2.8225115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28984508) q[0];
sx q[0];
rz(-0.3388437) q[0];
sx q[0];
rz(-1.7293365) q[0];
rz(0.051636592) q[1];
sx q[1];
rz(-2.5187571) q[1];
sx q[1];
rz(2.9488865) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4225005) q[0];
sx q[0];
rz(-1.0570881) q[0];
sx q[0];
rz(-1.5490319) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.66484837) q[2];
sx q[2];
rz(-2.2908205) q[2];
sx q[2];
rz(2.8262422) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.440408) q[1];
sx q[1];
rz(-1.1478979) q[1];
sx q[1];
rz(-2.6404523) q[1];
rz(-0.64151836) q[3];
sx q[3];
rz(-1.8810442) q[3];
sx q[3];
rz(1.9812685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0685588) q[2];
sx q[2];
rz(-1.2898338) q[2];
sx q[2];
rz(-1.7603091) q[2];
rz(-0.51501385) q[3];
sx q[3];
rz(-2.2541855) q[3];
sx q[3];
rz(1.7611586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.16159049) q[0];
sx q[0];
rz(-1.9965633) q[0];
sx q[0];
rz(2.7951151) q[0];
rz(2.1222291) q[1];
sx q[1];
rz(-2.4788224) q[1];
sx q[1];
rz(-1.6874708) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8603191) q[0];
sx q[0];
rz(-0.77206445) q[0];
sx q[0];
rz(2.4738929) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2242111) q[2];
sx q[2];
rz(-1.1251984) q[2];
sx q[2];
rz(0.080597045) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8365663) q[1];
sx q[1];
rz(-1.2946212) q[1];
sx q[1];
rz(2.1146875) q[1];
x q[2];
rz(-0.10054133) q[3];
sx q[3];
rz(-1.7243598) q[3];
sx q[3];
rz(-2.6996149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.60160294) q[2];
sx q[2];
rz(-1.3301962) q[2];
sx q[2];
rz(-0.23923624) q[2];
rz(2.7573977) q[3];
sx q[3];
rz(-2.4592063) q[3];
sx q[3];
rz(2.189883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9397028) q[0];
sx q[0];
rz(-1.6738142) q[0];
sx q[0];
rz(1.3772759) q[0];
rz(-1.9210303) q[1];
sx q[1];
rz(-0.86253291) q[1];
sx q[1];
rz(0.16622226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8140504) q[0];
sx q[0];
rz(-2.6084427) q[0];
sx q[0];
rz(-2.0943805) q[0];
rz(-pi) q[1];
rz(2.7923585) q[2];
sx q[2];
rz(-2.376261) q[2];
sx q[2];
rz(-2.1527596) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6580171) q[1];
sx q[1];
rz(-1.3596351) q[1];
sx q[1];
rz(-1.5394475) q[1];
rz(2.1020426) q[3];
sx q[3];
rz(-0.47787468) q[3];
sx q[3];
rz(0.22051624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.347747) q[2];
sx q[2];
rz(-0.85024992) q[2];
sx q[2];
rz(1.0469077) q[2];
rz(1.2547803) q[3];
sx q[3];
rz(-2.051765) q[3];
sx q[3];
rz(1.8449529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60085249) q[0];
sx q[0];
rz(-0.37373251) q[0];
sx q[0];
rz(2.0284213) q[0];
rz(0.81740776) q[1];
sx q[1];
rz(-1.4459041) q[1];
sx q[1];
rz(-0.36849749) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93732801) q[0];
sx q[0];
rz(-0.44631347) q[0];
sx q[0];
rz(-2.242779) q[0];
rz(1.4442367) q[2];
sx q[2];
rz(-0.83240763) q[2];
sx q[2];
rz(1.0711993) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8702183) q[1];
sx q[1];
rz(-1.8370359) q[1];
sx q[1];
rz(-0.72721796) q[1];
rz(1.0115252) q[3];
sx q[3];
rz(-0.91010053) q[3];
sx q[3];
rz(1.6370156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.14785279) q[2];
sx q[2];
rz(-2.5233614) q[2];
sx q[2];
rz(-1.3573793) q[2];
rz(2.3944858) q[3];
sx q[3];
rz(-2.1785469) q[3];
sx q[3];
rz(-2.4660306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52637446) q[0];
sx q[0];
rz(-0.45658657) q[0];
sx q[0];
rz(1.0427465) q[0];
rz(-2.3710947) q[1];
sx q[1];
rz(-2.1310525) q[1];
sx q[1];
rz(2.3811293) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2337991) q[0];
sx q[0];
rz(-1.2408419) q[0];
sx q[0];
rz(2.8889662) q[0];
rz(-pi) q[1];
rz(-1.7793525) q[2];
sx q[2];
rz(-0.18392662) q[2];
sx q[2];
rz(-0.49801258) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.059372455) q[1];
sx q[1];
rz(-0.43082985) q[1];
sx q[1];
rz(-0.18980356) q[1];
rz(-pi) q[2];
rz(2.4861927) q[3];
sx q[3];
rz(-1.8701347) q[3];
sx q[3];
rz(1.6884402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.20327917) q[2];
sx q[2];
rz(-0.61050296) q[2];
sx q[2];
rz(2.7321613) q[2];
rz(1.5832541) q[3];
sx q[3];
rz(-0.46054545) q[3];
sx q[3];
rz(-2.515365) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64774491) q[0];
sx q[0];
rz(-1.0418325) q[0];
sx q[0];
rz(-2.7000725) q[0];
rz(-1.5589177) q[1];
sx q[1];
rz(-1.1987004) q[1];
sx q[1];
rz(-0.62335062) q[1];
rz(-2.6263993) q[2];
sx q[2];
rz(-0.11820785) q[2];
sx q[2];
rz(-1.6109656) q[2];
rz(-1.4382451) q[3];
sx q[3];
rz(-1.3302531) q[3];
sx q[3];
rz(2.1147685) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
