OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.83752051) q[0];
sx q[0];
rz(-0.66523319) q[0];
sx q[0];
rz(1.2256149) q[0];
rz(2.4576814) q[1];
sx q[1];
rz(-0.39165762) q[1];
sx q[1];
rz(-2.439523) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0797146) q[0];
sx q[0];
rz(-1.6883381) q[0];
sx q[0];
rz(-2.4410309) q[0];
rz(1.443154) q[2];
sx q[2];
rz(-0.1662456) q[2];
sx q[2];
rz(-0.12285168) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4106284) q[1];
sx q[1];
rz(-2.0515039) q[1];
sx q[1];
rz(2.5334873) q[1];
x q[2];
rz(2.2625173) q[3];
sx q[3];
rz(-1.723515) q[3];
sx q[3];
rz(-0.16593753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2012653) q[2];
sx q[2];
rz(-1.6028812) q[2];
sx q[2];
rz(2.8931457) q[2];
rz(1.1343608) q[3];
sx q[3];
rz(-3.0076707) q[3];
sx q[3];
rz(-2.0094481) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0487173) q[0];
sx q[0];
rz(-0.55260125) q[0];
sx q[0];
rz(1.5485113) q[0];
rz(-2.6379207) q[1];
sx q[1];
rz(-1.2232989) q[1];
sx q[1];
rz(-0.37685397) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.037721264) q[0];
sx q[0];
rz(-1.590017) q[0];
sx q[0];
rz(0.47211418) q[0];
rz(-0.81266788) q[2];
sx q[2];
rz(-2.118131) q[2];
sx q[2];
rz(-3.0603882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.406817) q[1];
sx q[1];
rz(-1.3208913) q[1];
sx q[1];
rz(0.49244778) q[1];
rz(-pi) q[2];
rz(0.33263388) q[3];
sx q[3];
rz(-1.6432228) q[3];
sx q[3];
rz(2.4410332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5571931) q[2];
sx q[2];
rz(-2.445745) q[2];
sx q[2];
rz(-0.86563555) q[2];
rz(1.668476) q[3];
sx q[3];
rz(-0.98554635) q[3];
sx q[3];
rz(0.98751155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5675548) q[0];
sx q[0];
rz(-1.2873298) q[0];
sx q[0];
rz(2.0903184) q[0];
rz(1.6846664) q[1];
sx q[1];
rz(-1.3070062) q[1];
sx q[1];
rz(1.6251224) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54615669) q[0];
sx q[0];
rz(-1.2896898) q[0];
sx q[0];
rz(2.3670822) q[0];
rz(-pi) q[1];
rz(2.462492) q[2];
sx q[2];
rz(-2.1235941) q[2];
sx q[2];
rz(-2.2246085) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9697204) q[1];
sx q[1];
rz(-2.0341968) q[1];
sx q[1];
rz(-0.42831383) q[1];
x q[2];
rz(0.061013075) q[3];
sx q[3];
rz(-1.5170238) q[3];
sx q[3];
rz(-0.53053571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7463344) q[2];
sx q[2];
rz(-1.0552152) q[2];
sx q[2];
rz(0.22221097) q[2];
rz(0.5018417) q[3];
sx q[3];
rz(-1.7343438) q[3];
sx q[3];
rz(-2.3327904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2489081) q[0];
sx q[0];
rz(-0.48766708) q[0];
sx q[0];
rz(1.1647613) q[0];
rz(0.89272967) q[1];
sx q[1];
rz(-1.3803394) q[1];
sx q[1];
rz(-2.8597615) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.079275) q[0];
sx q[0];
rz(-1.6910416) q[0];
sx q[0];
rz(2.6577302) q[0];
rz(-pi) q[1];
x q[1];
rz(2.495419) q[2];
sx q[2];
rz(-2.3149688) q[2];
sx q[2];
rz(-1.7648362) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6161595) q[1];
sx q[1];
rz(-1.6405917) q[1];
sx q[1];
rz(-1.7341803) q[1];
rz(-pi) q[2];
rz(-3.115913) q[3];
sx q[3];
rz(-2.0168983) q[3];
sx q[3];
rz(-0.0080953117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.0622327) q[2];
sx q[2];
rz(-0.56228176) q[2];
sx q[2];
rz(0.53323659) q[2];
rz(-1.0821292) q[3];
sx q[3];
rz(-0.60492587) q[3];
sx q[3];
rz(-2.3281039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37394062) q[0];
sx q[0];
rz(-2.7426608) q[0];
sx q[0];
rz(1.8178513) q[0];
rz(2.3228877) q[1];
sx q[1];
rz(-0.74344126) q[1];
sx q[1];
rz(1.0680107) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49379738) q[0];
sx q[0];
rz(-1.6500705) q[0];
sx q[0];
rz(-0.1597516) q[0];
rz(-pi) q[1];
rz(2.4765268) q[2];
sx q[2];
rz(-2.1882727) q[2];
sx q[2];
rz(-1.1368874) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.9449759) q[1];
sx q[1];
rz(-0.91755897) q[1];
sx q[1];
rz(1.7004299) q[1];
rz(-1.0982008) q[3];
sx q[3];
rz(-1.7180746) q[3];
sx q[3];
rz(-0.46340273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9615122) q[2];
sx q[2];
rz(-2.027812) q[2];
sx q[2];
rz(-0.86165825) q[2];
rz(-2.1152451) q[3];
sx q[3];
rz(-1.15871) q[3];
sx q[3];
rz(-1.1177184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79477972) q[0];
sx q[0];
rz(-0.81951278) q[0];
sx q[0];
rz(0.89282194) q[0];
rz(-0.18094856) q[1];
sx q[1];
rz(-0.67829689) q[1];
sx q[1];
rz(0.13042626) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3149726) q[0];
sx q[0];
rz(-2.9493647) q[0];
sx q[0];
rz(-2.1703224) q[0];
rz(-pi) q[1];
rz(2.0753292) q[2];
sx q[2];
rz(-2.1158525) q[2];
sx q[2];
rz(-1.2912599) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.32864704) q[1];
sx q[1];
rz(-2.9025902) q[1];
sx q[1];
rz(-1.4471329) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.1993066) q[3];
sx q[3];
rz(-2.350507) q[3];
sx q[3];
rz(-1.4996604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9170407) q[2];
sx q[2];
rz(-1.9150534) q[2];
sx q[2];
rz(0.87635931) q[2];
rz(-1.2419491) q[3];
sx q[3];
rz(-2.2010937) q[3];
sx q[3];
rz(1.0103286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4001615) q[0];
sx q[0];
rz(-3.042996) q[0];
sx q[0];
rz(-2.7778991) q[0];
rz(-2.3997276) q[1];
sx q[1];
rz(-1.0262841) q[1];
sx q[1];
rz(0.40282869) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9972853) q[0];
sx q[0];
rz(-1.8107521) q[0];
sx q[0];
rz(1.994094) q[0];
rz(-pi) q[1];
rz(0.81554697) q[2];
sx q[2];
rz(-0.52482739) q[2];
sx q[2];
rz(1.8458402) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0707222) q[1];
sx q[1];
rz(-1.4938584) q[1];
sx q[1];
rz(1.6926195) q[1];
x q[2];
rz(0.63490156) q[3];
sx q[3];
rz(-1.7909808) q[3];
sx q[3];
rz(1.6600773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.3591298) q[2];
sx q[2];
rz(-0.368258) q[2];
sx q[2];
rz(1.5472319) q[2];
rz(-0.83526978) q[3];
sx q[3];
rz(-0.95649496) q[3];
sx q[3];
rz(-2.6529151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35992026) q[0];
sx q[0];
rz(-1.3363573) q[0];
sx q[0];
rz(2.6265327) q[0];
rz(1.4393648) q[1];
sx q[1];
rz(-1.8481588) q[1];
sx q[1];
rz(2.7424367) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7997417) q[0];
sx q[0];
rz(-0.13093311) q[0];
sx q[0];
rz(-0.013850613) q[0];
rz(-pi) q[1];
rz(-1.3005343) q[2];
sx q[2];
rz(-2.3944602) q[2];
sx q[2];
rz(3.0556553) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5134352) q[1];
sx q[1];
rz(-1.5062467) q[1];
sx q[1];
rz(-1.3099531) q[1];
rz(0.018980515) q[3];
sx q[3];
rz(-3.0759263) q[3];
sx q[3];
rz(1.1733023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1601552) q[2];
sx q[2];
rz(-1.8835386) q[2];
sx q[2];
rz(-1.0464279) q[2];
rz(0.6662755) q[3];
sx q[3];
rz(-2.8847238) q[3];
sx q[3];
rz(1.0506786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.94205725) q[0];
sx q[0];
rz(-0.10534795) q[0];
sx q[0];
rz(0.60607213) q[0];
rz(-0.82707682) q[1];
sx q[1];
rz(-1.3304293) q[1];
sx q[1];
rz(-1.8082632) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5904986) q[0];
sx q[0];
rz(-2.7685389) q[0];
sx q[0];
rz(2.8540552) q[0];
rz(-pi) q[1];
rz(-0.47918646) q[2];
sx q[2];
rz(-2.3385417) q[2];
sx q[2];
rz(0.89862862) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8136452) q[1];
sx q[1];
rz(-0.678181) q[1];
sx q[1];
rz(1.7575043) q[1];
rz(-pi) q[2];
rz(2.9700432) q[3];
sx q[3];
rz(-2.6772873) q[3];
sx q[3];
rz(2.7363079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0102319) q[2];
sx q[2];
rz(-0.92326814) q[2];
sx q[2];
rz(0.38491797) q[2];
rz(0.63079232) q[3];
sx q[3];
rz(-1.8600978) q[3];
sx q[3];
rz(-1.7990641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.887562) q[0];
sx q[0];
rz(-1.3986724) q[0];
sx q[0];
rz(1.4400462) q[0];
rz(2.4002659) q[1];
sx q[1];
rz(-1.7404375) q[1];
sx q[1];
rz(0.25873605) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7546623) q[0];
sx q[0];
rz(-1.5748236) q[0];
sx q[0];
rz(1.8167102) q[0];
x q[1];
rz(-2.0889241) q[2];
sx q[2];
rz(-1.8473704) q[2];
sx q[2];
rz(0.58108789) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1867793) q[1];
sx q[1];
rz(-0.80866058) q[1];
sx q[1];
rz(-1.0891799) q[1];
rz(-1.1734937) q[3];
sx q[3];
rz(-2.3081995) q[3];
sx q[3];
rz(-1.6034338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.72186333) q[2];
sx q[2];
rz(-2.0294919) q[2];
sx q[2];
rz(1.8224243) q[2];
rz(-2.3223274) q[3];
sx q[3];
rz(-2.0115439) q[3];
sx q[3];
rz(2.5949902) q[3];
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
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4010314) q[0];
sx q[0];
rz(-0.87420976) q[0];
sx q[0];
rz(2.5921205) q[0];
rz(1.7026547) q[1];
sx q[1];
rz(-2.4733652) q[1];
sx q[1];
rz(0.53818902) q[1];
rz(-2.2808711) q[2];
sx q[2];
rz(-1.5790944) q[2];
sx q[2];
rz(-0.38766833) q[2];
rz(1.9520252) q[3];
sx q[3];
rz(-0.59567957) q[3];
sx q[3];
rz(2.4761562) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
