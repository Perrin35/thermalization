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
rz(-0.67796081) q[0];
sx q[0];
rz(-0.27422658) q[0];
sx q[0];
rz(1.8507313) q[0];
rz(-0.21386799) q[1];
sx q[1];
rz(-2.8254421) q[1];
sx q[1];
rz(1.4996127) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9737273) q[0];
sx q[0];
rz(-1.8059547) q[0];
sx q[0];
rz(-0.27882378) q[0];
x q[1];
rz(1.2798645) q[2];
sx q[2];
rz(-2.3290344) q[2];
sx q[2];
rz(-1.1464553) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1823481) q[1];
sx q[1];
rz(-2.3492491) q[1];
sx q[1];
rz(-0.46087973) q[1];
x q[2];
rz(0.95176272) q[3];
sx q[3];
rz(-0.72975273) q[3];
sx q[3];
rz(2.8918134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9700254) q[2];
sx q[2];
rz(-1.0415404) q[2];
sx q[2];
rz(-2.2402666) q[2];
rz(0.46402913) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(0.06981167) q[3];
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
sx q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0932015) q[0];
sx q[0];
rz(-0.036660107) q[0];
sx q[0];
rz(-1.8638336) q[0];
rz(3.0473862) q[1];
sx q[1];
rz(-0.46368805) q[1];
sx q[1];
rz(1.6069848) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8600991) q[0];
sx q[0];
rz(-1.5380385) q[0];
sx q[0];
rz(2.1707235) q[0];
rz(-pi) q[1];
rz(1.2708475) q[2];
sx q[2];
rz(-2.9510164) q[2];
sx q[2];
rz(1.9924763) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.4581286) q[1];
sx q[1];
rz(-2.5741815) q[1];
sx q[1];
rz(0.50000425) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7369698) q[3];
sx q[3];
rz(-1.957486) q[3];
sx q[3];
rz(0.85961941) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.719912) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(-0.16033944) q[2];
rz(2.4028589) q[3];
sx q[3];
rz(-2.2952047) q[3];
sx q[3];
rz(-1.4275449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5147603) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(2.6318188) q[0];
rz(-1.431142) q[1];
sx q[1];
rz(-1.9606083) q[1];
sx q[1];
rz(0.74657718) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7187481) q[0];
sx q[0];
rz(-1.5908518) q[0];
sx q[0];
rz(-2.123318) q[0];
rz(-pi) q[1];
rz(0.96087281) q[2];
sx q[2];
rz(-0.93797937) q[2];
sx q[2];
rz(-0.52458143) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1521897) q[1];
sx q[1];
rz(-0.8968938) q[1];
sx q[1];
rz(2.6722355) q[1];
x q[2];
rz(2.7535149) q[3];
sx q[3];
rz(-0.92736926) q[3];
sx q[3];
rz(0.79659792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9868682) q[2];
sx q[2];
rz(-2.8747323) q[2];
sx q[2];
rz(-1.9179087) q[2];
rz(-1.0736505) q[3];
sx q[3];
rz(-1.4370388) q[3];
sx q[3];
rz(-2.2435718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8680442) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(0.37937382) q[0];
rz(0.1768449) q[1];
sx q[1];
rz(-1.481448) q[1];
sx q[1];
rz(0.79536974) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58461713) q[0];
sx q[0];
rz(-1.66357) q[0];
sx q[0];
rz(1.2014821) q[0];
rz(-1.2878809) q[2];
sx q[2];
rz(-0.77436111) q[2];
sx q[2];
rz(-1.1224358) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2242608) q[1];
sx q[1];
rz(-0.59228173) q[1];
sx q[1];
rz(1.0975973) q[1];
x q[2];
rz(0.12563184) q[3];
sx q[3];
rz(-1.5141308) q[3];
sx q[3];
rz(-1.1606596) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.4312326) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(-1.360652) q[2];
rz(-1.0294754) q[3];
sx q[3];
rz(-1.4808713) q[3];
sx q[3];
rz(1.8195389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65115702) q[0];
sx q[0];
rz(-1.5999726) q[0];
sx q[0];
rz(-1.5486451) q[0];
rz(-2.7410638) q[1];
sx q[1];
rz(-1.6533886) q[1];
sx q[1];
rz(-1.7318447) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.627305) q[0];
sx q[0];
rz(-2.3804733) q[0];
sx q[0];
rz(-2.850015) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.35946741) q[2];
sx q[2];
rz(-1.7181944) q[2];
sx q[2];
rz(1.1254252) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5286136) q[1];
sx q[1];
rz(-0.76957146) q[1];
sx q[1];
rz(-0.90958325) q[1];
x q[2];
rz(3.0733969) q[3];
sx q[3];
rz(-2.1917097) q[3];
sx q[3];
rz(-0.62486006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.8283525) q[2];
sx q[2];
rz(-1.4408377) q[2];
sx q[2];
rz(0.43133119) q[2];
rz(0.055015419) q[3];
sx q[3];
rz(-2.8300245) q[3];
sx q[3];
rz(-0.55606786) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
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
rz(0.7677652) q[0];
sx q[0];
rz(-2.8887833) q[0];
sx q[0];
rz(1.025169) q[0];
rz(2.688736) q[1];
sx q[1];
rz(-2.5621474) q[1];
sx q[1];
rz(2.2377009) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1548123) q[0];
sx q[0];
rz(-2.203699) q[0];
sx q[0];
rz(-2.7557719) q[0];
x q[1];
rz(0.63389961) q[2];
sx q[2];
rz(-2.1261151) q[2];
sx q[2];
rz(-0.429053) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.55906193) q[1];
sx q[1];
rz(-0.85147038) q[1];
sx q[1];
rz(-1.2079617) q[1];
rz(-pi) q[2];
rz(-2.3748906) q[3];
sx q[3];
rz(-0.23582102) q[3];
sx q[3];
rz(-0.23877777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.30770939) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(2.9252388) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-2.4560865) q[3];
sx q[3];
rz(-1.640813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0093339) q[0];
sx q[0];
rz(-1.1067156) q[0];
sx q[0];
rz(2.4427781) q[0];
rz(1.9578594) q[1];
sx q[1];
rz(-1.4212757) q[1];
sx q[1];
rz(-0.85174495) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4119371) q[0];
sx q[0];
rz(-1.0348399) q[0];
sx q[0];
rz(3.0974749) q[0];
rz(-3.1021031) q[2];
sx q[2];
rz(-2.563884) q[2];
sx q[2];
rz(-2.5414987) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.70515436) q[1];
sx q[1];
rz(-0.90696224) q[1];
sx q[1];
rz(-1.1052119) q[1];
rz(0.64166358) q[3];
sx q[3];
rz(-1.9681276) q[3];
sx q[3];
rz(-1.404055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.304004) q[2];
sx q[2];
rz(-1.3112661) q[2];
sx q[2];
rz(0.50789976) q[2];
rz(1.0287644) q[3];
sx q[3];
rz(-0.84380904) q[3];
sx q[3];
rz(2.540551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.69865882) q[0];
sx q[0];
rz(-1.3154987) q[0];
sx q[0];
rz(-3.1319295) q[0];
rz(-1.6161605) q[1];
sx q[1];
rz(-1.3776255) q[1];
sx q[1];
rz(2.756871) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1533536) q[0];
sx q[0];
rz(-1.8073449) q[0];
sx q[0];
rz(-2.9896215) q[0];
x q[1];
rz(-1.9571976) q[2];
sx q[2];
rz(-1.3829872) q[2];
sx q[2];
rz(-1.9822497) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.5009969) q[1];
sx q[1];
rz(-1.2944376) q[1];
sx q[1];
rz(2.9261057) q[1];
rz(-0.77120292) q[3];
sx q[3];
rz(-1.1386385) q[3];
sx q[3];
rz(1.3753396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.559451) q[2];
sx q[2];
rz(-0.95953512) q[2];
sx q[2];
rz(3.1257296) q[2];
rz(1.2265685) q[3];
sx q[3];
rz(-2.3358986) q[3];
sx q[3];
rz(-2.5198643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-1.7568307) q[0];
sx q[0];
rz(-2.4767196) q[0];
sx q[0];
rz(-1.0913947) q[0];
rz(-0.81820828) q[1];
sx q[1];
rz(-1.810377) q[1];
sx q[1];
rz(-0.6689201) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2125137) q[0];
sx q[0];
rz(-0.90796472) q[0];
sx q[0];
rz(0.79848358) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.0537655) q[2];
sx q[2];
rz(-2.3685799) q[2];
sx q[2];
rz(-0.20394606) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0099677) q[1];
sx q[1];
rz(-2.5796081) q[1];
sx q[1];
rz(-1.4373006) q[1];
rz(-pi) q[2];
rz(1.091033) q[3];
sx q[3];
rz(-2.5026813) q[3];
sx q[3];
rz(-1.0425488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5074671) q[2];
sx q[2];
rz(-2.6125245) q[2];
sx q[2];
rz(-0.96366209) q[2];
rz(-1.4108747) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(1.7760407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.1249579) q[0];
sx q[0];
rz(-2.5732915) q[0];
sx q[0];
rz(1.6868663) q[0];
rz(1.1085054) q[1];
sx q[1];
rz(-1.6051555) q[1];
sx q[1];
rz(-0.32807168) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.011400819) q[0];
sx q[0];
rz(-1.6284124) q[0];
sx q[0];
rz(-1.685623) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56888442) q[2];
sx q[2];
rz(-0.21368229) q[2];
sx q[2];
rz(-1.3485707) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.83392622) q[1];
sx q[1];
rz(-2.1824679) q[1];
sx q[1];
rz(-1.0990547) q[1];
x q[2];
rz(-1.0193908) q[3];
sx q[3];
rz(-2.9388722) q[3];
sx q[3];
rz(-2.3644476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.46128094) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(-0.4168365) q[2];
rz(0.69878116) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(0.39066395) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9417435) q[0];
sx q[0];
rz(-1.2171634) q[0];
sx q[0];
rz(-1.4456277) q[0];
rz(-2.4327714) q[1];
sx q[1];
rz(-0.91239057) q[1];
sx q[1];
rz(-0.053587996) q[1];
rz(1.3016635) q[2];
sx q[2];
rz(-2.4872645) q[2];
sx q[2];
rz(-0.80254868) q[2];
rz(0.82984701) q[3];
sx q[3];
rz(-2.0086858) q[3];
sx q[3];
rz(0.84926844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
