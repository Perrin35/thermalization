OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.29785922) q[0];
sx q[0];
rz(-2.5279186) q[0];
sx q[0];
rz(2.4224129) q[0];
rz(1.367388) q[1];
sx q[1];
rz(-0.24582882) q[1];
sx q[1];
rz(2.153102) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91627097) q[0];
sx q[0];
rz(-2.7164408) q[0];
sx q[0];
rz(-1.7122373) q[0];
x q[1];
rz(-0.97889401) q[2];
sx q[2];
rz(-1.4375293) q[2];
sx q[2];
rz(-0.29083458) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.5060252) q[1];
sx q[1];
rz(-0.68817455) q[1];
sx q[1];
rz(-2.3597005) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.037976102) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(1.2167041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6926379) q[2];
sx q[2];
rz(-1.2636893) q[2];
sx q[2];
rz(-0.56837481) q[2];
rz(-1.189399) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(-0.18251671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.73873591) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(-2.7785595) q[0];
rz(2.1733213) q[1];
sx q[1];
rz(-2.4666511) q[1];
sx q[1];
rz(1.2526858) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2371212) q[0];
sx q[0];
rz(-1.4112817) q[0];
sx q[0];
rz(-1.6617387) q[0];
rz(-pi) q[1];
rz(-2.230174) q[2];
sx q[2];
rz(-1.0936001) q[2];
sx q[2];
rz(2.0267682) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.15236552) q[1];
sx q[1];
rz(-0.75132912) q[1];
sx q[1];
rz(-2.3163124) q[1];
x q[2];
rz(-1.9545752) q[3];
sx q[3];
rz(-0.46758258) q[3];
sx q[3];
rz(-0.13320696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.12053717) q[2];
sx q[2];
rz(-1.3201822) q[2];
sx q[2];
rz(-2.7978314) q[2];
rz(-2.8072642) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(-2.6722369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47135982) q[0];
sx q[0];
rz(-0.25046644) q[0];
sx q[0];
rz(3.1345471) q[0];
rz(-2.7650611) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(-2.4287756) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3526488) q[0];
sx q[0];
rz(-1.3508571) q[0];
sx q[0];
rz(-0.21531944) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1586645) q[2];
sx q[2];
rz(-0.40908989) q[2];
sx q[2];
rz(1.2661753) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.30046001) q[1];
sx q[1];
rz(-1.7736048) q[1];
sx q[1];
rz(2.1122785) q[1];
x q[2];
rz(0.080428877) q[3];
sx q[3];
rz(-2.0518528) q[3];
sx q[3];
rz(-1.2642494) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(0.66398579) q[2];
rz(-1.239423) q[3];
sx q[3];
rz(-0.23012161) q[3];
sx q[3];
rz(0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5749213) q[0];
sx q[0];
rz(-3.1118588) q[0];
sx q[0];
rz(0.57408875) q[0];
rz(0.29218778) q[1];
sx q[1];
rz(-2.8657587) q[1];
sx q[1];
rz(2.2132197) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5260122) q[0];
sx q[0];
rz(-2.2768339) q[0];
sx q[0];
rz(-3.1403149) q[0];
rz(2.3525535) q[2];
sx q[2];
rz(-1.3256729) q[2];
sx q[2];
rz(-2.6546728) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.2239383) q[1];
sx q[1];
rz(-1.8173216) q[1];
sx q[1];
rz(2.0018105) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91289642) q[3];
sx q[3];
rz(-1.1607338) q[3];
sx q[3];
rz(1.0421154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0488247) q[2];
sx q[2];
rz(-2.6269045) q[2];
sx q[2];
rz(1.9862004) q[2];
rz(2.998735) q[3];
sx q[3];
rz(-2.6070049) q[3];
sx q[3];
rz(2.3891881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0578385) q[0];
sx q[0];
rz(-0.46828073) q[0];
sx q[0];
rz(0.074247867) q[0];
rz(-1.4986562) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(2.1898851) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041391011) q[0];
sx q[0];
rz(-2.3889184) q[0];
sx q[0];
rz(1.4647096) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1553467) q[2];
sx q[2];
rz(-0.84976053) q[2];
sx q[2];
rz(-0.69794387) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.530572) q[1];
sx q[1];
rz(-1.7360577) q[1];
sx q[1];
rz(2.1940439) q[1];
rz(-pi) q[2];
rz(-0.40162556) q[3];
sx q[3];
rz(-1.0394319) q[3];
sx q[3];
rz(-1.348996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.64530659) q[2];
sx q[2];
rz(-2.9149865) q[2];
sx q[2];
rz(-1.9892233) q[2];
rz(-1.0460098) q[3];
sx q[3];
rz(-2.4029616) q[3];
sx q[3];
rz(3.1402821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1482658) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(-2.6519725) q[0];
rz(-1.024225) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(1.6061868) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.09774694) q[0];
sx q[0];
rz(-1.4717719) q[0];
sx q[0];
rz(1.6188341) q[0];
rz(-0.77184446) q[2];
sx q[2];
rz(-0.8470042) q[2];
sx q[2];
rz(-2.6592902) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3052169) q[1];
sx q[1];
rz(-1.7292542) q[1];
sx q[1];
rz(1.8578908) q[1];
rz(-pi) q[2];
rz(-0.62932265) q[3];
sx q[3];
rz(-1.4368847) q[3];
sx q[3];
rz(2.8117992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.17807047) q[2];
sx q[2];
rz(-0.81729752) q[2];
sx q[2];
rz(-2.9658588) q[2];
rz(-2.0641616) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(-0.0065461672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0734633) q[0];
sx q[0];
rz(-2.926565) q[0];
sx q[0];
rz(-1.2699132) q[0];
rz(-0.1779671) q[1];
sx q[1];
rz(-0.83825076) q[1];
sx q[1];
rz(-1.3659182) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.022) q[0];
sx q[0];
rz(-0.79607841) q[0];
sx q[0];
rz(-2.6636366) q[0];
rz(-pi) q[1];
rz(1.3283417) q[2];
sx q[2];
rz(-1.889466) q[2];
sx q[2];
rz(-0.53675011) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7682225) q[1];
sx q[1];
rz(-1.7473514) q[1];
sx q[1];
rz(-2.8819487) q[1];
rz(2.6122983) q[3];
sx q[3];
rz(-1.4636453) q[3];
sx q[3];
rz(1.9107493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0697249) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(-2.243637) q[2];
rz(-2.9979624) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(0.34415054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9550069) q[0];
sx q[0];
rz(-2.2877559) q[0];
sx q[0];
rz(2.8198077) q[0];
rz(-2.2161662) q[1];
sx q[1];
rz(-1.4326452) q[1];
sx q[1];
rz(-0.61703533) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077438146) q[0];
sx q[0];
rz(-1.7700717) q[0];
sx q[0];
rz(-1.827924) q[0];
rz(0.033173325) q[2];
sx q[2];
rz(-1.7863986) q[2];
sx q[2];
rz(-2.6278091) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.033213381) q[1];
sx q[1];
rz(-1.9404267) q[1];
sx q[1];
rz(1.0934248) q[1];
rz(-pi) q[2];
rz(3.0934107) q[3];
sx q[3];
rz(-2.391624) q[3];
sx q[3];
rz(0.90660209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.59163219) q[2];
sx q[2];
rz(-2.154921) q[2];
sx q[2];
rz(3.0275596) q[2];
rz(0.36241254) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(1.9053649) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47700259) q[0];
sx q[0];
rz(-0.53628558) q[0];
sx q[0];
rz(-1.3569008) q[0];
rz(2.3214031) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(-1.483451) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26950715) q[0];
sx q[0];
rz(-1.5604696) q[0];
sx q[0];
rz(-1.6540065) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5513068) q[2];
sx q[2];
rz(-1.6620518) q[2];
sx q[2];
rz(1.1003189) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4087275) q[1];
sx q[1];
rz(-0.36204007) q[1];
sx q[1];
rz(-1.15508) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5914707) q[3];
sx q[3];
rz(-1.2391483) q[3];
sx q[3];
rz(1.5105997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0004398) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(1.8971987) q[2];
rz(0.42516285) q[3];
sx q[3];
rz(-2.3624698) q[3];
sx q[3];
rz(1.5104793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4193831) q[0];
sx q[0];
rz(-0.49383759) q[0];
sx q[0];
rz(0.37049946) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(-1.3409021) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1271034) q[0];
sx q[0];
rz(-1.9368441) q[0];
sx q[0];
rz(-2.8949379) q[0];
rz(0.71435931) q[2];
sx q[2];
rz(-1.6742801) q[2];
sx q[2];
rz(-1.6911182) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.51532981) q[1];
sx q[1];
rz(-1.9326107) q[1];
sx q[1];
rz(-2.8742909) q[1];
x q[2];
rz(-0.16593905) q[3];
sx q[3];
rz(-0.12291848) q[3];
sx q[3];
rz(-0.77792203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.6293388) q[2];
sx q[2];
rz(-0.24015716) q[2];
sx q[2];
rz(1.5226927) q[2];
rz(2.2807138) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(-3.1057152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9027949) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(-2.8181656) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(-1.9771489) q[2];
sx q[2];
rz(-1.5487557) q[2];
sx q[2];
rz(1.1974481) q[2];
rz(-1.0988416) q[3];
sx q[3];
rz(-1.2225716) q[3];
sx q[3];
rz(-0.93248758) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
