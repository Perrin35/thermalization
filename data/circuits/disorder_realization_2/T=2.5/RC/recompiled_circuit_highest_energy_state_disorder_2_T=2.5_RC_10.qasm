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
rz(2.9444115) q[0];
sx q[0];
rz(-1.0721167) q[0];
sx q[0];
rz(-1.5972692) q[0];
rz(0.3577258) q[1];
sx q[1];
rz(4.975714) q[1];
sx q[1];
rz(9.8006021) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.25426769) q[0];
sx q[0];
rz(-1.642881) q[0];
sx q[0];
rz(-2.3743126) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2959458) q[2];
sx q[2];
rz(-2.0811199) q[2];
sx q[2];
rz(-0.21462552) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5800505) q[1];
sx q[1];
rz(-2.1218178) q[1];
sx q[1];
rz(-2.2660794) q[1];
rz(-2.8055111) q[3];
sx q[3];
rz(-2.1770526) q[3];
sx q[3];
rz(1.1949902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.4481005) q[2];
sx q[2];
rz(-2.7585612) q[2];
sx q[2];
rz(2.7310272) q[2];
rz(-0.52452123) q[3];
sx q[3];
rz(-2.9295242) q[3];
sx q[3];
rz(-1.9922403) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9583424) q[0];
sx q[0];
rz(-0.14160937) q[0];
sx q[0];
rz(0.70567065) q[0];
rz(-1.8535829) q[1];
sx q[1];
rz(-2.5649004) q[1];
sx q[1];
rz(-2.0075331) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81728888) q[0];
sx q[0];
rz(-1.9729483) q[0];
sx q[0];
rz(-0.83950801) q[0];
rz(-pi) q[1];
rz(2.6666849) q[2];
sx q[2];
rz(-0.38198745) q[2];
sx q[2];
rz(1.8595921) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4422851) q[1];
sx q[1];
rz(-0.87934994) q[1];
sx q[1];
rz(-1.0710909) q[1];
rz(-2.5763578) q[3];
sx q[3];
rz(-1.7218523) q[3];
sx q[3];
rz(-0.24542576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.80295339) q[2];
sx q[2];
rz(-2.6889375) q[2];
sx q[2];
rz(-2.6111641) q[2];
rz(-2.8545634) q[3];
sx q[3];
rz(-2.65286) q[3];
sx q[3];
rz(2.4052525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0158143) q[0];
sx q[0];
rz(-1.1017355) q[0];
sx q[0];
rz(-2.0871479) q[0];
rz(-1.6619445) q[1];
sx q[1];
rz(-0.92250657) q[1];
sx q[1];
rz(-2.0747917) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.023533688) q[0];
sx q[0];
rz(-1.3830782) q[0];
sx q[0];
rz(-0.037506692) q[0];
x q[1];
rz(-2.4652723) q[2];
sx q[2];
rz(-1.6427077) q[2];
sx q[2];
rz(1.7184217) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.57745614) q[1];
sx q[1];
rz(-1.3423663) q[1];
sx q[1];
rz(2.0453358) q[1];
rz(-pi) q[2];
rz(2.8425806) q[3];
sx q[3];
rz(-0.87871534) q[3];
sx q[3];
rz(2.690154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.58831424) q[2];
sx q[2];
rz(-1.30043) q[2];
sx q[2];
rz(-1.1021357) q[2];
rz(2.6831902) q[3];
sx q[3];
rz(-1.6130684) q[3];
sx q[3];
rz(-0.63157356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-0.44963968) q[0];
sx q[0];
rz(-0.75943685) q[0];
sx q[0];
rz(2.6547883) q[0];
rz(1.5454769) q[1];
sx q[1];
rz(-2.7524452) q[1];
sx q[1];
rz(2.5130689) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5647873) q[0];
sx q[0];
rz(-2.724455) q[0];
sx q[0];
rz(1.9030626) q[0];
rz(-pi) q[1];
rz(-0.14735584) q[2];
sx q[2];
rz(-2.1389845) q[2];
sx q[2];
rz(2.1939366) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1284229) q[1];
sx q[1];
rz(-2.7047415) q[1];
sx q[1];
rz(2.0274202) q[1];
x q[2];
rz(3.0699203) q[3];
sx q[3];
rz(-2.1183156) q[3];
sx q[3];
rz(1.6184834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.5635166) q[2];
sx q[2];
rz(-3.0264644) q[2];
sx q[2];
rz(2.3797177) q[2];
rz(0.21841194) q[3];
sx q[3];
rz(-2.8438447) q[3];
sx q[3];
rz(1.0485605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8258284) q[0];
sx q[0];
rz(-2.789412) q[0];
sx q[0];
rz(-0.6045565) q[0];
rz(-1.8479337) q[1];
sx q[1];
rz(-0.81984055) q[1];
sx q[1];
rz(2.7972319) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6647756) q[0];
sx q[0];
rz(-1.4882384) q[0];
sx q[0];
rz(0.020391012) q[0];
rz(0.44565752) q[2];
sx q[2];
rz(-0.67699922) q[2];
sx q[2];
rz(2.0930365) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.02270123) q[1];
sx q[1];
rz(-1.1384607) q[1];
sx q[1];
rz(2.6328264) q[1];
rz(-pi) q[2];
rz(2.6743498) q[3];
sx q[3];
rz(-0.95047104) q[3];
sx q[3];
rz(-2.1324699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1206426) q[2];
sx q[2];
rz(-2.5072493) q[2];
sx q[2];
rz(2.9316588) q[2];
rz(1.1076934) q[3];
sx q[3];
rz(-1.6898797) q[3];
sx q[3];
rz(-0.2618739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8629859) q[0];
sx q[0];
rz(-2.5513273) q[0];
sx q[0];
rz(-3.0721831) q[0];
rz(0.033992652) q[1];
sx q[1];
rz(-2.4885204) q[1];
sx q[1];
rz(2.6896844) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2632217) q[0];
sx q[0];
rz(-0.42455205) q[0];
sx q[0];
rz(2.7431737) q[0];
rz(-pi) q[1];
rz(-2.4546409) q[2];
sx q[2];
rz(-2.070362) q[2];
sx q[2];
rz(-2.9925008) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0167061) q[1];
sx q[1];
rz(-1.5779902) q[1];
sx q[1];
rz(1.5924988) q[1];
x q[2];
rz(2.8699401) q[3];
sx q[3];
rz(-1.324904) q[3];
sx q[3];
rz(-3.0761443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7273093) q[2];
sx q[2];
rz(-0.80237979) q[2];
sx q[2];
rz(-1.9963973) q[2];
rz(0.9582054) q[3];
sx q[3];
rz(-1.9327791) q[3];
sx q[3];
rz(0.27829471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3713643) q[0];
sx q[0];
rz(-0.48786491) q[0];
sx q[0];
rz(0.86139739) q[0];
rz(2.1791012) q[1];
sx q[1];
rz(-2.1884514) q[1];
sx q[1];
rz(0.15693396) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.033534) q[0];
sx q[0];
rz(-0.80488759) q[0];
sx q[0];
rz(-1.0548841) q[0];
x q[1];
rz(0.10617039) q[2];
sx q[2];
rz(-1.9261126) q[2];
sx q[2];
rz(2.4715142) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.90060252) q[1];
sx q[1];
rz(-0.99731904) q[1];
sx q[1];
rz(2.1585224) q[1];
rz(-pi) q[2];
rz(-3.0928218) q[3];
sx q[3];
rz(-1.7836106) q[3];
sx q[3];
rz(1.770894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.018205) q[2];
sx q[2];
rz(-0.79424477) q[2];
sx q[2];
rz(1.0495074) q[2];
rz(-0.47979313) q[3];
sx q[3];
rz(-0.19194651) q[3];
sx q[3];
rz(2.979579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4341693) q[0];
sx q[0];
rz(-2.8987085) q[0];
sx q[0];
rz(-0.49651399) q[0];
rz(-1.6560582) q[1];
sx q[1];
rz(-0.86692202) q[1];
sx q[1];
rz(3.1190125) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3568087) q[0];
sx q[0];
rz(-1.7473146) q[0];
sx q[0];
rz(1.7872732) q[0];
x q[1];
rz(1.7143676) q[2];
sx q[2];
rz(-1.5036791) q[2];
sx q[2];
rz(2.3877833) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.034790967) q[1];
sx q[1];
rz(-2.3402434) q[1];
sx q[1];
rz(1.31557) q[1];
x q[2];
rz(2.1053505) q[3];
sx q[3];
rz(-1.4009985) q[3];
sx q[3];
rz(-1.2359197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.61985832) q[2];
sx q[2];
rz(-2.0010184) q[2];
sx q[2];
rz(0.015259585) q[2];
rz(-2.7742625) q[3];
sx q[3];
rz(-0.44706774) q[3];
sx q[3];
rz(2.7801311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5527282) q[0];
sx q[0];
rz(-2.9155154) q[0];
sx q[0];
rz(3.0638301) q[0];
rz(0.11847682) q[1];
sx q[1];
rz(-0.34093726) q[1];
sx q[1];
rz(2.4854787) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8828147) q[0];
sx q[0];
rz(-1.4442632) q[0];
sx q[0];
rz(-1.2650498) q[0];
rz(-pi) q[1];
rz(2.9703548) q[2];
sx q[2];
rz(-1.3067553) q[2];
sx q[2];
rz(-1.1080826) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.39086884) q[1];
sx q[1];
rz(-2.0489397) q[1];
sx q[1];
rz(-0.45420809) q[1];
rz(-pi) q[2];
rz(2.0875075) q[3];
sx q[3];
rz(-2.0626522) q[3];
sx q[3];
rz(2.3824177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.52350998) q[2];
sx q[2];
rz(-1.9205576) q[2];
sx q[2];
rz(0.54100424) q[2];
rz(-0.4506166) q[3];
sx q[3];
rz(-1.6559947) q[3];
sx q[3];
rz(-0.97644067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4377874) q[0];
sx q[0];
rz(-1.7923651) q[0];
sx q[0];
rz(3.034814) q[0];
rz(1.59185) q[1];
sx q[1];
rz(-0.50250643) q[1];
sx q[1];
rz(-1.7639311) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1049809) q[0];
sx q[0];
rz(-2.689765) q[0];
sx q[0];
rz(-1.4696448) q[0];
x q[1];
rz(-2.8286727) q[2];
sx q[2];
rz(-1.7831137) q[2];
sx q[2];
rz(2.2231399) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.71860524) q[1];
sx q[1];
rz(-1.8180625) q[1];
sx q[1];
rz(-0.24691564) q[1];
rz(-pi) q[2];
x q[2];
rz(1.000147) q[3];
sx q[3];
rz(-1.1091145) q[3];
sx q[3];
rz(-0.0014071597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0926853) q[2];
sx q[2];
rz(-1.4212298) q[2];
sx q[2];
rz(2.8468724) q[2];
rz(-2.5367391) q[3];
sx q[3];
rz(-1.9113144) q[3];
sx q[3];
rz(3.0321583) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43815628) q[0];
sx q[0];
rz(-1.3180757) q[0];
sx q[0];
rz(-0.99676589) q[0];
rz(0.017771487) q[1];
sx q[1];
rz(-1.8780864) q[1];
sx q[1];
rz(-1.3395739) q[1];
rz(-0.97311592) q[2];
sx q[2];
rz(-0.86116198) q[2];
sx q[2];
rz(1.8088928) q[2];
rz(-2.2733489) q[3];
sx q[3];
rz(-2.1001742) q[3];
sx q[3];
rz(-2.1471837) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
