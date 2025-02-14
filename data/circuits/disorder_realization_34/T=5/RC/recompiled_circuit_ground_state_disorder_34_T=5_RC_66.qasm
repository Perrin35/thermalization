OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4021969) q[0];
sx q[0];
rz(-0.88391179) q[0];
sx q[0];
rz(-2.28595) q[0];
rz(-0.17172509) q[1];
sx q[1];
rz(-0.11556927) q[1];
sx q[1];
rz(2.5791383) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2320975) q[0];
sx q[0];
rz(-1.2065071) q[0];
sx q[0];
rz(3.0845736) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2954584) q[2];
sx q[2];
rz(-0.52342969) q[2];
sx q[2];
rz(-1.1688423) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6537498) q[1];
sx q[1];
rz(-1.1636793) q[1];
sx q[1];
rz(-0.458325) q[1];
x q[2];
rz(0.39333087) q[3];
sx q[3];
rz(-1.9489904) q[3];
sx q[3];
rz(-0.26547576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.90613753) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(0.79929024) q[2];
rz(0.47131395) q[3];
sx q[3];
rz(-0.91336942) q[3];
sx q[3];
rz(-2.2057064) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039463194) q[0];
sx q[0];
rz(-1.3251745) q[0];
sx q[0];
rz(-1.348173) q[0];
rz(1.3735636) q[1];
sx q[1];
rz(-1.9919688) q[1];
sx q[1];
rz(1.1522393) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3593529) q[0];
sx q[0];
rz(-2.3639661) q[0];
sx q[0];
rz(-1.1471143) q[0];
x q[1];
rz(2.197108) q[2];
sx q[2];
rz(-2.1462893) q[2];
sx q[2];
rz(-1.4552417) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.79346953) q[1];
sx q[1];
rz(-1.4033699) q[1];
sx q[1];
rz(2.4504285) q[1];
rz(-2.3186734) q[3];
sx q[3];
rz(-0.69861087) q[3];
sx q[3];
rz(-2.5320092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.79933244) q[2];
sx q[2];
rz(-0.41877425) q[2];
sx q[2];
rz(0.93117923) q[2];
rz(3.0350507) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(-0.60025269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0837285) q[0];
sx q[0];
rz(-1.3902384) q[0];
sx q[0];
rz(-0.51783836) q[0];
rz(2.2528516) q[1];
sx q[1];
rz(-0.70777142) q[1];
sx q[1];
rz(2.6944366) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5260552) q[0];
sx q[0];
rz(-1.8319905) q[0];
sx q[0];
rz(-1.5047856) q[0];
x q[1];
rz(0.60751407) q[2];
sx q[2];
rz(-2.315763) q[2];
sx q[2];
rz(-1.5811282) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.67990002) q[1];
sx q[1];
rz(-0.60631207) q[1];
sx q[1];
rz(-1.8099422) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6154061) q[3];
sx q[3];
rz(-1.8136214) q[3];
sx q[3];
rz(-0.015004166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7670224) q[2];
sx q[2];
rz(-0.76408237) q[2];
sx q[2];
rz(0.53263295) q[2];
rz(-2.8940708) q[3];
sx q[3];
rz(-2.4041924) q[3];
sx q[3];
rz(2.1029162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-1.5175051) q[0];
sx q[0];
rz(-3.0563323) q[0];
sx q[0];
rz(-3.0349773) q[0];
rz(-0.34635776) q[1];
sx q[1];
rz(-0.84500161) q[1];
sx q[1];
rz(-1.1603629) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3702404) q[0];
sx q[0];
rz(-1.816907) q[0];
sx q[0];
rz(-1.3918124) q[0];
x q[1];
rz(2.3290655) q[2];
sx q[2];
rz(-0.90355325) q[2];
sx q[2];
rz(-3.0548801) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9306158) q[1];
sx q[1];
rz(-1.9257571) q[1];
sx q[1];
rz(2.8205754) q[1];
rz(-pi) q[2];
rz(2.6909573) q[3];
sx q[3];
rz(-2.051578) q[3];
sx q[3];
rz(-1.4522875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13566636) q[2];
sx q[2];
rz(-2.8432507) q[2];
sx q[2];
rz(-3.0653817) q[2];
rz(-0.58586079) q[3];
sx q[3];
rz(-1.1548837) q[3];
sx q[3];
rz(-1.3425672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0214486) q[0];
sx q[0];
rz(-1.3055389) q[0];
sx q[0];
rz(2.7914877) q[0];
rz(2.1856951) q[1];
sx q[1];
rz(-1.2799542) q[1];
sx q[1];
rz(-1.9116481) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0405925) q[0];
sx q[0];
rz(-0.44389566) q[0];
sx q[0];
rz(2.1568265) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7857871) q[2];
sx q[2];
rz(-2.9724398) q[2];
sx q[2];
rz(-0.40119888) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0240682) q[1];
sx q[1];
rz(-1.9000016) q[1];
sx q[1];
rz(-2.9419521) q[1];
rz(-pi) q[2];
rz(-0.1996207) q[3];
sx q[3];
rz(-0.99310447) q[3];
sx q[3];
rz(-2.9132089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.89796394) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(1.492307) q[2];
rz(0.40488511) q[3];
sx q[3];
rz(-2.3758774) q[3];
sx q[3];
rz(-0.83612061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.006007) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(2.5591922) q[0];
rz(2.4549585) q[1];
sx q[1];
rz(-1.735894) q[1];
sx q[1];
rz(-1.2581717) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1863886) q[0];
sx q[0];
rz(-2.0058647) q[0];
sx q[0];
rz(-0.13829903) q[0];
x q[1];
rz(2.1297087) q[2];
sx q[2];
rz(-0.61043533) q[2];
sx q[2];
rz(1.8545146) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88731474) q[1];
sx q[1];
rz(-1.6266903) q[1];
sx q[1];
rz(3.0718551) q[1];
x q[2];
rz(-1.889713) q[3];
sx q[3];
rz(-2.5237759) q[3];
sx q[3];
rz(1.8998673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.76576343) q[2];
sx q[2];
rz(-1.9932237) q[2];
sx q[2];
rz(-1.0180391) q[2];
rz(0.24614075) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(-1.5181946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70596424) q[0];
sx q[0];
rz(-2.3376597) q[0];
sx q[0];
rz(-0.58018082) q[0];
rz(2.9980581) q[1];
sx q[1];
rz(-0.48964557) q[1];
sx q[1];
rz(2.9339583) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.872179) q[0];
sx q[0];
rz(-0.84194833) q[0];
sx q[0];
rz(2.6711885) q[0];
rz(2.0775547) q[2];
sx q[2];
rz(-2.4993976) q[2];
sx q[2];
rz(-2.3395777) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5112111) q[1];
sx q[1];
rz(-0.89511725) q[1];
sx q[1];
rz(-0.33501321) q[1];
x q[2];
rz(2.7782441) q[3];
sx q[3];
rz(-1.4814113) q[3];
sx q[3];
rz(0.18401981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.84270728) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(0.60619727) q[2];
rz(0.68743622) q[3];
sx q[3];
rz(-1.1339374) q[3];
sx q[3];
rz(0.70639759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
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
rz(0.48924482) q[0];
sx q[0];
rz(-0.027996538) q[0];
sx q[0];
rz(-2.0835173) q[0];
rz(-0.10969133) q[1];
sx q[1];
rz(-2.0249764) q[1];
sx q[1];
rz(1.6995957) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82175151) q[0];
sx q[0];
rz(-1.1761001) q[0];
sx q[0];
rz(-1.234353) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50050723) q[2];
sx q[2];
rz(-1.8138514) q[2];
sx q[2];
rz(0.99832035) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.8311288) q[1];
sx q[1];
rz(-2.052124) q[1];
sx q[1];
rz(-0.83408611) q[1];
rz(-pi) q[2];
rz(-2.2820246) q[3];
sx q[3];
rz(-0.35214927) q[3];
sx q[3];
rz(-1.6363877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.66120061) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(-1.2083758) q[2];
rz(-0.66323534) q[3];
sx q[3];
rz(-1.6845208) q[3];
sx q[3];
rz(1.9251582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97063589) q[0];
sx q[0];
rz(-0.96681505) q[0];
sx q[0];
rz(-3.048625) q[0];
rz(1.8611106) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(-3.0063937) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7057892) q[0];
sx q[0];
rz(-3.0717141) q[0];
sx q[0];
rz(2.4326434) q[0];
rz(-pi) q[1];
rz(-0.2580169) q[2];
sx q[2];
rz(-1.4224367) q[2];
sx q[2];
rz(-0.49000636) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7635423) q[1];
sx q[1];
rz(-1.592118) q[1];
sx q[1];
rz(2.887602) q[1];
rz(1.1208833) q[3];
sx q[3];
rz(-2.9723047) q[3];
sx q[3];
rz(1.560488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.4325503) q[2];
sx q[2];
rz(-1.9229527) q[2];
sx q[2];
rz(-2.3700628) q[2];
rz(0.49154526) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(2.5468723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58525697) q[0];
sx q[0];
rz(-0.62194967) q[0];
sx q[0];
rz(2.0167895) q[0];
rz(-1.2987761) q[1];
sx q[1];
rz(-2.5262084) q[1];
sx q[1];
rz(2.3769456) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1021834) q[0];
sx q[0];
rz(-1.4847857) q[0];
sx q[0];
rz(3.0360704) q[0];
x q[1];
rz(2.2510193) q[2];
sx q[2];
rz(-1.169551) q[2];
sx q[2];
rz(0.97637343) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.56988) q[1];
sx q[1];
rz(-0.24001828) q[1];
sx q[1];
rz(1.0102788) q[1];
rz(-0.57101698) q[3];
sx q[3];
rz(-0.63435508) q[3];
sx q[3];
rz(0.77658021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3760959) q[2];
sx q[2];
rz(-1.8509879) q[2];
sx q[2];
rz(2.9343904) q[2];
rz(-2.1616705) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(0.19237147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1205263) q[0];
sx q[0];
rz(-1.4561894) q[0];
sx q[0];
rz(-0.86984632) q[0];
rz(0.5008685) q[1];
sx q[1];
rz(-2.905838) q[1];
sx q[1];
rz(0.9314608) q[1];
rz(1.6899213) q[2];
sx q[2];
rz(-1.7075734) q[2];
sx q[2];
rz(1.3592958) q[2];
rz(-1.408314) q[3];
sx q[3];
rz(-1.8058895) q[3];
sx q[3];
rz(-1.4598082) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
