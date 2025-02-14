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
rz(1.053831) q[0];
sx q[0];
rz(-2.5913936) q[0];
sx q[0];
rz(1.3681816) q[0];
rz(-0.47222459) q[1];
sx q[1];
rz(-3.0250186) q[1];
sx q[1];
rz(-0.20224686) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2496754) q[0];
sx q[0];
rz(-1.37171) q[0];
sx q[0];
rz(2.9249935) q[0];
rz(2.0933242) q[2];
sx q[2];
rz(-1.1236345) q[2];
sx q[2];
rz(1.9856404) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.579548) q[1];
sx q[1];
rz(-1.8534396) q[1];
sx q[1];
rz(-2.1971027) q[1];
rz(-1.6458166) q[3];
sx q[3];
rz(-2.6743691) q[3];
sx q[3];
rz(-0.4969767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.54174417) q[2];
sx q[2];
rz(-2.2969963) q[2];
sx q[2];
rz(2.2005626) q[2];
rz(-0.30432501) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(0.26722515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8180654) q[0];
sx q[0];
rz(-0.37779385) q[0];
sx q[0];
rz(-2.3622808) q[0];
rz(1.7794973) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(-2.0077226) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2529566) q[0];
sx q[0];
rz(-1.1373113) q[0];
sx q[0];
rz(-3.033328) q[0];
x q[1];
rz(-3.0236565) q[2];
sx q[2];
rz(-1.5124161) q[2];
sx q[2];
rz(-0.3922677) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.081959978) q[1];
sx q[1];
rz(-1.4631008) q[1];
sx q[1];
rz(-1.8376985) q[1];
rz(1.5725582) q[3];
sx q[3];
rz(-0.50927466) q[3];
sx q[3];
rz(2.5363896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.69089943) q[2];
sx q[2];
rz(-2.1096114) q[2];
sx q[2];
rz(-2.9449985) q[2];
rz(2.172566) q[3];
sx q[3];
rz(-1.6220379) q[3];
sx q[3];
rz(-1.2381747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75152385) q[0];
sx q[0];
rz(-2.5876434) q[0];
sx q[0];
rz(-0.73177904) q[0];
rz(-0.36830184) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(-0.36645737) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0822507) q[0];
sx q[0];
rz(-1.3069131) q[0];
sx q[0];
rz(1.7543955) q[0];
rz(0.24213893) q[2];
sx q[2];
rz(-1.3506827) q[2];
sx q[2];
rz(1.5508224) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4225821) q[1];
sx q[1];
rz(-1.9937464) q[1];
sx q[1];
rz(-2.1450536) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2499834) q[3];
sx q[3];
rz(-0.82538939) q[3];
sx q[3];
rz(2.2872383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1825819) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(-2.2491573) q[2];
rz(-2.6871032) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(0.23196001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7682122) q[0];
sx q[0];
rz(-2.9990271) q[0];
sx q[0];
rz(2.7349803) q[0];
rz(0.82798249) q[1];
sx q[1];
rz(-1.7384638) q[1];
sx q[1];
rz(2.3060395) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31553651) q[0];
sx q[0];
rz(-0.18504194) q[0];
sx q[0];
rz(1.3547784) q[0];
x q[1];
rz(-2.0876483) q[2];
sx q[2];
rz(-0.39363855) q[2];
sx q[2];
rz(0.37154276) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20488508) q[1];
sx q[1];
rz(-0.46006535) q[1];
sx q[1];
rz(-0.33554828) q[1];
x q[2];
rz(-1.1855679) q[3];
sx q[3];
rz(-1.6611757) q[3];
sx q[3];
rz(2.4848695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.7416209) q[2];
sx q[2];
rz(-2.8394832) q[2];
sx q[2];
rz(0.92140222) q[2];
rz(-1.7737927) q[3];
sx q[3];
rz(-2.1801345) q[3];
sx q[3];
rz(0.14149806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2530186) q[0];
sx q[0];
rz(-1.0971917) q[0];
sx q[0];
rz(-0.15884037) q[0];
rz(1.6244434) q[1];
sx q[1];
rz(-0.30888638) q[1];
sx q[1];
rz(1.812017) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81013524) q[0];
sx q[0];
rz(-1.5890501) q[0];
sx q[0];
rz(-1.5814596) q[0];
rz(-pi) q[1];
rz(2.1324124) q[2];
sx q[2];
rz(-1.3671285) q[2];
sx q[2];
rz(-0.18850279) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.6318887) q[1];
sx q[1];
rz(-1.77138) q[1];
sx q[1];
rz(1.9575809) q[1];
rz(-pi) q[2];
rz(-2.1383189) q[3];
sx q[3];
rz(-1.7470932) q[3];
sx q[3];
rz(2.3749173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.57881957) q[2];
sx q[2];
rz(-0.69391888) q[2];
sx q[2];
rz(-0.55207437) q[2];
rz(-2.3479346) q[3];
sx q[3];
rz(-0.59760439) q[3];
sx q[3];
rz(-0.98646599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3119222) q[0];
sx q[0];
rz(-0.86450082) q[0];
sx q[0];
rz(0.70575869) q[0];
rz(-3.0041079) q[1];
sx q[1];
rz(-0.64422137) q[1];
sx q[1];
rz(2.4286043) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1437837) q[0];
sx q[0];
rz(-1.1857872) q[0];
sx q[0];
rz(1.817817) q[0];
rz(0.6829459) q[2];
sx q[2];
rz(-0.6295528) q[2];
sx q[2];
rz(0.94972875) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.56604993) q[1];
sx q[1];
rz(-2.7311196) q[1];
sx q[1];
rz(-1.8156575) q[1];
x q[2];
rz(-0.1162545) q[3];
sx q[3];
rz(-1.861912) q[3];
sx q[3];
rz(-1.8965885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9299499) q[2];
sx q[2];
rz(-0.86923081) q[2];
sx q[2];
rz(-2.8575274) q[2];
rz(-0.12187135) q[3];
sx q[3];
rz(-2.8594696) q[3];
sx q[3];
rz(-1.6596644) q[3];
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
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.876494) q[0];
sx q[0];
rz(-0.6186741) q[0];
sx q[0];
rz(0.70736831) q[0];
rz(2.7012198) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(-1.3622989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4934568) q[0];
sx q[0];
rz(-1.3246944) q[0];
sx q[0];
rz(1.2854693) q[0];
rz(-0.84378924) q[2];
sx q[2];
rz(-1.2449045) q[2];
sx q[2];
rz(2.1383204) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.3063115) q[1];
sx q[1];
rz(-2.0871442) q[1];
sx q[1];
rz(0.55037127) q[1];
x q[2];
rz(-2.6853722) q[3];
sx q[3];
rz(-1.6743321) q[3];
sx q[3];
rz(0.52952784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.97544396) q[2];
sx q[2];
rz(-0.41836172) q[2];
sx q[2];
rz(2.5725906) q[2];
rz(-1.0373908) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(2.67498) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55050945) q[0];
sx q[0];
rz(-2.789848) q[0];
sx q[0];
rz(-1.2048703) q[0];
rz(-0.51271802) q[1];
sx q[1];
rz(-0.7690438) q[1];
sx q[1];
rz(0.18096322) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34535392) q[0];
sx q[0];
rz(-1.9441883) q[0];
sx q[0];
rz(2.7143584) q[0];
rz(-pi) q[1];
rz(0.36278533) q[2];
sx q[2];
rz(-2.7698667) q[2];
sx q[2];
rz(3.0759916) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1026336) q[1];
sx q[1];
rz(-2.5397083) q[1];
sx q[1];
rz(-1.9485056) q[1];
x q[2];
rz(1.4970786) q[3];
sx q[3];
rz(-1.6811007) q[3];
sx q[3];
rz(-3.0347787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.5391431) q[2];
sx q[2];
rz(-0.64843833) q[2];
sx q[2];
rz(0.096972018) q[2];
rz(-0.55475956) q[3];
sx q[3];
rz(-1.8192889) q[3];
sx q[3];
rz(-2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6612514) q[0];
sx q[0];
rz(-2.3619409) q[0];
sx q[0];
rz(-0.10777792) q[0];
rz(-2.5906471) q[1];
sx q[1];
rz(-0.44083732) q[1];
sx q[1];
rz(-1.617618) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0227528) q[0];
sx q[0];
rz(-1.4413537) q[0];
sx q[0];
rz(1.3368827) q[0];
rz(-pi) q[1];
x q[1];
rz(0.6438821) q[2];
sx q[2];
rz(-1.9877627) q[2];
sx q[2];
rz(-1.3326926) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2366166) q[1];
sx q[1];
rz(-1.1091241) q[1];
sx q[1];
rz(-0.048358976) q[1];
rz(-pi) q[2];
x q[2];
rz(0.097899072) q[3];
sx q[3];
rz(-2.6479841) q[3];
sx q[3];
rz(0.61033397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23065755) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(-1.1304193) q[2];
rz(0.23760992) q[3];
sx q[3];
rz(-0.41046023) q[3];
sx q[3];
rz(2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.0216574) q[0];
sx q[0];
rz(-0.1788685) q[0];
sx q[0];
rz(-2.2005431) q[0];
rz(2.7811116) q[1];
sx q[1];
rz(-1.4070114) q[1];
sx q[1];
rz(1.2233268) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1038641) q[0];
sx q[0];
rz(-0.99261221) q[0];
sx q[0];
rz(1.8072771) q[0];
x q[1];
rz(0.98127301) q[2];
sx q[2];
rz(-0.86893493) q[2];
sx q[2];
rz(1.6035994) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4810923) q[1];
sx q[1];
rz(-1.5005439) q[1];
sx q[1];
rz(0.41892799) q[1];
x q[2];
rz(1.2480601) q[3];
sx q[3];
rz(-0.56005037) q[3];
sx q[3];
rz(-1.7843475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90784043) q[2];
sx q[2];
rz(-1.3774104) q[2];
sx q[2];
rz(0.09093786) q[2];
rz(-0.20445538) q[3];
sx q[3];
rz(-2.3369868) q[3];
sx q[3];
rz(-2.0596152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37888708) q[0];
sx q[0];
rz(-1.5305516) q[0];
sx q[0];
rz(-1.332921) q[0];
rz(1.7145722) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(2.3007921) q[2];
sx q[2];
rz(-1.9304095) q[2];
sx q[2];
rz(-0.87113397) q[2];
rz(1.5132202) q[3];
sx q[3];
rz(-1.3575469) q[3];
sx q[3];
rz(-1.224106) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
