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
rz(2.565413) q[0];
sx q[0];
rz(-1.6771069) q[0];
sx q[0];
rz(-2.488774) q[0];
rz(2.6989812) q[1];
sx q[1];
rz(-1.1224597) q[1];
sx q[1];
rz(0.62243661) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5028337) q[0];
sx q[0];
rz(-1.2593102) q[0];
sx q[0];
rz(-3.019096) q[0];
x q[1];
rz(-1.9046225) q[2];
sx q[2];
rz(-2.1458652) q[2];
sx q[2];
rz(2.0890369) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.7622432) q[1];
sx q[1];
rz(-2.6062638) q[1];
sx q[1];
rz(1.9368656) q[1];
rz(-pi) q[2];
rz(-1.7951598) q[3];
sx q[3];
rz(-2.5208559) q[3];
sx q[3];
rz(-0.94514314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.5114674) q[2];
sx q[2];
rz(-1.7897391) q[2];
sx q[2];
rz(-2.5959065) q[2];
rz(-0.68583471) q[3];
sx q[3];
rz(-2.4617709) q[3];
sx q[3];
rz(2.1849476) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.086455258) q[0];
sx q[0];
rz(-0.70239037) q[0];
sx q[0];
rz(-1.0954683) q[0];
rz(2.144004) q[1];
sx q[1];
rz(-1.9065403) q[1];
sx q[1];
rz(1.9445317) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2491978) q[0];
sx q[0];
rz(-0.64606842) q[0];
sx q[0];
rz(1.9996253) q[0];
x q[1];
rz(2.1904693) q[2];
sx q[2];
rz(-2.24555) q[2];
sx q[2];
rz(-0.63979606) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.61651308) q[1];
sx q[1];
rz(-2.5933806) q[1];
sx q[1];
rz(1.0031149) q[1];
rz(-pi) q[2];
rz(-1.3838816) q[3];
sx q[3];
rz(-1.8946365) q[3];
sx q[3];
rz(-1.327654) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.41584388) q[2];
sx q[2];
rz(-1.7472605) q[2];
sx q[2];
rz(1.7572629) q[2];
rz(-1.3437126) q[3];
sx q[3];
rz(-1.5680771) q[3];
sx q[3];
rz(1.2725007) q[3];
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
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62190732) q[0];
sx q[0];
rz(-2.3540731) q[0];
sx q[0];
rz(2.7144879) q[0];
rz(1.2318132) q[1];
sx q[1];
rz(-1.6424664) q[1];
sx q[1];
rz(-0.61418358) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2285143) q[0];
sx q[0];
rz(-2.7102463) q[0];
sx q[0];
rz(-2.449663) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9236058) q[2];
sx q[2];
rz(-2.3828016) q[2];
sx q[2];
rz(-1.2365149) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8607448) q[1];
sx q[1];
rz(-3.0079746) q[1];
sx q[1];
rz(1.0722776) q[1];
x q[2];
rz(2.6355686) q[3];
sx q[3];
rz(-1.3576686) q[3];
sx q[3];
rz(-2.0863078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4517639) q[2];
sx q[2];
rz(-0.58388766) q[2];
sx q[2];
rz(0.27099398) q[2];
rz(-2.6594035) q[3];
sx q[3];
rz(-2.2995583) q[3];
sx q[3];
rz(1.2838001) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.86376205) q[0];
sx q[0];
rz(-1.3348802) q[0];
sx q[0];
rz(-0.41896391) q[0];
rz(2.9169182) q[1];
sx q[1];
rz(-1.9326262) q[1];
sx q[1];
rz(-0.077542543) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6200752) q[0];
sx q[0];
rz(-1.2780315) q[0];
sx q[0];
rz(-1.0241246) q[0];
rz(-pi) q[1];
rz(1.2214022) q[2];
sx q[2];
rz(-0.82942334) q[2];
sx q[2];
rz(1.0238613) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.71633881) q[1];
sx q[1];
rz(-1.1341061) q[1];
sx q[1];
rz(1.4136057) q[1];
rz(0.29234286) q[3];
sx q[3];
rz(-1.7708395) q[3];
sx q[3];
rz(-0.26250473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.055723995) q[2];
sx q[2];
rz(-2.8102977) q[2];
sx q[2];
rz(-1.2223988) q[2];
rz(-0.27082768) q[3];
sx q[3];
rz(-1.2374249) q[3];
sx q[3];
rz(0.45473155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5772575) q[0];
sx q[0];
rz(-2.1463647) q[0];
sx q[0];
rz(1.7876392) q[0];
rz(2.1063781) q[1];
sx q[1];
rz(-1.926492) q[1];
sx q[1];
rz(-1.6114906) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.852762) q[0];
sx q[0];
rz(-2.0977328) q[0];
sx q[0];
rz(-0.88165347) q[0];
x q[1];
rz(1.7351088) q[2];
sx q[2];
rz(-1.8949869) q[2];
sx q[2];
rz(-2.6099082) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.93188796) q[1];
sx q[1];
rz(-1.6242289) q[1];
sx q[1];
rz(-0.40219743) q[1];
x q[2];
rz(-3.0102481) q[3];
sx q[3];
rz(-2.7274611) q[3];
sx q[3];
rz(1.1826178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.80663854) q[2];
sx q[2];
rz(-2.1636476) q[2];
sx q[2];
rz(-3.0777001) q[2];
rz(2.4334) q[3];
sx q[3];
rz(-1.6876561) q[3];
sx q[3];
rz(1.3341058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42234364) q[0];
sx q[0];
rz(-0.023119211) q[0];
sx q[0];
rz(1.0759906) q[0];
rz(3.0883582) q[1];
sx q[1];
rz(-0.85317555) q[1];
sx q[1];
rz(2.7630973) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4814893) q[0];
sx q[0];
rz(-1.3162587) q[0];
sx q[0];
rz(1.7053563) q[0];
rz(-pi) q[1];
rz(-1.4827864) q[2];
sx q[2];
rz(-1.529379) q[2];
sx q[2];
rz(0.27829042) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.15684756) q[1];
sx q[1];
rz(-1.7498921) q[1];
sx q[1];
rz(-2.8791134) q[1];
rz(-pi) q[2];
x q[2];
rz(0.58121292) q[3];
sx q[3];
rz(-0.85678116) q[3];
sx q[3];
rz(0.17457822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9516912) q[2];
sx q[2];
rz(-1.4611763) q[2];
sx q[2];
rz(2.107673) q[2];
rz(-0.90384358) q[3];
sx q[3];
rz(-0.90141064) q[3];
sx q[3];
rz(3.097539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1164301) q[0];
sx q[0];
rz(-1.6950386) q[0];
sx q[0];
rz(-0.59301162) q[0];
rz(3.016839) q[1];
sx q[1];
rz(-1.1589103) q[1];
sx q[1];
rz(-0.4932901) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6007227) q[0];
sx q[0];
rz(-2.713795) q[0];
sx q[0];
rz(0.85036253) q[0];
rz(-0.10693018) q[2];
sx q[2];
rz(-0.22935805) q[2];
sx q[2];
rz(1.3170674) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64831454) q[1];
sx q[1];
rz(-1.9946792) q[1];
sx q[1];
rz(-2.5878169) q[1];
rz(-1.4629355) q[3];
sx q[3];
rz(-1.2965805) q[3];
sx q[3];
rz(1.3317684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7571681) q[2];
sx q[2];
rz(-1.9809456) q[2];
sx q[2];
rz(-1.9942795) q[2];
rz(2.3186963) q[3];
sx q[3];
rz(-1.1742679) q[3];
sx q[3];
rz(0.66250044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5230474) q[0];
sx q[0];
rz(-1.86919) q[0];
sx q[0];
rz(0.4185032) q[0];
rz(-3.0283527) q[1];
sx q[1];
rz(-1.6721882) q[1];
sx q[1];
rz(-2.0638827) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5961469) q[0];
sx q[0];
rz(-1.2823063) q[0];
sx q[0];
rz(-0.01878564) q[0];
rz(-pi) q[1];
rz(-1.244721) q[2];
sx q[2];
rz(-2.3387944) q[2];
sx q[2];
rz(1.7434415) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.019494836) q[1];
sx q[1];
rz(-1.5217402) q[1];
sx q[1];
rz(-8*pi/13) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6304794) q[3];
sx q[3];
rz(-2.2264997) q[3];
sx q[3];
rz(2.7969517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.30617994) q[2];
sx q[2];
rz(-0.70859185) q[2];
sx q[2];
rz(-1.6861247) q[2];
rz(2.9742187) q[3];
sx q[3];
rz(-0.51515976) q[3];
sx q[3];
rz(-1.0719871) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6587332) q[0];
sx q[0];
rz(-1.8287683) q[0];
sx q[0];
rz(-0.40153781) q[0];
rz(0.43831476) q[1];
sx q[1];
rz(-0.79151789) q[1];
sx q[1];
rz(-1.4567136) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7244723) q[0];
sx q[0];
rz(-0.03532413) q[0];
sx q[0];
rz(1.8457343) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3517604) q[2];
sx q[2];
rz(-0.92680762) q[2];
sx q[2];
rz(-2.7119659) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.40565421) q[1];
sx q[1];
rz(-0.88684139) q[1];
sx q[1];
rz(-2.4962462) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51299398) q[3];
sx q[3];
rz(-1.1502642) q[3];
sx q[3];
rz(0.094407439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.38176408) q[2];
sx q[2];
rz(-2.9324053) q[2];
sx q[2];
rz(-0.71167243) q[2];
rz(3.107374) q[3];
sx q[3];
rz(-2.4331369) q[3];
sx q[3];
rz(1.155352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3592247) q[0];
sx q[0];
rz(-1.6509667) q[0];
sx q[0];
rz(0.79886287) q[0];
rz(1.0460188) q[1];
sx q[1];
rz(-1.1963528) q[1];
sx q[1];
rz(0.23652133) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2940148) q[0];
sx q[0];
rz(-1.6315407) q[0];
sx q[0];
rz(-0.022367064) q[0];
rz(-pi) q[1];
rz(0.78153083) q[2];
sx q[2];
rz(-2.3252914) q[2];
sx q[2];
rz(0.81354173) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7812209) q[1];
sx q[1];
rz(-2.0016333) q[1];
sx q[1];
rz(-1.197351) q[1];
rz(-pi) q[2];
rz(-1.4528689) q[3];
sx q[3];
rz(-2.2783359) q[3];
sx q[3];
rz(0.70048571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.11655) q[2];
sx q[2];
rz(-2.5849085) q[2];
sx q[2];
rz(0.64209783) q[2];
rz(0.56946483) q[3];
sx q[3];
rz(-2.6537708) q[3];
sx q[3];
rz(0.08629442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(1.4702598) q[0];
sx q[0];
rz(-1.3339806) q[0];
sx q[0];
rz(0.95680923) q[0];
rz(1.9500465) q[1];
sx q[1];
rz(-1.5622495) q[1];
sx q[1];
rz(-1.539485) q[1];
rz(-1.446522) q[2];
sx q[2];
rz(-2.3588603) q[2];
sx q[2];
rz(-1.0609577) q[2];
rz(-1.1453876) q[3];
sx q[3];
rz(-1.8775826) q[3];
sx q[3];
rz(2.0043397) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
