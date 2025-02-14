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
rz(3.6917917) q[0];
sx q[0];
rz(10.79296) q[0];
rz(-0.47222459) q[1];
sx q[1];
rz(-3.0250186) q[1];
sx q[1];
rz(-0.20224686) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2496754) q[0];
sx q[0];
rz(-1.7698827) q[0];
sx q[0];
rz(-2.9249935) q[0];
x q[1];
rz(-2.6361385) q[2];
sx q[2];
rz(-1.1040282) q[2];
sx q[2];
rz(-0.17077489) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.56204462) q[1];
sx q[1];
rz(-1.8534396) q[1];
sx q[1];
rz(0.94448994) q[1];
rz(-3.1038001) q[3];
sx q[3];
rz(-1.1049912) q[3];
sx q[3];
rz(-2.56063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5998485) q[2];
sx q[2];
rz(-0.84459633) q[2];
sx q[2];
rz(-2.2005626) q[2];
rz(-2.8372676) q[3];
sx q[3];
rz(-0.86186886) q[3];
sx q[3];
rz(2.8743675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32352725) q[0];
sx q[0];
rz(-2.7637988) q[0];
sx q[0];
rz(0.7793119) q[0];
rz(-1.7794973) q[1];
sx q[1];
rz(-0.39924386) q[1];
sx q[1];
rz(-1.13387) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88863605) q[0];
sx q[0];
rz(-1.1373113) q[0];
sx q[0];
rz(-0.10826464) q[0];
rz(3.0236565) q[2];
sx q[2];
rz(-1.6291766) q[2];
sx q[2];
rz(2.749325) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8633921) q[1];
sx q[1];
rz(-0.28732936) q[1];
sx q[1];
rz(1.1817688) q[1];
rz(-2.0800703) q[3];
sx q[3];
rz(-1.5716553) q[3];
sx q[3];
rz(0.96405503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4506932) q[2];
sx q[2];
rz(-1.0319812) q[2];
sx q[2];
rz(2.9449985) q[2];
rz(-0.96902668) q[3];
sx q[3];
rz(-1.6220379) q[3];
sx q[3];
rz(-1.2381747) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75152385) q[0];
sx q[0];
rz(-0.5539493) q[0];
sx q[0];
rz(2.4098136) q[0];
rz(2.7732908) q[1];
sx q[1];
rz(-2.383547) q[1];
sx q[1];
rz(-0.36645737) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.059341934) q[0];
sx q[0];
rz(-1.8346795) q[0];
sx q[0];
rz(-1.7543955) q[0];
rz(-pi) q[1];
rz(-0.24213893) q[2];
sx q[2];
rz(-1.79091) q[2];
sx q[2];
rz(-1.5907703) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0302387) q[1];
sx q[1];
rz(-1.0525647) q[1];
sx q[1];
rz(2.6494762) q[1];
x q[2];
rz(0.89160925) q[3];
sx q[3];
rz(-0.82538939) q[3];
sx q[3];
rz(2.2872383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.1825819) q[2];
sx q[2];
rz(-2.3704447) q[2];
sx q[2];
rz(0.89243531) q[2];
rz(-2.6871032) q[3];
sx q[3];
rz(-2.317704) q[3];
sx q[3];
rz(0.23196001) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3733805) q[0];
sx q[0];
rz(-2.9990271) q[0];
sx q[0];
rz(2.7349803) q[0];
rz(-2.3136102) q[1];
sx q[1];
rz(-1.7384638) q[1];
sx q[1];
rz(-0.83555317) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8260561) q[0];
sx q[0];
rz(-2.9565507) q[0];
sx q[0];
rz(-1.3547784) q[0];
rz(1.0539444) q[2];
sx q[2];
rz(-0.39363855) q[2];
sx q[2];
rz(0.37154276) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.20488508) q[1];
sx q[1];
rz(-0.46006535) q[1];
sx q[1];
rz(-2.8060444) q[1];
rz(-pi) q[2];
rz(1.8074513) q[3];
sx q[3];
rz(-0.39517212) q[3];
sx q[3];
rz(-0.69505668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.39997175) q[2];
sx q[2];
rz(-0.30210945) q[2];
sx q[2];
rz(-0.92140222) q[2];
rz(1.7737927) q[3];
sx q[3];
rz(-0.96145815) q[3];
sx q[3];
rz(0.14149806) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2530186) q[0];
sx q[0];
rz(-2.044401) q[0];
sx q[0];
rz(2.9827523) q[0];
rz(1.5171492) q[1];
sx q[1];
rz(-2.8327063) q[1];
sx q[1];
rz(-1.3295757) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3314574) q[0];
sx q[0];
rz(-1.5525426) q[0];
sx q[0];
rz(-1.5601331) q[0];
rz(-pi) q[1];
rz(2.9022568) q[2];
sx q[2];
rz(-1.0221326) q[2];
sx q[2];
rz(1.2557097) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.5097039) q[1];
sx q[1];
rz(-1.3702127) q[1];
sx q[1];
rz(1.1840118) q[1];
rz(-pi) q[2];
rz(-0.20820372) q[3];
sx q[3];
rz(-2.1284687) q[3];
sx q[3];
rz(-0.69277908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.57881957) q[2];
sx q[2];
rz(-0.69391888) q[2];
sx q[2];
rz(0.55207437) q[2];
rz(2.3479346) q[3];
sx q[3];
rz(-2.5439883) q[3];
sx q[3];
rz(-0.98646599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3119222) q[0];
sx q[0];
rz(-0.86450082) q[0];
sx q[0];
rz(-0.70575869) q[0];
rz(-0.13748473) q[1];
sx q[1];
rz(-2.4973713) q[1];
sx q[1];
rz(-0.71298832) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1437837) q[0];
sx q[0];
rz(-1.1857872) q[0];
sx q[0];
rz(1.3237756) q[0];
x q[1];
rz(1.1399066) q[2];
sx q[2];
rz(-2.0451114) q[2];
sx q[2];
rz(2.9805019) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.5755427) q[1];
sx q[1];
rz(-2.7311196) q[1];
sx q[1];
rz(1.3259352) q[1];
x q[2];
rz(1.8637795) q[3];
sx q[3];
rz(-1.6821386) q[3];
sx q[3];
rz(2.7822943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2116427) q[2];
sx q[2];
rz(-2.2723618) q[2];
sx q[2];
rz(0.28406528) q[2];
rz(3.0197213) q[3];
sx q[3];
rz(-2.8594696) q[3];
sx q[3];
rz(1.4819283) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26509869) q[0];
sx q[0];
rz(-0.6186741) q[0];
sx q[0];
rz(-0.70736831) q[0];
rz(0.44037285) q[1];
sx q[1];
rz(-2.423954) q[1];
sx q[1];
rz(1.3622989) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4934568) q[0];
sx q[0];
rz(-1.8168983) q[0];
sx q[0];
rz(-1.2854693) q[0];
rz(-pi) q[1];
rz(-2.0411891) q[2];
sx q[2];
rz(-0.78436034) q[2];
sx q[2];
rz(-0.91299865) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.8352812) q[1];
sx q[1];
rz(-2.0871442) q[1];
sx q[1];
rz(0.55037127) q[1];
rz(-0.23162095) q[3];
sx q[3];
rz(-0.4670139) q[3];
sx q[3];
rz(1.8927595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.97544396) q[2];
sx q[2];
rz(-2.7232309) q[2];
sx q[2];
rz(-2.5725906) q[2];
rz(2.1042018) q[3];
sx q[3];
rz(-0.85809696) q[3];
sx q[3];
rz(-0.4666127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5910832) q[0];
sx q[0];
rz(-2.789848) q[0];
sx q[0];
rz(1.2048703) q[0];
rz(-0.51271802) q[1];
sx q[1];
rz(-2.3725489) q[1];
sx q[1];
rz(-0.18096322) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34535392) q[0];
sx q[0];
rz(-1.1974044) q[0];
sx q[0];
rz(0.42723421) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4333192) q[2];
sx q[2];
rz(-1.9172693) q[2];
sx q[2];
rz(2.8203143) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9258283) q[1];
sx q[1];
rz(-1.3604394) q[1];
sx q[1];
rz(1.0025566) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0309903) q[3];
sx q[3];
rz(-1.4975274) q[3];
sx q[3];
rz(-1.6694809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.60244954) q[2];
sx q[2];
rz(-2.4931543) q[2];
sx q[2];
rz(3.0446206) q[2];
rz(-2.5868331) q[3];
sx q[3];
rz(-1.8192889) q[3];
sx q[3];
rz(2.7831912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.6612514) q[0];
sx q[0];
rz(-0.77965176) q[0];
sx q[0];
rz(0.10777792) q[0];
rz(0.55094552) q[1];
sx q[1];
rz(-2.7007553) q[1];
sx q[1];
rz(-1.5239747) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11883987) q[0];
sx q[0];
rz(-1.700239) q[0];
sx q[0];
rz(1.8047099) q[0];
rz(1.0650159) q[2];
sx q[2];
rz(-2.1518101) q[2];
sx q[2];
rz(-2.6083824) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1283933) q[1];
sx q[1];
rz(-2.6775762) q[1];
sx q[1];
rz(1.4739407) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0436936) q[3];
sx q[3];
rz(-0.49360858) q[3];
sx q[3];
rz(-0.61033397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.23065755) q[2];
sx q[2];
rz(-0.50369889) q[2];
sx q[2];
rz(1.1304193) q[2];
rz(-0.23760992) q[3];
sx q[3];
rz(-0.41046023) q[3];
sx q[3];
rz(-2.3835278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0216574) q[0];
sx q[0];
rz(-2.9627242) q[0];
sx q[0];
rz(-2.2005431) q[0];
rz(-0.36048105) q[1];
sx q[1];
rz(-1.4070114) q[1];
sx q[1];
rz(-1.9182659) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1038641) q[0];
sx q[0];
rz(-2.1489804) q[0];
sx q[0];
rz(-1.8072771) q[0];
rz(-pi) q[1];
rz(-0.98127301) q[2];
sx q[2];
rz(-2.2726577) q[2];
sx q[2];
rz(-1.5379932) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.0666484) q[1];
sx q[1];
rz(-2.7171591) q[1];
sx q[1];
rz(-2.9702999) q[1];
x q[2];
rz(-1.8935325) q[3];
sx q[3];
rz(-0.56005037) q[3];
sx q[3];
rz(1.3572451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90784043) q[2];
sx q[2];
rz(-1.7641822) q[2];
sx q[2];
rz(3.0506548) q[2];
rz(2.9371373) q[3];
sx q[3];
rz(-0.8046059) q[3];
sx q[3];
rz(2.0596152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7627056) q[0];
sx q[0];
rz(-1.611041) q[0];
sx q[0];
rz(1.8086717) q[0];
rz(1.7145722) q[1];
sx q[1];
rz(-2.6418229) q[1];
sx q[1];
rz(-1.5508834) q[1];
rz(-1.0574404) q[2];
sx q[2];
rz(-2.3427137) q[2];
sx q[2];
rz(-2.8165934) q[2];
rz(-1.5132202) q[3];
sx q[3];
rz(-1.7840458) q[3];
sx q[3];
rz(1.9174867) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
