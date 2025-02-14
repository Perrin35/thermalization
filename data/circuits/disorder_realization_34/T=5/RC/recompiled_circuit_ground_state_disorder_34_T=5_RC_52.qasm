OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7393957) q[0];
sx q[0];
rz(4.0255044) q[0];
sx q[0];
rz(8.5691353) q[0];
rz(-0.17172509) q[1];
sx q[1];
rz(-0.11556927) q[1];
sx q[1];
rz(2.5791383) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2320975) q[0];
sx q[0];
rz(-1.9350855) q[0];
sx q[0];
rz(-0.057019071) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0777656) q[2];
sx q[2];
rz(-1.7071144) q[2];
sx q[2];
rz(2.9796114) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3823279) q[1];
sx q[1];
rz(-2.538343) q[1];
sx q[1];
rz(2.3690577) q[1];
rz(2.7482618) q[3];
sx q[3];
rz(-1.1926023) q[3];
sx q[3];
rz(2.8761169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.90613753) q[2];
sx q[2];
rz(-1.9635341) q[2];
sx q[2];
rz(0.79929024) q[2];
rz(-2.6702787) q[3];
sx q[3];
rz(-2.2282232) q[3];
sx q[3];
rz(2.2057064) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.039463194) q[0];
sx q[0];
rz(-1.3251745) q[0];
sx q[0];
rz(-1.348173) q[0];
rz(-1.3735636) q[1];
sx q[1];
rz(-1.1496239) q[1];
sx q[1];
rz(-1.9893533) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.099502787) q[0];
sx q[0];
rz(-1.2782017) q[0];
sx q[0];
rz(-2.3022404) q[0];
rz(-pi) q[1];
rz(2.197108) q[2];
sx q[2];
rz(-2.1462893) q[2];
sx q[2];
rz(-1.4552417) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5630398) q[1];
sx q[1];
rz(-0.7078979) q[1];
sx q[1];
rz(-2.8824174) q[1];
rz(-pi) q[2];
rz(-0.51898414) q[3];
sx q[3];
rz(-2.0618084) q[3];
sx q[3];
rz(-1.6512914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.3422602) q[2];
sx q[2];
rz(-2.7228184) q[2];
sx q[2];
rz(2.2104134) q[2];
rz(3.0350507) q[3];
sx q[3];
rz(-1.1139161) q[3];
sx q[3];
rz(2.54134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.057864144) q[0];
sx q[0];
rz(-1.3902384) q[0];
sx q[0];
rz(2.6237543) q[0];
rz(2.2528516) q[1];
sx q[1];
rz(-2.4338212) q[1];
sx q[1];
rz(0.4471561) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0797821) q[0];
sx q[0];
rz(-1.5070276) q[0];
sx q[0];
rz(-0.26173862) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5340786) q[2];
sx q[2];
rz(-0.82582966) q[2];
sx q[2];
rz(1.5811282) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0886317) q[1];
sx q[1];
rz(-1.7061894) q[1];
sx q[1];
rz(-0.9779344) q[1];
x q[2];
rz(-0.52618653) q[3];
sx q[3];
rz(-1.8136214) q[3];
sx q[3];
rz(0.015004166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.37457028) q[2];
sx q[2];
rz(-2.3775103) q[2];
sx q[2];
rz(-0.53263295) q[2];
rz(2.8940708) q[3];
sx q[3];
rz(-2.4041924) q[3];
sx q[3];
rz(1.0386764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.6240876) q[0];
sx q[0];
rz(-3.0563323) q[0];
sx q[0];
rz(-3.0349773) q[0];
rz(2.7952349) q[1];
sx q[1];
rz(-0.84500161) q[1];
sx q[1];
rz(-1.1603629) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7713523) q[0];
sx q[0];
rz(-1.3246857) q[0];
sx q[0];
rz(1.3918124) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82616634) q[2];
sx q[2];
rz(-1.0001422) q[2];
sx q[2];
rz(-2.014239) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9306158) q[1];
sx q[1];
rz(-1.2158356) q[1];
sx q[1];
rz(0.32101722) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.87507803) q[3];
sx q[3];
rz(-0.64662537) q[3];
sx q[3];
rz(-0.88133206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0059263) q[2];
sx q[2];
rz(-0.29834193) q[2];
sx q[2];
rz(3.0653817) q[2];
rz(-2.5557319) q[3];
sx q[3];
rz(-1.1548837) q[3];
sx q[3];
rz(1.3425672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12014408) q[0];
sx q[0];
rz(-1.8360538) q[0];
sx q[0];
rz(-2.7914877) q[0];
rz(2.1856951) q[1];
sx q[1];
rz(-1.2799542) q[1];
sx q[1];
rz(-1.9116481) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53287017) q[0];
sx q[0];
rz(-1.9367095) q[0];
sx q[0];
rz(-0.25718148) q[0];
rz(1.3558055) q[2];
sx q[2];
rz(-2.9724398) q[2];
sx q[2];
rz(-2.7403938) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.753641) q[1];
sx q[1];
rz(-1.759583) q[1];
sx q[1];
rz(1.2354047) q[1];
x q[2];
rz(-1.2754945) q[3];
sx q[3];
rz(-0.60747889) q[3];
sx q[3];
rz(-2.5584084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.89796394) q[2];
sx q[2];
rz(-2.882759) q[2];
sx q[2];
rz(1.492307) q[2];
rz(-2.7367075) q[3];
sx q[3];
rz(-2.3758774) q[3];
sx q[3];
rz(2.305472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.006007) q[0];
sx q[0];
rz(-2.0642991) q[0];
sx q[0];
rz(-2.5591922) q[0];
rz(-0.68663418) q[1];
sx q[1];
rz(-1.4056987) q[1];
sx q[1];
rz(1.2581717) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63623896) q[0];
sx q[0];
rz(-0.4551783) q[0];
sx q[0];
rz(-1.8591465) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1297087) q[2];
sx q[2];
rz(-2.5311573) q[2];
sx q[2];
rz(1.287078) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.68738378) q[1];
sx q[1];
rz(-1.6404248) q[1];
sx q[1];
rz(-1.5147665) q[1];
rz(-pi) q[2];
rz(0.97719394) q[3];
sx q[3];
rz(-1.3881637) q[3];
sx q[3];
rz(3.0754418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.3758292) q[2];
sx q[2];
rz(-1.9932237) q[2];
sx q[2];
rz(-2.1235535) q[2];
rz(-0.24614075) q[3];
sx q[3];
rz(-1.7459511) q[3];
sx q[3];
rz(-1.6233981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70596424) q[0];
sx q[0];
rz(-0.80393296) q[0];
sx q[0];
rz(2.5614118) q[0];
rz(2.9980581) q[1];
sx q[1];
rz(-2.6519471) q[1];
sx q[1];
rz(0.20763436) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.62791) q[0];
sx q[0];
rz(-1.2259036) q[0];
sx q[0];
rz(-2.3570127) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0775547) q[2];
sx q[2];
rz(-2.4993976) q[2];
sx q[2];
rz(-2.3395777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.0032867) q[1];
sx q[1];
rz(-2.3992743) q[1];
sx q[1];
rz(-1.9600541) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.36334857) q[3];
sx q[3];
rz(-1.6601813) q[3];
sx q[3];
rz(2.9575728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.84270728) q[2];
sx q[2];
rz(-1.2879813) q[2];
sx q[2];
rz(-0.60619727) q[2];
rz(2.4541564) q[3];
sx q[3];
rz(-1.1339374) q[3];
sx q[3];
rz(-0.70639759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6523478) q[0];
sx q[0];
rz(-3.1135961) q[0];
sx q[0];
rz(2.0835173) q[0];
rz(3.0319013) q[1];
sx q[1];
rz(-2.0249764) q[1];
sx q[1];
rz(-1.441997) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82175151) q[0];
sx q[0];
rz(-1.9654925) q[0];
sx q[0];
rz(-1.234353) q[0];
x q[1];
rz(-2.6410854) q[2];
sx q[2];
rz(-1.3277413) q[2];
sx q[2];
rz(2.1432723) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8311288) q[1];
sx q[1];
rz(-1.0894686) q[1];
sx q[1];
rz(-0.83408611) q[1];
rz(-pi) q[2];
rz(1.2992925) q[3];
sx q[3];
rz(-1.3437004) q[3];
sx q[3];
rz(0.61448594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.66120061) q[2];
sx q[2];
rz(-0.19442393) q[2];
sx q[2];
rz(-1.9332168) q[2];
rz(-0.66323534) q[3];
sx q[3];
rz(-1.4570718) q[3];
sx q[3];
rz(-1.9251582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97063589) q[0];
sx q[0];
rz(-0.96681505) q[0];
sx q[0];
rz(-3.048625) q[0];
rz(-1.8611106) q[1];
sx q[1];
rz(-2.4130776) q[1];
sx q[1];
rz(3.0063937) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8427348) q[0];
sx q[0];
rz(-1.6162685) q[0];
sx q[0];
rz(-3.088515) q[0];
x q[1];
rz(2.6117295) q[2];
sx q[2];
rz(-2.8447897) q[2];
sx q[2];
rz(2.5713845) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1982806) q[1];
sx q[1];
rz(-1.824728) q[1];
sx q[1];
rz(-1.5928245) q[1];
x q[2];
rz(1.1208833) q[3];
sx q[3];
rz(-0.16928798) q[3];
sx q[3];
rz(1.5811046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4325503) q[2];
sx q[2];
rz(-1.21864) q[2];
sx q[2];
rz(-2.3700628) q[2];
rz(-0.49154526) q[3];
sx q[3];
rz(-1.2794269) q[3];
sx q[3];
rz(0.5947203) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
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
rz(1.2987761) q[1];
sx q[1];
rz(-0.61538428) q[1];
sx q[1];
rz(2.3769456) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0394093) q[0];
sx q[0];
rz(-1.656807) q[0];
sx q[0];
rz(-3.0360704) q[0];
rz(-pi) q[1];
rz(-0.9773639) q[2];
sx q[2];
rz(-0.77319169) q[2];
sx q[2];
rz(2.9969281) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5950039) q[1];
sx q[1];
rz(-1.6975132) q[1];
sx q[1];
rz(1.3664043) q[1];
rz(-pi) q[2];
rz(-0.57101698) q[3];
sx q[3];
rz(-2.5072376) q[3];
sx q[3];
rz(2.3650124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7654968) q[2];
sx q[2];
rz(-1.2906047) q[2];
sx q[2];
rz(-0.20720227) q[2];
rz(-0.97992212) q[3];
sx q[3];
rz(-1.1332952) q[3];
sx q[3];
rz(2.9492212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1205263) q[0];
sx q[0];
rz(-1.4561894) q[0];
sx q[0];
rz(-0.86984632) q[0];
rz(2.6407241) q[1];
sx q[1];
rz(-0.23575467) q[1];
sx q[1];
rz(-2.2101319) q[1];
rz(2.4293368) q[2];
sx q[2];
rz(-2.9604572) q[2];
sx q[2];
rz(-1.0618718) q[2];
rz(1.408314) q[3];
sx q[3];
rz(-1.3357031) q[3];
sx q[3];
rz(1.6817844) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
