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
rz(1.3833157) q[0];
sx q[0];
rz(-1.436469) q[0];
sx q[0];
rz(0.96487784) q[0];
rz(2.3896253) q[1];
sx q[1];
rz(-2.7120092) q[1];
sx q[1];
rz(1.0493976) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0239068) q[0];
sx q[0];
rz(-1.2313594) q[0];
sx q[0];
rz(-0.95603671) q[0];
rz(-2.9994316) q[2];
sx q[2];
rz(-1.249525) q[2];
sx q[2];
rz(-2.4427593) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.4360848) q[1];
sx q[1];
rz(-0.75172808) q[1];
sx q[1];
rz(-0.99940325) q[1];
x q[2];
rz(-1.6635632) q[3];
sx q[3];
rz(-2.1481107) q[3];
sx q[3];
rz(-2.4924459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.45399484) q[2];
sx q[2];
rz(-1.6124085) q[2];
sx q[2];
rz(3.122984) q[2];
rz(-0.54443693) q[3];
sx q[3];
rz(-0.33014044) q[3];
sx q[3];
rz(-0.64793599) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675156) q[0];
sx q[0];
rz(-1.6870966) q[0];
sx q[0];
rz(-2.6065705) q[0];
rz(0.53994838) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(1.3998869) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9104509) q[0];
sx q[0];
rz(-2.5325091) q[0];
sx q[0];
rz(-0.85394359) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.56324701) q[2];
sx q[2];
rz(-2.0958825) q[2];
sx q[2];
rz(-1.4525177) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8230086) q[1];
sx q[1];
rz(-0.46097791) q[1];
sx q[1];
rz(-2.5557842) q[1];
rz(2.4055491) q[3];
sx q[3];
rz(-0.76220817) q[3];
sx q[3];
rz(-2.8413642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.0743559) q[2];
sx q[2];
rz(-1.3398193) q[2];
sx q[2];
rz(-2.2155217) q[2];
rz(2.1615084) q[3];
sx q[3];
rz(-2.1772549) q[3];
sx q[3];
rz(0.66942352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7922908) q[0];
sx q[0];
rz(-1.2783569) q[0];
sx q[0];
rz(-1.6287623) q[0];
rz(2.8577562) q[1];
sx q[1];
rz(-2.2161039) q[1];
sx q[1];
rz(-1.1121174) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.29942) q[0];
sx q[0];
rz(-1.5636235) q[0];
sx q[0];
rz(0.0070863574) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4623423) q[2];
sx q[2];
rz(-1.6466093) q[2];
sx q[2];
rz(0.52559847) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2123472) q[1];
sx q[1];
rz(-2.622859) q[1];
sx q[1];
rz(-0.38350819) q[1];
rz(-pi) q[2];
rz(-1.9341957) q[3];
sx q[3];
rz(-1.8135241) q[3];
sx q[3];
rz(-3.0128765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5230368) q[2];
sx q[2];
rz(-1.3448389) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(-2.2526422) q[3];
sx q[3];
rz(-0.44822732) q[3];
sx q[3];
rz(-0.86047188) q[3];
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
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3717644) q[0];
sx q[0];
rz(-2.9673321) q[0];
sx q[0];
rz(2.4523822) q[0];
rz(0.31461942) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(1.4124195) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.012216598) q[0];
sx q[0];
rz(-2.3123992) q[0];
sx q[0];
rz(2.0022805) q[0];
rz(-1.6666404) q[2];
sx q[2];
rz(-1.3587225) q[2];
sx q[2];
rz(-2.1532358) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8103579) q[1];
sx q[1];
rz(-2.7572933) q[1];
sx q[1];
rz(-2.9879301) q[1];
x q[2];
rz(-1.2127498) q[3];
sx q[3];
rz(-0.51135495) q[3];
sx q[3];
rz(-0.33530385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41669258) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(-0.55595428) q[2];
rz(1.8703095) q[3];
sx q[3];
rz(-1.8794182) q[3];
sx q[3];
rz(-2.8506193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81412643) q[0];
sx q[0];
rz(-3.0004582) q[0];
sx q[0];
rz(-2.5471174) q[0];
rz(0.56198436) q[1];
sx q[1];
rz(-2.2832506) q[1];
sx q[1];
rz(-2.1655703) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41244477) q[0];
sx q[0];
rz(-1.9322104) q[0];
sx q[0];
rz(-2.0301719) q[0];
rz(-0.75404928) q[2];
sx q[2];
rz(-1.591914) q[2];
sx q[2];
rz(1.6789951) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.70006338) q[1];
sx q[1];
rz(-2.2231327) q[1];
sx q[1];
rz(-1.9487914) q[1];
rz(0.30140437) q[3];
sx q[3];
rz(-1.3632953) q[3];
sx q[3];
rz(-1.0338155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.48656616) q[2];
sx q[2];
rz(-1.2532633) q[2];
sx q[2];
rz(-2.5308934) q[2];
rz(-0.30580172) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(-1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3139451) q[0];
sx q[0];
rz(-3.1231472) q[0];
sx q[0];
rz(-1.4661283) q[0];
rz(-0.77955359) q[1];
sx q[1];
rz(-1.9773217) q[1];
sx q[1];
rz(-2.4028042) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95979106) q[0];
sx q[0];
rz(-2.7914146) q[0];
sx q[0];
rz(-2.2018593) q[0];
x q[1];
rz(-3.1081057) q[2];
sx q[2];
rz(-2.8346363) q[2];
sx q[2];
rz(-2.592776) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.14911095) q[1];
sx q[1];
rz(-0.68435366) q[1];
sx q[1];
rz(-2.1696287) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1973279) q[3];
sx q[3];
rz(-1.2324782) q[3];
sx q[3];
rz(-0.23972971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.5256727) q[2];
sx q[2];
rz(-1.2391261) q[2];
sx q[2];
rz(0.8626779) q[2];
rz(0.45977965) q[3];
sx q[3];
rz(-1.6254057) q[3];
sx q[3];
rz(-1.1427243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.2343242) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(-1.020485) q[0];
rz(0.057295784) q[1];
sx q[1];
rz(-1.658193) q[1];
sx q[1];
rz(-0.78757706) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5334398) q[0];
sx q[0];
rz(-0.69204563) q[0];
sx q[0];
rz(-2.7638021) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8767005) q[2];
sx q[2];
rz(-0.26703003) q[2];
sx q[2];
rz(1.9445813) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.5912093) q[1];
sx q[1];
rz(-1.0799284) q[1];
sx q[1];
rz(-0.20259133) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1586886) q[3];
sx q[3];
rz(-1.5421151) q[3];
sx q[3];
rz(2.8010288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.82137498) q[2];
sx q[2];
rz(-1.3221952) q[2];
sx q[2];
rz(-2.5569432) q[2];
rz(0.65230495) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(-1.339213) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(-0.24340165) q[0];
sx q[0];
rz(-1.3404055) q[0];
sx q[0];
rz(-2.8133494) q[0];
rz(0.39168656) q[1];
sx q[1];
rz(-1.1513386) q[1];
sx q[1];
rz(-1.6627056) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1680964) q[0];
sx q[0];
rz(-0.35051051) q[0];
sx q[0];
rz(-0.71061937) q[0];
rz(-pi) q[1];
x q[1];
rz(0.72490643) q[2];
sx q[2];
rz(-2.4838243) q[2];
sx q[2];
rz(-2.7810514) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7110243) q[1];
sx q[1];
rz(-1.4757753) q[1];
sx q[1];
rz(-1.8501758) q[1];
x q[2];
rz(2.0746873) q[3];
sx q[3];
rz(-0.41620884) q[3];
sx q[3];
rz(1.6817301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.0216003) q[2];
sx q[2];
rz(-1.768521) q[2];
sx q[2];
rz(2.3528986) q[2];
rz(-0.24104077) q[3];
sx q[3];
rz(-2.2153416) q[3];
sx q[3];
rz(-2.327976) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4929844) q[0];
sx q[0];
rz(-2.6032175) q[0];
sx q[0];
rz(0.96187821) q[0];
rz(-2.9073763) q[1];
sx q[1];
rz(-1.1385671) q[1];
sx q[1];
rz(-0.73807565) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2797425) q[0];
sx q[0];
rz(-1.8075917) q[0];
sx q[0];
rz(-0.40646942) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0162651) q[2];
sx q[2];
rz(-2.1574852) q[2];
sx q[2];
rz(-2.4054804) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.16633655) q[1];
sx q[1];
rz(-1.8615906) q[1];
sx q[1];
rz(1.1329805) q[1];
rz(-pi) q[2];
rz(2.4752615) q[3];
sx q[3];
rz(-1.7937346) q[3];
sx q[3];
rz(-1.2957038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.48478475) q[2];
sx q[2];
rz(-1.022555) q[2];
sx q[2];
rz(-1.8812995) q[2];
rz(1.8079181) q[3];
sx q[3];
rz(-1.0013872) q[3];
sx q[3];
rz(2.8792152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(0.86724487) q[0];
rz(-2.1278837) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(2.0713461) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4152756) q[0];
sx q[0];
rz(-1.4196102) q[0];
sx q[0];
rz(-0.12355208) q[0];
rz(-pi) q[1];
rz(0.81999166) q[2];
sx q[2];
rz(-0.89897663) q[2];
sx q[2];
rz(-1.2206248) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.2667497) q[1];
sx q[1];
rz(-1.0581746) q[1];
sx q[1];
rz(2.0155409) q[1];
rz(-pi) q[2];
rz(-0.66859122) q[3];
sx q[3];
rz(-2.9644659) q[3];
sx q[3];
rz(-1.1531545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9762207) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(1.9099859) q[2];
rz(1.5008789) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-0.73178449) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83723849) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(2.7267743) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(1.5157221) q[2];
sx q[2];
rz(-2.016042) q[2];
sx q[2];
rz(2.945937) q[2];
rz(-2.1555156) q[3];
sx q[3];
rz(-0.90860962) q[3];
sx q[3];
rz(-1.2398401) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
