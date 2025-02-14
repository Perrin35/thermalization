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
rz(-2.1767148) q[0];
rz(-0.75196737) q[1];
sx q[1];
rz(-0.42958346) q[1];
sx q[1];
rz(2.092195) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8937738) q[0];
sx q[0];
rz(-2.4501094) q[0];
sx q[0];
rz(-1.0214424) q[0];
rz(-pi) q[1];
rz(-0.14216106) q[2];
sx q[2];
rz(-1.249525) q[2];
sx q[2];
rz(-0.69883332) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8559554) q[1];
sx q[1];
rz(-2.1826943) q[1];
sx q[1];
rz(-2.6735071) q[1];
rz(-pi) q[2];
rz(2.562306) q[3];
sx q[3];
rz(-1.6484954) q[3];
sx q[3];
rz(-2.169211) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.45399484) q[2];
sx q[2];
rz(-1.5291841) q[2];
sx q[2];
rz(-0.018608658) q[2];
rz(-2.5971557) q[3];
sx q[3];
rz(-2.8114522) q[3];
sx q[3];
rz(-0.64793599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3675156) q[0];
sx q[0];
rz(-1.4544961) q[0];
sx q[0];
rz(2.6065705) q[0];
rz(0.53994838) q[1];
sx q[1];
rz(-2.5060563) q[1];
sx q[1];
rz(1.3998869) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0947123) q[0];
sx q[0];
rz(-2.0167354) q[0];
sx q[0];
rz(0.42973862) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97008791) q[2];
sx q[2];
rz(-1.0905438) q[2];
sx q[2];
rz(2.9532972) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.68097673) q[1];
sx q[1];
rz(-1.1910805) q[1];
sx q[1];
rz(1.3028076) q[1];
rz(-pi) q[2];
rz(-1.0008282) q[3];
sx q[3];
rz(-2.1080351) q[3];
sx q[3];
rz(1.9443823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0743559) q[2];
sx q[2];
rz(-1.8017733) q[2];
sx q[2];
rz(-2.2155217) q[2];
rz(0.98008424) q[3];
sx q[3];
rz(-0.96433774) q[3];
sx q[3];
rz(0.66942352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7922908) q[0];
sx q[0];
rz(-1.8632357) q[0];
sx q[0];
rz(-1.5128304) q[0];
rz(-0.28383645) q[1];
sx q[1];
rz(-0.92548871) q[1];
sx q[1];
rz(-2.0294752) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6215206) q[0];
sx q[0];
rz(-3.1315098) q[0];
sx q[0];
rz(-0.79147379) q[0];
rz(-pi) q[1];
rz(0.076259344) q[2];
sx q[2];
rz(-1.6789376) q[2];
sx q[2];
rz(2.0881483) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.92924547) q[1];
sx q[1];
rz(-0.51873365) q[1];
sx q[1];
rz(-2.7580845) q[1];
x q[2];
rz(2.882631) q[3];
sx q[3];
rz(-1.2185214) q[3];
sx q[3];
rz(1.6083628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5230368) q[2];
sx q[2];
rz(-1.7967537) q[2];
sx q[2];
rz(1.1019361) q[2];
rz(0.88895041) q[3];
sx q[3];
rz(-0.44822732) q[3];
sx q[3];
rz(-0.86047188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76982826) q[0];
sx q[0];
rz(-0.17426057) q[0];
sx q[0];
rz(-2.4523822) q[0];
rz(2.8269732) q[1];
sx q[1];
rz(-0.92051053) q[1];
sx q[1];
rz(-1.4124195) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2570626) q[0];
sx q[0];
rz(-1.2572968) q[0];
sx q[0];
rz(-2.351981) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21302235) q[2];
sx q[2];
rz(-1.4771059) q[2];
sx q[2];
rz(2.5793864) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7594193) q[1];
sx q[1];
rz(-1.628211) q[1];
sx q[1];
rz(-2.7613954) q[1];
rz(1.0869157) q[3];
sx q[3];
rz(-1.3984507) q[3];
sx q[3];
rz(2.2215171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7249001) q[2];
sx q[2];
rz(-0.56513864) q[2];
sx q[2];
rz(2.5856384) q[2];
rz(1.2712831) q[3];
sx q[3];
rz(-1.8794182) q[3];
sx q[3];
rz(2.8506193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(2.3274662) q[0];
sx q[0];
rz(-0.1411345) q[0];
sx q[0];
rz(2.5471174) q[0];
rz(-2.5796083) q[1];
sx q[1];
rz(-0.85834208) q[1];
sx q[1];
rz(-0.97602239) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41244477) q[0];
sx q[0];
rz(-1.9322104) q[0];
sx q[0];
rz(1.1114208) q[0];
rz(-0.75404928) q[2];
sx q[2];
rz(-1.5496786) q[2];
sx q[2];
rz(1.4625975) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4415293) q[1];
sx q[1];
rz(-2.2231327) q[1];
sx q[1];
rz(-1.9487914) q[1];
x q[2];
rz(-2.8401883) q[3];
sx q[3];
rz(-1.7782974) q[3];
sx q[3];
rz(1.0338155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6550265) q[2];
sx q[2];
rz(-1.2532633) q[2];
sx q[2];
rz(-0.61069926) q[2];
rz(-2.8357909) q[3];
sx q[3];
rz(-0.96547258) q[3];
sx q[3];
rz(1.3293728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82764757) q[0];
sx q[0];
rz(-0.018445404) q[0];
sx q[0];
rz(-1.4661283) q[0];
rz(2.3620391) q[1];
sx q[1];
rz(-1.164271) q[1];
sx q[1];
rz(2.4028042) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95979106) q[0];
sx q[0];
rz(-2.7914146) q[0];
sx q[0];
rz(2.2018593) q[0];
x q[1];
rz(-1.5814085) q[2];
sx q[2];
rz(-1.877575) q[2];
sx q[2];
rz(-2.5576484) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.87118426) q[1];
sx q[1];
rz(-2.1201907) q[1];
sx q[1];
rz(0.43105468) q[1];
rz(-pi) q[2];
rz(-1.9442648) q[3];
sx q[3];
rz(-1.9091144) q[3];
sx q[3];
rz(0.23972971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5256727) q[2];
sx q[2];
rz(-1.2391261) q[2];
sx q[2];
rz(0.8626779) q[2];
rz(-0.45977965) q[3];
sx q[3];
rz(-1.516187) q[3];
sx q[3];
rz(-1.1427243) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2343242) q[0];
sx q[0];
rz(-0.37721226) q[0];
sx q[0];
rz(2.1211076) q[0];
rz(0.057295784) q[1];
sx q[1];
rz(-1.4833996) q[1];
sx q[1];
rz(0.78757706) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.60815281) q[0];
sx q[0];
rz(-2.449547) q[0];
sx q[0];
rz(0.37779053) q[0];
rz(0.25813132) q[2];
sx q[2];
rz(-1.5016593) q[2];
sx q[2];
rz(0.11786945) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0246525) q[1];
sx q[1];
rz(-1.3924011) q[1];
sx q[1];
rz(-1.0712888) q[1];
rz(-pi) q[2];
rz(3.1071289) q[3];
sx q[3];
rz(-2.1584145) q[3];
sx q[3];
rz(-1.8922488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.82137498) q[2];
sx q[2];
rz(-1.8193974) q[2];
sx q[2];
rz(0.58464948) q[2];
rz(-2.4892877) q[3];
sx q[3];
rz(-1.2290596) q[3];
sx q[3];
rz(1.8023796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898191) q[0];
sx q[0];
rz(-1.8011872) q[0];
sx q[0];
rz(2.8133494) q[0];
rz(0.39168656) q[1];
sx q[1];
rz(-1.990254) q[1];
sx q[1];
rz(1.6627056) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97349629) q[0];
sx q[0];
rz(-2.7910821) q[0];
sx q[0];
rz(0.71061937) q[0];
x q[1];
rz(-2.0441891) q[2];
sx q[2];
rz(-2.0461296) q[2];
sx q[2];
rz(1.2021827) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.7110243) q[1];
sx q[1];
rz(-1.4757753) q[1];
sx q[1];
rz(-1.2914168) q[1];
rz(-pi) q[2];
rz(-1.0669054) q[3];
sx q[3];
rz(-0.41620884) q[3];
sx q[3];
rz(1.6817301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0216003) q[2];
sx q[2];
rz(-1.3730717) q[2];
sx q[2];
rz(-0.78869406) q[2];
rz(2.9005519) q[3];
sx q[3];
rz(-0.92625109) q[3];
sx q[3];
rz(-0.81361667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.4929844) q[0];
sx q[0];
rz(-0.5383752) q[0];
sx q[0];
rz(-0.96187821) q[0];
rz(-2.9073763) q[1];
sx q[1];
rz(-2.0030256) q[1];
sx q[1];
rz(0.73807565) q[1];
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
rz(-0.57508399) q[2];
sx q[2];
rz(-0.7204537) q[2];
sx q[2];
rz(-1.4478113) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.35154147) q[3];
sx q[3];
rz(-0.69720399) q[3];
sx q[3];
rz(2.5923924) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.48478475) q[2];
sx q[2];
rz(-2.1190376) q[2];
sx q[2];
rz(1.8812995) q[2];
rz(-1.8079181) q[3];
sx q[3];
rz(-2.1402054) q[3];
sx q[3];
rz(-0.26237747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8986847) q[0];
sx q[0];
rz(-2.3299291) q[0];
sx q[0];
rz(0.86724487) q[0];
rz(1.0137089) q[1];
sx q[1];
rz(-1.8127245) q[1];
sx q[1];
rz(-1.0702466) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8673422) q[0];
sx q[0];
rz(-0.19495067) q[0];
sx q[0];
rz(-0.89063962) q[0];
rz(-pi) q[1];
rz(2.321601) q[2];
sx q[2];
rz(-0.89897663) q[2];
sx q[2];
rz(1.2206248) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.874843) q[1];
sx q[1];
rz(-2.0834181) q[1];
sx q[1];
rz(2.0155409) q[1];
rz(-2.4730014) q[3];
sx q[3];
rz(-0.1771268) q[3];
sx q[3];
rz(-1.1531545) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9762207) q[2];
sx q[2];
rz(-2.3207211) q[2];
sx q[2];
rz(1.2316068) q[2];
rz(1.5008789) q[3];
sx q[3];
rz(-2.9314163) q[3];
sx q[3];
rz(-0.73178449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83723849) q[0];
sx q[0];
rz(-1.9914892) q[0];
sx q[0];
rz(-1.7627841) q[0];
rz(-0.41481836) q[1];
sx q[1];
rz(-1.3777614) q[1];
sx q[1];
rz(-1.6090964) q[1];
rz(3.0267486) q[2];
sx q[2];
rz(-0.44841246) q[2];
sx q[2];
rz(2.8186225) q[2];
rz(0.75178643) q[3];
sx q[3];
rz(-1.1204168) q[3];
sx q[3];
rz(0.71747019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
