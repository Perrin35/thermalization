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
rz(0.87822479) q[0];
sx q[0];
rz(3.8642519) q[0];
sx q[0];
rz(10.470358) q[0];
rz(0.0027520952) q[1];
sx q[1];
rz(-1.4581008) q[1];
sx q[1];
rz(2.4207065) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84522879) q[0];
sx q[0];
rz(-2.2382444) q[0];
sx q[0];
rz(3.0240701) q[0];
rz(-pi) q[1];
rz(-0.051250967) q[2];
sx q[2];
rz(-1.7619216) q[2];
sx q[2];
rz(-0.65039907) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8869141) q[1];
sx q[1];
rz(-0.47770914) q[1];
sx q[1];
rz(2.5049247) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.74303307) q[3];
sx q[3];
rz(-1.475302) q[3];
sx q[3];
rz(-2.9955289) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2846994) q[2];
sx q[2];
rz(-1.928494) q[2];
sx q[2];
rz(-2.6846679) q[2];
rz(-0.65447718) q[3];
sx q[3];
rz(-0.5623397) q[3];
sx q[3];
rz(0.42403179) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65861312) q[0];
sx q[0];
rz(-2.787866) q[0];
sx q[0];
rz(2.5335627) q[0];
rz(-1.9425302) q[1];
sx q[1];
rz(-0.46747318) q[1];
sx q[1];
rz(0.79493585) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1300995) q[0];
sx q[0];
rz(-1.8602295) q[0];
sx q[0];
rz(-0.35057803) q[0];
rz(1.5995922) q[2];
sx q[2];
rz(-1.303728) q[2];
sx q[2];
rz(2.0350589) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7424189) q[1];
sx q[1];
rz(-1.8493506) q[1];
sx q[1];
rz(-0.919085) q[1];
x q[2];
rz(2.7057418) q[3];
sx q[3];
rz(-2.4727949) q[3];
sx q[3];
rz(-2.8613402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.46637154) q[2];
sx q[2];
rz(-1.0543062) q[2];
sx q[2];
rz(1.4804473) q[2];
rz(-0.36888567) q[3];
sx q[3];
rz(-1.2494913) q[3];
sx q[3];
rz(-0.7307542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22496255) q[0];
sx q[0];
rz(-2.1376305) q[0];
sx q[0];
rz(-0.68600255) q[0];
rz(-1.8983768) q[1];
sx q[1];
rz(-2.9153283) q[1];
sx q[1];
rz(-2.539198) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43322726) q[0];
sx q[0];
rz(-2.4658563) q[0];
sx q[0];
rz(2.2883718) q[0];
x q[1];
rz(-3.1092293) q[2];
sx q[2];
rz(-1.6113069) q[2];
sx q[2];
rz(1.4537741) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.408566) q[1];
sx q[1];
rz(-1.0045718) q[1];
sx q[1];
rz(3.0263961) q[1];
rz(-pi) q[2];
rz(2.0394348) q[3];
sx q[3];
rz(-1.8640362) q[3];
sx q[3];
rz(1.6654429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38597044) q[2];
sx q[2];
rz(-1.8381511) q[2];
sx q[2];
rz(1.2516359) q[2];
rz(-2.5899467) q[3];
sx q[3];
rz(-2.9988204) q[3];
sx q[3];
rz(0.84757203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(2.8400693) q[0];
sx q[0];
rz(-1.9020377) q[0];
sx q[0];
rz(0.10313343) q[0];
rz(2.3221807) q[1];
sx q[1];
rz(-2.2401919) q[1];
sx q[1];
rz(2.588396) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2870194) q[0];
sx q[0];
rz(-2.6691045) q[0];
sx q[0];
rz(2.2437139) q[0];
rz(-pi) q[1];
rz(-0.35049866) q[2];
sx q[2];
rz(-1.0523049) q[2];
sx q[2];
rz(-1.1120591) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.64606372) q[1];
sx q[1];
rz(-0.49139443) q[1];
sx q[1];
rz(1.2786521) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0023213) q[3];
sx q[3];
rz(-1.6628886) q[3];
sx q[3];
rz(-1.0669277) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4210522) q[2];
sx q[2];
rz(-0.74137509) q[2];
sx q[2];
rz(0.53483024) q[2];
rz(2.363291) q[3];
sx q[3];
rz(-2.4141267) q[3];
sx q[3];
rz(-2.9261869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68742037) q[0];
sx q[0];
rz(-1.146831) q[0];
sx q[0];
rz(-0.73690328) q[0];
rz(-3.0829644) q[1];
sx q[1];
rz(-2.3643989) q[1];
sx q[1];
rz(2.5545205) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5694052) q[0];
sx q[0];
rz(-2.5068034) q[0];
sx q[0];
rz(1.6247747) q[0];
rz(-pi) q[1];
rz(-1.3413853) q[2];
sx q[2];
rz(-1.5093813) q[2];
sx q[2];
rz(-3.1385076) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.2559511) q[1];
sx q[1];
rz(-1.5382861) q[1];
sx q[1];
rz(1.2684275) q[1];
rz(0.80758269) q[3];
sx q[3];
rz(-2.7135757) q[3];
sx q[3];
rz(1.310809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.29207841) q[2];
sx q[2];
rz(-2.2167315) q[2];
sx q[2];
rz(-2.4936131) q[2];
rz(0.61602151) q[3];
sx q[3];
rz(-1.5024622) q[3];
sx q[3];
rz(-3.093921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.7306526) q[0];
sx q[0];
rz(-2.3766282) q[0];
sx q[0];
rz(2.1867645) q[0];
rz(3.094589) q[1];
sx q[1];
rz(-2.4689597) q[1];
sx q[1];
rz(1.0161374) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8076626) q[0];
sx q[0];
rz(-0.74650812) q[0];
sx q[0];
rz(0.73005116) q[0];
rz(2.5437433) q[2];
sx q[2];
rz(-0.81839222) q[2];
sx q[2];
rz(-2.1482836) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.816427) q[1];
sx q[1];
rz(-1.0652958) q[1];
sx q[1];
rz(-1.5542404) q[1];
rz(-pi) q[2];
rz(1.5136857) q[3];
sx q[3];
rz(-2.1942003) q[3];
sx q[3];
rz(-0.36950612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.59812349) q[2];
sx q[2];
rz(-1.2657961) q[2];
sx q[2];
rz(-0.050696105) q[2];
rz(3.0905881) q[3];
sx q[3];
rz(-1.3011322) q[3];
sx q[3];
rz(0.72371975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-2.8766668) q[0];
sx q[0];
rz(-0.51487881) q[0];
sx q[0];
rz(1.4798973) q[0];
rz(-2.0199846) q[1];
sx q[1];
rz(-1.2622967) q[1];
sx q[1];
rz(-2.6720537) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8136217) q[0];
sx q[0];
rz(-1.6780954) q[0];
sx q[0];
rz(0.046949887) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2079427) q[2];
sx q[2];
rz(-0.41703014) q[2];
sx q[2];
rz(2.1977958) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4457629) q[1];
sx q[1];
rz(-1.1053021) q[1];
sx q[1];
rz(0.19449921) q[1];
rz(-2.2121053) q[3];
sx q[3];
rz(-0.90228122) q[3];
sx q[3];
rz(-2.6098982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40662128) q[2];
sx q[2];
rz(-1.5903641) q[2];
sx q[2];
rz(2.0242019) q[2];
rz(-2.0495074) q[3];
sx q[3];
rz(-1.4039682) q[3];
sx q[3];
rz(-2.7058069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0206873) q[0];
sx q[0];
rz(-2.0285323) q[0];
sx q[0];
rz(-0.06509617) q[0];
rz(-0.33196017) q[1];
sx q[1];
rz(-1.3061701) q[1];
sx q[1];
rz(-1.8665705) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37133163) q[0];
sx q[0];
rz(-1.341371) q[0];
sx q[0];
rz(3.1317932) q[0];
x q[1];
rz(-2.0135856) q[2];
sx q[2];
rz(-1.4472359) q[2];
sx q[2];
rz(-2.3905001) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4091585) q[1];
sx q[1];
rz(-2.2526764) q[1];
sx q[1];
rz(-1.5049008) q[1];
x q[2];
rz(-0.17586768) q[3];
sx q[3];
rz(-2.3308836) q[3];
sx q[3];
rz(0.46739235) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.28978213) q[2];
sx q[2];
rz(-1.449457) q[2];
sx q[2];
rz(-2.1688478) q[2];
rz(-2.4408477) q[3];
sx q[3];
rz(-0.86207977) q[3];
sx q[3];
rz(0.8249445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7026611) q[0];
sx q[0];
rz(-0.40533608) q[0];
sx q[0];
rz(-2.7324556) q[0];
rz(1.9955955) q[1];
sx q[1];
rz(-1.0640249) q[1];
sx q[1];
rz(1.8786059) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1564558) q[0];
sx q[0];
rz(-0.0081774769) q[0];
sx q[0];
rz(-1.5550434) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1767597) q[2];
sx q[2];
rz(-1.3140643) q[2];
sx q[2];
rz(-2.9213411) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.9085051) q[1];
sx q[1];
rz(-2.246664) q[1];
sx q[1];
rz(0.29260537) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0060482) q[3];
sx q[3];
rz(-0.57430797) q[3];
sx q[3];
rz(1.3804466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.24141773) q[2];
sx q[2];
rz(-1.2215542) q[2];
sx q[2];
rz(-0.22132604) q[2];
rz(2.6545702) q[3];
sx q[3];
rz(-2.2723787) q[3];
sx q[3];
rz(-1.937449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0814334) q[0];
sx q[0];
rz(-2.8964323) q[0];
sx q[0];
rz(1.2231476) q[0];
rz(-0.080991117) q[1];
sx q[1];
rz(-1.5206189) q[1];
sx q[1];
rz(-1.5222668) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8232097) q[0];
sx q[0];
rz(-1.902474) q[0];
sx q[0];
rz(-1.7451502) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1606101) q[2];
sx q[2];
rz(-2.207617) q[2];
sx q[2];
rz(0.75904075) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4829068) q[1];
sx q[1];
rz(-1.7665837) q[1];
sx q[1];
rz(1.8808603) q[1];
rz(-1.7459938) q[3];
sx q[3];
rz(-2.0591551) q[3];
sx q[3];
rz(1.7217896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.12485518) q[2];
sx q[2];
rz(-1.7502461) q[2];
sx q[2];
rz(2.4511231) q[2];
rz(-2.8705583) q[3];
sx q[3];
rz(-0.46509585) q[3];
sx q[3];
rz(-1.5439699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4616213) q[0];
sx q[0];
rz(-2.2053056) q[0];
sx q[0];
rz(2.9619138) q[0];
rz(0.24196729) q[1];
sx q[1];
rz(-2.2503743) q[1];
sx q[1];
rz(1.888884) q[1];
rz(-2.719328) q[2];
sx q[2];
rz(-1.761407) q[2];
sx q[2];
rz(1.882886) q[2];
rz(1.9492502) q[3];
sx q[3];
rz(-0.71810953) q[3];
sx q[3];
rz(0.22267846) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
