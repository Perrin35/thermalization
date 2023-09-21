OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5390227) q[0];
sx q[0];
rz(-2.5780926) q[0];
sx q[0];
rz(-0.45698693) q[0];
rz(-0.62178388) q[1];
sx q[1];
rz(-0.68067247) q[1];
sx q[1];
rz(1.2759804) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5029992) q[0];
sx q[0];
rz(-2.8832957) q[0];
sx q[0];
rz(-2.5845085) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0511191) q[2];
sx q[2];
rz(-1.4960559) q[2];
sx q[2];
rz(-0.59404101) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1162122) q[1];
sx q[1];
rz(-2.930021) q[1];
sx q[1];
rz(-1.8616574) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.071322) q[3];
sx q[3];
rz(-0.71159092) q[3];
sx q[3];
rz(2.4430032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1093381) q[2];
sx q[2];
rz(-2.2685752) q[2];
sx q[2];
rz(-0.16513744) q[2];
rz(1.6312381) q[3];
sx q[3];
rz(-0.32838467) q[3];
sx q[3];
rz(0.74222773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5052658) q[0];
sx q[0];
rz(-2.9716182) q[0];
sx q[0];
rz(-0.045036137) q[0];
rz(0.33915195) q[1];
sx q[1];
rz(-1.3749326) q[1];
sx q[1];
rz(-2.3388458) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4883603) q[0];
sx q[0];
rz(-2.687504) q[0];
sx q[0];
rz(-1.0351719) q[0];
rz(-1.0674092) q[2];
sx q[2];
rz(-1.6620087) q[2];
sx q[2];
rz(-0.78546055) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2903366) q[1];
sx q[1];
rz(-2.1888121) q[1];
sx q[1];
rz(3.0783155) q[1];
rz(1.4071776) q[3];
sx q[3];
rz(-2.46582) q[3];
sx q[3];
rz(0.27392745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7654483) q[2];
sx q[2];
rz(-0.62952289) q[2];
sx q[2];
rz(-1.4260028) q[2];
rz(0.60570335) q[3];
sx q[3];
rz(-2.1798445) q[3];
sx q[3];
rz(-3.0595996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1894492) q[0];
sx q[0];
rz(-1.209963) q[0];
sx q[0];
rz(0.28513518) q[0];
rz(0.51672283) q[1];
sx q[1];
rz(-1.6500094) q[1];
sx q[1];
rz(1.8018988) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0403255) q[0];
sx q[0];
rz(-2.2650532) q[0];
sx q[0];
rz(-2.5800152) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2013024) q[2];
sx q[2];
rz(-0.72417799) q[2];
sx q[2];
rz(0.61831123) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.50946402) q[1];
sx q[1];
rz(-2.8948935) q[1];
sx q[1];
rz(0.86172608) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4524868) q[3];
sx q[3];
rz(-1.2561241) q[3];
sx q[3];
rz(-0.93022197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.67119917) q[2];
sx q[2];
rz(-2.2097887) q[2];
sx q[2];
rz(1.408067) q[2];
rz(2.9584598) q[3];
sx q[3];
rz(-1.0147084) q[3];
sx q[3];
rz(-0.071921913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.4801487) q[0];
sx q[0];
rz(-2.792795) q[0];
sx q[0];
rz(-2.9273422) q[0];
rz(2.7125773) q[1];
sx q[1];
rz(-1.1209826) q[1];
sx q[1];
rz(-0.66657153) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.727107) q[0];
sx q[0];
rz(-1.7761199) q[0];
sx q[0];
rz(-1.3798825) q[0];
rz(-pi) q[1];
rz(0.71422691) q[2];
sx q[2];
rz(-1.7923317) q[2];
sx q[2];
rz(-0.54745882) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.8909292) q[1];
sx q[1];
rz(-0.86859413) q[1];
sx q[1];
rz(-0.82171792) q[1];
rz(1.4227082) q[3];
sx q[3];
rz(-1.7184966) q[3];
sx q[3];
rz(1.5642017) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.8551222) q[2];
sx q[2];
rz(-0.21229599) q[2];
sx q[2];
rz(0.66037035) q[2];
rz(1.1874229) q[3];
sx q[3];
rz(-1.9027477) q[3];
sx q[3];
rz(2.58113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51263556) q[0];
sx q[0];
rz(-1.9218788) q[0];
sx q[0];
rz(-0.83706013) q[0];
rz(-2.191026) q[1];
sx q[1];
rz(-0.66851139) q[1];
sx q[1];
rz(1.1594835) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449074) q[0];
sx q[0];
rz(-2.3527745) q[0];
sx q[0];
rz(0.44992723) q[0];
x q[1];
rz(2.9858066) q[2];
sx q[2];
rz(-2.163379) q[2];
sx q[2];
rz(-2.8197545) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78427181) q[1];
sx q[1];
rz(-1.9958152) q[1];
sx q[1];
rz(-0.97370633) q[1];
rz(2.106835) q[3];
sx q[3];
rz(-2.1376107) q[3];
sx q[3];
rz(2.00622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.435047) q[2];
sx q[2];
rz(-1.1112301) q[2];
sx q[2];
rz(1.2716028) q[2];
rz(2.0909677) q[3];
sx q[3];
rz(-1.8836861) q[3];
sx q[3];
rz(0.064855382) q[3];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5669252) q[0];
sx q[0];
rz(-0.25952473) q[0];
sx q[0];
rz(1.1151474) q[0];
rz(-3.0976345) q[1];
sx q[1];
rz(-1.7936488) q[1];
sx q[1];
rz(0.10087068) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4911762) q[0];
sx q[0];
rz(-0.6498148) q[0];
sx q[0];
rz(-2.8894436) q[0];
rz(-2.6249044) q[2];
sx q[2];
rz(-2.5942205) q[2];
sx q[2];
rz(-1.9043497) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.069177376) q[1];
sx q[1];
rz(-2.8848856) q[1];
sx q[1];
rz(-2.2290345) q[1];
rz(0.50562596) q[3];
sx q[3];
rz(-2.6856832) q[3];
sx q[3];
rz(2.0924007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.7373401) q[2];
sx q[2];
rz(-1.6254144) q[2];
sx q[2];
rz(-0.92528382) q[2];
rz(2.960079) q[3];
sx q[3];
rz(-2.4003024) q[3];
sx q[3];
rz(-1.8484176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3826542) q[0];
sx q[0];
rz(-0.42025867) q[0];
sx q[0];
rz(2.1784901) q[0];
rz(-1.5513647) q[1];
sx q[1];
rz(-1.0179049) q[1];
sx q[1];
rz(0.68774736) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4884268) q[0];
sx q[0];
rz(-1.352172) q[0];
sx q[0];
rz(1.2290918) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8887927) q[2];
sx q[2];
rz(-1.7858292) q[2];
sx q[2];
rz(3.0766578) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.93393275) q[1];
sx q[1];
rz(-1.7001517) q[1];
sx q[1];
rz(-1.594286) q[1];
x q[2];
rz(1.0024768) q[3];
sx q[3];
rz(-1.5842614) q[3];
sx q[3];
rz(0.79893597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3380276) q[2];
sx q[2];
rz(-1.7269208) q[2];
sx q[2];
rz(-0.93758279) q[2];
rz(0.98085105) q[3];
sx q[3];
rz(-0.3347446) q[3];
sx q[3];
rz(0.78398314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5279609) q[0];
sx q[0];
rz(-2.1870446) q[0];
sx q[0];
rz(-3.0623073) q[0];
rz(1.1212564) q[1];
sx q[1];
rz(-0.93833485) q[1];
sx q[1];
rz(1.942873) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51533651) q[0];
sx q[0];
rz(-1.5021828) q[0];
sx q[0];
rz(-2.7068044) q[0];
x q[1];
rz(3.0555658) q[2];
sx q[2];
rz(-1.2089529) q[2];
sx q[2];
rz(1.3081683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3938013) q[1];
sx q[1];
rz(-2.0732905) q[1];
sx q[1];
rz(1.5138813) q[1];
rz(-pi) q[2];
rz(-1.8317322) q[3];
sx q[3];
rz(-2.0203777) q[3];
sx q[3];
rz(2.159312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.85598677) q[2];
sx q[2];
rz(-2.9957643) q[2];
sx q[2];
rz(-0.78906995) q[2];
rz(-0.35813913) q[3];
sx q[3];
rz(-1.5538244) q[3];
sx q[3];
rz(1.2808799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57551861) q[0];
sx q[0];
rz(-1.386336) q[0];
sx q[0];
rz(-2.5532706) q[0];
rz(-0.026780216) q[1];
sx q[1];
rz(-0.90590042) q[1];
sx q[1];
rz(2.2081597) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8844922) q[0];
sx q[0];
rz(-1.4643837) q[0];
sx q[0];
rz(-1.6305001) q[0];
x q[1];
rz(-0.33414267) q[2];
sx q[2];
rz(-1.2103989) q[2];
sx q[2];
rz(-2.6956218) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.583608) q[1];
sx q[1];
rz(-0.68478157) q[1];
sx q[1];
rz(3.0834553) q[1];
rz(-0.65255717) q[3];
sx q[3];
rz(-1.8592632) q[3];
sx q[3];
rz(1.0126589) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0788706) q[2];
sx q[2];
rz(-2.3214919) q[2];
sx q[2];
rz(0.36994568) q[2];
rz(2.2423819) q[3];
sx q[3];
rz(-1.4827385) q[3];
sx q[3];
rz(-0.95329154) q[3];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9545492) q[0];
sx q[0];
rz(-1.3267013) q[0];
sx q[0];
rz(-0.6533587) q[0];
rz(1.5006784) q[1];
sx q[1];
rz(-1.6710284) q[1];
sx q[1];
rz(-2.9577589) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1856954) q[0];
sx q[0];
rz(-2.9736608) q[0];
sx q[0];
rz(0.76783617) q[0];
rz(-pi) q[1];
rz(2.3887563) q[2];
sx q[2];
rz(-2.2194153) q[2];
sx q[2];
rz(-0.93968771) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.93563423) q[1];
sx q[1];
rz(-1.7365672) q[1];
sx q[1];
rz(2.2531855) q[1];
rz(-3.1179908) q[3];
sx q[3];
rz(-0.58348237) q[3];
sx q[3];
rz(2.6445146) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2227778) q[2];
sx q[2];
rz(-0.99514014) q[2];
sx q[2];
rz(0.20239057) q[2];
rz(-1.0732132) q[3];
sx q[3];
rz(-1.9066633) q[3];
sx q[3];
rz(-1.7741268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4595173) q[0];
sx q[0];
rz(-0.85953241) q[0];
sx q[0];
rz(-0.55534242) q[0];
rz(-2.7783685) q[1];
sx q[1];
rz(-1.3365311) q[1];
sx q[1];
rz(2.8844759) q[1];
rz(-2.6828962) q[2];
sx q[2];
rz(-1.1970604) q[2];
sx q[2];
rz(1.6470439) q[2];
rz(-1.8300874) q[3];
sx q[3];
rz(-2.5355831) q[3];
sx q[3];
rz(-2.2021962) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
