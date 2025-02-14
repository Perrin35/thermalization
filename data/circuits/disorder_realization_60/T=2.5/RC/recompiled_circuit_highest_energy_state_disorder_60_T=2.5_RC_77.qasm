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
rz(-2.5118339) q[0];
sx q[0];
rz(-0.30183733) q[0];
sx q[0];
rz(-1.8344185) q[0];
rz(-0.52840003) q[1];
sx q[1];
rz(-0.83146787) q[1];
sx q[1];
rz(-2.6822579) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89090251) q[0];
sx q[0];
rz(-2.8807682) q[0];
sx q[0];
rz(-1.7623259) q[0];
rz(-1.2940965) q[2];
sx q[2];
rz(-1.3802229) q[2];
sx q[2];
rz(-1.1016359) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.6428549) q[1];
sx q[1];
rz(-0.6371405) q[1];
sx q[1];
rz(2.8481507) q[1];
x q[2];
rz(-1.5993871) q[3];
sx q[3];
rz(-0.66900476) q[3];
sx q[3];
rz(2.1137848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.1186195) q[2];
sx q[2];
rz(-1.8034673) q[2];
sx q[2];
rz(2.5845134) q[2];
rz(0.45595512) q[3];
sx q[3];
rz(-0.7619226) q[3];
sx q[3];
rz(2.5383811) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2073681) q[0];
sx q[0];
rz(-2.4226483) q[0];
sx q[0];
rz(-1.3907322) q[0];
rz(1.3181744) q[1];
sx q[1];
rz(-1.7134106) q[1];
sx q[1];
rz(0.21600977) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7129546) q[0];
sx q[0];
rz(-2.1508576) q[0];
sx q[0];
rz(-0.38016502) q[0];
rz(2.7223396) q[2];
sx q[2];
rz(-1.6800095) q[2];
sx q[2];
rz(2.2855503) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.60571721) q[1];
sx q[1];
rz(-0.89369666) q[1];
sx q[1];
rz(-2.2435921) q[1];
rz(-pi) q[2];
rz(1.9253821) q[3];
sx q[3];
rz(-0.56963339) q[3];
sx q[3];
rz(2.5550752) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4231437) q[2];
sx q[2];
rz(-0.37280145) q[2];
sx q[2];
rz(0.049840363) q[2];
rz(-0.54346624) q[3];
sx q[3];
rz(-1.1246357) q[3];
sx q[3];
rz(1.0928833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8922888) q[0];
sx q[0];
rz(-1.1019305) q[0];
sx q[0];
rz(2.1732543) q[0];
rz(3.1210506) q[1];
sx q[1];
rz(-1.3803218) q[1];
sx q[1];
rz(-1.4302018) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.641159) q[0];
sx q[0];
rz(-0.64167385) q[0];
sx q[0];
rz(3.1042751) q[0];
x q[1];
rz(1.2015351) q[2];
sx q[2];
rz(-2.103269) q[2];
sx q[2];
rz(0.42064253) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.30114445) q[1];
sx q[1];
rz(-1.9663875) q[1];
sx q[1];
rz(-1.7823094) q[1];
rz(-1.537451) q[3];
sx q[3];
rz(-1.9474267) q[3];
sx q[3];
rz(-1.4550524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.7057544) q[2];
sx q[2];
rz(-2.4987554) q[2];
sx q[2];
rz(1.359681) q[2];
rz(-0.5082353) q[3];
sx q[3];
rz(-1.619092) q[3];
sx q[3];
rz(-1.9967509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.085676) q[0];
sx q[0];
rz(-0.29656947) q[0];
sx q[0];
rz(2.1348409) q[0];
rz(2.309917) q[1];
sx q[1];
rz(-2.3995903) q[1];
sx q[1];
rz(-2.0337909) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2058973) q[0];
sx q[0];
rz(-2.1855121) q[0];
sx q[0];
rz(0.45386916) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.37563373) q[2];
sx q[2];
rz(-1.3357329) q[2];
sx q[2];
rz(1.3974853) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2017758) q[1];
sx q[1];
rz(-0.30758938) q[1];
sx q[1];
rz(-1.5053476) q[1];
x q[2];
rz(-0.54385984) q[3];
sx q[3];
rz(-0.50854581) q[3];
sx q[3];
rz(-0.16112215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6882249) q[2];
sx q[2];
rz(-2.0757165) q[2];
sx q[2];
rz(2.8284956) q[2];
rz(-2.744216) q[3];
sx q[3];
rz(-1.4704967) q[3];
sx q[3];
rz(2.9562922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6322286) q[0];
sx q[0];
rz(-2.2780184) q[0];
sx q[0];
rz(-2.2160227) q[0];
rz(0.46982345) q[1];
sx q[1];
rz(-1.8269822) q[1];
sx q[1];
rz(-2.125461) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7604664) q[0];
sx q[0];
rz(-2.010049) q[0];
sx q[0];
rz(-2.683987) q[0];
rz(2.2158898) q[2];
sx q[2];
rz(-2.7113554) q[2];
sx q[2];
rz(1.0412585) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.092246902) q[1];
sx q[1];
rz(-1.7723485) q[1];
sx q[1];
rz(1.3231897) q[1];
rz(-1.6658832) q[3];
sx q[3];
rz(-0.32470266) q[3];
sx q[3];
rz(-1.269817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4887345) q[2];
sx q[2];
rz(-2.8159339) q[2];
sx q[2];
rz(0.48750901) q[2];
rz(-2.2440535) q[3];
sx q[3];
rz(-1.2582015) q[3];
sx q[3];
rz(1.0645617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94051903) q[0];
sx q[0];
rz(-0.93795332) q[0];
sx q[0];
rz(0.67365375) q[0];
rz(-1.9081839) q[1];
sx q[1];
rz(-2.1234832) q[1];
sx q[1];
rz(0.76603755) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7097283) q[0];
sx q[0];
rz(-1.0610631) q[0];
sx q[0];
rz(-0.30987375) q[0];
rz(-pi) q[1];
rz(0.38389194) q[2];
sx q[2];
rz(-2.9523179) q[2];
sx q[2];
rz(-1.6757019) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1704857) q[1];
sx q[1];
rz(-1.0909214) q[1];
sx q[1];
rz(2.9779153) q[1];
rz(-pi) q[2];
rz(-1.4333731) q[3];
sx q[3];
rz(-1.8738462) q[3];
sx q[3];
rz(-2.0772484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.172714) q[2];
sx q[2];
rz(-1.4977027) q[2];
sx q[2];
rz(-2.3806351) q[2];
rz(0.40192762) q[3];
sx q[3];
rz(-0.26508489) q[3];
sx q[3];
rz(-2.5529805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48259476) q[0];
sx q[0];
rz(-2.3880279) q[0];
sx q[0];
rz(1.6931417) q[0];
rz(-0.50085577) q[1];
sx q[1];
rz(-1.294699) q[1];
sx q[1];
rz(-0.81333152) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2091917) q[0];
sx q[0];
rz(-1.4228357) q[0];
sx q[0];
rz(-0.39011057) q[0];
rz(-1.3766889) q[2];
sx q[2];
rz(-2.5849205) q[2];
sx q[2];
rz(-2.106032) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8706559) q[1];
sx q[1];
rz(-1.7084048) q[1];
sx q[1];
rz(0.44160053) q[1];
rz(-0.058815033) q[3];
sx q[3];
rz(-1.4487189) q[3];
sx q[3];
rz(-0.59497817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8616051) q[2];
sx q[2];
rz(-2.5364752) q[2];
sx q[2];
rz(-2.81847) q[2];
rz(1.7396287) q[3];
sx q[3];
rz(-1.2819382) q[3];
sx q[3];
rz(1.3824979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-2.062227) q[0];
sx q[0];
rz(-2.0923738) q[0];
sx q[0];
rz(-2.3754689) q[0];
rz(-2.3221723) q[1];
sx q[1];
rz(-2.1111646) q[1];
sx q[1];
rz(1.5310418) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3643506) q[0];
sx q[0];
rz(-1.3421699) q[0];
sx q[0];
rz(-1.6407439) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.28240164) q[2];
sx q[2];
rz(-0.67495433) q[2];
sx q[2];
rz(0.16004983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.8360525) q[1];
sx q[1];
rz(-0.96287268) q[1];
sx q[1];
rz(-2.4120283) q[1];
rz(-pi) q[2];
rz(-0.29115486) q[3];
sx q[3];
rz(-1.2586354) q[3];
sx q[3];
rz(2.1812584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3766342) q[2];
sx q[2];
rz(-0.19516334) q[2];
sx q[2];
rz(0.39806077) q[2];
rz(-2.5352488) q[3];
sx q[3];
rz(-1.5747993) q[3];
sx q[3];
rz(2.2280367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23676087) q[0];
sx q[0];
rz(-0.30192152) q[0];
sx q[0];
rz(2.6171369) q[0];
rz(0.66871387) q[1];
sx q[1];
rz(-1.5293744) q[1];
sx q[1];
rz(1.3800157) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76261824) q[0];
sx q[0];
rz(-1.5272104) q[0];
sx q[0];
rz(2.144078) q[0];
x q[1];
rz(-1.0247714) q[2];
sx q[2];
rz(-0.54440343) q[2];
sx q[2];
rz(-2.4449206) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6898247) q[1];
sx q[1];
rz(-2.5344026) q[1];
sx q[1];
rz(-0.78801207) q[1];
rz(-2.5502167) q[3];
sx q[3];
rz(-2.5246361) q[3];
sx q[3];
rz(2.5643831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1606719) q[2];
sx q[2];
rz(-2.394684) q[2];
sx q[2];
rz(-2.6778632) q[2];
rz(1.4240228) q[3];
sx q[3];
rz(-2.0537036) q[3];
sx q[3];
rz(-3.1080642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56228191) q[0];
sx q[0];
rz(-2.4233241) q[0];
sx q[0];
rz(0.41807362) q[0];
rz(-2.383291) q[1];
sx q[1];
rz(-2.1577991) q[1];
sx q[1];
rz(-1.8565348) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67659527) q[0];
sx q[0];
rz(-0.9228068) q[0];
sx q[0];
rz(1.4761488) q[0];
rz(-pi) q[1];
rz(-2.7964716) q[2];
sx q[2];
rz(-2.1639544) q[2];
sx q[2];
rz(2.5760004) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9662719) q[1];
sx q[1];
rz(-2.113225) q[1];
sx q[1];
rz(-1.346051) q[1];
rz(-pi) q[2];
x q[2];
rz(0.22104903) q[3];
sx q[3];
rz(-1.2411102) q[3];
sx q[3];
rz(-0.21176618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.9717504) q[2];
sx q[2];
rz(-1.0047793) q[2];
sx q[2];
rz(2.8954835) q[2];
rz(-1.2717815) q[3];
sx q[3];
rz(-1.5747986) q[3];
sx q[3];
rz(-1.3625712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
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
rz(-1.9825738) q[0];
sx q[0];
rz(-1.4763426) q[0];
sx q[0];
rz(2.939298) q[0];
rz(-2.5166439) q[1];
sx q[1];
rz(-0.85660558) q[1];
sx q[1];
rz(0.4013335) q[1];
rz(-1.779378) q[2];
sx q[2];
rz(-1.7558395) q[2];
sx q[2];
rz(0.21370733) q[2];
rz(1.2558595) q[3];
sx q[3];
rz(-1.3797912) q[3];
sx q[3];
rz(2.2906863) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
