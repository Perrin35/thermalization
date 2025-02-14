OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6310298) q[0];
sx q[0];
rz(3.6462311) q[0];
sx q[0];
rz(12.991821) q[0];
rz(-2.5692441) q[1];
sx q[1];
rz(-1.0971789) q[1];
sx q[1];
rz(1.9414577) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5050738) q[0];
sx q[0];
rz(-2.8905792) q[0];
sx q[0];
rz(-3.1025629) q[0];
rz(-pi) q[1];
rz(-0.27165551) q[2];
sx q[2];
rz(-1.2363529) q[2];
sx q[2];
rz(1.3525427) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8737826) q[1];
sx q[1];
rz(-1.7846196) q[1];
sx q[1];
rz(-2.7604483) q[1];
rz(-pi) q[2];
rz(2.4704504) q[3];
sx q[3];
rz(-1.6832388) q[3];
sx q[3];
rz(-0.46026106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9878865) q[2];
sx q[2];
rz(-1.989863) q[2];
sx q[2];
rz(-0.64725867) q[2];
rz(-2.9636532) q[3];
sx q[3];
rz(-1.8837594) q[3];
sx q[3];
rz(1.9378763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7385638) q[0];
sx q[0];
rz(-0.084349923) q[0];
sx q[0];
rz(2.6179598) q[0];
rz(1.2298443) q[1];
sx q[1];
rz(-2.3975027) q[1];
sx q[1];
rz(-2.8935166) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46892525) q[0];
sx q[0];
rz(-1.4644196) q[0];
sx q[0];
rz(-0.74161462) q[0];
rz(1.1311943) q[2];
sx q[2];
rz(-1.1710376) q[2];
sx q[2];
rz(0.93610379) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4011127) q[1];
sx q[1];
rz(-1.7164299) q[1];
sx q[1];
rz(-3.1023953) q[1];
rz(-pi) q[2];
rz(-2.4753597) q[3];
sx q[3];
rz(-2.0012451) q[3];
sx q[3];
rz(-0.43678771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5292042) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(-2.6017792) q[2];
rz(0.16048935) q[3];
sx q[3];
rz(-1.548998) q[3];
sx q[3];
rz(0.81015712) q[3];
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
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983343) q[0];
sx q[0];
rz(-0.67688268) q[0];
sx q[0];
rz(0.13370378) q[0];
rz(-1.5191822) q[1];
sx q[1];
rz(-1.2956053) q[1];
sx q[1];
rz(-2.349283) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83082047) q[0];
sx q[0];
rz(-1.7541274) q[0];
sx q[0];
rz(-1.2123327) q[0];
x q[1];
rz(-2.9970005) q[2];
sx q[2];
rz(-1.4404313) q[2];
sx q[2];
rz(-3.0693288) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.6349259) q[1];
sx q[1];
rz(-0.95639765) q[1];
sx q[1];
rz(0.8393112) q[1];
rz(-pi) q[2];
rz(-0.00040690502) q[3];
sx q[3];
rz(-0.30012977) q[3];
sx q[3];
rz(2.2426734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.85492674) q[2];
sx q[2];
rz(-0.67015219) q[2];
sx q[2];
rz(0.83038846) q[2];
rz(1.3487799) q[3];
sx q[3];
rz(-2.1068137) q[3];
sx q[3];
rz(-2.5428298) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778359) q[0];
sx q[0];
rz(-2.1649375) q[0];
sx q[0];
rz(0.42309716) q[0];
rz(-1.1571723) q[1];
sx q[1];
rz(-1.6559699) q[1];
sx q[1];
rz(0.62209904) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3836648) q[0];
sx q[0];
rz(-2.1970941) q[0];
sx q[0];
rz(-0.48547283) q[0];
x q[1];
rz(-2.9329002) q[2];
sx q[2];
rz(-1.1596173) q[2];
sx q[2];
rz(1.5087939) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8129261) q[1];
sx q[1];
rz(-1.5541967) q[1];
sx q[1];
rz(1.9783201) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1142523) q[3];
sx q[3];
rz(-0.84544824) q[3];
sx q[3];
rz(0.10408653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.28079924) q[2];
sx q[2];
rz(-1.4192899) q[2];
sx q[2];
rz(0.61817509) q[2];
rz(1.482796) q[3];
sx q[3];
rz(-1.3525454) q[3];
sx q[3];
rz(1.5084722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9789155) q[0];
sx q[0];
rz(-1.3293043) q[0];
sx q[0];
rz(-0.83258122) q[0];
rz(-0.32514462) q[1];
sx q[1];
rz(-0.99601662) q[1];
sx q[1];
rz(0.064373374) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58080855) q[0];
sx q[0];
rz(-0.34750313) q[0];
sx q[0];
rz(-2.7827713) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.60466296) q[2];
sx q[2];
rz(-0.69398601) q[2];
sx q[2];
rz(-0.99550216) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.0400914) q[1];
sx q[1];
rz(-1.027193) q[1];
sx q[1];
rz(-2.5341847) q[1];
rz(2.1983002) q[3];
sx q[3];
rz(-1.8648913) q[3];
sx q[3];
rz(-2.7315745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9084106) q[2];
sx q[2];
rz(-0.57047129) q[2];
sx q[2];
rz(-1.849966) q[2];
rz(0.49606797) q[3];
sx q[3];
rz(-0.78947624) q[3];
sx q[3];
rz(1.1424278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18601501) q[0];
sx q[0];
rz(-0.29009524) q[0];
sx q[0];
rz(-2.7225323) q[0];
rz(-0.31173197) q[1];
sx q[1];
rz(-1.5985951) q[1];
sx q[1];
rz(2.9980803) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9547894) q[0];
sx q[0];
rz(-1.1214897) q[0];
sx q[0];
rz(1.9687998) q[0];
rz(0.70563282) q[2];
sx q[2];
rz(-2.7876283) q[2];
sx q[2];
rz(1.7662545) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.4048481) q[1];
sx q[1];
rz(-0.84538922) q[1];
sx q[1];
rz(2.0015697) q[1];
rz(0.26007248) q[3];
sx q[3];
rz(-1.7100157) q[3];
sx q[3];
rz(1.8782488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.36025563) q[2];
sx q[2];
rz(-0.90136734) q[2];
sx q[2];
rz(-2.0085013) q[2];
rz(-1.1059149) q[3];
sx q[3];
rz(-1.4191041) q[3];
sx q[3];
rz(-0.47479409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.356242) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(1.6957977) q[0];
rz(0.71714199) q[1];
sx q[1];
rz(-1.278806) q[1];
sx q[1];
rz(0.76146567) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26869795) q[0];
sx q[0];
rz(-0.94572645) q[0];
sx q[0];
rz(1.3891267) q[0];
rz(-pi) q[1];
rz(2.6638159) q[2];
sx q[2];
rz(-0.91467664) q[2];
sx q[2];
rz(-2.853924) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.0822507) q[1];
sx q[1];
rz(-2.5656156) q[1];
sx q[1];
rz(1.4800998) q[1];
x q[2];
rz(-2.981212) q[3];
sx q[3];
rz(-1.7564991) q[3];
sx q[3];
rz(-1.6689614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.36934272) q[2];
sx q[2];
rz(-2.6476634) q[2];
sx q[2];
rz(0.93144766) q[2];
rz(0.61819589) q[3];
sx q[3];
rz(-1.6341011) q[3];
sx q[3];
rz(0.96562323) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6922927) q[0];
sx q[0];
rz(-1.7789142) q[0];
sx q[0];
rz(-1.3344673) q[0];
rz(-0.5131228) q[1];
sx q[1];
rz(-0.98315364) q[1];
sx q[1];
rz(1.2293053) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9216825) q[0];
sx q[0];
rz(-1.1510885) q[0];
sx q[0];
rz(0.76292636) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9077577) q[2];
sx q[2];
rz(-0.2744199) q[2];
sx q[2];
rz(-1.929259) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8905764) q[1];
sx q[1];
rz(-0.15726798) q[1];
sx q[1];
rz(3.0096439) q[1];
rz(-pi) q[2];
rz(2.4793998) q[3];
sx q[3];
rz(-1.5569038) q[3];
sx q[3];
rz(-2.3628949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8339771) q[2];
sx q[2];
rz(-0.5270842) q[2];
sx q[2];
rz(0.21044593) q[2];
rz(-3.0757507) q[3];
sx q[3];
rz(-2.6688771) q[3];
sx q[3];
rz(-2.3426447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54302067) q[0];
sx q[0];
rz(-2.4871171) q[0];
sx q[0];
rz(-1.2731113) q[0];
rz(-1.8674564) q[1];
sx q[1];
rz(-2.4576063) q[1];
sx q[1];
rz(2.2575016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9724204) q[0];
sx q[0];
rz(-2.5247249) q[0];
sx q[0];
rz(1.5703809) q[0];
rz(1.8294677) q[2];
sx q[2];
rz(-1.0610559) q[2];
sx q[2];
rz(0.81196751) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.99779656) q[1];
sx q[1];
rz(-1.9854768) q[1];
sx q[1];
rz(2.8383377) q[1];
rz(-2.8330634) q[3];
sx q[3];
rz(-1.3236227) q[3];
sx q[3];
rz(-0.17251523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.990443) q[2];
sx q[2];
rz(-2.1572025) q[2];
sx q[2];
rz(1.3746064) q[2];
rz(0.19045842) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
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
rz(2.9780438) q[0];
sx q[0];
rz(-1.4530285) q[0];
sx q[0];
rz(1.8769886) q[0];
rz(3.0679852) q[1];
sx q[1];
rz(-1.0818447) q[1];
sx q[1];
rz(-2.5216865) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6524175) q[0];
sx q[0];
rz(-0.76120725) q[0];
sx q[0];
rz(-1.8106145) q[0];
x q[1];
rz(2.5307199) q[2];
sx q[2];
rz(-1.1430642) q[2];
sx q[2];
rz(1.7868702) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.1312841) q[1];
sx q[1];
rz(-2.8493053) q[1];
sx q[1];
rz(0.75914219) q[1];
x q[2];
rz(-1.8100753) q[3];
sx q[3];
rz(-2.0258697) q[3];
sx q[3];
rz(-2.6309155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8668883) q[2];
sx q[2];
rz(-0.70211774) q[2];
sx q[2];
rz(0.14599027) q[2];
rz(2.919096) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(-0.72993025) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9324026) q[0];
sx q[0];
rz(-1.9575735) q[0];
sx q[0];
rz(0.14458543) q[0];
rz(1.2296386) q[1];
sx q[1];
rz(-1.790779) q[1];
sx q[1];
rz(1.889224) q[1];
rz(0.59496224) q[2];
sx q[2];
rz(-1.7528201) q[2];
sx q[2];
rz(2.5217944) q[2];
rz(-1.4208117) q[3];
sx q[3];
rz(-0.21258988) q[3];
sx q[3];
rz(-1.7439738) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
