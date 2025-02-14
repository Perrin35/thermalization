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
rz(2.1386327) q[0];
sx q[0];
rz(-0.47337368) q[0];
sx q[0];
rz(-1.0869429) q[0];
rz(-0.76397693) q[1];
sx q[1];
rz(2.6978701) q[1];
sx q[1];
rz(13.397747) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2646541) q[0];
sx q[0];
rz(-1.042141) q[0];
sx q[0];
rz(2.1318046) q[0];
x q[1];
rz(-1.126312) q[2];
sx q[2];
rz(-2.630027) q[2];
sx q[2];
rz(-1.7351406) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8352141) q[1];
sx q[1];
rz(-1.6740546) q[1];
sx q[1];
rz(1.4500965) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0696297) q[3];
sx q[3];
rz(-1.5441893) q[3];
sx q[3];
rz(-0.047544971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.6143371) q[2];
sx q[2];
rz(-1.8680806) q[2];
sx q[2];
rz(0.41198507) q[2];
rz(-0.24474239) q[3];
sx q[3];
rz(-0.63260806) q[3];
sx q[3];
rz(-2.8542724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7016542) q[0];
sx q[0];
rz(-0.97173062) q[0];
sx q[0];
rz(2.7213851) q[0];
rz(-0.48928753) q[1];
sx q[1];
rz(-0.91541618) q[1];
sx q[1];
rz(1.9583826) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1294043) q[0];
sx q[0];
rz(-0.79372915) q[0];
sx q[0];
rz(0.94214006) q[0];
x q[1];
rz(-0.80090307) q[2];
sx q[2];
rz(-0.43673968) q[2];
sx q[2];
rz(2.2855121) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.015731363) q[1];
sx q[1];
rz(-2.1441048) q[1];
sx q[1];
rz(-1.1124951) q[1];
rz(-pi) q[2];
rz(3.0239931) q[3];
sx q[3];
rz(-1.3983852) q[3];
sx q[3];
rz(-0.93061479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6861787) q[2];
sx q[2];
rz(-1.4718461) q[2];
sx q[2];
rz(-3.0012644) q[2];
rz(0.27648196) q[3];
sx q[3];
rz(-0.77767196) q[3];
sx q[3];
rz(2.8592143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33271933) q[0];
sx q[0];
rz(-2.7820899) q[0];
sx q[0];
rz(1.3746388) q[0];
rz(-0.91284347) q[1];
sx q[1];
rz(-0.54250598) q[1];
sx q[1];
rz(-1.0309781) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.08374005) q[0];
sx q[0];
rz(-0.3763323) q[0];
sx q[0];
rz(2.0876838) q[0];
x q[1];
rz(0.39856492) q[2];
sx q[2];
rz(-2.0953278) q[2];
sx q[2];
rz(0.87510671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1304969) q[1];
sx q[1];
rz(-1.5807932) q[1];
sx q[1];
rz(1.1586055) q[1];
x q[2];
rz(-2.616373) q[3];
sx q[3];
rz(-2.7708443) q[3];
sx q[3];
rz(1.4300089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.012152) q[2];
sx q[2];
rz(-1.0544216) q[2];
sx q[2];
rz(-2.5035109) q[2];
rz(0.10073999) q[3];
sx q[3];
rz(-1.0120069) q[3];
sx q[3];
rz(-2.6428599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3070372) q[0];
sx q[0];
rz(-1.5191673) q[0];
sx q[0];
rz(-1.3822973) q[0];
rz(2.8221829) q[1];
sx q[1];
rz(-1.6326135) q[1];
sx q[1];
rz(0.75739783) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35505193) q[0];
sx q[0];
rz(-3.0664223) q[0];
sx q[0];
rz(-3.0274903) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4057233) q[2];
sx q[2];
rz(-2.9506548) q[2];
sx q[2];
rz(1.5519993) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0590225) q[1];
sx q[1];
rz(-0.76434128) q[1];
sx q[1];
rz(-2.2405626) q[1];
x q[2];
rz(2.4457127) q[3];
sx q[3];
rz(-2.4081488) q[3];
sx q[3];
rz(2.2242198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.77183977) q[2];
sx q[2];
rz(-0.87551337) q[2];
sx q[2];
rz(-0.80779752) q[2];
rz(1.7317023) q[3];
sx q[3];
rz(-1.8818703) q[3];
sx q[3];
rz(2.4894449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0659502) q[0];
sx q[0];
rz(-2.763651) q[0];
sx q[0];
rz(2.156303) q[0];
rz(1.4957042) q[1];
sx q[1];
rz(-2.0031877) q[1];
sx q[1];
rz(0.92934242) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.028469298) q[0];
sx q[0];
rz(-1.3036393) q[0];
sx q[0];
rz(2.5384253) q[0];
x q[1];
rz(2.5337436) q[2];
sx q[2];
rz(-1.2298541) q[2];
sx q[2];
rz(2.8729977) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4080491) q[1];
sx q[1];
rz(-2.0966623) q[1];
sx q[1];
rz(0.2307363) q[1];
x q[2];
rz(-3.1272917) q[3];
sx q[3];
rz(-2.1319444) q[3];
sx q[3];
rz(2.1545269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0890961) q[2];
sx q[2];
rz(-1.7868944) q[2];
sx q[2];
rz(-1.0864786) q[2];
rz(-0.9440445) q[3];
sx q[3];
rz(-3.0510674) q[3];
sx q[3];
rz(1.121678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1679967) q[0];
sx q[0];
rz(-0.43742988) q[0];
sx q[0];
rz(2.9214389) q[0];
rz(1.5036229) q[1];
sx q[1];
rz(-1.3238246) q[1];
sx q[1];
rz(0.55353177) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39237994) q[0];
sx q[0];
rz(-0.22870453) q[0];
sx q[0];
rz(-0.52992757) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36350162) q[2];
sx q[2];
rz(-1.7365626) q[2];
sx q[2];
rz(1.7004418) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5157317) q[1];
sx q[1];
rz(-1.5222094) q[1];
sx q[1];
rz(2.4481985) q[1];
rz(-pi) q[2];
rz(0.68562724) q[3];
sx q[3];
rz(-1.3925288) q[3];
sx q[3];
rz(-2.7741711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0963875) q[2];
sx q[2];
rz(-1.6738946) q[2];
sx q[2];
rz(-0.20827797) q[2];
rz(2.0221209) q[3];
sx q[3];
rz(-1.2866311) q[3];
sx q[3];
rz(2.313405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49523062) q[0];
sx q[0];
rz(-1.3883075) q[0];
sx q[0];
rz(-2.0679423) q[0];
rz(3.011009) q[1];
sx q[1];
rz(-1.5740296) q[1];
sx q[1];
rz(-2.1317587) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45859858) q[0];
sx q[0];
rz(-1.7560335) q[0];
sx q[0];
rz(-0.40249975) q[0];
rz(-pi) q[1];
rz(2.9931941) q[2];
sx q[2];
rz(-1.6434533) q[2];
sx q[2];
rz(-2.1138482) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.1853789) q[1];
sx q[1];
rz(-1.662743) q[1];
sx q[1];
rz(0.57422178) q[1];
rz(2.6920986) q[3];
sx q[3];
rz(-2.426034) q[3];
sx q[3];
rz(-0.64330949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4449571) q[2];
sx q[2];
rz(-2.770165) q[2];
sx q[2];
rz(0.32336393) q[2];
rz(-2.4671386) q[3];
sx q[3];
rz(-0.89265299) q[3];
sx q[3];
rz(-3.118012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12482878) q[0];
sx q[0];
rz(-0.49377307) q[0];
sx q[0];
rz(2.0523409) q[0];
rz(-0.59843868) q[1];
sx q[1];
rz(-1.2804223) q[1];
sx q[1];
rz(0.79900297) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3865144) q[0];
sx q[0];
rz(-1.8505972) q[0];
sx q[0];
rz(-0.31679074) q[0];
rz(-pi) q[1];
rz(1.4510462) q[2];
sx q[2];
rz(-2.471912) q[2];
sx q[2];
rz(-1.0543752) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5889783) q[1];
sx q[1];
rz(-1.2689021) q[1];
sx q[1];
rz(1.7266858) q[1];
rz(-pi) q[2];
rz(-0.022074583) q[3];
sx q[3];
rz(-2.8997053) q[3];
sx q[3];
rz(2.8245408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.66990718) q[2];
sx q[2];
rz(-0.50614637) q[2];
sx q[2];
rz(1.252582) q[2];
rz(-2.0910828) q[3];
sx q[3];
rz(-0.59058467) q[3];
sx q[3];
rz(-0.93241507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2893386) q[0];
sx q[0];
rz(-1.7409356) q[0];
sx q[0];
rz(-2.6259212) q[0];
rz(1.4562666) q[1];
sx q[1];
rz(-1.3412424) q[1];
sx q[1];
rz(-2.8175443) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1103678) q[0];
sx q[0];
rz(-2.661781) q[0];
sx q[0];
rz(-0.57439248) q[0];
rz(-2.134974) q[2];
sx q[2];
rz(-0.7340275) q[2];
sx q[2];
rz(-0.013040868) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0258603) q[1];
sx q[1];
rz(-2.6792791) q[1];
sx q[1];
rz(1.2532534) q[1];
rz(1.0131272) q[3];
sx q[3];
rz(-1.5032167) q[3];
sx q[3];
rz(-2.9014719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1590165) q[2];
sx q[2];
rz(-2.0197208) q[2];
sx q[2];
rz(2.4328361) q[2];
rz(0.7835663) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(-1.5172575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5358955) q[0];
sx q[0];
rz(-2.504183) q[0];
sx q[0];
rz(2.9851483) q[0];
rz(-0.63977301) q[1];
sx q[1];
rz(-0.6548869) q[1];
sx q[1];
rz(-2.1515062) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7759686) q[0];
sx q[0];
rz(-0.98593119) q[0];
sx q[0];
rz(-2.3192899) q[0];
rz(-pi) q[1];
rz(-1.9205356) q[2];
sx q[2];
rz(-1.4863401) q[2];
sx q[2];
rz(1.3049558) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.2459979) q[1];
sx q[1];
rz(-0.57064547) q[1];
sx q[1];
rz(-0.24773189) q[1];
rz(-2.5699432) q[3];
sx q[3];
rz(-2.902973) q[3];
sx q[3];
rz(2.1941136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5618374) q[2];
sx q[2];
rz(-1.729915) q[2];
sx q[2];
rz(-2.5282395) q[2];
rz(-2.8776045) q[3];
sx q[3];
rz(-1.8050906) q[3];
sx q[3];
rz(-2.4700375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.092125208) q[0];
sx q[0];
rz(-1.5277852) q[0];
sx q[0];
rz(1.1396136) q[0];
rz(-1.4641948) q[1];
sx q[1];
rz(-0.72036998) q[1];
sx q[1];
rz(1.5710685) q[1];
rz(1.9665267) q[2];
sx q[2];
rz(-1.1968975) q[2];
sx q[2];
rz(0.71002985) q[2];
rz(0.055394854) q[3];
sx q[3];
rz(-0.61943409) q[3];
sx q[3];
rz(-2.7243765) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
