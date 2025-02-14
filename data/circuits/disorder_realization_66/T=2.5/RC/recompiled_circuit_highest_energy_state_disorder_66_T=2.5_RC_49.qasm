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
rz(2.0546497) q[0];
rz(2.3776157) q[1];
sx q[1];
rz(-2.6978701) q[1];
sx q[1];
rz(0.8313764) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87693857) q[0];
sx q[0];
rz(-1.042141) q[0];
sx q[0];
rz(-1.009788) q[0];
x q[1];
rz(2.0399206) q[2];
sx q[2];
rz(-1.3587102) q[2];
sx q[2];
rz(2.9123127) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.30637851) q[1];
sx q[1];
rz(-1.6740546) q[1];
sx q[1];
rz(-1.4500965) q[1];
rz(-0.071962996) q[3];
sx q[3];
rz(-1.5441893) q[3];
sx q[3];
rz(0.047544971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5272556) q[2];
sx q[2];
rz(-1.8680806) q[2];
sx q[2];
rz(2.7296076) q[2];
rz(-2.8968503) q[3];
sx q[3];
rz(-2.5089846) q[3];
sx q[3];
rz(-2.8542724) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43993846) q[0];
sx q[0];
rz(-0.97173062) q[0];
sx q[0];
rz(0.42020759) q[0];
rz(2.6523051) q[1];
sx q[1];
rz(-0.91541618) q[1];
sx q[1];
rz(1.9583826) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32589697) q[0];
sx q[0];
rz(-0.95615471) q[0];
sx q[0];
rz(2.6026898) q[0];
rz(-0.31416201) q[2];
sx q[2];
rz(-1.8793686) q[2];
sx q[2];
rz(-1.4663638) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4192761) q[1];
sx q[1];
rz(-2.4240383) q[1];
sx q[1];
rz(2.5408387) q[1];
rz(-pi) q[2];
rz(-1.39721) q[3];
sx q[3];
rz(-1.4549482) q[3];
sx q[3];
rz(0.66044745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6861787) q[2];
sx q[2];
rz(-1.4718461) q[2];
sx q[2];
rz(-0.14032826) q[2];
rz(-2.8651107) q[3];
sx q[3];
rz(-0.77767196) q[3];
sx q[3];
rz(2.8592143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33271933) q[0];
sx q[0];
rz(-0.35950279) q[0];
sx q[0];
rz(1.3746388) q[0];
rz(-2.2287492) q[1];
sx q[1];
rz(-2.5990867) q[1];
sx q[1];
rz(2.1106145) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1682273) q[0];
sx q[0];
rz(-1.388167) q[0];
sx q[0];
rz(1.9017066) q[0];
x q[1];
rz(1.0101843) q[2];
sx q[2];
rz(-1.9133781) q[2];
sx q[2];
rz(-0.48784384) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1304969) q[1];
sx q[1];
rz(-1.5607995) q[1];
sx q[1];
rz(-1.1586055) q[1];
rz(-pi) q[2];
rz(-2.8171478) q[3];
sx q[3];
rz(-1.7534755) q[3];
sx q[3];
rz(0.35445359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.012152) q[2];
sx q[2];
rz(-1.0544216) q[2];
sx q[2];
rz(-0.63808179) q[2];
rz(0.10073999) q[3];
sx q[3];
rz(-1.0120069) q[3];
sx q[3];
rz(-2.6428599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3070372) q[0];
sx q[0];
rz(-1.6224253) q[0];
sx q[0];
rz(1.7592953) q[0];
rz(-2.8221829) q[1];
sx q[1];
rz(-1.5089792) q[1];
sx q[1];
rz(0.75739783) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7865407) q[0];
sx q[0];
rz(-0.075170366) q[0];
sx q[0];
rz(3.0274903) q[0];
rz(2.9992835) q[2];
sx q[2];
rz(-1.6985296) q[2];
sx q[2];
rz(2.3960428) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.74906236) q[1];
sx q[1];
rz(-2.1442701) q[1];
sx q[1];
rz(0.5368781) q[1];
rz(-2.5365165) q[3];
sx q[3];
rz(-2.0143442) q[3];
sx q[3];
rz(-0.098066948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3697529) q[2];
sx q[2];
rz(-0.87551337) q[2];
sx q[2];
rz(-0.80779752) q[2];
rz(1.4098903) q[3];
sx q[3];
rz(-1.8818703) q[3];
sx q[3];
rz(-2.4894449) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0756425) q[0];
sx q[0];
rz(-0.37794161) q[0];
sx q[0];
rz(2.156303) q[0];
rz(-1.4957042) q[1];
sx q[1];
rz(-1.138405) q[1];
sx q[1];
rz(-2.2122502) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2332704) q[0];
sx q[0];
rz(-0.65289557) q[0];
sx q[0];
rz(2.6920431) q[0];
x q[1];
rz(2.5337436) q[2];
sx q[2];
rz(-1.9117386) q[2];
sx q[2];
rz(0.268595) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4080491) q[1];
sx q[1];
rz(-2.0966623) q[1];
sx q[1];
rz(-0.2307363) q[1];
rz(-pi) q[2];
rz(-3.1272917) q[3];
sx q[3];
rz(-2.1319444) q[3];
sx q[3];
rz(2.1545269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0890961) q[2];
sx q[2];
rz(-1.3546983) q[2];
sx q[2];
rz(1.0864786) q[2];
rz(2.1975482) q[3];
sx q[3];
rz(-3.0510674) q[3];
sx q[3];
rz(1.121678) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1679967) q[0];
sx q[0];
rz(-2.7041628) q[0];
sx q[0];
rz(-2.9214389) q[0];
rz(1.5036229) q[1];
sx q[1];
rz(-1.3238246) q[1];
sx q[1];
rz(0.55353177) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7492127) q[0];
sx q[0];
rz(-0.22870453) q[0];
sx q[0];
rz(-2.6116651) q[0];
rz(-pi) q[1];
rz(0.36350162) q[2];
sx q[2];
rz(-1.40503) q[2];
sx q[2];
rz(1.4411508) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.014723881) q[1];
sx q[1];
rz(-2.2632102) q[1];
sx q[1];
rz(1.5076553) q[1];
rz(-pi) q[2];
x q[2];
rz(1.342085) q[3];
sx q[3];
rz(-0.89806496) q[3];
sx q[3];
rz(-1.3474238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0963875) q[2];
sx q[2];
rz(-1.4676981) q[2];
sx q[2];
rz(-2.9333147) q[2];
rz(-1.1194718) q[3];
sx q[3];
rz(-1.8549615) q[3];
sx q[3];
rz(0.82818762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49523062) q[0];
sx q[0];
rz(-1.7532852) q[0];
sx q[0];
rz(-2.0679423) q[0];
rz(-0.13058361) q[1];
sx q[1];
rz(-1.5675631) q[1];
sx q[1];
rz(2.1317587) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9511418) q[0];
sx q[0];
rz(-1.1755687) q[0];
sx q[0];
rz(-1.3698831) q[0];
rz(-pi) q[1];
rz(1.4973348) q[2];
sx q[2];
rz(-1.7188004) q[2];
sx q[2];
rz(-2.5876887) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.9562137) q[1];
sx q[1];
rz(-1.662743) q[1];
sx q[1];
rz(0.57422178) q[1];
x q[2];
rz(1.9319264) q[3];
sx q[3];
rz(-2.2029366) q[3];
sx q[3];
rz(1.9295285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.69663557) q[2];
sx q[2];
rz(-0.37142763) q[2];
sx q[2];
rz(-0.32336393) q[2];
rz(0.67445406) q[3];
sx q[3];
rz(-2.2489397) q[3];
sx q[3];
rz(3.118012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12482878) q[0];
sx q[0];
rz(-2.6478196) q[0];
sx q[0];
rz(-1.0892518) q[0];
rz(-2.543154) q[1];
sx q[1];
rz(-1.8611703) q[1];
sx q[1];
rz(0.79900297) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.625811) q[0];
sx q[0];
rz(-2.7220753) q[0];
sx q[0];
rz(-0.74504344) q[0];
rz(-pi) q[1];
rz(-1.6905464) q[2];
sx q[2];
rz(-0.66968067) q[2];
sx q[2];
rz(-2.0872175) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5889783) q[1];
sx q[1];
rz(-1.2689021) q[1];
sx q[1];
rz(-1.7266858) q[1];
rz(-pi) q[2];
x q[2];
rz(2.899762) q[3];
sx q[3];
rz(-1.5655091) q[3];
sx q[3];
rz(1.9092803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.66990718) q[2];
sx q[2];
rz(-2.6354463) q[2];
sx q[2];
rz(-1.8890107) q[2];
rz(-1.0505098) q[3];
sx q[3];
rz(-2.551008) q[3];
sx q[3];
rz(-0.93241507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-2.2893386) q[0];
sx q[0];
rz(-1.7409356) q[0];
sx q[0];
rz(0.5156714) q[0];
rz(-1.6853261) q[1];
sx q[1];
rz(-1.3412424) q[1];
sx q[1];
rz(-2.8175443) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7406612) q[0];
sx q[0];
rz(-1.9687511) q[0];
sx q[0];
rz(1.8463355) q[0];
rz(-2.6921047) q[2];
sx q[2];
rz(-0.96908334) q[2];
sx q[2];
rz(0.692761) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9726213) q[1];
sx q[1];
rz(-1.7105127) q[1];
sx q[1];
rz(1.1286382) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0619929) q[3];
sx q[3];
rz(-2.1270424) q[3];
sx q[3];
rz(1.85301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9825762) q[2];
sx q[2];
rz(-1.1218718) q[2];
sx q[2];
rz(-0.70875657) q[2];
rz(-2.3580264) q[3];
sx q[3];
rz(-1.9400027) q[3];
sx q[3];
rz(1.6243352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5358955) q[0];
sx q[0];
rz(-0.63740969) q[0];
sx q[0];
rz(-0.15644431) q[0];
rz(2.5018196) q[1];
sx q[1];
rz(-2.4867058) q[1];
sx q[1];
rz(-0.99008647) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7759686) q[0];
sx q[0];
rz(-2.1556615) q[0];
sx q[0];
rz(-2.3192899) q[0];
rz(1.813011) q[2];
sx q[2];
rz(-0.35938423) q[2];
sx q[2];
rz(0.49312544) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.53452021) q[1];
sx q[1];
rz(-1.4379518) q[1];
sx q[1];
rz(-0.55664363) q[1];
rz(-pi) q[2];
rz(1.7016497) q[3];
sx q[3];
rz(-1.7709117) q[3];
sx q[3];
rz(0.36267422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.5797552) q[2];
sx q[2];
rz(-1.4116776) q[2];
sx q[2];
rz(2.5282395) q[2];
rz(-0.26398811) q[3];
sx q[3];
rz(-1.3365021) q[3];
sx q[3];
rz(-2.4700375) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0494674) q[0];
sx q[0];
rz(-1.5277852) q[0];
sx q[0];
rz(1.1396136) q[0];
rz(1.4641948) q[1];
sx q[1];
rz(-2.4212227) q[1];
sx q[1];
rz(-1.5705241) q[1];
rz(2.3650305) q[2];
sx q[2];
rz(-0.53755098) q[2];
sx q[2];
rz(2.9992486) q[2];
rz(2.522884) q[3];
sx q[3];
rz(-1.6029463) q[3];
sx q[3];
rz(-1.1986986) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
