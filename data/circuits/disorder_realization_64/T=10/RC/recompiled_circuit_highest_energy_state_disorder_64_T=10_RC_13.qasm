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
rz(2.0848701) q[0];
sx q[0];
rz(-1.8615847) q[0];
sx q[0];
rz(-2.4155389) q[0];
rz(-2.5441406) q[1];
sx q[1];
rz(-0.51550454) q[1];
sx q[1];
rz(-1.0322303) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5846378) q[0];
sx q[0];
rz(-0.9171066) q[0];
sx q[0];
rz(1.4435221) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.039041877) q[2];
sx q[2];
rz(-1.0532951) q[2];
sx q[2];
rz(1.7797178) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0192248) q[1];
sx q[1];
rz(-1.6926973) q[1];
sx q[1];
rz(2.8278964) q[1];
rz(-pi) q[2];
rz(-2.0809523) q[3];
sx q[3];
rz(-0.83612305) q[3];
sx q[3];
rz(-1.6250798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.44077474) q[2];
sx q[2];
rz(-0.60443193) q[2];
sx q[2];
rz(0.087892858) q[2];
rz(-1.6739738) q[3];
sx q[3];
rz(-1.847495) q[3];
sx q[3];
rz(0.99883336) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0198332) q[0];
sx q[0];
rz(-1.3061981) q[0];
sx q[0];
rz(-1.649296) q[0];
rz(-2.3536033) q[1];
sx q[1];
rz(-1.183527) q[1];
sx q[1];
rz(-0.96510395) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3639587) q[0];
sx q[0];
rz(-0.84375644) q[0];
sx q[0];
rz(1.4476089) q[0];
rz(0.73678686) q[2];
sx q[2];
rz(-1.4649142) q[2];
sx q[2];
rz(-2.4491765) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.9763042) q[1];
sx q[1];
rz(-1.3852264) q[1];
sx q[1];
rz(1.1052422) q[1];
rz(-pi) q[2];
rz(2.9645637) q[3];
sx q[3];
rz(-0.47613371) q[3];
sx q[3];
rz(-0.04410049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1284156) q[2];
sx q[2];
rz(-2.8030277) q[2];
sx q[2];
rz(1.817912) q[2];
rz(2.2549818) q[3];
sx q[3];
rz(-1.9804136) q[3];
sx q[3];
rz(-2.9429341) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24432467) q[0];
sx q[0];
rz(-2.1710945) q[0];
sx q[0];
rz(0.6955198) q[0];
rz(-0.8356525) q[1];
sx q[1];
rz(-1.4202159) q[1];
sx q[1];
rz(0.64170352) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1276961) q[0];
sx q[0];
rz(-2.3900716) q[0];
sx q[0];
rz(-1.7836545) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4303248) q[2];
sx q[2];
rz(-1.4431653) q[2];
sx q[2];
rz(0.39943275) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1759626) q[1];
sx q[1];
rz(-1.2568736) q[1];
sx q[1];
rz(-2.0759275) q[1];
rz(2.1459177) q[3];
sx q[3];
rz(-2.174587) q[3];
sx q[3];
rz(0.010494516) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.16331638) q[2];
sx q[2];
rz(-1.1288319) q[2];
sx q[2];
rz(1.3345435) q[2];
rz(1.7415907) q[3];
sx q[3];
rz(-1.6440697) q[3];
sx q[3];
rz(-1.1388206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1432994) q[0];
sx q[0];
rz(-2.2704953) q[0];
sx q[0];
rz(2.368108) q[0];
rz(-0.086291226) q[1];
sx q[1];
rz(-0.50420612) q[1];
sx q[1];
rz(-2.6533244) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.047866743) q[0];
sx q[0];
rz(-2.6742088) q[0];
sx q[0];
rz(-2.9899238) q[0];
x q[1];
rz(-0.62520844) q[2];
sx q[2];
rz(-2.3397988) q[2];
sx q[2];
rz(1.9156574) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.5871856) q[1];
sx q[1];
rz(-0.76495612) q[1];
sx q[1];
rz(1.2423963) q[1];
x q[2];
rz(1.1639488) q[3];
sx q[3];
rz(-1.196612) q[3];
sx q[3];
rz(-1.0519303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6733072) q[2];
sx q[2];
rz(-2.0279334) q[2];
sx q[2];
rz(2.9735273) q[2];
rz(2.7189861) q[3];
sx q[3];
rz(-0.97413078) q[3];
sx q[3];
rz(-2.9882123) q[3];
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
rz(2.3733658) q[0];
sx q[0];
rz(-0.90665561) q[0];
sx q[0];
rz(-1.1905097) q[0];
rz(0.52781421) q[1];
sx q[1];
rz(-1.0820791) q[1];
sx q[1];
rz(-3.1324918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86316696) q[0];
sx q[0];
rz(-2.3519807) q[0];
sx q[0];
rz(-0.65152438) q[0];
rz(-pi) q[1];
rz(-1.6318237) q[2];
sx q[2];
rz(-1.6313071) q[2];
sx q[2];
rz(-2.4946545) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2796265) q[1];
sx q[1];
rz(-1.1512966) q[1];
sx q[1];
rz(0.90200938) q[1];
rz(-2.6713773) q[3];
sx q[3];
rz(-2.145916) q[3];
sx q[3];
rz(2.0755538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8268383) q[2];
sx q[2];
rz(-2.645731) q[2];
sx q[2];
rz(-0.66519386) q[2];
rz(2.1151755) q[3];
sx q[3];
rz(-1.0746936) q[3];
sx q[3];
rz(2.3608666) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5196395) q[0];
sx q[0];
rz(-0.53549796) q[0];
sx q[0];
rz(2.7435379) q[0];
rz(-1.6710501) q[1];
sx q[1];
rz(-0.87659756) q[1];
sx q[1];
rz(-0.7801396) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1127284) q[0];
sx q[0];
rz(-2.5327589) q[0];
sx q[0];
rz(-0.91797773) q[0];
rz(-pi) q[1];
x q[1];
rz(0.24274893) q[2];
sx q[2];
rz(-1.374259) q[2];
sx q[2];
rz(-0.32699531) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.4698339) q[1];
sx q[1];
rz(-2.2398754) q[1];
sx q[1];
rz(-1.9356739) q[1];
x q[2];
rz(0.31731242) q[3];
sx q[3];
rz(-0.91083357) q[3];
sx q[3];
rz(2.8803006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.857343) q[2];
sx q[2];
rz(-2.0390859) q[2];
sx q[2];
rz(0.16172376) q[2];
rz(0.11180793) q[3];
sx q[3];
rz(-0.72607741) q[3];
sx q[3];
rz(2.4017754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1424471) q[0];
sx q[0];
rz(-1.4302) q[0];
sx q[0];
rz(-2.4430742) q[0];
rz(-2.0493719) q[1];
sx q[1];
rz(-0.75463808) q[1];
sx q[1];
rz(3.0879367) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92019161) q[0];
sx q[0];
rz(-2.2760253) q[0];
sx q[0];
rz(0.34366519) q[0];
x q[1];
rz(-2.0094107) q[2];
sx q[2];
rz(-1.918186) q[2];
sx q[2];
rz(-1.5127104) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.084890993) q[1];
sx q[1];
rz(-2.6680601) q[1];
sx q[1];
rz(0.11736203) q[1];
rz(-pi) q[2];
rz(2.2174683) q[3];
sx q[3];
rz(-1.4507626) q[3];
sx q[3];
rz(-1.8450575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.0003537) q[2];
sx q[2];
rz(-0.69811368) q[2];
sx q[2];
rz(0.21256438) q[2];
rz(-2.8138748) q[3];
sx q[3];
rz(-2.1543584) q[3];
sx q[3];
rz(-1.7879965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89606365) q[0];
sx q[0];
rz(-1.117417) q[0];
sx q[0];
rz(-2.8118706) q[0];
rz(-2.3700736) q[1];
sx q[1];
rz(-0.78116575) q[1];
sx q[1];
rz(0.82569295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59300789) q[0];
sx q[0];
rz(-2.4720044) q[0];
sx q[0];
rz(3.0010745) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0631376) q[2];
sx q[2];
rz(-2.1657073) q[2];
sx q[2];
rz(1.9878146) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.030071478) q[1];
sx q[1];
rz(-2.3276405) q[1];
sx q[1];
rz(-1.7205283) q[1];
rz(-pi) q[2];
rz(-1.3729893) q[3];
sx q[3];
rz(-1.8231492) q[3];
sx q[3];
rz(2.7276602) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.19994911) q[2];
sx q[2];
rz(-1.7295094) q[2];
sx q[2];
rz(-1.4818209) q[2];
rz(3.0485349) q[3];
sx q[3];
rz(-1.8424282) q[3];
sx q[3];
rz(2.602128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2411497) q[0];
sx q[0];
rz(-2.7261782) q[0];
sx q[0];
rz(-2.7442617) q[0];
rz(0.98995248) q[1];
sx q[1];
rz(-1.3497458) q[1];
sx q[1];
rz(1.0362524) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62972126) q[0];
sx q[0];
rz(-2.7366834) q[0];
sx q[0];
rz(-3.1121503) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.36340275) q[2];
sx q[2];
rz(-0.37523233) q[2];
sx q[2];
rz(-2.4979532) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81779382) q[1];
sx q[1];
rz(-1.7897871) q[1];
sx q[1];
rz(-1.6621672) q[1];
x q[2];
rz(-2.1734222) q[3];
sx q[3];
rz(-0.60634469) q[3];
sx q[3];
rz(1.9023638) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.022126023) q[2];
sx q[2];
rz(-2.4194722) q[2];
sx q[2];
rz(-1.3163346) q[2];
rz(-2.8890166) q[3];
sx q[3];
rz(-1.794869) q[3];
sx q[3];
rz(2.778497) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13755688) q[0];
sx q[0];
rz(-0.79020774) q[0];
sx q[0];
rz(-0.23605119) q[0];
rz(3.0421742) q[1];
sx q[1];
rz(-2.3853018) q[1];
sx q[1];
rz(1.258446) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61533046) q[0];
sx q[0];
rz(-1.0122504) q[0];
sx q[0];
rz(-1.4790776) q[0];
x q[1];
rz(-0.61014097) q[2];
sx q[2];
rz(-2.379866) q[2];
sx q[2];
rz(-1.6118647) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.7459877) q[1];
sx q[1];
rz(-1.6853991) q[1];
sx q[1];
rz(2.458861) q[1];
x q[2];
rz(1.997273) q[3];
sx q[3];
rz(-1.1955373) q[3];
sx q[3];
rz(-2.8464908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.9749757) q[2];
sx q[2];
rz(-0.48305837) q[2];
sx q[2];
rz(-2.100259) q[2];
rz(-2.5552022) q[3];
sx q[3];
rz(-0.47724637) q[3];
sx q[3];
rz(-2.3311116) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36622421) q[0];
sx q[0];
rz(-1.0564221) q[0];
sx q[0];
rz(-1.139241) q[0];
rz(2.1495023) q[1];
sx q[1];
rz(-1.5403668) q[1];
sx q[1];
rz(-1.5425727) q[1];
rz(-0.31789657) q[2];
sx q[2];
rz(-1.7902725) q[2];
sx q[2];
rz(-0.27223311) q[2];
rz(-0.59320368) q[3];
sx q[3];
rz(-1.9402258) q[3];
sx q[3];
rz(2.6711834) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
