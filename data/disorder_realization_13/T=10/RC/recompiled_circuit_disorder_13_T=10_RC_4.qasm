OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.17833248) q[0];
sx q[0];
rz(-1.4890716) q[0];
sx q[0];
rz(2.2464377) q[0];
rz(2.826638) q[1];
sx q[1];
rz(-1.0840253) q[1];
sx q[1];
rz(-1.4562343) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0460912) q[0];
sx q[0];
rz(-1.4399733) q[0];
sx q[0];
rz(0.019898947) q[0];
rz(2.1452791) q[2];
sx q[2];
rz(-0.62553863) q[2];
sx q[2];
rz(1.1686981) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.40741062) q[1];
sx q[1];
rz(-1.4200746) q[1];
sx q[1];
rz(-2.461344) q[1];
rz(-pi) q[2];
rz(-2.9772894) q[3];
sx q[3];
rz(-2.8176753) q[3];
sx q[3];
rz(-1.2795554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3866117) q[2];
sx q[2];
rz(-1.745801) q[2];
sx q[2];
rz(-0.45271978) q[2];
rz(2.9833941) q[3];
sx q[3];
rz(-0.69163624) q[3];
sx q[3];
rz(0.89481568) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2125856) q[0];
sx q[0];
rz(-2.043262) q[0];
sx q[0];
rz(-1.1520977) q[0];
rz(-1.903803) q[1];
sx q[1];
rz(-1.5367616) q[1];
sx q[1];
rz(-0.47098413) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5841056) q[0];
sx q[0];
rz(-0.4852681) q[0];
sx q[0];
rz(-1.6886061) q[0];
rz(-pi) q[1];
rz(-1.2330301) q[2];
sx q[2];
rz(-2.517189) q[2];
sx q[2];
rz(2.117346) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2742548) q[1];
sx q[1];
rz(-1.6254289) q[1];
sx q[1];
rz(-2.0591303) q[1];
rz(1.4906293) q[3];
sx q[3];
rz(-0.87682322) q[3];
sx q[3];
rz(0.44192867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0835138) q[2];
sx q[2];
rz(-2.5364272) q[2];
sx q[2];
rz(-0.2557959) q[2];
rz(-1.6563709) q[3];
sx q[3];
rz(-1.1896313) q[3];
sx q[3];
rz(-0.16168693) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26329041) q[0];
sx q[0];
rz(-1.1061763) q[0];
sx q[0];
rz(2.8702452) q[0];
rz(2.4052606) q[1];
sx q[1];
rz(-1.5356179) q[1];
sx q[1];
rz(2.7022865) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8205748) q[0];
sx q[0];
rz(-1.4903755) q[0];
sx q[0];
rz(3.1119425) q[0];
x q[1];
rz(1.5191684) q[2];
sx q[2];
rz(-1.2106967) q[2];
sx q[2];
rz(-2.1215631) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1767133) q[1];
sx q[1];
rz(-2.5039154) q[1];
sx q[1];
rz(1.8646851) q[1];
rz(-pi) q[2];
rz(0.4337173) q[3];
sx q[3];
rz(-2.0842413) q[3];
sx q[3];
rz(-2.123326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1674041) q[2];
sx q[2];
rz(-1.5891275) q[2];
sx q[2];
rz(-2.9411194) q[2];
rz(0.75508562) q[3];
sx q[3];
rz(-2.1217767) q[3];
sx q[3];
rz(-0.42373207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.077483594) q[0];
sx q[0];
rz(-1.5923201) q[0];
sx q[0];
rz(-1.2444929) q[0];
rz(2.3311133) q[1];
sx q[1];
rz(-1.3296209) q[1];
sx q[1];
rz(0.8746075) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32828242) q[0];
sx q[0];
rz(-0.62424849) q[0];
sx q[0];
rz(1.8391795) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7941197) q[2];
sx q[2];
rz(-0.407019) q[2];
sx q[2];
rz(2.4328872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.306327) q[1];
sx q[1];
rz(-1.7976465) q[1];
sx q[1];
rz(0.38362417) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0365385) q[3];
sx q[3];
rz(-1.4964364) q[3];
sx q[3];
rz(0.33945938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.13742927) q[2];
sx q[2];
rz(-2.2518297) q[2];
sx q[2];
rz(1.7144263) q[2];
rz(3.0754722) q[3];
sx q[3];
rz(-0.36589208) q[3];
sx q[3];
rz(1.6920413) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3604597) q[0];
sx q[0];
rz(-2.9678678) q[0];
sx q[0];
rz(-2.5710035) q[0];
rz(-0.55496201) q[1];
sx q[1];
rz(-2.404232) q[1];
sx q[1];
rz(0.74329174) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63440454) q[0];
sx q[0];
rz(-2.2162063) q[0];
sx q[0];
rz(-0.27642823) q[0];
rz(-pi) q[1];
rz(-2.6402316) q[2];
sx q[2];
rz(-1.6622346) q[2];
sx q[2];
rz(-2.3102592) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35768269) q[1];
sx q[1];
rz(-0.3158814) q[1];
sx q[1];
rz(0.79343474) q[1];
rz(0.80446135) q[3];
sx q[3];
rz(-1.3692229) q[3];
sx q[3];
rz(-1.3988914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1935929) q[2];
sx q[2];
rz(-1.3698545) q[2];
sx q[2];
rz(-0.26838475) q[2];
rz(2.0949481) q[3];
sx q[3];
rz(-2.7768551) q[3];
sx q[3];
rz(2.823765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.627581) q[0];
sx q[0];
rz(-1.528897) q[0];
sx q[0];
rz(0.50672379) q[0];
rz(-0.2535893) q[1];
sx q[1];
rz(-1.8702303) q[1];
sx q[1];
rz(2.0862897) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9297766) q[0];
sx q[0];
rz(-1.7535216) q[0];
sx q[0];
rz(-2.1972448) q[0];
rz(-pi) q[1];
rz(-2.2248473) q[2];
sx q[2];
rz(-1.142821) q[2];
sx q[2];
rz(0.59631729) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.036990449) q[1];
sx q[1];
rz(-1.9676419) q[1];
sx q[1];
rz(-2.2657822) q[1];
rz(-2.8860693) q[3];
sx q[3];
rz(-2.3821085) q[3];
sx q[3];
rz(1.7373191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.48173299) q[2];
sx q[2];
rz(-1.0446171) q[2];
sx q[2];
rz(2.2591023) q[2];
rz(-2.4957538) q[3];
sx q[3];
rz(-1.1487938) q[3];
sx q[3];
rz(-1.8576436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(-1.5053453) q[0];
sx q[0];
rz(-1.9286276) q[0];
sx q[0];
rz(2.3690467) q[0];
rz(1.729471) q[1];
sx q[1];
rz(-1.9344784) q[1];
sx q[1];
rz(-2.5678182) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6360146) q[0];
sx q[0];
rz(-1.6806707) q[0];
sx q[0];
rz(-1.3789603) q[0];
rz(-3.119486) q[2];
sx q[2];
rz(-1.3706285) q[2];
sx q[2];
rz(-2.9526763) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3473914) q[1];
sx q[1];
rz(-0.83354356) q[1];
sx q[1];
rz(-2.3458523) q[1];
rz(1.0049099) q[3];
sx q[3];
rz(-0.27540576) q[3];
sx q[3];
rz(-2.3435081) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.039375719) q[2];
sx q[2];
rz(-2.6874459) q[2];
sx q[2];
rz(0.77073628) q[2];
rz(-0.43631521) q[3];
sx q[3];
rz(-1.8728914) q[3];
sx q[3];
rz(-1.8541981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8687826) q[0];
sx q[0];
rz(-1.9406809) q[0];
sx q[0];
rz(-1.1707206) q[0];
rz(-0.51013485) q[1];
sx q[1];
rz(-1.3565823) q[1];
sx q[1];
rz(1.8458813) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0062795) q[0];
sx q[0];
rz(-1.8989925) q[0];
sx q[0];
rz(-2.6562064) q[0];
rz(-2.2175118) q[2];
sx q[2];
rz(-0.43281049) q[2];
sx q[2];
rz(0.88976394) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1295373) q[1];
sx q[1];
rz(-0.66857282) q[1];
sx q[1];
rz(-1.2052016) q[1];
rz(-pi) q[2];
rz(1.6298953) q[3];
sx q[3];
rz(-2.700138) q[3];
sx q[3];
rz(-2.3578701) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.43508139) q[2];
sx q[2];
rz(-1.4979829) q[2];
sx q[2];
rz(0.56813017) q[2];
rz(-0.98012296) q[3];
sx q[3];
rz(-1.0390493) q[3];
sx q[3];
rz(2.9476416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4762964) q[0];
sx q[0];
rz(-2.3275573) q[0];
sx q[0];
rz(-0.67767674) q[0];
rz(-2.9455345) q[1];
sx q[1];
rz(-2.129107) q[1];
sx q[1];
rz(2.303404) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3121376) q[0];
sx q[0];
rz(-1.2078309) q[0];
sx q[0];
rz(1.3296933) q[0];
rz(-1.499275) q[2];
sx q[2];
rz(-2.3841249) q[2];
sx q[2];
rz(-1.0345936) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.2182525) q[1];
sx q[1];
rz(-0.76949161) q[1];
sx q[1];
rz(-1.283265) q[1];
rz(-1.4872876) q[3];
sx q[3];
rz(-1.828308) q[3];
sx q[3];
rz(2.6058692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.80205408) q[2];
sx q[2];
rz(-0.69028091) q[2];
sx q[2];
rz(-1.8748803) q[2];
rz(-2.1789815) q[3];
sx q[3];
rz(-1.5701141) q[3];
sx q[3];
rz(0.094749711) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72120136) q[0];
sx q[0];
rz(-1.9045916) q[0];
sx q[0];
rz(0.23751968) q[0];
rz(2.1233842) q[1];
sx q[1];
rz(-0.84914452) q[1];
sx q[1];
rz(-2.9097897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34323877) q[0];
sx q[0];
rz(-0.51734561) q[0];
sx q[0];
rz(-1.7645287) q[0];
rz(-2.899029) q[2];
sx q[2];
rz(-1.3360268) q[2];
sx q[2];
rz(-2.501542) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.787549) q[1];
sx q[1];
rz(-1.7122867) q[1];
sx q[1];
rz(-0.267412) q[1];
x q[2];
rz(-1.3123355) q[3];
sx q[3];
rz(-2.4371394) q[3];
sx q[3];
rz(-0.36650141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1101749) q[2];
sx q[2];
rz(-1.2538223) q[2];
sx q[2];
rz(-1.0882264) q[2];
rz(2.7534289) q[3];
sx q[3];
rz(-2.4813014) q[3];
sx q[3];
rz(0.80374074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5647472) q[0];
sx q[0];
rz(-1.7871465) q[0];
sx q[0];
rz(-0.47252895) q[0];
rz(2.172773) q[1];
sx q[1];
rz(-2.4333654) q[1];
sx q[1];
rz(-2.416837) q[1];
rz(-0.3248365) q[2];
sx q[2];
rz(-1.799121) q[2];
sx q[2];
rz(0.15098235) q[2];
rz(-0.10235056) q[3];
sx q[3];
rz(-0.1889189) q[3];
sx q[3];
rz(1.2253996) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
