OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.1563675) q[0];
sx q[0];
rz(1.8591783) q[0];
sx q[0];
rz(9.5556762) q[0];
rz(0.52783293) q[1];
sx q[1];
rz(3.5117709) q[1];
sx q[1];
rz(12.389046) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89761868) q[0];
sx q[0];
rz(-1.8332687) q[0];
sx q[0];
rz(2.3489683) q[0];
rz(-pi) q[1];
rz(-2.913397) q[2];
sx q[2];
rz(-1.8642117) q[2];
sx q[2];
rz(2.5876837) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.51715358) q[1];
sx q[1];
rz(-0.53688795) q[1];
sx q[1];
rz(1.0978903) q[1];
rz(-pi) q[2];
rz(-2.1471094) q[3];
sx q[3];
rz(-1.8015993) q[3];
sx q[3];
rz(-2.1891914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52239546) q[2];
sx q[2];
rz(-2.0530687) q[2];
sx q[2];
rz(-1.8001451) q[2];
rz(1.7893192) q[3];
sx q[3];
rz(-1.5315346) q[3];
sx q[3];
rz(1.8488098) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.751048) q[0];
sx q[0];
rz(-1.9123257) q[0];
sx q[0];
rz(2.6265662) q[0];
rz(0.36188778) q[1];
sx q[1];
rz(-1.7444976) q[1];
sx q[1];
rz(0.99641689) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5866885) q[0];
sx q[0];
rz(-1.3350272) q[0];
sx q[0];
rz(0.20637189) q[0];
x q[1];
rz(-2.9662232) q[2];
sx q[2];
rz(-1.1423472) q[2];
sx q[2];
rz(1.0340978) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.540058) q[1];
sx q[1];
rz(-0.52943474) q[1];
sx q[1];
rz(-1.7537746) q[1];
rz(-pi) q[2];
rz(2.6985136) q[3];
sx q[3];
rz(-1.0765809) q[3];
sx q[3];
rz(2.6576603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.4126052) q[2];
sx q[2];
rz(-0.94634405) q[2];
sx q[2];
rz(-1.6607355) q[2];
rz(0.49731538) q[3];
sx q[3];
rz(-2.1931084) q[3];
sx q[3];
rz(0.724154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5363252) q[0];
sx q[0];
rz(-1.8999506) q[0];
sx q[0];
rz(1.190825) q[0];
rz(-2.7438927) q[1];
sx q[1];
rz(-1.560248) q[1];
sx q[1];
rz(-0.89888987) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5918698) q[0];
sx q[0];
rz(-1.6373349) q[0];
sx q[0];
rz(-1.9107242) q[0];
x q[1];
rz(1.6562881) q[2];
sx q[2];
rz(-1.7608661) q[2];
sx q[2];
rz(-2.9051733) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0278416) q[1];
sx q[1];
rz(-1.5261302) q[1];
sx q[1];
rz(1.0381446) q[1];
rz(-1.8679138) q[3];
sx q[3];
rz(-0.53830244) q[3];
sx q[3];
rz(2.3915714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0209823) q[2];
sx q[2];
rz(-1.9526498) q[2];
sx q[2];
rz(-2.9494542) q[2];
rz(-2.4118679) q[3];
sx q[3];
rz(-0.14484043) q[3];
sx q[3];
rz(2.5016968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9925053) q[0];
sx q[0];
rz(-2.3893116) q[0];
sx q[0];
rz(1.3080904) q[0];
rz(-2.2970842) q[1];
sx q[1];
rz(-1.4627855) q[1];
sx q[1];
rz(0.62215296) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7956411) q[0];
sx q[0];
rz(-1.9824195) q[0];
sx q[0];
rz(-0.53070416) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8573158) q[2];
sx q[2];
rz(-0.19294365) q[2];
sx q[2];
rz(-2.3920453) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.51650652) q[1];
sx q[1];
rz(-1.7343905) q[1];
sx q[1];
rz(-1.518599) q[1];
x q[2];
rz(2.8243869) q[3];
sx q[3];
rz(-1.6807846) q[3];
sx q[3];
rz(3.0134137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.017612351) q[2];
sx q[2];
rz(-1.7698741) q[2];
sx q[2];
rz(2.2464216) q[2];
rz(-0.34919843) q[3];
sx q[3];
rz(-2.5694191) q[3];
sx q[3];
rz(0.23381843) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8643841) q[0];
sx q[0];
rz(-1.9338436) q[0];
sx q[0];
rz(-0.54339093) q[0];
rz(3.0282989) q[1];
sx q[1];
rz(-1.3985059) q[1];
sx q[1];
rz(-1.8494122) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6627412) q[0];
sx q[0];
rz(-2.2355192) q[0];
sx q[0];
rz(0.49481884) q[0];
rz(-pi) q[1];
rz(1.0903484) q[2];
sx q[2];
rz(-2.4329429) q[2];
sx q[2];
rz(-0.50856579) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1317392) q[1];
sx q[1];
rz(-0.92127548) q[1];
sx q[1];
rz(2.3765537) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4025147) q[3];
sx q[3];
rz(-1.0024286) q[3];
sx q[3];
rz(0.06423244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.31305227) q[2];
sx q[2];
rz(-1.7834168) q[2];
sx q[2];
rz(0.48946112) q[2];
rz(-2.475907) q[3];
sx q[3];
rz(-0.61167115) q[3];
sx q[3];
rz(-2.3556975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8996443) q[0];
sx q[0];
rz(-2.2388832) q[0];
sx q[0];
rz(-0.06074252) q[0];
rz(1.863106) q[1];
sx q[1];
rz(-2.2272031) q[1];
sx q[1];
rz(-1.0260822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36878186) q[0];
sx q[0];
rz(-2.368481) q[0];
sx q[0];
rz(2.2863273) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8158978) q[2];
sx q[2];
rz(-2.1148588) q[2];
sx q[2];
rz(-1.0808672) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.4602604) q[1];
sx q[1];
rz(-1.7074013) q[1];
sx q[1];
rz(0.25306074) q[1];
rz(-pi) q[2];
rz(0.52512759) q[3];
sx q[3];
rz(-2.5487102) q[3];
sx q[3];
rz(1.8754995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8817899) q[2];
sx q[2];
rz(-1.5851574) q[2];
sx q[2];
rz(-0.078710236) q[2];
rz(1.6591266) q[3];
sx q[3];
rz(-0.31646287) q[3];
sx q[3];
rz(-0.44221529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-0.92597961) q[0];
sx q[0];
rz(-2.7613566) q[0];
sx q[0];
rz(1.6967787) q[0];
rz(2.3346057) q[1];
sx q[1];
rz(-1.746256) q[1];
sx q[1];
rz(-3.1296465) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58633198) q[0];
sx q[0];
rz(-1.5975448) q[0];
sx q[0];
rz(-0.93179758) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9482466) q[2];
sx q[2];
rz(-2.6259661) q[2];
sx q[2];
rz(2.8597997) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0199737) q[1];
sx q[1];
rz(-1.498292) q[1];
sx q[1];
rz(-2.1825779) q[1];
x q[2];
rz(-1.1720522) q[3];
sx q[3];
rz(-1.4280025) q[3];
sx q[3];
rz(2.2384584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.5936467) q[2];
sx q[2];
rz(-2.839489) q[2];
sx q[2];
rz(-0.83615237) q[2];
rz(-3.0892843) q[3];
sx q[3];
rz(-1.6736504) q[3];
sx q[3];
rz(-0.93594319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9713822) q[0];
sx q[0];
rz(-2.2844071) q[0];
sx q[0];
rz(0.49825391) q[0];
rz(-2.1360629) q[1];
sx q[1];
rz(-1.6906831) q[1];
sx q[1];
rz(3.0301869) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94994369) q[0];
sx q[0];
rz(-0.43131098) q[0];
sx q[0];
rz(-1.8032719) q[0];
x q[1];
rz(-2.7132052) q[2];
sx q[2];
rz(-1.7991788) q[2];
sx q[2];
rz(-1.0699748) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-3.1131683) q[1];
sx q[1];
rz(-2.3481124) q[1];
sx q[1];
rz(-1.9774578) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4440246) q[3];
sx q[3];
rz(-1.2283192) q[3];
sx q[3];
rz(2.3449183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.19693836) q[2];
sx q[2];
rz(-1.3613746) q[2];
sx q[2];
rz(1.4062175) q[2];
rz(-1.1574636) q[3];
sx q[3];
rz(-2.1486053) q[3];
sx q[3];
rz(-1.5520613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.120753) q[0];
sx q[0];
rz(-2.4036305) q[0];
sx q[0];
rz(-1.7027759) q[0];
rz(2.2976047) q[1];
sx q[1];
rz(-1.3071209) q[1];
sx q[1];
rz(-2.9356975) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5000847) q[0];
sx q[0];
rz(-2.5091268) q[0];
sx q[0];
rz(3.0728805) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0768422) q[2];
sx q[2];
rz(-2.4528385) q[2];
sx q[2];
rz(-1.233686) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3033324) q[1];
sx q[1];
rz(-2.9853959) q[1];
sx q[1];
rz(-0.52429838) q[1];
rz(-0.4959373) q[3];
sx q[3];
rz(-1.933799) q[3];
sx q[3];
rz(-2.7285226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0042808) q[2];
sx q[2];
rz(-2.0840804) q[2];
sx q[2];
rz(1.6005969) q[2];
rz(2.6565523) q[3];
sx q[3];
rz(-1.78616) q[3];
sx q[3];
rz(0.11390991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6523022) q[0];
sx q[0];
rz(-3.089383) q[0];
sx q[0];
rz(-1.8452277) q[0];
rz(-2.0908053) q[1];
sx q[1];
rz(-1.6810345) q[1];
sx q[1];
rz(-2.5349862) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3171995) q[0];
sx q[0];
rz(-1.6034295) q[0];
sx q[0];
rz(1.7343108) q[0];
rz(0.20175378) q[2];
sx q[2];
rz(-0.56900185) q[2];
sx q[2];
rz(-3.0425827) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7672528) q[1];
sx q[1];
rz(-1.2657968) q[1];
sx q[1];
rz(-1.6348331) q[1];
rz(1.4539155) q[3];
sx q[3];
rz(-1.8532527) q[3];
sx q[3];
rz(0.64843897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4198833) q[2];
sx q[2];
rz(-1.1315283) q[2];
sx q[2];
rz(2.4427872) q[2];
rz(-1.6803668) q[3];
sx q[3];
rz(-2.1278087) q[3];
sx q[3];
rz(-2.6644126) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89194311) q[0];
sx q[0];
rz(-1.5105381) q[0];
sx q[0];
rz(1.426209) q[0];
rz(-2.7632948) q[1];
sx q[1];
rz(-2.5143647) q[1];
sx q[1];
rz(-0.73077269) q[1];
rz(-0.30834352) q[2];
sx q[2];
rz(-2.0258122) q[2];
sx q[2];
rz(-3.0292055) q[2];
rz(-1.6017687) q[3];
sx q[3];
rz(-1.0145368) q[3];
sx q[3];
rz(-1.9817286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
