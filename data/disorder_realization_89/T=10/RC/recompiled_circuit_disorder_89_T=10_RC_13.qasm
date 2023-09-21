OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3857631) q[0];
sx q[0];
rz(-1.7321777) q[0];
sx q[0];
rz(-2.8470319) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(2.6309738) q[1];
sx q[1];
rz(9.0030158) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.784214) q[0];
sx q[0];
rz(-1.5970018) q[0];
sx q[0];
rz(1.5062917) q[0];
rz(1.5316891) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(1.6959015) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0114853) q[1];
sx q[1];
rz(-1.752578) q[1];
sx q[1];
rz(0.45244042) q[1];
rz(-pi) q[2];
rz(-2.9362039) q[3];
sx q[3];
rz(-0.60384149) q[3];
sx q[3];
rz(-0.83128923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0191779) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(-2.1315234) q[2];
rz(1.646237) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(-8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(0.85900599) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(-1.4650311) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61629907) q[0];
sx q[0];
rz(-1.787961) q[0];
sx q[0];
rz(-2.8263682) q[0];
rz(-1.9081406) q[2];
sx q[2];
rz(-2.0806899) q[2];
sx q[2];
rz(-2.9166729) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3611991) q[1];
sx q[1];
rz(-2.2269899) q[1];
sx q[1];
rz(-1.111163) q[1];
rz(1.2831849) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(2.9420497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7039965) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(-1.0393418) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(2.5487652) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52477437) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(0.16361374) q[0];
rz(0.71584654) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(1.7680426) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4931902) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(0.52669749) q[0];
x q[1];
rz(0.85767304) q[2];
sx q[2];
rz(-1.8048986) q[2];
sx q[2];
rz(0.26299325) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8858582) q[1];
sx q[1];
rz(-2.246937) q[1];
sx q[1];
rz(-0.83791344) q[1];
x q[2];
rz(-0.22294873) q[3];
sx q[3];
rz(-1.3461337) q[3];
sx q[3];
rz(-2.6509681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(1.019657) q[2];
rz(-0.9807469) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402128) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(2.3024978) q[0];
rz(1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(2.531321) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2903039) q[0];
sx q[0];
rz(-1.5595058) q[0];
sx q[0];
rz(1.6800866) q[0];
x q[1];
rz(-0.82579124) q[2];
sx q[2];
rz(-0.50160393) q[2];
sx q[2];
rz(2.290291) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8295633) q[1];
sx q[1];
rz(-0.89013541) q[1];
sx q[1];
rz(-2.8042253) q[1];
rz(-0.88704349) q[3];
sx q[3];
rz(-0.700637) q[3];
sx q[3];
rz(0.62473434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(-2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(-1.5326112) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.8283432) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4603235) q[0];
sx q[0];
rz(-1.2460421) q[0];
sx q[0];
rz(0.89694174) q[0];
rz(-2.1336742) q[2];
sx q[2];
rz(-3.0745227) q[2];
sx q[2];
rz(1.1941393) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9338867) q[1];
sx q[1];
rz(-1.943214) q[1];
sx q[1];
rz(-1.9718584) q[1];
rz(-pi) q[2];
rz(-2.26236) q[3];
sx q[3];
rz(-2.0428265) q[3];
sx q[3];
rz(1.2513245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.1398754) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(0.66429794) q[2];
rz(-2.4364566) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(-0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649081) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(3.0867807) q[0];
rz(2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(0.48318133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38835634) q[0];
sx q[0];
rz(-1.6325132) q[0];
sx q[0];
rz(1.7524936) q[0];
rz(-pi) q[1];
rz(-2.0334843) q[2];
sx q[2];
rz(-1.6001284) q[2];
sx q[2];
rz(0.5043475) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9890081) q[1];
sx q[1];
rz(-2.2629988) q[1];
sx q[1];
rz(-2.2425376) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11404927) q[3];
sx q[3];
rz(-1.9012791) q[3];
sx q[3];
rz(1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(0.085263578) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(-0.15032642) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36713704) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(0.95463395) q[0];
rz(-2.5700991) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(0.25209299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.945767) q[0];
sx q[0];
rz(-0.44647631) q[0];
sx q[0];
rz(-2.8004942) q[0];
rz(-2.3977631) q[2];
sx q[2];
rz(-2.2093536) q[2];
sx q[2];
rz(-3.0628052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.39602173) q[1];
sx q[1];
rz(-0.77116291) q[1];
sx q[1];
rz(-1.039617) q[1];
rz(-pi) q[2];
rz(1.1105359) q[3];
sx q[3];
rz(-1.2292976) q[3];
sx q[3];
rz(-2.5550492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(-2.2371116) q[2];
rz(-1.0036184) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-0.65902567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7082108) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(1.4138387) q[0];
rz(-2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.7699014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1641952) q[0];
sx q[0];
rz(-1.482555) q[0];
sx q[0];
rz(2.1910358) q[0];
rz(-pi) q[1];
rz(1.781109) q[2];
sx q[2];
rz(-0.82092199) q[2];
sx q[2];
rz(2.6471777) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15118571) q[1];
sx q[1];
rz(-1.8670261) q[1];
sx q[1];
rz(2.5469261) q[1];
rz(-pi) q[2];
rz(-0.38512226) q[3];
sx q[3];
rz(-0.38673863) q[3];
sx q[3];
rz(-0.56798565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(2.7436942) q[2];
rz(2.6509616) q[3];
sx q[3];
rz(-1.5964973) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(0.21815498) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(2.871002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.035707) q[0];
sx q[0];
rz(-0.31292715) q[0];
sx q[0];
rz(0.77625113) q[0];
x q[1];
rz(-1.2241237) q[2];
sx q[2];
rz(-1.7520095) q[2];
sx q[2];
rz(0.62435645) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2968263) q[1];
sx q[1];
rz(-1.4335732) q[1];
sx q[1];
rz(1.5674577) q[1];
rz(2.050428) q[3];
sx q[3];
rz(-1.661146) q[3];
sx q[3];
rz(-2.3144212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.33971912) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.489893) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40437317) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-0.15429601) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-0.74238366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5091632) q[0];
sx q[0];
rz(-2.2262648) q[0];
sx q[0];
rz(-2.3497033) q[0];
rz(-pi) q[1];
rz(-0.97787751) q[2];
sx q[2];
rz(-1.9354543) q[2];
sx q[2];
rz(-2.8713771) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87634516) q[1];
sx q[1];
rz(-2.7956388) q[1];
sx q[1];
rz(-0.032851263) q[1];
rz(-pi) q[2];
rz(1.5417682) q[3];
sx q[3];
rz(-1.4283984) q[3];
sx q[3];
rz(-2.8673429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8538889) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(-1.5819736) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(-0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(-1.3416946) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(-2.3201597) q[2];
sx q[2];
rz(-1.8100304) q[2];
sx q[2];
rz(2.9678154) q[2];
rz(2.3711754) q[3];
sx q[3];
rz(-2.9478248) q[3];
sx q[3];
rz(-0.41268681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
