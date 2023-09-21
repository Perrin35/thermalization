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
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5428107) q[0];
sx q[0];
rz(-0.069617696) q[0];
sx q[0];
rz(1.1845864) q[0];
x q[1];
rz(1.6099036) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(-1.6959015) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.20323682) q[1];
sx q[1];
rz(-0.48523808) q[1];
sx q[1];
rz(-2.7435702) q[1];
rz(-pi) q[2];
rz(1.4310322) q[3];
sx q[3];
rz(-2.1602109) q[3];
sx q[3];
rz(-2.0624269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1224147) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(-1.646237) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(-8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0533957) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(-0.37047085) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(-1.4650311) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61629907) q[0];
sx q[0];
rz(-1.3536317) q[0];
sx q[0];
rz(-0.3152245) q[0];
rz(1.2334521) q[2];
sx q[2];
rz(-1.0609027) q[2];
sx q[2];
rz(-0.22491977) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.78039353) q[1];
sx q[1];
rz(-2.2269899) q[1];
sx q[1];
rz(-2.0304297) q[1];
rz(-1.8584077) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(2.9420497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(-1.0393418) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(-2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6168183) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-0.16361374) q[0];
rz(-2.4257461) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(-1.7680426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.93689504) q[0];
sx q[0];
rz(-1.0442797) q[0];
sx q[0];
rz(-1.599647) q[0];
rz(0.30544124) q[2];
sx q[2];
rz(-2.2605596) q[2];
sx q[2];
rz(-1.5058215) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.19792412) q[1];
sx q[1];
rz(-1.0218044) q[1];
sx q[1];
rz(-2.3180069) q[1];
x q[2];
rz(-0.80188607) q[3];
sx q[3];
rz(-0.31517866) q[3];
sx q[3];
rz(-1.8568045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(-1.019657) q[2];
rz(-0.9807469) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402128) q[0];
sx q[0];
rz(-2.4432683) q[0];
sx q[0];
rz(-0.83909488) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2903039) q[0];
sx q[0];
rz(-1.5820869) q[0];
sx q[0];
rz(1.6800866) q[0];
rz(-pi) q[1];
rz(-0.82579124) q[2];
sx q[2];
rz(-2.6399887) q[2];
sx q[2];
rz(0.85130168) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1000881) q[1];
sx q[1];
rz(-1.3106292) q[1];
sx q[1];
rz(2.2799904) q[1];
x q[2];
rz(2.652076) q[3];
sx q[3];
rz(-1.0474516) q[3];
sx q[3];
rz(2.948992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.443976) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(-2.234263) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(-1.2835519) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5089371) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.5326112) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.8746904) q[1];
sx q[1];
rz(-1.8283432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6812692) q[0];
sx q[0];
rz(-1.8955505) q[0];
sx q[0];
rz(-0.89694174) q[0];
rz(3.1057642) q[2];
sx q[2];
rz(-1.6275068) q[2];
sx q[2];
rz(-0.63024516) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6254127) q[1];
sx q[1];
rz(-1.9429632) q[1];
sx q[1];
rz(-2.7402997) q[1];
rz(-2.26236) q[3];
sx q[3];
rz(-2.0428265) q[3];
sx q[3];
rz(-1.8902682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.0017172) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(0.66429794) q[2];
rz(-2.4364566) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076684549) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(0.054811906) q[0];
rz(1.0143771) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(-0.48318133) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38835634) q[0];
sx q[0];
rz(-1.5090794) q[0];
sx q[0];
rz(-1.7524936) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5051571) q[2];
sx q[2];
rz(-2.6780431) q[2];
sx q[2];
rz(1.0077196) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3997765) q[1];
sx q[1];
rz(-2.2175334) q[1];
sx q[1];
rz(2.497655) q[1];
x q[2];
rz(-1.9032848) q[3];
sx q[3];
rz(-1.4629435) q[3];
sx q[3];
rz(2.7319542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(-3.0563291) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(-2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(0.95463395) q[0];
rz(-0.57149354) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(-0.25209299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3207021) q[0];
sx q[0];
rz(-1.989869) q[0];
sx q[0];
rz(1.4120031) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3977631) q[2];
sx q[2];
rz(-2.2093536) q[2];
sx q[2];
rz(3.0628052) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5735057) q[1];
sx q[1];
rz(-1.2099669) q[1];
sx q[1];
rz(0.87330694) q[1];
rz(-pi) q[2];
rz(2.7639286) q[3];
sx q[3];
rz(-1.1389684) q[3];
sx q[3];
rz(-2.3218732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.38348848) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(2.2371116) q[2];
rz(-2.1379743) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-2.482567) q[3];
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
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7082108) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.3716912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46946457) q[0];
sx q[0];
rz(-0.95333316) q[0];
sx q[0];
rz(-3.0332964) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3604836) q[2];
sx q[2];
rz(-2.3206707) q[2];
sx q[2];
rz(2.6471777) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0121507) q[1];
sx q[1];
rz(-2.4852942) q[1];
sx q[1];
rz(2.6427569) q[1];
rz(-pi) q[2];
rz(-1.4189818) q[3];
sx q[3];
rz(-1.9278798) q[3];
sx q[3];
rz(-2.9861772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3147605) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(-2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21815498) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(0.23432215) q[0];
rz(-1.6802457) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-2.871002) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71352495) q[0];
sx q[0];
rz(-1.3534091) q[0];
sx q[0];
rz(0.2268976) q[0];
rz(1.9174689) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(2.5172362) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.2724243) q[1];
sx q[1];
rz(-0.13726343) q[1];
sx q[1];
rz(3.1174201) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0398265) q[3];
sx q[3];
rz(-2.0483077) q[3];
sx q[3];
rz(-2.4448642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.33971912) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(1.489893) q[2];
rz(2.4387032) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-2.9872966) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(2.399209) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48001227) q[0];
sx q[0];
rz(-2.1614657) q[0];
sx q[0];
rz(2.3175879) q[0];
rz(-pi) q[1];
rz(2.7102091) q[2];
sx q[2];
rz(-2.1200392) q[2];
sx q[2];
rz(2.0768349) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.66354499) q[1];
sx q[1];
rz(-1.5819342) q[1];
sx q[1];
rz(-0.34578172) q[1];
rz(-pi) q[2];
rz(1.5998245) q[3];
sx q[3];
rz(-1.4283984) q[3];
sx q[3];
rz(-0.27424973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5596191) q[2];
rz(-0.21073267) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(-2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8515274) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(1.7998981) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(-2.3201597) q[2];
sx q[2];
rz(-1.8100304) q[2];
sx q[2];
rz(2.9678154) q[2];
rz(-3.001694) q[3];
sx q[3];
rz(-1.4362873) q[3];
sx q[3];
rz(0.39713058) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];