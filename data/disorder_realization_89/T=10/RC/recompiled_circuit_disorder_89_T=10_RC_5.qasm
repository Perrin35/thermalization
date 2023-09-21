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
rz(1.9298676) q[0];
sx q[0];
rz(-1.6352788) q[0];
sx q[0];
rz(3.1153326) q[0];
x q[1];
rz(-1.5316891) q[2];
sx q[2];
rz(-2.2570059) q[2];
sx q[2];
rz(-1.4456911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0114853) q[1];
sx q[1];
rz(-1.3890146) q[1];
sx q[1];
rz(-0.45244042) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7105605) q[3];
sx q[3];
rz(-0.98138175) q[3];
sx q[3];
rz(1.0791657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1224147) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(1.0100693) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(-8*pi/11) q[3];
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
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(0.85900599) q[0];
rz(0.37047085) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(-1.6765615) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252936) q[0];
sx q[0];
rz(-1.3536317) q[0];
sx q[0];
rz(0.3152245) q[0];
x q[1];
rz(0.53441647) q[2];
sx q[2];
rz(-0.60306163) q[2];
sx q[2];
rz(0.39819983) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78039353) q[1];
sx q[1];
rz(-2.2269899) q[1];
sx q[1];
rz(-2.0304297) q[1];
rz(1.2831849) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(-0.19954296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-2.583288) q[2];
rz(-1.0393418) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(-2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(0.16361374) q[0];
rz(-0.71584654) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64840245) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(2.6148952) q[0];
rz(2.2839196) q[2];
sx q[2];
rz(-1.8048986) q[2];
sx q[2];
rz(-0.26299325) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2557345) q[1];
sx q[1];
rz(-2.246937) q[1];
sx q[1];
rz(-2.3036792) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9186439) q[3];
sx q[3];
rz(-1.795459) q[3];
sx q[3];
rz(-2.6509681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402128) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(-2.3024978) q[0];
rz(-1.6038731) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85128879) q[0];
sx q[0];
rz(-1.5595058) q[0];
sx q[0];
rz(-1.4615061) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3158014) q[2];
sx q[2];
rz(-2.6399887) q[2];
sx q[2];
rz(0.85130168) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.8205386) q[1];
sx q[1];
rz(-2.39403) q[1];
sx q[1];
rz(1.9588406) q[1];
rz(2.652076) q[3];
sx q[3];
rz(-1.0474516) q[3];
sx q[3];
rz(2.948992) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.443976) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(0.237341) q[2];
rz(0.90732968) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(-1.6089815) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.8283432) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8719296) q[0];
sx q[0];
rz(-2.4047244) q[0];
sx q[0];
rz(-2.065573) q[0];
rz(1.6275431) q[2];
sx q[2];
rz(-1.5350255) q[2];
sx q[2];
rz(2.2030731) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.207706) q[1];
sx q[1];
rz(-1.1983786) q[1];
sx q[1];
rz(-1.1697342) q[1];
rz(-2.2458514) q[3];
sx q[3];
rz(-0.81479077) q[3];
sx q[3];
rz(2.9588587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.1398754) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(2.4772947) q[2];
rz(2.4364566) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.076684549) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(0.054811906) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(2.6584113) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9478215) q[0];
sx q[0];
rz(-1.7521439) q[0];
sx q[0];
rz(-0.062747196) q[0];
rz(-pi) q[1];
rz(-2.0334843) q[2];
sx q[2];
rz(-1.5414642) q[2];
sx q[2];
rz(2.6372452) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.15258458) q[1];
sx q[1];
rz(-2.2629988) q[1];
sx q[1];
rz(0.89905507) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.11404927) q[3];
sx q[3];
rz(-1.2403135) q[3];
sx q[3];
rz(-1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744556) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(-2.5700991) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(-0.25209299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82089059) q[0];
sx q[0];
rz(-1.989869) q[0];
sx q[0];
rz(1.7295895) q[0];
x q[1];
rz(-0.7438296) q[2];
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
sx q[0];
rz(-pi/2) q[0];
rz(-2.8512307) q[1];
sx q[1];
rz(-0.92612672) q[1];
sx q[1];
rz(0.45749299) q[1];
rz(-2.7639286) q[3];
sx q[3];
rz(-2.0026243) q[3];
sx q[3];
rz(-2.3218732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7581042) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43338183) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(0.66043234) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(1.7699014) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46946457) q[0];
sx q[0];
rz(-2.1882595) q[0];
sx q[0];
rz(-0.10829627) q[0];
x q[1];
rz(-0.76099446) q[2];
sx q[2];
rz(-1.7241663) q[2];
sx q[2];
rz(2.209687) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.614537) q[1];
sx q[1];
rz(-2.136288) q[1];
sx q[1];
rz(-1.9238227) q[1];
rz(-1.7226108) q[3];
sx q[3];
rz(-1.9278798) q[3];
sx q[3];
rz(-0.15541542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(-0.49063101) q[3];
sx q[3];
rz(-1.5964973) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9234377) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-2.9072705) q[0];
rz(1.6802457) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(2.871002) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9070248) q[0];
sx q[0];
rz(-1.7922635) q[0];
sx q[0];
rz(1.347876) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0762392) q[2];
sx q[2];
rz(-2.7521172) q[2];
sx q[2];
rz(-1.4091834) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8447664) q[1];
sx q[1];
rz(-1.7080194) q[1];
sx q[1];
rz(-1.5674577) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0911646) q[3];
sx q[3];
rz(-1.4804466) q[3];
sx q[3];
rz(2.3144212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.33971912) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(-1.489893) q[2];
rz(-0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7372195) q[0];
sx q[0];
rz(-1.568856) q[0];
sx q[0];
rz(-2.9872966) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(2.399209) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5091632) q[0];
sx q[0];
rz(-2.2262648) q[0];
sx q[0];
rz(2.3497033) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97147271) q[2];
sx q[2];
rz(-2.4571672) q[2];
sx q[2];
rz(1.7873834) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.2652475) q[1];
sx q[1];
rz(-0.34595385) q[1];
sx q[1];
rz(-3.1087414) q[1];
x q[2];
rz(2.9991355) q[3];
sx q[3];
rz(-1.5995306) q[3];
sx q[3];
rz(1.8491668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(0.21073267) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29006526) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(1.7998981) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(-2.8200091) q[2];
sx q[2];
rz(-2.2939773) q[2];
sx q[2];
rz(1.1800223) q[2];
rz(-1.4349764) q[3];
sx q[3];
rz(-1.7094231) q[3];
sx q[3];
rz(1.9490449) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
