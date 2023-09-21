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
rz(0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298676) q[0];
sx q[0];
rz(-1.5063138) q[0];
sx q[0];
rz(-3.1153326) q[0];
rz(-0.68658457) q[2];
sx q[2];
rz(-1.6010487) q[2];
sx q[2];
rz(0.1498915) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4946343) q[1];
sx q[1];
rz(-2.0152433) q[1];
sx q[1];
rz(-1.3691982) q[1];
rz(-pi) q[2];
x q[2];
rz(0.59395091) q[3];
sx q[3];
rz(-1.4547326) q[3];
sx q[3];
rz(-2.5719197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1224147) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(-1.0100693) q[2];
rz(-1.4953556) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(-3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0533957) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(0.85900599) q[0];
rz(-0.37047085) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(-1.6765615) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252936) q[0];
sx q[0];
rz(-1.787961) q[0];
sx q[0];
rz(2.8263682) q[0];
rz(-pi) q[1];
rz(-0.53441647) q[2];
sx q[2];
rz(-0.60306163) q[2];
sx q[2];
rz(2.7433928) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3611991) q[1];
sx q[1];
rz(-2.2269899) q[1];
sx q[1];
rz(1.111163) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2831849) q[3];
sx q[3];
rz(-2.0809485) q[3];
sx q[3];
rz(2.9420497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.7039965) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52477437) q[0];
sx q[0];
rz(-0.95239788) q[0];
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
rz(-0.87953075) q[0];
sx q[0];
rz(-2.6143605) q[0];
sx q[0];
rz(-0.049588163) q[0];
rz(1.221237) q[2];
sx q[2];
rz(-0.74410838) q[2];
sx q[2];
rz(-2.0958054) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9436685) q[1];
sx q[1];
rz(-2.1197882) q[1];
sx q[1];
rz(2.3180069) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9186439) q[3];
sx q[3];
rz(-1.795459) q[3];
sx q[3];
rz(2.6509681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0107161) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(-1.019657) q[2];
rz(-2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90137988) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(0.83909488) q[0];
rz(1.5377195) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.423324) q[0];
sx q[0];
rz(-1.6800796) q[0];
sx q[0];
rz(0.011358326) q[0];
x q[1];
rz(2.3158014) q[2];
sx q[2];
rz(-2.6399887) q[2];
sx q[2];
rz(-2.290291) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.8295633) q[1];
sx q[1];
rz(-2.2514572) q[1];
sx q[1];
rz(2.8042253) q[1];
rz(-2.1498333) q[3];
sx q[3];
rz(-1.990253) q[3];
sx q[3];
rz(-1.5031682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-2.9042517) q[2];
rz(0.90732968) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(-1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6326555) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(-1.5326112) q[0];
rz(-0.73348796) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.8283432) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6812692) q[0];
sx q[0];
rz(-1.8955505) q[0];
sx q[0];
rz(-0.89694174) q[0];
rz(-1.6275431) q[2];
sx q[2];
rz(-1.6065671) q[2];
sx q[2];
rz(2.2030731) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.9338867) q[1];
sx q[1];
rz(-1.943214) q[1];
sx q[1];
rz(1.1697342) q[1];
rz(-pi) q[2];
rz(0.89574121) q[3];
sx q[3];
rz(-2.3268019) q[3];
sx q[3];
rz(-2.9588587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0017172) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(0.66429794) q[2];
rz(0.70513606) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649081) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(-3.0867807) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(-2.6584113) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85855243) q[0];
sx q[0];
rz(-0.1917834) q[0];
sx q[0];
rz(1.900308) q[0];
rz(-1.6364355) q[2];
sx q[2];
rz(-0.46354957) q[2];
sx q[2];
rz(-2.133873) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.74181611) q[1];
sx q[1];
rz(-2.2175334) q[1];
sx q[1];
rz(2.497655) q[1];
rz(0.11404927) q[3];
sx q[3];
rz(-1.2403135) q[3];
sx q[3];
rz(1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(-3.0563291) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(0.15032642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744556) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(-2.1869587) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(2.8894997) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4567586) q[0];
sx q[0];
rz(-1.7157468) q[0];
sx q[0];
rz(-2.7177939) q[0];
x q[1];
rz(2.3102975) q[2];
sx q[2];
rz(-2.202946) q[2];
sx q[2];
rz(-1.0747773) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.29036197) q[1];
sx q[1];
rz(-2.2154659) q[1];
sx q[1];
rz(-0.45749299) q[1];
rz(-pi) q[2];
rz(1.1105359) q[3];
sx q[3];
rz(-1.2292976) q[3];
sx q[3];
rz(0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7581042) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(-2.2371116) q[2];
rz(2.1379743) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7082108) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(0.66043234) q[1];
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
rz(2.6721281) q[0];
sx q[0];
rz(-2.1882595) q[0];
sx q[0];
rz(0.10829627) q[0];
rz(-2.3805982) q[2];
sx q[2];
rz(-1.7241663) q[2];
sx q[2];
rz(0.93190565) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.15118571) q[1];
sx q[1];
rz(-1.2745665) q[1];
sx q[1];
rz(0.59466655) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.780704) q[3];
sx q[3];
rz(-1.4286255) q[3];
sx q[3];
rz(1.7796381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(2.7436942) q[2];
rz(-2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21815498) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(-1.6802457) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-2.871002) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71352495) q[0];
sx q[0];
rz(-1.3534091) q[0];
sx q[0];
rz(2.9146951) q[0];
x q[1];
rz(0.19240304) q[2];
sx q[2];
rz(-1.2300328) q[2];
sx q[2];
rz(-0.88142384) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2724243) q[1];
sx q[1];
rz(-3.0043292) q[1];
sx q[1];
rz(0.024172524) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.050428) q[3];
sx q[3];
rz(-1.661146) q[3];
sx q[3];
rz(2.3144212) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8018735) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(-1.489893) q[2];
rz(-0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(0.15429601) q[0];
rz(2.1986296) q[1];
sx q[1];
rz(-2.0097201) q[1];
sx q[1];
rz(-0.74238366) q[1];
rz(pi/2) q[2];
sx q[2];
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
rz(0.43138357) q[2];
sx q[2];
rz(-1.0215534) q[2];
sx q[2];
rz(2.0768349) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4780477) q[1];
sx q[1];
rz(-1.5596584) q[1];
sx q[1];
rz(2.7958109) q[1];
rz(-pi) q[2];
rz(-0.14245716) q[3];
sx q[3];
rz(-1.5420621) q[3];
sx q[3];
rz(1.2924259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5819736) q[2];
rz(2.93086) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(0.5334841) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(1.7998981) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(2.3201597) q[2];
sx q[2];
rz(-1.3315622) q[2];
sx q[2];
rz(-0.17377725) q[2];
rz(3.001694) q[3];
sx q[3];
rz(-1.7053053) q[3];
sx q[3];
rz(-2.7444621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];