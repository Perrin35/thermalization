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
rz(-2.8566868) q[1];
sx q[1];
rz(-2.6309738) q[1];
sx q[1];
rz(2.7198305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.784214) q[0];
sx q[0];
rz(-1.5970018) q[0];
sx q[0];
rz(1.5062917) q[0];
rz(-pi) q[1];
rz(0.047702567) q[2];
sx q[2];
rz(-0.68714266) q[2];
sx q[2];
rz(1.7575761) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.20323682) q[1];
sx q[1];
rz(-0.48523808) q[1];
sx q[1];
rz(2.7435702) q[1];
rz(-pi) q[2];
rz(-0.59395091) q[3];
sx q[3];
rz(-1.4547326) q[3];
sx q[3];
rz(-0.56967294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.1224147) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-2.1315234) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(3*pi/11) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0533957) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(2.2825867) q[0];
rz(-0.37047085) q[1];
sx q[1];
rz(-1.5444688) q[1];
sx q[1];
rz(-1.6765615) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1169491) q[0];
sx q[0];
rz(-1.878371) q[0];
sx q[0];
rz(1.3427539) q[0];
rz(-2.6071762) q[2];
sx q[2];
rz(-2.538531) q[2];
sx q[2];
rz(-0.39819983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.78039353) q[1];
sx q[1];
rz(-0.91460278) q[1];
sx q[1];
rz(-1.111163) q[1];
rz(0.52821918) q[3];
sx q[3];
rz(-1.8209407) q[3];
sx q[3];
rz(-1.5147097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(1.0393418) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(-0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52477437) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-2.9779789) q[0];
rz(0.71584654) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93689504) q[0];
sx q[0];
rz(-2.097313) q[0];
sx q[0];
rz(-1.599647) q[0];
rz(-1.9203556) q[2];
sx q[2];
rz(-0.74410838) q[2];
sx q[2];
rz(1.0457872) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.8858582) q[1];
sx q[1];
rz(-0.89465562) q[1];
sx q[1];
rz(-0.83791344) q[1];
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
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0107161) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(-1.019657) q[2];
rz(2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2402128) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(-0.83909488) q[0];
rz(1.5377195) q[1];
sx q[1];
rz(-1.6042177) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82204098) q[0];
sx q[0];
rz(-3.0317231) q[0];
sx q[0];
rz(-1.4676453) q[0];
rz(-0.82579124) q[2];
sx q[2];
rz(-0.50160393) q[2];
sx q[2];
rz(-0.85130168) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.1000881) q[1];
sx q[1];
rz(-1.8309635) q[1];
sx q[1];
rz(2.2799904) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2545492) q[3];
sx q[3];
rz(-0.700637) q[3];
sx q[3];
rz(2.5168583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.443976) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(0.237341) q[2];
rz(2.234263) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(1.2835519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(2.2446509) q[0];
x q[1];
rz(3.1057642) q[2];
sx q[2];
rz(-1.5140859) q[2];
sx q[2];
rz(0.63024516) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.65391958) q[1];
sx q[1];
rz(-0.54033414) q[1];
sx q[1];
rz(-2.3565156) q[1];
rz(0.58532183) q[3];
sx q[3];
rz(-0.96671852) q[3];
sx q[3];
rz(-2.4620591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0017172) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(0.66429794) q[2];
rz(-0.70513606) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(-3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
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
rz(2.2830402) q[0];
sx q[0];
rz(-2.9498093) q[0];
sx q[0];
rz(1.2412846) q[0];
x q[1];
rz(-3.1088164) q[2];
sx q[2];
rz(-2.0332697) q[2];
sx q[2];
rz(-1.0810766) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8878057) q[1];
sx q[1];
rz(-1.0711545) q[1];
sx q[1];
rz(-2.3274724) q[1];
rz(-0.11404927) q[3];
sx q[3];
rz(-1.2403135) q[3];
sx q[3];
rz(-1.198311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(-0.085263578) q[2];
rz(-1.9412458) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(0.15032642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7744556) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(-2.1869587) q[0];
rz(2.5700991) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(-0.25209299) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.945767) q[0];
sx q[0];
rz(-2.6951163) q[0];
sx q[0];
rz(-2.8004942) q[0];
x q[1];
rz(-2.3102975) q[2];
sx q[2];
rz(-2.202946) q[2];
sx q[2];
rz(-2.0668154) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.29036197) q[1];
sx q[1];
rz(-0.92612672) q[1];
sx q[1];
rz(-0.45749299) q[1];
rz(-1.1105359) q[3];
sx q[3];
rz(-1.9122951) q[3];
sx q[3];
rz(0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(-0.90448109) q[2];
rz(-1.0036184) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(0.65902567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7082108) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(1.7277539) q[0];
rz(2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(1.7699014) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6721281) q[0];
sx q[0];
rz(-2.1882595) q[0];
sx q[0];
rz(0.10829627) q[0];
rz(-pi) q[1];
rz(1.781109) q[2];
sx q[2];
rz(-0.82092199) q[2];
sx q[2];
rz(-0.49441499) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.129442) q[1];
sx q[1];
rz(-2.4852942) q[1];
sx q[1];
rz(-2.6427569) q[1];
x q[2];
rz(1.4189818) q[3];
sx q[3];
rz(-1.9278798) q[3];
sx q[3];
rz(-0.15541542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21815498) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(2.9072705) q[0];
rz(1.6802457) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(2.871002) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.035707) q[0];
sx q[0];
rz(-2.8286655) q[0];
sx q[0];
rz(-0.77625113) q[0];
rz(1.9174689) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(-0.62435645) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2724243) q[1];
sx q[1];
rz(-3.0043292) q[1];
sx q[1];
rz(3.1174201) q[1];
rz(-pi) q[2];
rz(3.0398265) q[3];
sx q[3];
rz(-1.093285) q[3];
sx q[3];
rz(0.69672841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.33971912) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(1.6516997) q[2];
rz(-0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(-1.9809013) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.568856) q[0];
sx q[0];
rz(-0.15429601) q[0];
rz(-0.94296304) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-2.399209) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61481793) q[0];
sx q[0];
rz(-2.1702538) q[0];
sx q[0];
rz(2.4012698) q[0];
x q[1];
rz(-0.97787751) q[2];
sx q[2];
rz(-1.9354543) q[2];
sx q[2];
rz(0.2702156) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2652475) q[1];
sx q[1];
rz(-0.34595385) q[1];
sx q[1];
rz(0.032851263) q[1];
rz(-pi) q[2];
rz(-0.19974444) q[3];
sx q[3];
rz(-2.9962857) q[3];
sx q[3];
rz(-2.6655281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8538889) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(-1.5596191) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29006526) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(1.7998981) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(0.82143299) q[2];
sx q[2];
rz(-1.8100304) q[2];
sx q[2];
rz(2.9678154) q[2];
rz(-0.7704173) q[3];
sx q[3];
rz(-2.9478248) q[3];
sx q[3];
rz(-0.41268681) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
