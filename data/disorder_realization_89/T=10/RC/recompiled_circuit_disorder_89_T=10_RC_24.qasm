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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2117251) q[0];
sx q[0];
rz(-1.6352788) q[0];
sx q[0];
rz(0.026260016) q[0];
rz(0.68658457) q[2];
sx q[2];
rz(-1.6010487) q[2];
sx q[2];
rz(2.9917012) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.4946343) q[1];
sx q[1];
rz(-1.1263493) q[1];
sx q[1];
rz(-1.7723945) q[1];
rz(-pi) q[2];
rz(0.20538879) q[3];
sx q[3];
rz(-2.5377512) q[3];
sx q[3];
rz(-2.3103034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1224147) q[2];
sx q[2];
rz(-0.77029595) q[2];
sx q[2];
rz(1.0100693) q[2];
rz(-1.646237) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(-8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(-0.37047085) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(-1.4650311) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37061781) q[0];
sx q[0];
rz(-2.7608702) q[0];
sx q[0];
rz(-2.5230663) q[0];
x q[1];
rz(-1.2334521) q[2];
sx q[2];
rz(-2.0806899) q[2];
sx q[2];
rz(-0.22491977) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3611991) q[1];
sx q[1];
rz(-2.2269899) q[1];
sx q[1];
rz(2.0304297) q[1];
x q[2];
rz(2.6133735) q[3];
sx q[3];
rz(-1.8209407) q[3];
sx q[3];
rz(-1.6268829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7039965) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(-0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(2.4257461) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(-1.7680426) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64840245) q[0];
sx q[0];
rz(-1.545854) q[0];
sx q[0];
rz(-2.6148952) q[0];
rz(2.8361514) q[2];
sx q[2];
rz(-2.2605596) q[2];
sx q[2];
rz(-1.6357712) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.92257567) q[1];
sx q[1];
rz(-2.1891928) q[1];
sx q[1];
rz(-2.4465357) q[1];
rz(-pi) q[2];
rz(-0.80188607) q[3];
sx q[3];
rz(-0.31517866) q[3];
sx q[3];
rz(1.2847881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(2.1219357) q[2];
rz(2.1608458) q[3];
sx q[3];
rz(-0.5766944) q[3];
sx q[3];
rz(-0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137988) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(-0.83909488) q[0];
rz(-1.5377195) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71826868) q[0];
sx q[0];
rz(-1.4615131) q[0];
sx q[0];
rz(-3.1302343) q[0];
rz(1.9539815) q[2];
sx q[2];
rz(-1.2387347) q[2];
sx q[2];
rz(-3.1021169) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.8205386) q[1];
sx q[1];
rz(-2.39403) q[1];
sx q[1];
rz(1.9588406) q[1];
rz(2.2545492) q[3];
sx q[3];
rz(-2.4409557) q[3];
sx q[3];
rz(-0.62473434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6976167) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(-0.90732968) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(1.5326112) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(-1.3132494) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7816313) q[0];
sx q[0];
rz(-2.2035723) q[0];
sx q[0];
rz(0.40681337) q[0];
rz(-1.6275431) q[2];
sx q[2];
rz(-1.5350255) q[2];
sx q[2];
rz(-2.2030731) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6254127) q[1];
sx q[1];
rz(-1.1986294) q[1];
sx q[1];
rz(-2.7402997) q[1];
x q[2];
rz(0.89574121) q[3];
sx q[3];
rz(-2.3268019) q[3];
sx q[3];
rz(-2.9588587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1398754) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(-2.4772947) q[2];
rz(0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(0.10372182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649081) q[0];
sx q[0];
rz(-0.10184558) q[0];
sx q[0];
rz(3.0867807) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(-0.48318133) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9478215) q[0];
sx q[0];
rz(-1.7521439) q[0];
sx q[0];
rz(-3.0788455) q[0];
rz(-pi) q[1];
rz(-1.1081084) q[2];
sx q[2];
rz(-1.6001284) q[2];
sx q[2];
rz(2.6372452) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3997765) q[1];
sx q[1];
rz(-2.2175334) q[1];
sx q[1];
rz(-2.497655) q[1];
rz(-pi) q[2];
rz(-3.0275434) q[3];
sx q[3];
rz(-1.9012791) q[3];
sx q[3];
rz(1.9432817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(-0.085263578) q[2];
rz(-1.9412458) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7744556) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(-0.95463395) q[0];
rz(2.5700991) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(0.25209299) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.945767) q[0];
sx q[0];
rz(-2.6951163) q[0];
sx q[0];
rz(-2.8004942) q[0];
x q[1];
rz(0.7438296) q[2];
sx q[2];
rz(-0.93223909) q[2];
sx q[2];
rz(3.0628052) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.568087) q[1];
sx q[1];
rz(-1.2099669) q[1];
sx q[1];
rz(-0.87330694) q[1];
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
rz(-pi/2) q[1];
rz(-2.7581042) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(-0.90448109) q[2];
rz(2.1379743) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-0.65902567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
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
rz(0.43338183) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(-1.4138387) q[0];
rz(-2.4811603) q[1];
sx q[1];
rz(-1.753189) q[1];
sx q[1];
rz(1.7699014) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8577514) q[0];
sx q[0];
rz(-0.62566602) q[0];
sx q[0];
rz(1.7218504) q[0];
x q[1];
rz(1.3604836) q[2];
sx q[2];
rz(-2.3206707) q[2];
sx q[2];
rz(2.6471777) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.129442) q[1];
sx q[1];
rz(-0.65629849) q[1];
sx q[1];
rz(-2.6427569) q[1];
rz(-pi) q[2];
rz(1.7226108) q[3];
sx q[3];
rz(-1.2137128) q[3];
sx q[3];
rz(2.9861772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-2.7436942) q[2];
rz(-0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.8060961) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21815498) q[0];
sx q[0];
rz(-1.9240009) q[0];
sx q[0];
rz(-2.9072705) q[0];
rz(1.6802457) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(-0.27059069) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4280677) q[0];
sx q[0];
rz(-1.3534091) q[0];
sx q[0];
rz(-0.2268976) q[0];
x q[1];
rz(2.9491896) q[2];
sx q[2];
rz(-1.9115598) q[2];
sx q[2];
rz(-0.88142384) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.8671659) q[1];
sx q[1];
rz(-1.5741036) q[1];
sx q[1];
rz(-0.13722384) q[1];
rz(-pi) q[2];
rz(-1.0911646) q[3];
sx q[3];
rz(-1.4804466) q[3];
sx q[3];
rz(-0.8271715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.8018735) q[2];
sx q[2];
rz(-0.084771307) q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40437317) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-2.9872966) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-2.0097201) q[1];
sx q[1];
rz(0.74238366) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48001227) q[0];
sx q[0];
rz(-0.98012692) q[0];
sx q[0];
rz(2.3175879) q[0];
rz(-pi) q[1];
rz(-2.1637151) q[2];
sx q[2];
rz(-1.2061384) q[2];
sx q[2];
rz(-2.8713771) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2652475) q[1];
sx q[1];
rz(-2.7956388) q[1];
sx q[1];
rz(-3.1087414) q[1];
rz(-pi) q[2];
rz(0.14245716) q[3];
sx q[3];
rz(-1.5995306) q[3];
sx q[3];
rz(-1.8491668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(-1.5596191) q[2];
rz(2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
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
rz(-1.7998981) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(0.32158357) q[2];
sx q[2];
rz(-2.2939773) q[2];
sx q[2];
rz(1.1800223) q[2];
rz(1.7066163) q[3];
sx q[3];
rz(-1.7094231) q[3];
sx q[3];
rz(1.9490449) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
