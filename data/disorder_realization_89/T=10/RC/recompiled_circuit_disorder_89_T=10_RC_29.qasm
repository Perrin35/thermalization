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
rz(0.4217622) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298676) q[0];
sx q[0];
rz(-1.6352788) q[0];
sx q[0];
rz(-0.026260016) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5316891) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(1.4456911) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1301073) q[1];
sx q[1];
rz(-1.3890146) q[1];
sx q[1];
rz(0.45244042) q[1];
x q[2];
rz(-0.59395091) q[3];
sx q[3];
rz(-1.4547326) q[3];
sx q[3];
rz(-0.56967294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(1.0100693) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(2.2825867) q[0];
rz(-0.37047085) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.6765615) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0246436) q[0];
sx q[0];
rz(-1.2632217) q[0];
sx q[0];
rz(1.3427539) q[0];
rz(-pi) q[1];
rz(1.9081406) q[2];
sx q[2];
rz(-2.0806899) q[2];
sx q[2];
rz(-0.22491977) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.78039353) q[1];
sx q[1];
rz(-0.91460278) q[1];
sx q[1];
rz(-1.111163) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2831849) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(-0.19954296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4375962) q[2];
sx q[2];
rz(-1.9251172) q[2];
sx q[2];
rz(0.55830467) q[2];
rz(2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(2.5487652) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-2.1891948) q[0];
sx q[0];
rz(-2.9779789) q[0];
rz(-0.71584654) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(pi/2) q[2];
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
x q[1];
rz(-2.2839196) q[2];
sx q[2];
rz(-1.8048986) q[2];
sx q[2];
rz(-2.8785994) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9436685) q[1];
sx q[1];
rz(-1.0218044) q[1];
sx q[1];
rz(2.3180069) q[1];
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
rz(-2.5961582) q[2];
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
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90137988) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(-0.83909488) q[0];
rz(1.6038731) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85128879) q[0];
sx q[0];
rz(-1.5820869) q[0];
sx q[0];
rz(1.6800866) q[0];
rz(-pi) q[1];
rz(2.7856366) q[2];
sx q[2];
rz(-1.2095371) q[2];
sx q[2];
rz(-1.661983) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0415045) q[1];
sx q[1];
rz(-1.3106292) q[1];
sx q[1];
rz(0.86160223) q[1];
rz(-2.2545492) q[3];
sx q[3];
rz(-0.700637) q[3];
sx q[3];
rz(-0.62473434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.6976167) q[2];
sx q[2];
rz(-0.51091754) q[2];
sx q[2];
rz(2.9042517) q[2];
rz(2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6326555) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(-1.5326112) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.8283432) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7816313) q[0];
sx q[0];
rz(-2.2035723) q[0];
sx q[0];
rz(2.7347793) q[0];
rz(-1.6275431) q[2];
sx q[2];
rz(-1.6065671) q[2];
sx q[2];
rz(-0.93851954) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.65391958) q[1];
sx q[1];
rz(-2.6012585) q[1];
sx q[1];
rz(-2.3565156) q[1];
x q[2];
rz(0.89574121) q[3];
sx q[3];
rz(-0.81479077) q[3];
sx q[3];
rz(2.9588587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1398754) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(-0.66429794) q[2];
rz(2.4364566) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(-3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.076684549) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(-0.054811906) q[0];
rz(2.1272155) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(2.6584113) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2830402) q[0];
sx q[0];
rz(-0.1917834) q[0];
sx q[0];
rz(-1.2412846) q[0];
x q[1];
rz(0.032776253) q[2];
sx q[2];
rz(-1.1083229) q[2];
sx q[2];
rz(1.0810766) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.74181611) q[1];
sx q[1];
rz(-0.92405926) q[1];
sx q[1];
rz(0.6439376) q[1];
x q[2];
rz(1.2383078) q[3];
sx q[3];
rz(-1.4629435) q[3];
sx q[3];
rz(-0.40963848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(-1.9479472) q[1];
sx q[1];
rz(0.25209299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68483401) q[0];
sx q[0];
rz(-1.7157468) q[0];
sx q[0];
rz(-0.42379871) q[0];
rz(-2.3102975) q[2];
sx q[2];
rz(-2.202946) q[2];
sx q[2];
rz(-2.0668154) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.29036197) q[1];
sx q[1];
rz(-0.92612672) q[1];
sx q[1];
rz(-0.45749299) q[1];
x q[2];
rz(2.2456456) q[3];
sx q[3];
rz(-0.5657256) q[3];
sx q[3];
rz(1.5783527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(2.2371116) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-2.482567) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(-0.66043234) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(1.7699014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9773974) q[0];
sx q[0];
rz(-1.6590377) q[0];
sx q[0];
rz(2.1910358) q[0];
rz(-pi) q[1];
rz(-2.9210864) q[2];
sx q[2];
rz(-2.3683511) q[2];
sx q[2];
rz(2.3436433) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9904069) q[1];
sx q[1];
rz(-1.8670261) q[1];
sx q[1];
rz(0.59466655) q[1];
rz(-1.7226108) q[3];
sx q[3];
rz(-1.2137128) q[3];
sx q[3];
rz(-2.9861772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3147605) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-0.39789847) q[2];
rz(-2.6509616) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9234377) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(0.23432215) q[0];
rz(1.6802457) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(0.27059069) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10588564) q[0];
sx q[0];
rz(-2.8286655) q[0];
sx q[0];
rz(2.3653415) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9491896) q[2];
sx q[2];
rz(-1.2300328) q[2];
sx q[2];
rz(-2.2601688) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.27442676) q[1];
sx q[1];
rz(-1.5741036) q[1];
sx q[1];
rz(3.0043688) q[1];
rz(-pi) q[2];
x q[2];
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
rz(1.1606914) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-2.9872966) q[0];
rz(0.94296304) q[1];
sx q[1];
rz(-1.1318726) q[1];
sx q[1];
rz(-0.74238366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61481793) q[0];
sx q[0];
rz(-0.97133884) q[0];
sx q[0];
rz(-2.4012698) q[0];
rz(2.1701199) q[2];
sx q[2];
rz(-2.4571672) q[2];
sx q[2];
rz(-1.7873834) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.91126373) q[1];
sx q[1];
rz(-1.9165557) q[1];
sx q[1];
rz(1.5826349) q[1];
x q[2];
rz(-0.14245716) q[3];
sx q[3];
rz(-1.5995306) q[3];
sx q[3];
rz(-1.2924259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2877038) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5596191) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-0.78453523) q[3];
sx q[3];
rz(2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
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
rz(-0.32158357) q[2];
sx q[2];
rz(-0.84761534) q[2];
sx q[2];
rz(-1.9615704) q[2];
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
