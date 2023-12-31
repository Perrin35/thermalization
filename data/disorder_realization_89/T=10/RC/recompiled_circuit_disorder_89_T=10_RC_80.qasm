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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3573787) q[0];
sx q[0];
rz(-1.5970018) q[0];
sx q[0];
rz(-1.5062917) q[0];
rz(-pi) q[1];
rz(-2.4550081) q[2];
sx q[2];
rz(-1.540544) q[2];
sx q[2];
rz(0.1498915) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1301073) q[1];
sx q[1];
rz(-1.3890146) q[1];
sx q[1];
rz(-2.6891522) q[1];
x q[2];
rz(2.5476417) q[3];
sx q[3];
rz(-1.4547326) q[3];
sx q[3];
rz(2.5719197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(2.1315234) q[2];
rz(-1.646237) q[3];
sx q[3];
rz(-1.3388747) q[3];
sx q[3];
rz(3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.0533957) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(-2.2825867) q[0];
rz(0.37047085) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.4650311) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5252936) q[0];
sx q[0];
rz(-1.787961) q[0];
sx q[0];
rz(2.8263682) q[0];
rz(2.6066166) q[2];
sx q[2];
rz(-1.8638532) q[2];
sx q[2];
rz(1.5154293) q[2];
rz(-pi) q[3];
x q[3];
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
rz(-1.2831849) q[3];
sx q[3];
rz(-1.0606442) q[3];
sx q[3];
rz(0.19954296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-2.583288) q[2];
rz(1.0393418) q[3];
sx q[3];
rz(-1.9266409) q[3];
sx q[3];
rz(2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52477437) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(2.4257461) q[1];
sx q[1];
rz(-2.2952081) q[1];
sx q[1];
rz(1.3735501) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2620619) q[0];
sx q[0];
rz(-0.52723215) q[0];
sx q[0];
rz(0.049588163) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85767304) q[2];
sx q[2];
rz(-1.8048986) q[2];
sx q[2];
rz(2.8785994) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.19792412) q[1];
sx q[1];
rz(-2.1197882) q[1];
sx q[1];
rz(-0.82358574) q[1];
rz(-2.9186439) q[3];
sx q[3];
rz(-1.3461337) q[3];
sx q[3];
rz(-0.49062452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
sx q[2];
rz(1.019657) q[2];
rz(0.9807469) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(2.5286634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
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
rz(-1.537375) q[1];
sx q[1];
rz(-2.531321) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2903039) q[0];
sx q[0];
rz(-1.5820869) q[0];
sx q[0];
rz(-1.4615061) q[0];
rz(0.82579124) q[2];
sx q[2];
rz(-2.6399887) q[2];
sx q[2];
rz(2.290291) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1000881) q[1];
sx q[1];
rz(-1.8309635) q[1];
sx q[1];
rz(2.2799904) q[1];
rz(-pi) q[2];
rz(-0.48951666) q[3];
sx q[3];
rz(-2.094141) q[3];
sx q[3];
rz(0.19260064) q[3];
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
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6326555) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(1.5326112) q[0];
rz(-2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.3132494) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4603235) q[0];
sx q[0];
rz(-1.8955505) q[0];
sx q[0];
rz(2.2446509) q[0];
rz(-3.1057642) q[2];
sx q[2];
rz(-1.5140859) q[2];
sx q[2];
rz(2.5113475) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.4876731) q[1];
sx q[1];
rz(-2.6012585) q[1];
sx q[1];
rz(-2.3565156) q[1];
rz(2.2458514) q[3];
sx q[3];
rz(-0.81479077) q[3];
sx q[3];
rz(-2.9588587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0017172) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(-0.66429794) q[2];
rz(0.70513606) q[3];
sx q[3];
rz(-0.64283979) q[3];
sx q[3];
rz(3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649081) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(3.0867807) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(-0.48318133) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7532363) q[0];
sx q[0];
rz(-1.5090794) q[0];
sx q[0];
rz(1.7524936) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.032776253) q[2];
sx q[2];
rz(-1.1083229) q[2];
sx q[2];
rz(-1.0810766) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9890081) q[1];
sx q[1];
rz(-0.87859381) q[1];
sx q[1];
rz(2.2425376) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8910965) q[3];
sx q[3];
rz(-0.34892504) q[3];
sx q[3];
rz(0.85897173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(1.9412458) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(-2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(2.7744556) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(2.1869587) q[0];
rz(-0.57149354) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(-0.25209299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4567586) q[0];
sx q[0];
rz(-1.7157468) q[0];
sx q[0];
rz(-0.42379871) q[0];
rz(-pi) q[1];
rz(-0.78105314) q[2];
sx q[2];
rz(-0.99596802) q[2];
sx q[2];
rz(2.1512254) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.5735057) q[1];
sx q[1];
rz(-1.2099669) q[1];
sx q[1];
rz(0.87330694) q[1];
x q[2];
rz(-0.37766405) q[3];
sx q[3];
rz(-1.1389684) q[3];
sx q[3];
rz(0.8197195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7581042) q[2];
sx q[2];
rz(-2.1540116) q[2];
sx q[2];
rz(0.90448109) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-0.86528722) q[3];
sx q[3];
rz(-2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(0.66043234) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(-1.7699014) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28384128) q[0];
sx q[0];
rz(-0.62566602) q[0];
sx q[0];
rz(1.7218504) q[0];
rz(-pi) q[1];
rz(0.22050627) q[2];
sx q[2];
rz(-2.3683511) q[2];
sx q[2];
rz(2.3436433) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.15118571) q[1];
sx q[1];
rz(-1.2745665) q[1];
sx q[1];
rz(0.59466655) q[1];
rz(-pi) q[2];
rz(-2.7564704) q[3];
sx q[3];
rz(-0.38673863) q[3];
sx q[3];
rz(-2.573607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.82683212) q[2];
sx q[2];
rz(-0.62425745) q[2];
sx q[2];
rz(0.39789847) q[2];
rz(0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.9234377) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(2.9072705) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(0.27059069) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.035707) q[0];
sx q[0];
rz(-2.8286655) q[0];
sx q[0];
rz(-2.3653415) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2241237) q[2];
sx q[2];
rz(-1.3895831) q[2];
sx q[2];
rz(0.62435645) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2724243) q[1];
sx q[1];
rz(-0.13726343) q[1];
sx q[1];
rz(-0.024172524) q[1];
x q[2];
rz(1.7646592) q[3];
sx q[3];
rz(-2.6541775) q[3];
sx q[3];
rz(-0.91538115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8018735) q[2];
sx q[2];
rz(-3.0568213) q[2];
sx q[2];
rz(1.489893) q[2];
rz(0.70288944) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(1.9809013) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
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
rz(-1.1318726) q[1];
sx q[1];
rz(-2.399209) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5267747) q[0];
sx q[0];
rz(-2.1702538) q[0];
sx q[0];
rz(-2.4012698) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1701199) q[2];
sx q[2];
rz(-0.68442548) q[2];
sx q[2];
rz(-1.3542092) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.66354499) q[1];
sx q[1];
rz(-1.5596584) q[1];
sx q[1];
rz(2.7958109) q[1];
rz(-pi) q[2];
rz(-1.5417682) q[3];
sx q[3];
rz(-1.7131943) q[3];
sx q[3];
rz(0.27424973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2877038) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5596191) q[2];
rz(-2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(0.5334841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8515274) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(-1.7998981) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-0.32158357) q[2];
sx q[2];
rz(-0.84761534) q[2];
sx q[2];
rz(-1.9615704) q[2];
rz(-1.7066163) q[3];
sx q[3];
rz(-1.4321696) q[3];
sx q[3];
rz(-1.1925478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
