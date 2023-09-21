OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.75582957) q[0];
sx q[0];
rz(-1.4094149) q[0];
sx q[0];
rz(-0.29456079) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(-0.51061881) q[1];
sx q[1];
rz(-2.7198305) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298676) q[0];
sx q[0];
rz(-1.6352788) q[0];
sx q[0];
rz(-3.1153326) q[0];
rz(-1.6099036) q[2];
sx q[2];
rz(-0.88458672) q[2];
sx q[2];
rz(-1.4456911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1301073) q[1];
sx q[1];
rz(-1.752578) q[1];
sx q[1];
rz(-0.45244042) q[1];
rz(-pi) q[2];
rz(-1.7105605) q[3];
sx q[3];
rz(-2.1602109) q[3];
sx q[3];
rz(-2.0624269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.1224147) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-2.1315234) q[2];
rz(-1.646237) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(8*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.2258376) q[0];
sx q[0];
rz(2.2825867) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.4650311) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0246436) q[0];
sx q[0];
rz(-1.878371) q[0];
sx q[0];
rz(-1.7988388) q[0];
x q[1];
rz(-2.6066166) q[2];
sx q[2];
rz(-1.2777395) q[2];
sx q[2];
rz(1.5154293) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4620004) q[1];
sx q[1];
rz(-0.781171) q[1];
sx q[1];
rz(2.6189234) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.8584077) q[3];
sx q[3];
rz(-2.0809485) q[3];
sx q[3];
rz(-2.9420497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4375962) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-2.583288) q[2];
rz(-2.1022508) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(0.59282747) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
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
rz(2.4257461) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(1.7680426) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64840245) q[0];
sx q[0];
rz(-1.5957386) q[0];
sx q[0];
rz(-2.6148952) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.221237) q[2];
sx q[2];
rz(-0.74410838) q[2];
sx q[2];
rz(2.0958054) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92257567) q[1];
sx q[1];
rz(-2.1891928) q[1];
sx q[1];
rz(2.4465357) q[1];
rz(-2.3397066) q[3];
sx q[3];
rz(-2.826414) q[3];
sx q[3];
rz(-1.8568045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.13087656) q[2];
sx q[2];
rz(-0.54543442) q[2];
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
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2402128) q[0];
sx q[0];
rz(-0.69832435) q[0];
sx q[0];
rz(0.83909488) q[0];
rz(-1.5377195) q[1];
sx q[1];
rz(-1.537375) q[1];
sx q[1];
rz(0.61027169) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3195517) q[0];
sx q[0];
rz(-0.10986957) q[0];
sx q[0];
rz(-1.4676453) q[0];
rz(2.3158014) q[2];
sx q[2];
rz(-2.6399887) q[2];
sx q[2];
rz(-2.290291) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1000881) q[1];
sx q[1];
rz(-1.8309635) q[1];
sx q[1];
rz(-2.2799904) q[1];
rz(-pi) q[2];
rz(-0.9917594) q[3];
sx q[3];
rz(-1.990253) q[3];
sx q[3];
rz(-1.6384244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6976167) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(-0.237341) q[2];
rz(-0.90732968) q[3];
sx q[3];
rz(-1.5483587) q[3];
sx q[3];
rz(-1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6326555) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(1.6089815) q[0];
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
rz(-1.6812692) q[0];
sx q[0];
rz(-1.8955505) q[0];
sx q[0];
rz(-2.2446509) q[0];
x q[1];
rz(1.5140495) q[2];
sx q[2];
rz(-1.5350255) q[2];
sx q[2];
rz(0.93851954) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.65391958) q[1];
sx q[1];
rz(-2.6012585) q[1];
sx q[1];
rz(-0.78507702) q[1];
rz(-pi) q[2];
rz(0.87923268) q[3];
sx q[3];
rz(-2.0428265) q[3];
sx q[3];
rz(-1.8902682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.1398754) q[2];
sx q[2];
rz(-1.8687318) q[2];
sx q[2];
rz(2.4772947) q[2];
rz(0.70513606) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(-3.0378708) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0649081) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(3.0867807) q[0];
rz(-2.1272155) q[1];
sx q[1];
rz(-1.4784112) q[1];
sx q[1];
rz(0.48318133) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38835634) q[0];
sx q[0];
rz(-1.5090794) q[0];
sx q[0];
rz(1.389099) q[0];
rz(-pi) q[1];
rz(1.5051571) q[2];
sx q[2];
rz(-0.46354957) q[2];
sx q[2];
rz(-2.133873) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.15258458) q[1];
sx q[1];
rz(-2.2629988) q[1];
sx q[1];
rz(2.2425376) q[1];
rz(-0.11404927) q[3];
sx q[3];
rz(-1.2403135) q[3];
sx q[3];
rz(1.9432817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9873535) q[2];
sx q[2];
rz(-2.8581212) q[2];
sx q[2];
rz(3.0563291) q[2];
rz(-1.9412458) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(-2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(0.95463395) q[0];
rz(2.5700991) q[1];
sx q[1];
rz(-1.1936455) q[1];
sx q[1];
rz(-0.25209299) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68483401) q[0];
sx q[0];
rz(-1.7157468) q[0];
sx q[0];
rz(-2.7177939) q[0];
rz(-pi) q[1];
rz(2.3977631) q[2];
sx q[2];
rz(-0.93223909) q[2];
sx q[2];
rz(-3.0628052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.568087) q[1];
sx q[1];
rz(-1.9316257) q[1];
sx q[1];
rz(2.2682857) q[1];
rz(-pi) q[2];
rz(-0.37766405) q[3];
sx q[3];
rz(-1.1389684) q[3];
sx q[3];
rz(0.8197195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7581042) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(-2.2371116) q[2];
rz(-2.1379743) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(-0.65902567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
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
rz(-0.43338183) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(1.4138387) q[0];
rz(2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(1.7699014) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1641952) q[0];
sx q[0];
rz(-1.482555) q[0];
sx q[0];
rz(-2.1910358) q[0];
rz(-pi) q[1];
rz(-2.9210864) q[2];
sx q[2];
rz(-2.3683511) q[2];
sx q[2];
rz(2.3436433) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9904069) q[1];
sx q[1];
rz(-1.8670261) q[1];
sx q[1];
rz(-2.5469261) q[1];
rz(-pi) q[2];
rz(-0.36088862) q[3];
sx q[3];
rz(-1.7129671) q[3];
sx q[3];
rz(1.7796381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82683212) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(2.7436942) q[2];
rz(0.49063101) q[3];
sx q[3];
rz(-1.5450954) q[3];
sx q[3];
rz(-1.3354966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(2.9234377) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(-0.23432215) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-2.318858) q[1];
sx q[1];
rz(0.27059069) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9070248) q[0];
sx q[0];
rz(-1.3493291) q[0];
sx q[0];
rz(-1.347876) q[0];
rz(-pi) q[1];
rz(1.0762392) q[2];
sx q[2];
rz(-0.38947546) q[2];
sx q[2];
rz(1.7324093) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.2968263) q[1];
sx q[1];
rz(-1.7080194) q[1];
sx q[1];
rz(1.574135) q[1];
rz(-0.10176615) q[3];
sx q[3];
rz(-1.093285) q[3];
sx q[3];
rz(-2.4448642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.33971912) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.6516997) q[2];
rz(-0.70288944) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(-1.1606914) q[3];
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
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7372195) q[0];
sx q[0];
rz(-1.5727366) q[0];
sx q[0];
rz(-0.15429601) q[0];
rz(-2.1986296) q[1];
sx q[1];
rz(-2.0097201) q[1];
sx q[1];
rz(0.74238366) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5091632) q[0];
sx q[0];
rz(-0.91532781) q[0];
sx q[0];
rz(0.79188939) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7102091) q[2];
sx q[2];
rz(-1.0215534) q[2];
sx q[2];
rz(1.0647578) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.2652475) q[1];
sx q[1];
rz(-2.7956388) q[1];
sx q[1];
rz(0.032851263) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5417682) q[3];
sx q[3];
rz(-1.7131943) q[3];
sx q[3];
rz(0.27424973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2877038) q[2];
sx q[2];
rz(-2.0159857) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(0.21073267) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(-2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(2.8515274) q[0];
sx q[0];
rz(-2.1061438) q[0];
sx q[0];
rz(0.60594546) q[0];
rz(-1.7998981) q[1];
sx q[1];
rz(-1.1732027) q[1];
sx q[1];
rz(-1.8684594) q[1];
rz(-1.9146391) q[2];
sx q[2];
rz(-2.362102) q[2];
sx q[2];
rz(-1.4951928) q[2];
rz(-0.13989862) q[3];
sx q[3];
rz(-1.7053053) q[3];
sx q[3];
rz(-2.7444621) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
