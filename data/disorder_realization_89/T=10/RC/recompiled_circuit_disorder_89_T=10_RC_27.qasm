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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5428107) q[0];
sx q[0];
rz(-0.069617696) q[0];
sx q[0];
rz(-1.1845864) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4550081) q[2];
sx q[2];
rz(-1.6010487) q[2];
sx q[2];
rz(-0.1498915) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.9383558) q[1];
sx q[1];
rz(-2.6563546) q[1];
sx q[1];
rz(-0.39802246) q[1];
x q[2];
rz(0.59395091) q[3];
sx q[3];
rz(-1.4547326) q[3];
sx q[3];
rz(0.56967294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(2.1315234) q[2];
rz(1.646237) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0881969) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(0.85900599) q[0];
rz(-2.7711218) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.4650311) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7709748) q[0];
sx q[0];
rz(-2.7608702) q[0];
sx q[0];
rz(-0.61852635) q[0];
x q[1];
rz(2.6071762) q[2];
sx q[2];
rz(-2.538531) q[2];
sx q[2];
rz(-2.7433928) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.49711455) q[1];
sx q[1];
rz(-1.9299572) q[1];
sx q[1];
rz(-2.4317846) q[1];
rz(-0.52821918) q[3];
sx q[3];
rz(-1.3206519) q[3];
sx q[3];
rz(-1.5147097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4375962) q[2];
sx q[2];
rz(-1.2164755) q[2];
sx q[2];
rz(-0.55830467) q[2];
rz(1.0393418) q[3];
sx q[3];
rz(-1.2149518) q[3];
sx q[3];
rz(-2.5487652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
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
rz(1.3735501) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4931902) q[0];
sx q[0];
rz(-1.5957386) q[0];
sx q[0];
rz(2.6148952) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.30544124) q[2];
sx q[2];
rz(-0.88103308) q[2];
sx q[2];
rz(-1.5058215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2557345) q[1];
sx q[1];
rz(-0.89465562) q[1];
sx q[1];
rz(-0.83791344) q[1];
rz(-1.340629) q[3];
sx q[3];
rz(-1.7880511) q[3];
sx q[3];
rz(1.0297071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13087656) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(-2.1219357) q[2];
rz(0.9807469) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(-0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
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
rz(-1.537375) q[1];
sx q[1];
rz(-0.61027169) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82204098) q[0];
sx q[0];
rz(-0.10986957) q[0];
sx q[0];
rz(1.6739474) q[0];
rz(-1.9539815) q[2];
sx q[2];
rz(-1.9028579) q[2];
sx q[2];
rz(0.039475723) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3210541) q[1];
sx q[1];
rz(-2.39403) q[1];
sx q[1];
rz(-1.9588406) q[1];
x q[2];
rz(0.9917594) q[3];
sx q[3];
rz(-1.990253) q[3];
sx q[3];
rz(1.6384244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6326555) q[0];
sx q[0];
rz(-3.0314358) q[0];
sx q[0];
rz(1.6089815) q[0];
rz(2.4081047) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(1.8283432) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8719296) q[0];
sx q[0];
rz(-0.73686826) q[0];
sx q[0];
rz(-2.065573) q[0];
rz(-1.0079185) q[2];
sx q[2];
rz(-3.0745227) q[2];
sx q[2];
rz(1.9474533) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.65391958) q[1];
sx q[1];
rz(-0.54033414) q[1];
sx q[1];
rz(-2.3565156) q[1];
rz(-0.89574121) q[3];
sx q[3];
rz(-2.3268019) q[3];
sx q[3];
rz(-0.18273396) q[3];
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
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649081) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(0.054811906) q[0];
rz(-1.0143771) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(0.48318133) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7532363) q[0];
sx q[0];
rz(-1.5090794) q[0];
sx q[0];
rz(1.7524936) q[0];
rz(-pi) q[1];
rz(3.1088164) q[2];
sx q[2];
rz(-2.0332697) q[2];
sx q[2];
rz(1.0810766) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8878057) q[1];
sx q[1];
rz(-1.0711545) q[1];
sx q[1];
rz(-2.3274724) q[1];
x q[2];
rz(1.2383078) q[3];
sx q[3];
rz(-1.4629435) q[3];
sx q[3];
rz(-0.40963848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.15423916) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(0.085263578) q[2];
rz(-1.2003468) q[3];
sx q[3];
rz(-1.5162568) q[3];
sx q[3];
rz(2.9912662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36713704) q[0];
sx q[0];
rz(-3.0052003) q[0];
sx q[0];
rz(0.95463395) q[0];
rz(-2.5700991) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(-0.25209299) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1958256) q[0];
sx q[0];
rz(-2.6951163) q[0];
sx q[0];
rz(0.3410985) q[0];
rz(-0.78105314) q[2];
sx q[2];
rz(-0.99596802) q[2];
sx q[2];
rz(-0.99036723) q[2];
rz(-pi) q[3];
x q[3];
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
rz(-0.87330694) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2456456) q[3];
sx q[3];
rz(-0.5657256) q[3];
sx q[3];
rz(-1.56324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.38348848) q[2];
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
sx q[1];
rz(-pi) q[2];
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
rz(2.7082108) q[0];
sx q[0];
rz(-0.0032341783) q[0];
sx q[0];
rz(-1.7277539) q[0];
rz(-2.4811603) q[1];
sx q[1];
rz(-1.3884037) q[1];
sx q[1];
rz(1.3716912) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1641952) q[0];
sx q[0];
rz(-1.482555) q[0];
sx q[0];
rz(2.1910358) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.781109) q[2];
sx q[2];
rz(-2.3206707) q[2];
sx q[2];
rz(2.6471777) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5270556) q[1];
sx q[1];
rz(-1.0053047) q[1];
sx q[1];
rz(1.9238227) q[1];
rz(-pi) q[2];
rz(-2.7564704) q[3];
sx q[3];
rz(-0.38673863) q[3];
sx q[3];
rz(-2.573607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82683212) q[2];
sx q[2];
rz(-2.5173352) q[2];
sx q[2];
rz(-0.39789847) q[2];
rz(0.49063101) q[3];
sx q[3];
rz(-1.5964973) q[3];
sx q[3];
rz(-1.8060961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21815498) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(2.9072705) q[0];
rz(-1.461347) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(-2.871002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4280677) q[0];
sx q[0];
rz(-1.3534091) q[0];
sx q[0];
rz(2.9146951) q[0];
x q[1];
rz(-2.9491896) q[2];
sx q[2];
rz(-1.9115598) q[2];
sx q[2];
rz(0.88142384) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8671659) q[1];
sx q[1];
rz(-1.567489) q[1];
sx q[1];
rz(0.13722384) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0911646) q[3];
sx q[3];
rz(-1.4804466) q[3];
sx q[3];
rz(0.8271715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.33971912) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.489893) q[2];
rz(-2.4387032) q[3];
sx q[3];
rz(-1.1542164) q[3];
sx q[3];
rz(1.1606914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.61481793) q[0];
sx q[0];
rz(-2.1702538) q[0];
sx q[0];
rz(0.74032289) q[0];
x q[1];
rz(-2.1701199) q[2];
sx q[2];
rz(-0.68442548) q[2];
sx q[2];
rz(1.3542092) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.87634516) q[1];
sx q[1];
rz(-0.34595385) q[1];
sx q[1];
rz(0.032851263) q[1];
x q[2];
rz(0.19974444) q[3];
sx q[3];
rz(-0.14530694) q[3];
sx q[3];
rz(-2.6655281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8538889) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(-1.5596191) q[2];
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
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.29006526) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(-1.3416946) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(-0.32158357) q[2];
sx q[2];
rz(-0.84761534) q[2];
sx q[2];
rz(-1.9615704) q[2];
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