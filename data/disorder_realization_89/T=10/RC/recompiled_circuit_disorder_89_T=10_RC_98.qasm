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
rz(4.8737704) q[0];
sx q[0];
rz(9.1302172) q[0];
rz(0.28490588) q[1];
sx q[1];
rz(2.6309738) q[1];
sx q[1];
rz(9.0030158) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2117251) q[0];
sx q[0];
rz(-1.5063138) q[0];
sx q[0];
rz(-3.1153326) q[0];
rz(-pi) q[1];
rz(-2.4550081) q[2];
sx q[2];
rz(-1.6010487) q[2];
sx q[2];
rz(-0.1498915) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9383558) q[1];
sx q[1];
rz(-2.6563546) q[1];
sx q[1];
rz(-0.39802246) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5476417) q[3];
sx q[3];
rz(-1.4547326) q[3];
sx q[3];
rz(-0.56967294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.0191779) q[2];
sx q[2];
rz(-2.3712967) q[2];
sx q[2];
rz(-2.1315234) q[2];
rz(1.4953556) q[3];
sx q[3];
rz(-1.8027179) q[3];
sx q[3];
rz(-3*pi/11) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(1.0881969) q[0];
sx q[0];
rz(-1.915755) q[0];
sx q[0];
rz(-0.85900599) q[0];
rz(2.7711218) q[1];
sx q[1];
rz(-1.5971239) q[1];
sx q[1];
rz(1.6765615) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0246436) q[0];
sx q[0];
rz(-1.2632217) q[0];
sx q[0];
rz(1.3427539) q[0];
rz(0.53497603) q[2];
sx q[2];
rz(-1.8638532) q[2];
sx q[2];
rz(-1.5154293) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.78039353) q[1];
sx q[1];
rz(-0.91460278) q[1];
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
sx q[1];
rz(-pi/2) q[1];
rz(-1.7039965) q[2];
sx q[2];
rz(-1.9251172) q[2];
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
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6168183) q[0];
sx q[0];
rz(-0.95239788) q[0];
sx q[0];
rz(2.9779789) q[0];
rz(-2.4257461) q[1];
sx q[1];
rz(-0.84638458) q[1];
sx q[1];
rz(-1.7680426) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2620619) q[0];
sx q[0];
rz(-0.52723215) q[0];
sx q[0];
rz(-3.0920045) q[0];
x q[1];
rz(0.30544124) q[2];
sx q[2];
rz(-2.2605596) q[2];
sx q[2];
rz(1.6357712) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8858582) q[1];
sx q[1];
rz(-2.246937) q[1];
sx q[1];
rz(-2.3036792) q[1];
x q[2];
rz(0.80188607) q[3];
sx q[3];
rz(-0.31517866) q[3];
sx q[3];
rz(1.8568045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.0107161) q[2];
sx q[2];
rz(-2.5961582) q[2];
sx q[2];
rz(-1.019657) q[2];
rz(0.9807469) q[3];
sx q[3];
rz(-2.5648983) q[3];
sx q[3];
rz(-0.61292928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.423324) q[0];
sx q[0];
rz(-1.4615131) q[0];
sx q[0];
rz(-0.011358326) q[0];
rz(-pi) q[1];
rz(-0.82579124) q[2];
sx q[2];
rz(-0.50160393) q[2];
sx q[2];
rz(2.290291) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1000881) q[1];
sx q[1];
rz(-1.3106292) q[1];
sx q[1];
rz(0.86160223) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1498333) q[3];
sx q[3];
rz(-1.1513396) q[3];
sx q[3];
rz(-1.5031682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.443976) q[2];
sx q[2];
rz(-2.6306751) q[2];
sx q[2];
rz(0.237341) q[2];
rz(-2.234263) q[3];
sx q[3];
rz(-1.5932339) q[3];
sx q[3];
rz(-1.8580407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5089371) q[0];
sx q[0];
rz(-0.11015686) q[0];
sx q[0];
rz(-1.5326112) q[0];
rz(-0.73348796) q[1];
sx q[1];
rz(-1.2669022) q[1];
sx q[1];
rz(-1.3132494) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26966306) q[0];
sx q[0];
rz(-2.4047244) q[0];
sx q[0];
rz(-1.0760197) q[0];
x q[1];
rz(-1.6275431) q[2];
sx q[2];
rz(-1.6065671) q[2];
sx q[2];
rz(-0.93851954) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5161799) q[1];
sx q[1];
rz(-1.1986294) q[1];
sx q[1];
rz(-2.7402997) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5562708) q[3];
sx q[3];
rz(-0.96671852) q[3];
sx q[3];
rz(2.4620591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0017172) q[2];
sx q[2];
rz(-1.2728609) q[2];
sx q[2];
rz(-2.4772947) q[2];
rz(-2.4364566) q[3];
sx q[3];
rz(-2.4987529) q[3];
sx q[3];
rz(-3.0378708) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0649081) q[0];
sx q[0];
rz(-3.0397471) q[0];
sx q[0];
rz(-3.0867807) q[0];
rz(2.1272155) q[1];
sx q[1];
rz(-1.6631815) q[1];
sx q[1];
rz(-2.6584113) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38835634) q[0];
sx q[0];
rz(-1.5090794) q[0];
sx q[0];
rz(-1.389099) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5051571) q[2];
sx q[2];
rz(-0.46354957) q[2];
sx q[2];
rz(-1.0077196) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.253787) q[1];
sx q[1];
rz(-2.0704381) q[1];
sx q[1];
rz(-0.81412022) q[1];
rz(0.11404927) q[3];
sx q[3];
rz(-1.2403135) q[3];
sx q[3];
rz(-1.9432817) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9873535) q[2];
sx q[2];
rz(-0.2834715) q[2];
sx q[2];
rz(0.085263578) q[2];
rz(1.9412458) q[3];
sx q[3];
rz(-1.6253358) q[3];
sx q[3];
rz(0.15032642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36713704) q[0];
sx q[0];
rz(-0.13639233) q[0];
sx q[0];
rz(-2.1869587) q[0];
rz(0.57149354) q[1];
sx q[1];
rz(-1.9479472) q[1];
sx q[1];
rz(2.8894997) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.945767) q[0];
sx q[0];
rz(-0.44647631) q[0];
sx q[0];
rz(-0.3410985) q[0];
rz(0.78105314) q[2];
sx q[2];
rz(-0.99596802) q[2];
sx q[2];
rz(-2.1512254) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7455709) q[1];
sx q[1];
rz(-2.3704297) q[1];
sx q[1];
rz(1.039617) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1105359) q[3];
sx q[3];
rz(-1.9122951) q[3];
sx q[3];
rz(-0.5865435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.38348848) q[2];
sx q[2];
rz(-0.98758101) q[2];
sx q[2];
rz(0.90448109) q[2];
rz(1.0036184) q[3];
sx q[3];
rz(-2.2763054) q[3];
sx q[3];
rz(2.482567) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.43338183) q[0];
sx q[0];
rz(-3.1383585) q[0];
sx q[0];
rz(1.7277539) q[0];
rz(0.66043234) q[1];
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
rz(-2.8577514) q[0];
sx q[0];
rz(-2.5159266) q[0];
sx q[0];
rz(1.7218504) q[0];
rz(-pi) q[1];
rz(2.9210864) q[2];
sx q[2];
rz(-2.3683511) q[2];
sx q[2];
rz(-2.3436433) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.9904069) q[1];
sx q[1];
rz(-1.2745665) q[1];
sx q[1];
rz(-0.59466655) q[1];
rz(-2.7564704) q[3];
sx q[3];
rz(-0.38673863) q[3];
sx q[3];
rz(0.56798565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3147605) q[2];
sx q[2];
rz(-0.62425745) q[2];
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
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9234377) q[0];
sx q[0];
rz(-1.2175918) q[0];
sx q[0];
rz(-2.9072705) q[0];
rz(1.461347) q[1];
sx q[1];
rz(-0.82273465) q[1];
sx q[1];
rz(2.871002) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9070248) q[0];
sx q[0];
rz(-1.3493291) q[0];
sx q[0];
rz(1.7937167) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0762392) q[2];
sx q[2];
rz(-0.38947546) q[2];
sx q[2];
rz(1.4091834) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8691683) q[1];
sx q[1];
rz(-3.0043292) q[1];
sx q[1];
rz(-3.1174201) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0398265) q[3];
sx q[3];
rz(-2.0483077) q[3];
sx q[3];
rz(-0.69672841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.8018735) q[2];
sx q[2];
rz(-0.084771307) q[2];
sx q[2];
rz(1.489893) q[2];
rz(-2.4387032) q[3];
sx q[3];
rz(-1.9873762) q[3];
sx q[3];
rz(1.9809013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
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
rz(2.399209) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6615804) q[0];
sx q[0];
rz(-2.1614657) q[0];
sx q[0];
rz(-0.82400479) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7102091) q[2];
sx q[2];
rz(-2.1200392) q[2];
sx q[2];
rz(-2.0768349) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4780477) q[1];
sx q[1];
rz(-1.5596584) q[1];
sx q[1];
rz(-2.7958109) q[1];
x q[2];
rz(-2.9991355) q[3];
sx q[3];
rz(-1.5995306) q[3];
sx q[3];
rz(1.2924259) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2877038) q[2];
sx q[2];
rz(-1.1256069) q[2];
sx q[2];
rz(1.5819736) q[2];
rz(2.93086) q[3];
sx q[3];
rz(-2.3570574) q[3];
sx q[3];
rz(2.6081086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8515274) q[0];
sx q[0];
rz(-1.0354488) q[0];
sx q[0];
rz(-2.5356472) q[0];
rz(-1.3416946) q[1];
sx q[1];
rz(-1.9683899) q[1];
sx q[1];
rz(1.2731332) q[1];
rz(0.82143299) q[2];
sx q[2];
rz(-1.8100304) q[2];
sx q[2];
rz(2.9678154) q[2];
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
