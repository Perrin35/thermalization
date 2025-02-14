OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.27920029) q[0];
sx q[0];
rz(-1.2503799) q[0];
sx q[0];
rz(3.1412636) q[0];
rz(1.1286796) q[1];
sx q[1];
rz(5.1976701) q[1];
sx q[1];
rz(12.357098) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34680172) q[0];
sx q[0];
rz(-0.27898327) q[0];
sx q[0];
rz(-2.1643815) q[0];
x q[1];
rz(-1.7369274) q[2];
sx q[2];
rz(-2.1239933) q[2];
sx q[2];
rz(-0.90546135) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5857201) q[1];
sx q[1];
rz(-1.2798847) q[1];
sx q[1];
rz(-0.32679273) q[1];
rz(-0.23926034) q[3];
sx q[3];
rz(-1.8094899) q[3];
sx q[3];
rz(0.25588122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2613232) q[2];
sx q[2];
rz(-1.7797194) q[2];
sx q[2];
rz(0.23639354) q[2];
rz(-2.7586625) q[3];
sx q[3];
rz(-2.0536486) q[3];
sx q[3];
rz(-1.6898164) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3799389) q[0];
sx q[0];
rz(-2.8528657) q[0];
sx q[0];
rz(-2.5987103) q[0];
rz(0.44034964) q[1];
sx q[1];
rz(-1.1232702) q[1];
sx q[1];
rz(-1.1276668) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8790487) q[0];
sx q[0];
rz(-2.194756) q[0];
sx q[0];
rz(-1.1406374) q[0];
rz(-pi) q[1];
rz(2.2067546) q[2];
sx q[2];
rz(-1.4439772) q[2];
sx q[2];
rz(0.93709125) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4043413) q[1];
sx q[1];
rz(-1.7187498) q[1];
sx q[1];
rz(-0.80929359) q[1];
x q[2];
rz(-1.9502968) q[3];
sx q[3];
rz(-2.2436525) q[3];
sx q[3];
rz(1.6516453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.27309624) q[2];
sx q[2];
rz(-0.79600969) q[2];
sx q[2];
rz(-1.2924755) q[2];
rz(0.83186197) q[3];
sx q[3];
rz(-1.4757194) q[3];
sx q[3];
rz(-0.44124916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.454994) q[0];
sx q[0];
rz(-0.43231493) q[0];
sx q[0];
rz(1.5469714) q[0];
rz(0.0066241344) q[1];
sx q[1];
rz(-0.5358271) q[1];
sx q[1];
rz(-0.63820249) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2062195) q[0];
sx q[0];
rz(-1.6068335) q[0];
sx q[0];
rz(-1.0956148) q[0];
x q[1];
rz(-2.3715921) q[2];
sx q[2];
rz(-1.8714617) q[2];
sx q[2];
rz(-1.9274003) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7591175) q[1];
sx q[1];
rz(-2.6231172) q[1];
sx q[1];
rz(-2.3692959) q[1];
rz(-pi) q[2];
rz(-2.1177111) q[3];
sx q[3];
rz(-1.3860834) q[3];
sx q[3];
rz(0.46177542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.39614761) q[2];
sx q[2];
rz(-3.1148995) q[2];
sx q[2];
rz(0.8566345) q[2];
rz(-2.4062697) q[3];
sx q[3];
rz(-1.7263128) q[3];
sx q[3];
rz(1.0480688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51744866) q[0];
sx q[0];
rz(-0.051001661) q[0];
sx q[0];
rz(2.2080102) q[0];
rz(0.75992641) q[1];
sx q[1];
rz(-2.3257207) q[1];
sx q[1];
rz(1.5504206) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4226993) q[0];
sx q[0];
rz(-0.29445364) q[0];
sx q[0];
rz(-2.8317949) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9924852) q[2];
sx q[2];
rz(-0.44118525) q[2];
sx q[2];
rz(-0.18277482) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67956808) q[1];
sx q[1];
rz(-1.3408519) q[1];
sx q[1];
rz(0.70337062) q[1];
x q[2];
rz(-2.5752034) q[3];
sx q[3];
rz(-2.8479214) q[3];
sx q[3];
rz(-1.3704513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4178702) q[2];
sx q[2];
rz(-1.7524065) q[2];
sx q[2];
rz(-2.1088481) q[2];
rz(2.6241153) q[3];
sx q[3];
rz(-1.7101219) q[3];
sx q[3];
rz(-1.8364068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0378157) q[0];
sx q[0];
rz(-3.0396437) q[0];
sx q[0];
rz(-1.284449) q[0];
rz(-2.2880554) q[1];
sx q[1];
rz(-1.6910911) q[1];
sx q[1];
rz(-2.6654713) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.170661) q[0];
sx q[0];
rz(-1.8722152) q[0];
sx q[0];
rz(-1.0440676) q[0];
rz(1.2849152) q[2];
sx q[2];
rz(-1.21884) q[2];
sx q[2];
rz(2.1536323) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0619701) q[1];
sx q[1];
rz(-1.274612) q[1];
sx q[1];
rz(-2.4825063) q[1];
rz(-pi) q[2];
rz(-0.072612343) q[3];
sx q[3];
rz(-2.3092306) q[3];
sx q[3];
rz(1.2439963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7242754) q[2];
sx q[2];
rz(-0.83160916) q[2];
sx q[2];
rz(1.5403436) q[2];
rz(-0.075751461) q[3];
sx q[3];
rz(-1.3906393) q[3];
sx q[3];
rz(-0.72995228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.041357111) q[0];
sx q[0];
rz(-2.5900109) q[0];
sx q[0];
rz(-1.7919737) q[0];
rz(-1.2760108) q[1];
sx q[1];
rz(-2.3581235) q[1];
sx q[1];
rz(-1.0412201) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2305377) q[0];
sx q[0];
rz(-1.3528224) q[0];
sx q[0];
rz(-0.75454069) q[0];
rz(-1.5346933) q[2];
sx q[2];
rz(-2.3182333) q[2];
sx q[2];
rz(2.2359816) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0372842) q[1];
sx q[1];
rz(-1.3194024) q[1];
sx q[1];
rz(-0.3839375) q[1];
rz(-pi) q[2];
rz(-2.842343) q[3];
sx q[3];
rz(-0.78093442) q[3];
sx q[3];
rz(-2.2921103) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2714046) q[2];
sx q[2];
rz(-1.9726334) q[2];
sx q[2];
rz(0.99883396) q[2];
rz(0.53423229) q[3];
sx q[3];
rz(-1.9637008) q[3];
sx q[3];
rz(1.6629201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6967195) q[0];
sx q[0];
rz(-0.98170009) q[0];
sx q[0];
rz(-1.556742) q[0];
rz(2.059767) q[1];
sx q[1];
rz(-1.6422681) q[1];
sx q[1];
rz(0.36344847) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34521056) q[0];
sx q[0];
rz(-1.78336) q[0];
sx q[0];
rz(1.5770771) q[0];
rz(0.3360904) q[2];
sx q[2];
rz(-2.828183) q[2];
sx q[2];
rz(-1.4614507) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4070426) q[1];
sx q[1];
rz(-0.871823) q[1];
sx q[1];
rz(2.1457325) q[1];
rz(2.0601252) q[3];
sx q[3];
rz(-2.8301468) q[3];
sx q[3];
rz(-1.8925199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.58793679) q[2];
sx q[2];
rz(-2.9755972) q[2];
sx q[2];
rz(2.0641649) q[2];
rz(-1.5466127) q[3];
sx q[3];
rz(-1.3225821) q[3];
sx q[3];
rz(-0.67302978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8107635) q[0];
sx q[0];
rz(-1.6442693) q[0];
sx q[0];
rz(0.28054917) q[0];
rz(0.55889091) q[1];
sx q[1];
rz(-2.2160896) q[1];
sx q[1];
rz(1.9299318) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86211156) q[0];
sx q[0];
rz(-2.0923847) q[0];
sx q[0];
rz(0.83569367) q[0];
rz(-pi) q[1];
rz(1.3468387) q[2];
sx q[2];
rz(-2.0528762) q[2];
sx q[2];
rz(-0.48165762) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7560517) q[1];
sx q[1];
rz(-1.5935399) q[1];
sx q[1];
rz(2.3616516) q[1];
rz(-pi) q[2];
rz(-1.3727851) q[3];
sx q[3];
rz(-1.565838) q[3];
sx q[3];
rz(1.9123506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9676548) q[2];
sx q[2];
rz(-1.8123947) q[2];
sx q[2];
rz(2.3363414) q[2];
rz(-2.0179181) q[3];
sx q[3];
rz(-2.0011963) q[3];
sx q[3];
rz(0.32127389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6174018) q[0];
sx q[0];
rz(-0.62365714) q[0];
sx q[0];
rz(-0.50530857) q[0];
rz(0.71120039) q[1];
sx q[1];
rz(-2.2732747) q[1];
sx q[1];
rz(1.1127081) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8998535) q[0];
sx q[0];
rz(-2.02075) q[0];
sx q[0];
rz(-1.8983272) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64404393) q[2];
sx q[2];
rz(-0.13453787) q[2];
sx q[2];
rz(2.3969039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.45468293) q[1];
sx q[1];
rz(-2.0078604) q[1];
sx q[1];
rz(3.0019762) q[1];
rz(-pi) q[2];
rz(2.1264344) q[3];
sx q[3];
rz(-1.5969513) q[3];
sx q[3];
rz(-1.2873611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.50735551) q[2];
sx q[2];
rz(-1.2398182) q[2];
sx q[2];
rz(0.235454) q[2];
rz(-2.7638451) q[3];
sx q[3];
rz(-0.38769671) q[3];
sx q[3];
rz(3.0915251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78558755) q[0];
sx q[0];
rz(-1.8536785) q[0];
sx q[0];
rz(0.054314286) q[0];
rz(0.76447019) q[1];
sx q[1];
rz(-0.46934325) q[1];
sx q[1];
rz(-2.4321709) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55290907) q[0];
sx q[0];
rz(-1.3257265) q[0];
sx q[0];
rz(0.10677307) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4357752) q[2];
sx q[2];
rz(-1.3121989) q[2];
sx q[2];
rz(1.0279581) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6290789) q[1];
sx q[1];
rz(-0.98328749) q[1];
sx q[1];
rz(1.7085938) q[1];
rz(2.798779) q[3];
sx q[3];
rz(-1.4488932) q[3];
sx q[3];
rz(-1.3212622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28598049) q[2];
sx q[2];
rz(-1.6667655) q[2];
sx q[2];
rz(-0.093185514) q[2];
rz(-0.5992254) q[3];
sx q[3];
rz(-0.56120509) q[3];
sx q[3];
rz(-2.3563201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(1.7080606) q[0];
sx q[0];
rz(-2.1716433) q[0];
sx q[0];
rz(-0.94950983) q[0];
rz(0.064619725) q[1];
sx q[1];
rz(-0.043793543) q[1];
sx q[1];
rz(0.66088062) q[1];
rz(-0.074421616) q[2];
sx q[2];
rz(-1.0205129) q[2];
sx q[2];
rz(-0.27866904) q[2];
rz(0.23537469) q[3];
sx q[3];
rz(-2.4625851) q[3];
sx q[3];
rz(-2.4877297) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
