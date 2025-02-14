OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.2607245) q[0];
sx q[0];
rz(-0.27624929) q[0];
sx q[0];
rz(-2.2978388) q[0];
rz(1.002797) q[1];
sx q[1];
rz(-2.306814) q[1];
sx q[1];
rz(0.53258449) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40520378) q[0];
sx q[0];
rz(-1.4215934) q[0];
sx q[0];
rz(2.2264541) q[0];
rz(-pi) q[1];
rz(1.7108828) q[2];
sx q[2];
rz(-1.0647961) q[2];
sx q[2];
rz(-1.9078474) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.66825672) q[1];
sx q[1];
rz(-1.7219667) q[1];
sx q[1];
rz(-2.2676047) q[1];
rz(2.769773) q[3];
sx q[3];
rz(-1.9247247) q[3];
sx q[3];
rz(1.33385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3089932) q[2];
sx q[2];
rz(-1.4048046) q[2];
sx q[2];
rz(-3.0187606) q[2];
rz(0.041736688) q[3];
sx q[3];
rz(-1.5580956) q[3];
sx q[3];
rz(-0.48833716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9915344) q[0];
sx q[0];
rz(-0.3287065) q[0];
sx q[0];
rz(-2.026189) q[0];
rz(0.54488048) q[1];
sx q[1];
rz(-1.3976401) q[1];
sx q[1];
rz(-0.56745183) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3504336) q[0];
sx q[0];
rz(-1.6164403) q[0];
sx q[0];
rz(-1.6045536) q[0];
rz(-0.77623244) q[2];
sx q[2];
rz(-1.5481595) q[2];
sx q[2];
rz(2.5843888) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6524446) q[1];
sx q[1];
rz(-2.5622383) q[1];
sx q[1];
rz(-2.9577319) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6593245) q[3];
sx q[3];
rz(-1.8390788) q[3];
sx q[3];
rz(2.0577459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51022092) q[2];
sx q[2];
rz(-3.0192182) q[2];
sx q[2];
rz(-3.0369634) q[2];
rz(2.6369324) q[3];
sx q[3];
rz(-1.3480836) q[3];
sx q[3];
rz(-2.4095355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1402682) q[0];
sx q[0];
rz(-0.19135419) q[0];
sx q[0];
rz(-1.6562756) q[0];
rz(-2.5775919) q[1];
sx q[1];
rz(-1.0537078) q[1];
sx q[1];
rz(2.3174813) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7566273) q[0];
sx q[0];
rz(-0.55776309) q[0];
sx q[0];
rz(1.2560448) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97177695) q[2];
sx q[2];
rz(-2.2702771) q[2];
sx q[2];
rz(0.60764193) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.66131683) q[1];
sx q[1];
rz(-1.3696028) q[1];
sx q[1];
rz(2.7072705) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.070057292) q[3];
sx q[3];
rz(-1.3590711) q[3];
sx q[3];
rz(-2.5222561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.04909) q[2];
sx q[2];
rz(-2.8503032) q[2];
sx q[2];
rz(1.606344) q[2];
rz(-0.90318471) q[3];
sx q[3];
rz(-1.817037) q[3];
sx q[3];
rz(2.7170392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9400738) q[0];
sx q[0];
rz(-2.1551082) q[0];
sx q[0];
rz(-0.41341138) q[0];
rz(0.19551936) q[1];
sx q[1];
rz(-2.1045411) q[1];
sx q[1];
rz(-1.4541385) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91823309) q[0];
sx q[0];
rz(-2.4392895) q[0];
sx q[0];
rz(1.7683517) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6624751) q[2];
sx q[2];
rz(-2.7450941) q[2];
sx q[2];
rz(0.90603196) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.153943) q[1];
sx q[1];
rz(-2.5763303) q[1];
sx q[1];
rz(2.115519) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8995684) q[3];
sx q[3];
rz(-0.70902642) q[3];
sx q[3];
rz(-0.55870648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.89895189) q[2];
sx q[2];
rz(-0.61598888) q[2];
sx q[2];
rz(-1.2658489) q[2];
rz(-0.85593456) q[3];
sx q[3];
rz(-1.7549691) q[3];
sx q[3];
rz(1.9267193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6619381) q[0];
sx q[0];
rz(-2.9582773) q[0];
sx q[0];
rz(1.4187752) q[0];
rz(1.4310369) q[1];
sx q[1];
rz(-1.5242256) q[1];
sx q[1];
rz(-2.4180744) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91416178) q[0];
sx q[0];
rz(-2.1484366) q[0];
sx q[0];
rz(-1.3595057) q[0];
rz(-pi) q[1];
x q[1];
rz(0.75690143) q[2];
sx q[2];
rz(-0.98148966) q[2];
sx q[2];
rz(2.2841931) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0192467) q[1];
sx q[1];
rz(-2.2775688) q[1];
sx q[1];
rz(-3.138423) q[1];
rz(-pi) q[2];
rz(-2.789322) q[3];
sx q[3];
rz(-1.5776424) q[3];
sx q[3];
rz(-2.7628243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.7001223) q[2];
sx q[2];
rz(-2.1258326) q[2];
sx q[2];
rz(2.7777242) q[2];
rz(-1.3077334) q[3];
sx q[3];
rz(-1.4738844) q[3];
sx q[3];
rz(-1.7893121) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.42809197) q[0];
sx q[0];
rz(-2.2245753) q[0];
sx q[0];
rz(-0.8999024) q[0];
rz(-1.046754) q[1];
sx q[1];
rz(-2.7622107) q[1];
sx q[1];
rz(-2.5894763) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13487694) q[0];
sx q[0];
rz(-1.2326041) q[0];
sx q[0];
rz(1.6007559) q[0];
rz(-pi) q[1];
rz(-0.61814536) q[2];
sx q[2];
rz(-0.80354881) q[2];
sx q[2];
rz(-1.8185563) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6450198) q[1];
sx q[1];
rz(-1.7420118) q[1];
sx q[1];
rz(-0.091793072) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5512929) q[3];
sx q[3];
rz(-1.5029782) q[3];
sx q[3];
rz(-2.6026311) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3209352) q[2];
sx q[2];
rz(-2.4652017) q[2];
sx q[2];
rz(-2.9429842) q[2];
rz(2.0123539) q[3];
sx q[3];
rz(-1.6695453) q[3];
sx q[3];
rz(2.3614597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8010913) q[0];
sx q[0];
rz(-1.3022364) q[0];
sx q[0];
rz(-2.4275725) q[0];
rz(1.9188312) q[1];
sx q[1];
rz(-0.87632767) q[1];
sx q[1];
rz(-2.8661935) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3472766) q[0];
sx q[0];
rz(-1.1755318) q[0];
sx q[0];
rz(-2.9121132) q[0];
rz(2.8267639) q[2];
sx q[2];
rz(-0.98199516) q[2];
sx q[2];
rz(-0.77532979) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8620257) q[1];
sx q[1];
rz(-1.0890402) q[1];
sx q[1];
rz(-0.53558366) q[1];
rz(-3.1263859) q[3];
sx q[3];
rz(-2.4162724) q[3];
sx q[3];
rz(-2.050658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.3068646) q[2];
sx q[2];
rz(-1.4707668) q[2];
sx q[2];
rz(-1.4677706) q[2];
rz(-1.5214527) q[3];
sx q[3];
rz(-0.68632564) q[3];
sx q[3];
rz(-2.2036688) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15243212) q[0];
sx q[0];
rz(-2.2280362) q[0];
sx q[0];
rz(-2.1754225) q[0];
rz(-0.78978157) q[1];
sx q[1];
rz(-0.25895324) q[1];
sx q[1];
rz(2.1678179) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.98449303) q[0];
sx q[0];
rz(-2.8678082) q[0];
sx q[0];
rz(-2.0272001) q[0];
rz(-pi) q[1];
rz(-1.9386693) q[2];
sx q[2];
rz(-1.1689567) q[2];
sx q[2];
rz(1.7669058) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5681947) q[1];
sx q[1];
rz(-2.1713082) q[1];
sx q[1];
rz(0.37094122) q[1];
rz(-pi) q[2];
rz(1.6867181) q[3];
sx q[3];
rz(-0.478607) q[3];
sx q[3];
rz(-0.49540621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3070273) q[2];
sx q[2];
rz(-1.4296738) q[2];
sx q[2];
rz(-2.5359421) q[2];
rz(1.0904795) q[3];
sx q[3];
rz(-2.8993789) q[3];
sx q[3];
rz(-3.0428913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.403648) q[0];
sx q[0];
rz(-3.0927959) q[0];
sx q[0];
rz(-2.5700289) q[0];
rz(-0.38772186) q[1];
sx q[1];
rz(-0.78920904) q[1];
sx q[1];
rz(2.2472084) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8799202) q[0];
sx q[0];
rz(-2.399647) q[0];
sx q[0];
rz(-2.3319753) q[0];
rz(1.0036976) q[2];
sx q[2];
rz(-1.6693448) q[2];
sx q[2];
rz(-2.4273317) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.5928985) q[1];
sx q[1];
rz(-1.1709524) q[1];
sx q[1];
rz(0.26229761) q[1];
rz(-pi) q[2];
rz(-2.6329805) q[3];
sx q[3];
rz(-1.5023352) q[3];
sx q[3];
rz(0.48229846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.9992708) q[2];
sx q[2];
rz(-2.2597376) q[2];
sx q[2];
rz(-2.2692661) q[2];
rz(3.0669323) q[3];
sx q[3];
rz(-1.229769) q[3];
sx q[3];
rz(-0.57620302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9999303) q[0];
sx q[0];
rz(-2.7003728) q[0];
sx q[0];
rz(-0.73750752) q[0];
rz(2.4420338) q[1];
sx q[1];
rz(-2.12205) q[1];
sx q[1];
rz(-0.84651822) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2071211) q[0];
sx q[0];
rz(-1.9686598) q[0];
sx q[0];
rz(0.69402309) q[0];
x q[1];
rz(3.0902683) q[2];
sx q[2];
rz(-1.996417) q[2];
sx q[2];
rz(-2.1164621) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.032585) q[1];
sx q[1];
rz(-1.4334588) q[1];
sx q[1];
rz(2.0020383) q[1];
x q[2];
rz(-2.722214) q[3];
sx q[3];
rz(-2.7079795) q[3];
sx q[3];
rz(2.8948642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.49125853) q[2];
sx q[2];
rz(-1.1639872) q[2];
sx q[2];
rz(2.7872046) q[2];
rz(-2.3645511) q[3];
sx q[3];
rz(-0.6207501) q[3];
sx q[3];
rz(-2.0508155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95018321) q[0];
sx q[0];
rz(-1.2744899) q[0];
sx q[0];
rz(0.50344678) q[0];
rz(1.3632111) q[1];
sx q[1];
rz(-2.6111205) q[1];
sx q[1];
rz(3.0725239) q[1];
rz(1.6949881) q[2];
sx q[2];
rz(-2.4249981) q[2];
sx q[2];
rz(-0.80719765) q[2];
rz(-0.34240889) q[3];
sx q[3];
rz(-1.1384321) q[3];
sx q[3];
rz(-2.1428718) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
