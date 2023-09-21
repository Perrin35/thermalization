OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(-1.7295184) q[0];
rz(0.15481678) q[1];
sx q[1];
rz(-2.545949) q[1];
sx q[1];
rz(1.6593978) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94443653) q[0];
sx q[0];
rz(-1.7261788) q[0];
sx q[0];
rz(2.7128501) q[0];
rz(-pi) q[1];
rz(-1.4380768) q[2];
sx q[2];
rz(-1.8258397) q[2];
sx q[2];
rz(-2.1357352) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.23956242) q[1];
sx q[1];
rz(-0.67443496) q[1];
sx q[1];
rz(-0.85427888) q[1];
rz(-0.10098884) q[3];
sx q[3];
rz(-2.1312993) q[3];
sx q[3];
rz(-1.8141754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(2.1872897) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.8538063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-0.99825478) q[0];
sx q[0];
rz(-1.7049494) q[0];
sx q[0];
rz(0.026219333) q[0];
rz(1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-0.96347934) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.022097691) q[0];
sx q[0];
rz(-0.95675981) q[0];
sx q[0];
rz(-3.1387781) q[0];
x q[1];
rz(-2.8112667) q[2];
sx q[2];
rz(-2.1147554) q[2];
sx q[2];
rz(-0.75737539) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.0986833) q[1];
sx q[1];
rz(-0.60815647) q[1];
sx q[1];
rz(-2.152918) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1727932) q[3];
sx q[3];
rz(-0.41799823) q[3];
sx q[3];
rz(0.19817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6271237) q[2];
sx q[2];
rz(-2.0141979) q[2];
sx q[2];
rz(-3.0070686) q[2];
rz(-0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(-0.9427332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(-2.3441558) q[0];
rz(1.047661) q[1];
sx q[1];
rz(-0.14973775) q[1];
sx q[1];
rz(0.55999666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3968351) q[0];
sx q[0];
rz(-2.7042537) q[0];
sx q[0];
rz(-2.0646981) q[0];
x q[1];
rz(2.838344) q[2];
sx q[2];
rz(-1.5911284) q[2];
sx q[2];
rz(-3.0603527) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.8289889) q[1];
sx q[1];
rz(-1.3723515) q[1];
sx q[1];
rz(0.36638422) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0542164) q[3];
sx q[3];
rz(-1.9994352) q[3];
sx q[3];
rz(2.7408858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3893163) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(0.17253549) q[2];
rz(2.1595188) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-1.0579695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(0.28451434) q[0];
rz(2.8248887) q[1];
sx q[1];
rz(-2.7088294) q[1];
sx q[1];
rz(-1.8428615) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38628681) q[0];
sx q[0];
rz(-2.5884429) q[0];
sx q[0];
rz(-1.132071) q[0];
x q[1];
rz(0.0012871731) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(-0.019891642) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1620996) q[1];
sx q[1];
rz(-0.59826189) q[1];
sx q[1];
rz(-1.23566) q[1];
rz(1.8989765) q[3];
sx q[3];
rz(-1.6771852) q[3];
sx q[3];
rz(-1.002841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.6056885) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(0.38468012) q[2];
rz(0.7540594) q[3];
sx q[3];
rz(-2.0850756) q[3];
sx q[3];
rz(-1.6872905) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6032747) q[0];
sx q[0];
rz(-0.98709995) q[0];
sx q[0];
rz(1.3866562) q[0];
rz(0.23100135) q[1];
sx q[1];
rz(-1.341154) q[1];
sx q[1];
rz(0.2968266) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291527) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(0.57106437) q[0];
rz(-pi) q[1];
rz(-2.9752521) q[2];
sx q[2];
rz(-2.3918249) q[2];
sx q[2];
rz(1.9094085) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2424803) q[1];
sx q[1];
rz(-0.63226262) q[1];
sx q[1];
rz(-2.2805023) q[1];
x q[2];
rz(-2.4279168) q[3];
sx q[3];
rz(-1.6380966) q[3];
sx q[3];
rz(2.0590559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(0.39247593) q[2];
rz(1.1522419) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(2.8241482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.0157938) q[0];
sx q[0];
rz(-1.5690465) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(1.3279351) q[1];
sx q[1];
rz(-1.2633879) q[1];
sx q[1];
rz(-2.5352535) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4088926) q[0];
sx q[0];
rz(-0.048763976) q[0];
sx q[0];
rz(-0.34838895) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.41646429) q[2];
sx q[2];
rz(-2.205924) q[2];
sx q[2];
rz(0.093402775) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.025758) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(-1.238766) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0665928) q[3];
sx q[3];
rz(-1.787775) q[3];
sx q[3];
rz(1.0523588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-1.0417465) q[2];
sx q[2];
rz(-1.139337) q[2];
rz(1.6566488) q[3];
sx q[3];
rz(-1.9610201) q[3];
sx q[3];
rz(3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6095603) q[0];
sx q[0];
rz(-0.74815265) q[0];
sx q[0];
rz(2.6334921) q[0];
rz(1.5787026) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(-2.3513444) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19827851) q[0];
sx q[0];
rz(-2.4662848) q[0];
sx q[0];
rz(2.6576463) q[0];
rz(-1.5590645) q[2];
sx q[2];
rz(-1.7844229) q[2];
sx q[2];
rz(0.27660433) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8040647) q[1];
sx q[1];
rz(-1.84066) q[1];
sx q[1];
rz(-0.40562628) q[1];
rz(-pi) q[2];
x q[2];
rz(0.4831794) q[3];
sx q[3];
rz(-2.4348767) q[3];
sx q[3];
rz(-2.8578575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.053085176) q[2];
sx q[2];
rz(-2.6997456) q[2];
sx q[2];
rz(-1.4132168) q[2];
rz(-1.4767856) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(3.0800381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7168032) q[0];
sx q[0];
rz(-0.030310832) q[0];
sx q[0];
rz(2.0943663) q[0];
rz(2.5324902) q[1];
sx q[1];
rz(-1.7276238) q[1];
sx q[1];
rz(-1.3887127) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6435218) q[0];
sx q[0];
rz(-2.5812491) q[0];
sx q[0];
rz(1.5146921) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65865626) q[2];
sx q[2];
rz(-2.5052862) q[2];
sx q[2];
rz(-0.051740019) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.3988406) q[1];
sx q[1];
rz(-1.5821777) q[1];
sx q[1];
rz(-3.0497754) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96447585) q[3];
sx q[3];
rz(-1.8165605) q[3];
sx q[3];
rz(-0.59734674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(-2.2593373) q[2];
rz(-1.7404209) q[3];
sx q[3];
rz(-1.975235) q[3];
sx q[3];
rz(1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27154487) q[0];
sx q[0];
rz(-2.7247868) q[0];
sx q[0];
rz(-1.4260938) q[0];
rz(0.081461279) q[1];
sx q[1];
rz(-1.1625682) q[1];
sx q[1];
rz(-2.5833599) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8919864) q[0];
sx q[0];
rz(-0.59392801) q[0];
sx q[0];
rz(-2.3114165) q[0];
rz(1.105537) q[2];
sx q[2];
rz(-0.12013398) q[2];
sx q[2];
rz(2.5533954) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.25591125) q[1];
sx q[1];
rz(-0.7464039) q[1];
sx q[1];
rz(-1.6649151) q[1];
x q[2];
rz(-2.8006323) q[3];
sx q[3];
rz(-0.51135671) q[3];
sx q[3];
rz(-2.2380059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8490863) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(0.212184) q[2];
rz(0.21197453) q[3];
sx q[3];
rz(-2.4583355) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6367209) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(1.5378392) q[0];
rz(-2.3161855) q[1];
sx q[1];
rz(-0.67276612) q[1];
sx q[1];
rz(2.6182981) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494203) q[0];
sx q[0];
rz(-1.6329375) q[0];
sx q[0];
rz(0.6092351) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6388355) q[2];
sx q[2];
rz(-2.5685852) q[2];
sx q[2];
rz(-2.9266561) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.0303505) q[1];
sx q[1];
rz(-0.48644629) q[1];
sx q[1];
rz(-0.72279795) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3862615) q[3];
sx q[3];
rz(-0.33077251) q[3];
sx q[3];
rz(1.8473234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(0.26930299) q[2];
rz(2.6473911) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(-1.0860898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8158648) q[0];
sx q[0];
rz(-1.5300735) q[0];
sx q[0];
rz(-1.6515401) q[0];
rz(1.4670463) q[1];
sx q[1];
rz(-2.84927) q[1];
sx q[1];
rz(-1.897859) q[1];
rz(3.1303828) q[2];
sx q[2];
rz(-1.3080454) q[2];
sx q[2];
rz(-1.0946454) q[2];
rz(-2.9071517) q[3];
sx q[3];
rz(-0.81080484) q[3];
sx q[3];
rz(-0.57639359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];