OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.55968094) q[0];
sx q[0];
rz(-1.0868602) q[0];
sx q[0];
rz(-1.342919) q[0];
rz(1.5496594) q[1];
sx q[1];
rz(-0.078443371) q[1];
sx q[1];
rz(2.4989541) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4896048) q[0];
sx q[0];
rz(-1.542924) q[0];
sx q[0];
rz(-2.6023988) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50813714) q[2];
sx q[2];
rz(-1.4027486) q[2];
sx q[2];
rz(0.44067581) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3712284) q[1];
sx q[1];
rz(-2.9712354) q[1];
sx q[1];
rz(0.78070663) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5343127) q[3];
sx q[3];
rz(-1.3368703) q[3];
sx q[3];
rz(2.8521188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47444433) q[2];
sx q[2];
rz(-2.4953304) q[2];
sx q[2];
rz(-1.2791963) q[2];
rz(2.4228418) q[3];
sx q[3];
rz(-1.5703392) q[3];
sx q[3];
rz(0.32354245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46368018) q[0];
sx q[0];
rz(-1.7872515) q[0];
sx q[0];
rz(-2.1610778) q[0];
rz(-2.9837043) q[1];
sx q[1];
rz(-1.6442559) q[1];
sx q[1];
rz(-2.352879) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6533587) q[0];
sx q[0];
rz(-1.6846859) q[0];
sx q[0];
rz(-1.6949523) q[0];
rz(-pi) q[1];
rz(-1.5390225) q[2];
sx q[2];
rz(-1.8268158) q[2];
sx q[2];
rz(0.94871828) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.84856725) q[1];
sx q[1];
rz(-1.8752521) q[1];
sx q[1];
rz(-2.5679563) q[1];
rz(-2.4103568) q[3];
sx q[3];
rz(-0.9160708) q[3];
sx q[3];
rz(0.54192858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.7754037) q[2];
sx q[2];
rz(-1.0752233) q[2];
sx q[2];
rz(0.186084) q[2];
rz(-2.4880593) q[3];
sx q[3];
rz(-1.7553522) q[3];
sx q[3];
rz(0.11793605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2114975) q[0];
sx q[0];
rz(-2.1989172) q[0];
sx q[0];
rz(0.63013664) q[0];
rz(-3.0139626) q[1];
sx q[1];
rz(-2.4928513) q[1];
sx q[1];
rz(-2.4198467) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9741309) q[0];
sx q[0];
rz(-0.92739096) q[0];
sx q[0];
rz(1.2533623) q[0];
rz(-0.69408728) q[2];
sx q[2];
rz(-1.9352479) q[2];
sx q[2];
rz(1.8592161) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5210515) q[1];
sx q[1];
rz(-0.44577814) q[1];
sx q[1];
rz(-3.0522703) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.33629041) q[3];
sx q[3];
rz(-0.97067562) q[3];
sx q[3];
rz(2.1093413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7599941) q[2];
sx q[2];
rz(-2.1888032) q[2];
sx q[2];
rz(-2.2325113) q[2];
rz(0.51820731) q[3];
sx q[3];
rz(-1.0639233) q[3];
sx q[3];
rz(-2.853493) q[3];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5575314) q[0];
sx q[0];
rz(-0.73636213) q[0];
sx q[0];
rz(3.0990565) q[0];
rz(2.361239) q[1];
sx q[1];
rz(-2.6413554) q[1];
sx q[1];
rz(-2.3775878) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4064179) q[0];
sx q[0];
rz(-1.58332) q[0];
sx q[0];
rz(2.2494715) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4970384) q[2];
sx q[2];
rz(-1.6719712) q[2];
sx q[2];
rz(-2.8855756) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.9909775) q[1];
sx q[1];
rz(-2.5533278) q[1];
sx q[1];
rz(0.88612835) q[1];
rz(0.66391151) q[3];
sx q[3];
rz(-1.1949364) q[3];
sx q[3];
rz(-1.0192878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.8445231) q[2];
sx q[2];
rz(-1.917118) q[2];
sx q[2];
rz(-2.6085473) q[2];
rz(-2.879203) q[3];
sx q[3];
rz(-2.5049987) q[3];
sx q[3];
rz(-2.6791402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2247291) q[0];
sx q[0];
rz(-0.30039766) q[0];
sx q[0];
rz(1.2438783) q[0];
rz(0.22661701) q[1];
sx q[1];
rz(-2.367327) q[1];
sx q[1];
rz(0.40333834) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.078743155) q[0];
sx q[0];
rz(-1.3855977) q[0];
sx q[0];
rz(1.7395822) q[0];
x q[1];
rz(1.2071768) q[2];
sx q[2];
rz(-0.91797963) q[2];
sx q[2];
rz(1.2448685) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.4407318) q[1];
sx q[1];
rz(-1.7248389) q[1];
sx q[1];
rz(0.53149077) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5087318) q[3];
sx q[3];
rz(-0.78213464) q[3];
sx q[3];
rz(-0.49304214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8205745) q[2];
sx q[2];
rz(-1.0514739) q[2];
sx q[2];
rz(0.78249758) q[2];
rz(2.0292422) q[3];
sx q[3];
rz(-0.4370884) q[3];
sx q[3];
rz(-2.1330244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7111506) q[0];
sx q[0];
rz(-0.40200457) q[0];
sx q[0];
rz(-0.45561403) q[0];
rz(0.09952155) q[1];
sx q[1];
rz(-2.0030622) q[1];
sx q[1];
rz(3.1351556) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2086439) q[0];
sx q[0];
rz(-0.83983487) q[0];
sx q[0];
rz(-1.2293925) q[0];
x q[1];
rz(-0.10613425) q[2];
sx q[2];
rz(-1.1960293) q[2];
sx q[2];
rz(2.8449164) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9812614) q[1];
sx q[1];
rz(-2.6604974) q[1];
sx q[1];
rz(1.0465924) q[1];
x q[2];
rz(0.75140679) q[3];
sx q[3];
rz(-2.2304428) q[3];
sx q[3];
rz(-2.4305962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7727938) q[2];
sx q[2];
rz(-1.6364748) q[2];
sx q[2];
rz(2.2951365) q[2];
rz(-2.1438697) q[3];
sx q[3];
rz(-0.80934757) q[3];
sx q[3];
rz(-1.458582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0434175) q[0];
sx q[0];
rz(-0.88413969) q[0];
sx q[0];
rz(0.055710677) q[0];
rz(-2.3588691) q[1];
sx q[1];
rz(-1.1306154) q[1];
sx q[1];
rz(1.1605211) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8502064) q[0];
sx q[0];
rz(-2.0178447) q[0];
sx q[0];
rz(-0.38711754) q[0];
x q[1];
rz(-0.47126982) q[2];
sx q[2];
rz(-2.0435213) q[2];
sx q[2];
rz(1.8007235) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.9142368) q[1];
sx q[1];
rz(-2.1319234) q[1];
sx q[1];
rz(2.6941872) q[1];
x q[2];
rz(-2.3485687) q[3];
sx q[3];
rz(-2.6245955) q[3];
sx q[3];
rz(-1.6263863) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.065585) q[2];
sx q[2];
rz(-2.2202754) q[2];
sx q[2];
rz(0.78545061) q[2];
rz(0.75585946) q[3];
sx q[3];
rz(-1.2952341) q[3];
sx q[3];
rz(2.7261962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3128368) q[0];
sx q[0];
rz(-0.090933196) q[0];
sx q[0];
rz(1.998741) q[0];
rz(-1.8354592) q[1];
sx q[1];
rz(-1.9530692) q[1];
sx q[1];
rz(2.725504) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.04836719) q[0];
sx q[0];
rz(-2.0443633) q[0];
sx q[0];
rz(0.22305365) q[0];
x q[1];
rz(1.847319) q[2];
sx q[2];
rz(-2.6575436) q[2];
sx q[2];
rz(0.36177847) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5059698) q[1];
sx q[1];
rz(-1.1715679) q[1];
sx q[1];
rz(-1.7446306) q[1];
x q[2];
rz(1.0981406) q[3];
sx q[3];
rz(-1.4720159) q[3];
sx q[3];
rz(-3.0702101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1064421) q[2];
sx q[2];
rz(-1.687655) q[2];
sx q[2];
rz(-2.8184334) q[2];
rz(-0.20251003) q[3];
sx q[3];
rz(-1.860362) q[3];
sx q[3];
rz(2.5642853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37339661) q[0];
sx q[0];
rz(-0.62921262) q[0];
sx q[0];
rz(-1.1449822) q[0];
rz(-1.1960944) q[1];
sx q[1];
rz(-2.9856666) q[1];
sx q[1];
rz(-2.6224565) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11408344) q[0];
sx q[0];
rz(-0.4013831) q[0];
sx q[0];
rz(-1.2252349) q[0];
rz(2.9197308) q[2];
sx q[2];
rz(-1.8010745) q[2];
sx q[2];
rz(-0.78267539) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.12431006) q[1];
sx q[1];
rz(-2.3553684) q[1];
sx q[1];
rz(0.12676523) q[1];
rz(-0.2378283) q[3];
sx q[3];
rz(-1.3999709) q[3];
sx q[3];
rz(-1.2549787) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8273932) q[2];
sx q[2];
rz(-1.2131571) q[2];
sx q[2];
rz(-0.040977565) q[2];
rz(-0.86769062) q[3];
sx q[3];
rz(-2.4882081) q[3];
sx q[3];
rz(-2.6303671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39559078) q[0];
sx q[0];
rz(-0.87930167) q[0];
sx q[0];
rz(1.5270773) q[0];
rz(1.4279667) q[1];
sx q[1];
rz(-2.7797996) q[1];
sx q[1];
rz(1.6428927) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7258543) q[0];
sx q[0];
rz(-0.07809528) q[0];
sx q[0];
rz(-1.8321091) q[0];
rz(-pi) q[1];
x q[1];
rz(0.047949009) q[2];
sx q[2];
rz(-1.7257294) q[2];
sx q[2];
rz(2.5028554) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.0000671) q[1];
sx q[1];
rz(-2.2728517) q[1];
sx q[1];
rz(2.6932004) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4478217) q[3];
sx q[3];
rz(-1.3696559) q[3];
sx q[3];
rz(-0.5826544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3830118) q[2];
sx q[2];
rz(-1.0772394) q[2];
sx q[2];
rz(0.14653462) q[2];
rz(2.3274029) q[3];
sx q[3];
rz(-0.79629961) q[3];
sx q[3];
rz(-0.41671419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9183337) q[0];
sx q[0];
rz(-1.8948566) q[0];
sx q[0];
rz(-2.1444453) q[0];
rz(-1.8756443) q[1];
sx q[1];
rz(-2.1389778) q[1];
sx q[1];
rz(-1.9139342) q[1];
rz(-3.0284991) q[2];
sx q[2];
rz(-1.332915) q[2];
sx q[2];
rz(-1.1168196) q[2];
rz(-0.90445789) q[3];
sx q[3];
rz(-0.96257985) q[3];
sx q[3];
rz(-1.0812159) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];