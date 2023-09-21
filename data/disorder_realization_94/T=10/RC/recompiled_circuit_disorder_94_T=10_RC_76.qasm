OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.2774529) q[0];
sx q[0];
rz(-1.5885408) q[0];
sx q[0];
rz(1.5074402) q[0];
rz(1.5965257) q[1];
sx q[1];
rz(-0.59626055) q[1];
sx q[1];
rz(-2.526386) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9278487) q[0];
sx q[0];
rz(-0.97457492) q[0];
sx q[0];
rz(-0.56791373) q[0];
x q[1];
rz(0.46703672) q[2];
sx q[2];
rz(-0.39614284) q[2];
sx q[2];
rz(3.0774088) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.93278904) q[1];
sx q[1];
rz(-0.99901576) q[1];
sx q[1];
rz(-2.6245481) q[1];
rz(-pi) q[2];
rz(1.7513566) q[3];
sx q[3];
rz(-1.2689586) q[3];
sx q[3];
rz(2.3328822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.4404099) q[2];
sx q[2];
rz(-1.5298693) q[2];
sx q[2];
rz(-0.33828503) q[2];
rz(-1.7017378) q[3];
sx q[3];
rz(-2.2262636) q[3];
sx q[3];
rz(-2.2556944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97025362) q[0];
sx q[0];
rz(-2.4304424) q[0];
sx q[0];
rz(-3.1112444) q[0];
rz(-0.066210315) q[1];
sx q[1];
rz(-0.98774424) q[1];
sx q[1];
rz(1.5240086) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.679927) q[0];
sx q[0];
rz(-1.3711509) q[0];
sx q[0];
rz(-0.0017077831) q[0];
x q[1];
rz(1.5288058) q[2];
sx q[2];
rz(-2.6868372) q[2];
sx q[2];
rz(3.0595879) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.46088947) q[1];
sx q[1];
rz(-1.6398805) q[1];
sx q[1];
rz(2.1417888) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0321235) q[3];
sx q[3];
rz(-2.9328049) q[3];
sx q[3];
rz(1.8148282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7559738) q[2];
sx q[2];
rz(-2.1531343) q[2];
sx q[2];
rz(-1.1478109) q[2];
rz(-1.3267481) q[3];
sx q[3];
rz(-1.3245405) q[3];
sx q[3];
rz(-2.9045048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1266992) q[0];
sx q[0];
rz(-2.6601057) q[0];
sx q[0];
rz(0.31578627) q[0];
rz(-0.93859998) q[1];
sx q[1];
rz(-1.4626075) q[1];
sx q[1];
rz(-2.8895203) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.789327) q[0];
sx q[0];
rz(-2.0354712) q[0];
sx q[0];
rz(-0.69791039) q[0];
rz(-pi) q[1];
rz(2.2601068) q[2];
sx q[2];
rz(-0.93886095) q[2];
sx q[2];
rz(1.2607247) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1356126) q[1];
sx q[1];
rz(-0.95893919) q[1];
sx q[1];
rz(2.9115885) q[1];
rz(-0.73022233) q[3];
sx q[3];
rz(-2.632004) q[3];
sx q[3];
rz(2.2024221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0198274) q[2];
sx q[2];
rz(-0.44744197) q[2];
sx q[2];
rz(3.1075409) q[2];
rz(-0.017459067) q[3];
sx q[3];
rz(-1.7826467) q[3];
sx q[3];
rz(1.0954558) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62717342) q[0];
sx q[0];
rz(-1.5486516) q[0];
sx q[0];
rz(1.6148286) q[0];
rz(-1.0871672) q[1];
sx q[1];
rz(-0.68030578) q[1];
sx q[1];
rz(0.70708752) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3084761) q[0];
sx q[0];
rz(-2.0154698) q[0];
sx q[0];
rz(-2.0067257) q[0];
x q[1];
rz(2.884042) q[2];
sx q[2];
rz(-2.6817245) q[2];
sx q[2];
rz(2.3515153) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.1781588) q[1];
sx q[1];
rz(-1.6903094) q[1];
sx q[1];
rz(-0.91576373) q[1];
rz(-pi) q[2];
rz(-1.8791734) q[3];
sx q[3];
rz(-1.3652507) q[3];
sx q[3];
rz(1.8431078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.2924071) q[2];
sx q[2];
rz(-1.1881928) q[2];
sx q[2];
rz(-0.46009955) q[2];
rz(-1.397331) q[3];
sx q[3];
rz(-1.5528691) q[3];
sx q[3];
rz(2.8989255) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3209155) q[0];
sx q[0];
rz(-0.21235947) q[0];
sx q[0];
rz(1.3943577) q[0];
rz(2.0460515) q[1];
sx q[1];
rz(-1.54116) q[1];
sx q[1];
rz(-2.8869693) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4655612) q[0];
sx q[0];
rz(-1.4916294) q[0];
sx q[0];
rz(0.013750793) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.99322929) q[2];
sx q[2];
rz(-1.5830056) q[2];
sx q[2];
rz(2.4075367) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.188272) q[1];
sx q[1];
rz(-0.70536648) q[1];
sx q[1];
rz(-2.7642247) q[1];
rz(-pi) q[2];
rz(2.0650495) q[3];
sx q[3];
rz(-0.77270618) q[3];
sx q[3];
rz(-0.62437526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.022481) q[2];
sx q[2];
rz(-0.20038651) q[2];
sx q[2];
rz(1.3767892) q[2];
rz(-1.6453751) q[3];
sx q[3];
rz(-1.6330556) q[3];
sx q[3];
rz(2.0549324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
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
rz(-2.6901533) q[0];
sx q[0];
rz(-0.7398766) q[0];
sx q[0];
rz(-2.8421463) q[0];
rz(2.1014138) q[1];
sx q[1];
rz(-1.4458011) q[1];
sx q[1];
rz(-0.20656955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8812013) q[0];
sx q[0];
rz(-0.83575373) q[0];
sx q[0];
rz(0.3221237) q[0];
rz(-pi) q[1];
rz(2.0986404) q[2];
sx q[2];
rz(-1.9163418) q[2];
sx q[2];
rz(-2.8258459) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.088591136) q[1];
sx q[1];
rz(-2.2676761) q[1];
sx q[1];
rz(2.9104396) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8387186) q[3];
sx q[3];
rz(-2.2040963) q[3];
sx q[3];
rz(2.5686222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.841659) q[2];
sx q[2];
rz(-1.724023) q[2];
sx q[2];
rz(-2.858813) q[2];
rz(-0.81280604) q[3];
sx q[3];
rz(-2.7225284) q[3];
sx q[3];
rz(-0.16684428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92641002) q[0];
sx q[0];
rz(-1.0535425) q[0];
sx q[0];
rz(2.7600631) q[0];
rz(-0.58386699) q[1];
sx q[1];
rz(-0.54324141) q[1];
sx q[1];
rz(-1.3279703) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8570003) q[0];
sx q[0];
rz(-0.45495957) q[0];
sx q[0];
rz(1.3083463) q[0];
rz(2.6793924) q[2];
sx q[2];
rz(-1.000058) q[2];
sx q[2];
rz(0.58909033) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.92762891) q[1];
sx q[1];
rz(-2.2447526) q[1];
sx q[1];
rz(-0.2143292) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.51289576) q[3];
sx q[3];
rz(-1.9478056) q[3];
sx q[3];
rz(0.60037724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5232089) q[2];
sx q[2];
rz(-1.0655468) q[2];
sx q[2];
rz(1.3605114) q[2];
rz(1.7112188) q[3];
sx q[3];
rz(-2.1332707) q[3];
sx q[3];
rz(0.84806228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1469864) q[0];
sx q[0];
rz(-1.1598347) q[0];
sx q[0];
rz(2.9597136) q[0];
rz(-2.6673642) q[1];
sx q[1];
rz(-2.1209746) q[1];
sx q[1];
rz(-2.1906733) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3527746) q[0];
sx q[0];
rz(-2.7633861) q[0];
sx q[0];
rz(-2.2703855) q[0];
rz(-1.3457001) q[2];
sx q[2];
rz(-2.8689119) q[2];
sx q[2];
rz(2.8358012) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7656895) q[1];
sx q[1];
rz(-1.6357058) q[1];
sx q[1];
rz(2.9325571) q[1];
rz(-1.9740281) q[3];
sx q[3];
rz(-1.2207165) q[3];
sx q[3];
rz(1.10266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.2237079) q[2];
sx q[2];
rz(-2.6612838) q[2];
sx q[2];
rz(-3.0656832) q[2];
rz(0.54801303) q[3];
sx q[3];
rz(-1.3242105) q[3];
sx q[3];
rz(2.5089335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52371812) q[0];
sx q[0];
rz(-2.0848367) q[0];
sx q[0];
rz(-1.7653718) q[0];
rz(0.41704047) q[1];
sx q[1];
rz(-1.4191671) q[1];
sx q[1];
rz(-2.4818647) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2840246) q[0];
sx q[0];
rz(-2.095247) q[0];
sx q[0];
rz(0.88733034) q[0];
rz(2.3599239) q[2];
sx q[2];
rz(-1.4147007) q[2];
sx q[2];
rz(2.8761656) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.33471853) q[1];
sx q[1];
rz(-1.2804619) q[1];
sx q[1];
rz(0.90805407) q[1];
rz(-3.0845853) q[3];
sx q[3];
rz(-1.032864) q[3];
sx q[3];
rz(-1.2008592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.518121) q[2];
sx q[2];
rz(-0.76449624) q[2];
sx q[2];
rz(1.0260322) q[2];
rz(-2.9927411) q[3];
sx q[3];
rz(-1.0326577) q[3];
sx q[3];
rz(0.025645105) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.898107) q[0];
sx q[0];
rz(-2.4004816) q[0];
sx q[0];
rz(1.3056668) q[0];
rz(1.1765515) q[1];
sx q[1];
rz(-1.2780317) q[1];
sx q[1];
rz(2.1059039) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.068533) q[0];
sx q[0];
rz(-1.2464628) q[0];
sx q[0];
rz(-2.7102094) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5752441) q[2];
sx q[2];
rz(-2.1841335) q[2];
sx q[2];
rz(-1.7042421) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.90696883) q[1];
sx q[1];
rz(-1.3303489) q[1];
sx q[1];
rz(-0.30827) q[1];
x q[2];
rz(-2.3922937) q[3];
sx q[3];
rz(-1.7710847) q[3];
sx q[3];
rz(-1.1009969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.5130561) q[2];
sx q[2];
rz(-0.81130242) q[2];
sx q[2];
rz(1.9899842) q[2];
rz(-2.628905) q[3];
sx q[3];
rz(-2.0435464) q[3];
sx q[3];
rz(0.36469665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37968996) q[0];
sx q[0];
rz(-1.7871478) q[0];
sx q[0];
rz(-2.3085069) q[0];
rz(-1.5079386) q[1];
sx q[1];
rz(-2.5588551) q[1];
sx q[1];
rz(-0.48164639) q[1];
rz(0.41082906) q[2];
sx q[2];
rz(-1.6098235) q[2];
sx q[2];
rz(3.1237684) q[2];
rz(-3.1291943) q[3];
sx q[3];
rz(-2.5231902) q[3];
sx q[3];
rz(-1.4004422) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
