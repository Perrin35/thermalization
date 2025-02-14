OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7856287) q[0];
sx q[0];
rz(-0.2956737) q[0];
sx q[0];
rz(0.41391882) q[0];
rz(0.26894459) q[1];
sx q[1];
rz(-0.36965951) q[1];
sx q[1];
rz(-1.2764021) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.39525014) q[0];
sx q[0];
rz(-1.6621886) q[0];
sx q[0];
rz(1.4797828) q[0];
x q[1];
rz(-2.2391255) q[2];
sx q[2];
rz(-1.5215708) q[2];
sx q[2];
rz(-2.7572981) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1943097) q[1];
sx q[1];
rz(-2.3230246) q[1];
sx q[1];
rz(1.1771997) q[1];
rz(-0.26764457) q[3];
sx q[3];
rz(-2.5306314) q[3];
sx q[3];
rz(-2.5087422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74497491) q[2];
sx q[2];
rz(-1.2367542) q[2];
sx q[2];
rz(-1.940894) q[2];
rz(2.4602304) q[3];
sx q[3];
rz(-1.5228289) q[3];
sx q[3];
rz(-0.38667005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.50614) q[0];
sx q[0];
rz(-2.8400087) q[0];
sx q[0];
rz(-2.8801081) q[0];
rz(-1.5226978) q[1];
sx q[1];
rz(-0.50917429) q[1];
sx q[1];
rz(1.2759298) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8573924) q[0];
sx q[0];
rz(-2.5821857) q[0];
sx q[0];
rz(-2.9593234) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4894756) q[2];
sx q[2];
rz(-2.0154833) q[2];
sx q[2];
rz(0.43018489) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.2535227) q[1];
sx q[1];
rz(-1.8948529) q[1];
sx q[1];
rz(0.56323207) q[1];
rz(-0.034364465) q[3];
sx q[3];
rz(-1.4219163) q[3];
sx q[3];
rz(-1.9158165) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.2168938) q[2];
sx q[2];
rz(-1.0328707) q[2];
sx q[2];
rz(0.72803289) q[2];
rz(-1.7714455) q[3];
sx q[3];
rz(-2.5244505) q[3];
sx q[3];
rz(-0.035577687) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.034123357) q[0];
sx q[0];
rz(-2.0296445) q[0];
sx q[0];
rz(2.4566101) q[0];
rz(2.0383539) q[1];
sx q[1];
rz(-2.5220242) q[1];
sx q[1];
rz(-1.0122976) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82498589) q[0];
sx q[0];
rz(-1.1825145) q[0];
sx q[0];
rz(2.7690918) q[0];
rz(-0.84350297) q[2];
sx q[2];
rz(-0.25878227) q[2];
sx q[2];
rz(2.7121787) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.60437459) q[1];
sx q[1];
rz(-0.24362016) q[1];
sx q[1];
rz(-0.55018665) q[1];
x q[2];
rz(-1.0676882) q[3];
sx q[3];
rz(-0.87260428) q[3];
sx q[3];
rz(-0.13867913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.95969069) q[2];
sx q[2];
rz(-1.5896553) q[2];
sx q[2];
rz(-2.336179) q[2];
rz(0.81625932) q[3];
sx q[3];
rz(-0.62556848) q[3];
sx q[3];
rz(3.0474385) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1212946) q[0];
sx q[0];
rz(-0.90885201) q[0];
sx q[0];
rz(-0.62776172) q[0];
rz(0.10920564) q[1];
sx q[1];
rz(-0.86995482) q[1];
sx q[1];
rz(2.3039718) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.254346) q[0];
sx q[0];
rz(-1.5043048) q[0];
sx q[0];
rz(-2.6496135) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97561809) q[2];
sx q[2];
rz(-1.383923) q[2];
sx q[2];
rz(-0.65480168) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7594436) q[1];
sx q[1];
rz(-1.4370222) q[1];
sx q[1];
rz(0.70258883) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.37552278) q[3];
sx q[3];
rz(-1.2831472) q[3];
sx q[3];
rz(-2.875553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8236905) q[2];
sx q[2];
rz(-2.6573942) q[2];
sx q[2];
rz(-2.79706) q[2];
rz(-1.4065929) q[3];
sx q[3];
rz(-2.78648) q[3];
sx q[3];
rz(3.0936892) q[3];
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
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.542881) q[0];
sx q[0];
rz(-0.6364091) q[0];
sx q[0];
rz(2.6357546) q[0];
rz(2.089031) q[1];
sx q[1];
rz(-2.274175) q[1];
sx q[1];
rz(0.94598407) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81954581) q[0];
sx q[0];
rz(-1.844054) q[0];
sx q[0];
rz(1.3083434) q[0];
rz(-pi) q[1];
x q[1];
rz(1.714811) q[2];
sx q[2];
rz(-1.9249232) q[2];
sx q[2];
rz(1.9812552) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.2767503) q[1];
sx q[1];
rz(-2.4227243) q[1];
sx q[1];
rz(-2.5548773) q[1];
x q[2];
rz(1.4220174) q[3];
sx q[3];
rz(-1.5031723) q[3];
sx q[3];
rz(-2.5927582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.08789173) q[2];
sx q[2];
rz(-2.044007) q[2];
sx q[2];
rz(1.5778479) q[2];
rz(1.0846694) q[3];
sx q[3];
rz(-0.8756777) q[3];
sx q[3];
rz(-0.55115551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56187335) q[0];
sx q[0];
rz(-0.60573524) q[0];
sx q[0];
rz(-0.24000034) q[0];
rz(-0.9217681) q[1];
sx q[1];
rz(-1.0488989) q[1];
sx q[1];
rz(-1.9711432) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83607996) q[0];
sx q[0];
rz(-1.5275947) q[0];
sx q[0];
rz(0.57010285) q[0];
rz(-2.7851849) q[2];
sx q[2];
rz(-1.6898167) q[2];
sx q[2];
rz(-1.0072034) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.45615087) q[1];
sx q[1];
rz(-0.69760453) q[1];
sx q[1];
rz(0.82932034) q[1];
rz(-pi) q[2];
rz(-0.097436029) q[3];
sx q[3];
rz(-1.1717005) q[3];
sx q[3];
rz(2.4533085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3971098) q[2];
sx q[2];
rz(-0.8195256) q[2];
sx q[2];
rz(-0.58103713) q[2];
rz(2.2033612) q[3];
sx q[3];
rz(-2.5751028) q[3];
sx q[3];
rz(-2.4375622) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0051603) q[0];
sx q[0];
rz(-2.4649824) q[0];
sx q[0];
rz(0.50049472) q[0];
rz(-0.23807921) q[1];
sx q[1];
rz(-1.6903189) q[1];
sx q[1];
rz(-0.64661017) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80386418) q[0];
sx q[0];
rz(-3.1067305) q[0];
sx q[0];
rz(3.0906223) q[0];
x q[1];
rz(-0.033954577) q[2];
sx q[2];
rz(-1.4563378) q[2];
sx q[2];
rz(-1.3634699) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.889786) q[1];
sx q[1];
rz(-1.4507035) q[1];
sx q[1];
rz(-2.0577862) q[1];
rz(-pi) q[2];
rz(-0.061731652) q[3];
sx q[3];
rz(-1.2128856) q[3];
sx q[3];
rz(-0.45321143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2639192) q[2];
sx q[2];
rz(-0.71920243) q[2];
sx q[2];
rz(-0.70362299) q[2];
rz(1.8536812) q[3];
sx q[3];
rz(-1.0602919) q[3];
sx q[3];
rz(-0.50584403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8877761) q[0];
sx q[0];
rz(-2.0075338) q[0];
sx q[0];
rz(1.8164841) q[0];
rz(-2.3690986) q[1];
sx q[1];
rz(-2.019181) q[1];
sx q[1];
rz(-0.17556369) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85925519) q[0];
sx q[0];
rz(-1.757375) q[0];
sx q[0];
rz(0.76754359) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5042261) q[2];
sx q[2];
rz(-0.91203472) q[2];
sx q[2];
rz(0.8085685) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.404322) q[1];
sx q[1];
rz(-2.0312738) q[1];
sx q[1];
rz(-0.15562017) q[1];
x q[2];
rz(-0.71638338) q[3];
sx q[3];
rz(-1.4759403) q[3];
sx q[3];
rz(1.5945292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.10302155) q[2];
sx q[2];
rz(-1.0744289) q[2];
sx q[2];
rz(-1.0938905) q[2];
rz(1.966018) q[3];
sx q[3];
rz(-0.75785494) q[3];
sx q[3];
rz(0.28483835) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71838251) q[0];
sx q[0];
rz(-3.0960313) q[0];
sx q[0];
rz(-1.8157995) q[0];
rz(2.6153053) q[1];
sx q[1];
rz(-0.86295366) q[1];
sx q[1];
rz(-0.78904271) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51368749) q[0];
sx q[0];
rz(-0.18842489) q[0];
sx q[0];
rz(1.0303251) q[0];
rz(-pi) q[1];
rz(-2.7858235) q[2];
sx q[2];
rz(-1.3311183) q[2];
sx q[2];
rz(-3.0955448) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.583786) q[1];
sx q[1];
rz(-2.4135488) q[1];
sx q[1];
rz(-1.0398204) q[1];
x q[2];
rz(-2.2317722) q[3];
sx q[3];
rz(-2.2299521) q[3];
sx q[3];
rz(2.5499663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9290756) q[2];
sx q[2];
rz(-0.93364659) q[2];
sx q[2];
rz(0.90544256) q[2];
rz(-2.2255911) q[3];
sx q[3];
rz(-1.4429049) q[3];
sx q[3];
rz(1.0441095) q[3];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4395897) q[0];
sx q[0];
rz(-2.3083394) q[0];
sx q[0];
rz(-0.92779094) q[0];
rz(1.0935498) q[1];
sx q[1];
rz(-0.55892006) q[1];
sx q[1];
rz(0.13339001) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0722712) q[0];
sx q[0];
rz(-1.7786553) q[0];
sx q[0];
rz(2.8448807) q[0];
rz(-pi) q[1];
rz(-2.9355644) q[2];
sx q[2];
rz(-1.0090172) q[2];
sx q[2];
rz(0.59625193) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.16006) q[1];
sx q[1];
rz(-1.0602177) q[1];
sx q[1];
rz(-2.4909944) q[1];
rz(-pi) q[2];
rz(-1.978558) q[3];
sx q[3];
rz(-2.6490903) q[3];
sx q[3];
rz(-1.6655079) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6371969) q[2];
sx q[2];
rz(-1.7375676) q[2];
sx q[2];
rz(-1.8782328) q[2];
rz(-2.0718306) q[3];
sx q[3];
rz(-2.1061149) q[3];
sx q[3];
rz(3.0647762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
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
rz(2.1509811) q[0];
sx q[0];
rz(-0.69775109) q[0];
sx q[0];
rz(1.3569111) q[0];
rz(-1.6928584) q[1];
sx q[1];
rz(-2.3285463) q[1];
sx q[1];
rz(2.9978233) q[1];
rz(-2.8303296) q[2];
sx q[2];
rz(-1.0672206) q[2];
sx q[2];
rz(-0.5819566) q[2];
rz(2.1830758) q[3];
sx q[3];
rz(-0.41194852) q[3];
sx q[3];
rz(1.3923399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
