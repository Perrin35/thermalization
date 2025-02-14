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
rz(-0.4001652) q[0];
sx q[0];
rz(-0.53181177) q[0];
sx q[0];
rz(-0.80172932) q[0];
rz(3.7705295) q[1];
sx q[1];
rz(4.9744199) q[1];
sx q[1];
rz(9.0355102) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9297732) q[0];
sx q[0];
rz(-0.7302098) q[0];
sx q[0];
rz(-0.41359652) q[0];
x q[1];
rz(-0.24495895) q[2];
sx q[2];
rz(-2.9548822) q[2];
sx q[2];
rz(2.1706131) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.8402581) q[1];
sx q[1];
rz(-2.0724899) q[1];
sx q[1];
rz(-0.69712104) q[1];
rz(-pi) q[2];
rz(-1.6153159) q[3];
sx q[3];
rz(-0.72784892) q[3];
sx q[3];
rz(2.0522709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4929216) q[2];
sx q[2];
rz(-2.6914983) q[2];
sx q[2];
rz(2.8802803) q[2];
rz(-2.785545) q[3];
sx q[3];
rz(-1.1084403) q[3];
sx q[3];
rz(-2.0480919) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68591958) q[0];
sx q[0];
rz(-2.5643667) q[0];
sx q[0];
rz(-0.24396954) q[0];
rz(-0.4295373) q[1];
sx q[1];
rz(-1.4725087) q[1];
sx q[1];
rz(-0.74091774) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2527494) q[0];
sx q[0];
rz(-2.6359632) q[0];
sx q[0];
rz(-2.3028785) q[0];
rz(-pi) q[1];
rz(-2.7602729) q[2];
sx q[2];
rz(-1.6090074) q[2];
sx q[2];
rz(-0.080121839) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.59830211) q[1];
sx q[1];
rz(-1.0183698) q[1];
sx q[1];
rz(1.8364947) q[1];
rz(-2.5899825) q[3];
sx q[3];
rz(-0.79539652) q[3];
sx q[3];
rz(1.43184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4247158) q[2];
sx q[2];
rz(-1.6191142) q[2];
sx q[2];
rz(-0.87872046) q[2];
rz(2.1667513) q[3];
sx q[3];
rz(-1.8941555) q[3];
sx q[3];
rz(1.6802906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9413302) q[0];
sx q[0];
rz(-0.55568475) q[0];
sx q[0];
rz(-0.19361198) q[0];
rz(0.10784736) q[1];
sx q[1];
rz(-2.4272608) q[1];
sx q[1];
rz(2.329619) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0064751) q[0];
sx q[0];
rz(-2.0251946) q[0];
sx q[0];
rz(-3.019089) q[0];
x q[1];
rz(1.0103364) q[2];
sx q[2];
rz(-2.3046846) q[2];
sx q[2];
rz(2.3488655) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.0206611) q[1];
sx q[1];
rz(-1.963208) q[1];
sx q[1];
rz(-0.25041469) q[1];
rz(1.1319086) q[3];
sx q[3];
rz(-1.697007) q[3];
sx q[3];
rz(0.47078339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0783483) q[2];
sx q[2];
rz(-1.2816659) q[2];
sx q[2];
rz(3.1380624) q[2];
rz(2.6053536) q[3];
sx q[3];
rz(-0.25501525) q[3];
sx q[3];
rz(-0.53853881) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8435159) q[0];
sx q[0];
rz(-1.902782) q[0];
sx q[0];
rz(-2.3115944) q[0];
rz(-1.7371381) q[1];
sx q[1];
rz(-0.91157118) q[1];
sx q[1];
rz(-0.27127582) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41554444) q[0];
sx q[0];
rz(-1.3726639) q[0];
sx q[0];
rz(-1.3130441) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.098255946) q[2];
sx q[2];
rz(-0.66731221) q[2];
sx q[2];
rz(-1.7010393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7978366) q[1];
sx q[1];
rz(-0.35825729) q[1];
sx q[1];
rz(3.0628662) q[1];
rz(0.15437281) q[3];
sx q[3];
rz(-0.8027146) q[3];
sx q[3];
rz(-1.5824535) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6259049) q[2];
sx q[2];
rz(-2.7582464) q[2];
sx q[2];
rz(-0.77451998) q[2];
rz(-2.1382616) q[3];
sx q[3];
rz(-2.7480875) q[3];
sx q[3];
rz(-0.5459319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9970053) q[0];
sx q[0];
rz(-1.0969176) q[0];
sx q[0];
rz(1.0928094) q[0];
rz(1.3358759) q[1];
sx q[1];
rz(-2.3922258) q[1];
sx q[1];
rz(-2.5978079) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0896801) q[0];
sx q[0];
rz(-0.26482115) q[0];
sx q[0];
rz(1.7811243) q[0];
rz(1.4663962) q[2];
sx q[2];
rz(-2.246665) q[2];
sx q[2];
rz(2.3817348) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1256225) q[1];
sx q[1];
rz(-1.984375) q[1];
sx q[1];
rz(-0.57012635) q[1];
x q[2];
rz(-1.5690992) q[3];
sx q[3];
rz(-0.6203531) q[3];
sx q[3];
rz(-1.9242632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.233923) q[2];
sx q[2];
rz(-0.69207865) q[2];
sx q[2];
rz(2.862759) q[2];
rz(-0.28576609) q[3];
sx q[3];
rz(-0.91364342) q[3];
sx q[3];
rz(2.7260776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.435442) q[0];
sx q[0];
rz(-0.64318648) q[0];
sx q[0];
rz(-1.0253133) q[0];
rz(-3.0030491) q[1];
sx q[1];
rz(-1.5899315) q[1];
sx q[1];
rz(2.9655546) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6556102) q[0];
sx q[0];
rz(-1.8261251) q[0];
sx q[0];
rz(-0.53389733) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.3761212) q[2];
sx q[2];
rz(-1.2174213) q[2];
sx q[2];
rz(0.88549685) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.306539) q[1];
sx q[1];
rz(-1.3471148) q[1];
sx q[1];
rz(2.064631) q[1];
rz(1.9344091) q[3];
sx q[3];
rz(-2.6904958) q[3];
sx q[3];
rz(2.0157984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.90870086) q[2];
sx q[2];
rz(-0.42753926) q[2];
sx q[2];
rz(-2.3517081) q[2];
rz(-2.8120153) q[3];
sx q[3];
rz(-1.7695844) q[3];
sx q[3];
rz(2.4901539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71094197) q[0];
sx q[0];
rz(-2.9281404) q[0];
sx q[0];
rz(-0.99367225) q[0];
rz(2.1395394) q[1];
sx q[1];
rz(-0.99169815) q[1];
sx q[1];
rz(-1.4659945) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77540175) q[0];
sx q[0];
rz(-2.4537773) q[0];
sx q[0];
rz(-2.8596876) q[0];
rz(-pi) q[1];
rz(1.1118717) q[2];
sx q[2];
rz(-1.4676158) q[2];
sx q[2];
rz(2.9261719) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1217782) q[1];
sx q[1];
rz(-1.7516859) q[1];
sx q[1];
rz(-1.719285) q[1];
rz(-pi) q[2];
rz(-0.44963533) q[3];
sx q[3];
rz(-2.2581266) q[3];
sx q[3];
rz(-0.80865142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.34622908) q[2];
sx q[2];
rz(-0.84285223) q[2];
sx q[2];
rz(-2.8818434) q[2];
rz(-2.3527457) q[3];
sx q[3];
rz(-2.2987821) q[3];
sx q[3];
rz(1.3938676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9144834) q[0];
sx q[0];
rz(-1.8552584) q[0];
sx q[0];
rz(-0.9911384) q[0];
rz(-1.8046851) q[1];
sx q[1];
rz(-2.3146345) q[1];
sx q[1];
rz(-1.7327259) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0175655) q[0];
sx q[0];
rz(-2.6691737) q[0];
sx q[0];
rz(2.8851465) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0524247) q[2];
sx q[2];
rz(-1.0341511) q[2];
sx q[2];
rz(2.9111957) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.046890586) q[1];
sx q[1];
rz(-2.6841503) q[1];
sx q[1];
rz(-0.66384683) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4662631) q[3];
sx q[3];
rz(-1.7033782) q[3];
sx q[3];
rz(1.7426566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1651429) q[2];
sx q[2];
rz(-1.6565448) q[2];
sx q[2];
rz(2.5531947) q[2];
rz(-0.31383651) q[3];
sx q[3];
rz(-2.2321759) q[3];
sx q[3];
rz(-0.61217827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2345851) q[0];
sx q[0];
rz(-2.0186277) q[0];
sx q[0];
rz(-3.0871952) q[0];
rz(-2.3811978) q[1];
sx q[1];
rz(-2.1957928) q[1];
sx q[1];
rz(-1.0992345) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71610036) q[0];
sx q[0];
rz(-1.4045225) q[0];
sx q[0];
rz(-0.99207741) q[0];
rz(-pi) q[1];
rz(-0.77581232) q[2];
sx q[2];
rz(-0.96539984) q[2];
sx q[2];
rz(-2.2697732) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.526872) q[1];
sx q[1];
rz(-0.81382591) q[1];
sx q[1];
rz(-0.24095638) q[1];
rz(-pi) q[2];
rz(-2.0012747) q[3];
sx q[3];
rz(-3.0196689) q[3];
sx q[3];
rz(2.6170128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.2796563) q[2];
sx q[2];
rz(-1.6013689) q[2];
sx q[2];
rz(0.81735617) q[2];
rz(-3.127408) q[3];
sx q[3];
rz(-2.7458906) q[3];
sx q[3];
rz(-2.0524249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8297183) q[0];
sx q[0];
rz(-2.0901966) q[0];
sx q[0];
rz(-2.9201087) q[0];
rz(2.7102176) q[1];
sx q[1];
rz(-0.90310496) q[1];
sx q[1];
rz(2.121714) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9572302) q[0];
sx q[0];
rz(-1.5505322) q[0];
sx q[0];
rz(0.046837687) q[0];
rz(-0.20163147) q[2];
sx q[2];
rz(-1.4678737) q[2];
sx q[2];
rz(2.0992744) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4012667) q[1];
sx q[1];
rz(-1.5366842) q[1];
sx q[1];
rz(1.9032065) q[1];
rz(-pi) q[2];
rz(-1.4388891) q[3];
sx q[3];
rz(-1.7057837) q[3];
sx q[3];
rz(1.5832981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7027616) q[2];
sx q[2];
rz(-1.2063382) q[2];
sx q[2];
rz(-1.2552525) q[2];
rz(-3.1158279) q[3];
sx q[3];
rz(-1.8119101) q[3];
sx q[3];
rz(-2.6222031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7230566) q[0];
sx q[0];
rz(-2.3522455) q[0];
sx q[0];
rz(1.1242207) q[0];
rz(-0.33456805) q[1];
sx q[1];
rz(-2.6523013) q[1];
sx q[1];
rz(2.1025067) q[1];
rz(-0.46229565) q[2];
sx q[2];
rz(-1.3318599) q[2];
sx q[2];
rz(-0.75564043) q[2];
rz(-2.6565537) q[3];
sx q[3];
rz(-2.1342604) q[3];
sx q[3];
rz(0.66935183) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
