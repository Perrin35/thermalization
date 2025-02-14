OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98439944) q[0];
sx q[0];
rz(-1.1097044) q[0];
sx q[0];
rz(-2.1362526) q[0];
rz(0.73468626) q[1];
sx q[1];
rz(3.4975657) q[1];
sx q[1];
rz(11.674292) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6280373) q[0];
sx q[0];
rz(-1.6694067) q[0];
sx q[0];
rz(-0.74752083) q[0];
rz(-pi) q[1];
rz(0.96718915) q[2];
sx q[2];
rz(-1.4481067) q[2];
sx q[2];
rz(-1.3792147) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8651651) q[1];
sx q[1];
rz(-0.49404198) q[1];
sx q[1];
rz(2.9345291) q[1];
x q[2];
rz(2.4446176) q[3];
sx q[3];
rz(-2.4066952) q[3];
sx q[3];
rz(-1.01222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72508183) q[2];
sx q[2];
rz(-1.1434914) q[2];
sx q[2];
rz(1.7007281) q[2];
rz(-2.1848988) q[3];
sx q[3];
rz(-2.0243093) q[3];
sx q[3];
rz(-2.8982437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(2.6956534) q[0];
sx q[0];
rz(-2.74701) q[0];
sx q[0];
rz(-2.0126427) q[0];
rz(2.8961862) q[1];
sx q[1];
rz(-1.3134198) q[1];
sx q[1];
rz(0.66171563) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.322987) q[0];
sx q[0];
rz(-1.2862236) q[0];
sx q[0];
rz(-0.63708441) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7595046) q[2];
sx q[2];
rz(-1.3616865) q[2];
sx q[2];
rz(-1.8764302) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.405141) q[1];
sx q[1];
rz(-1.5390522) q[1];
sx q[1];
rz(-0.27145465) q[1];
rz(-pi) q[2];
rz(-1.2564959) q[3];
sx q[3];
rz(-0.30499015) q[3];
sx q[3];
rz(-1.631032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5976065) q[2];
sx q[2];
rz(-2.643955) q[2];
sx q[2];
rz(2.0987161) q[2];
rz(1.6889702) q[3];
sx q[3];
rz(-0.85113227) q[3];
sx q[3];
rz(-1.1037306) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76688981) q[0];
sx q[0];
rz(-0.68793982) q[0];
sx q[0];
rz(-1.3091298) q[0];
rz(-2.6922928) q[1];
sx q[1];
rz(-2.4520912) q[1];
sx q[1];
rz(1.4264533) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4065789) q[0];
sx q[0];
rz(-1.3258223) q[0];
sx q[0];
rz(-1.2981775) q[0];
x q[1];
rz(1.9932991) q[2];
sx q[2];
rz(-1.587217) q[2];
sx q[2];
rz(-0.33994477) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.215308) q[1];
sx q[1];
rz(-2.4325772) q[1];
sx q[1];
rz(-3.0224817) q[1];
x q[2];
rz(1.8168114) q[3];
sx q[3];
rz(-1.6520471) q[3];
sx q[3];
rz(1.2556835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7766777) q[2];
sx q[2];
rz(-1.9675576) q[2];
sx q[2];
rz(0.066369973) q[2];
rz(-2.5868609) q[3];
sx q[3];
rz(-1.6603671) q[3];
sx q[3];
rz(-1.5167282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9471112) q[0];
sx q[0];
rz(-0.28488657) q[0];
sx q[0];
rz(-0.77199212) q[0];
rz(-2.4783065) q[1];
sx q[1];
rz(-2.3978077) q[1];
sx q[1];
rz(-3.0139121) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0136702) q[0];
sx q[0];
rz(-1.575673) q[0];
sx q[0];
rz(2.4749066) q[0];
rz(-pi) q[1];
rz(-0.56945412) q[2];
sx q[2];
rz(-0.63210154) q[2];
sx q[2];
rz(-2.5394627) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3290594) q[1];
sx q[1];
rz(-0.61186463) q[1];
sx q[1];
rz(2.3298323) q[1];
rz(-pi) q[2];
rz(0.33416602) q[3];
sx q[3];
rz(-0.45229707) q[3];
sx q[3];
rz(2.3826007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.53106236) q[2];
sx q[2];
rz(-0.9610815) q[2];
sx q[2];
rz(0.5985716) q[2];
rz(-2.5564204) q[3];
sx q[3];
rz(-1.7681311) q[3];
sx q[3];
rz(3.0857871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35036206) q[0];
sx q[0];
rz(-1.6329916) q[0];
sx q[0];
rz(2.1685261) q[0];
rz(-0.19634253) q[1];
sx q[1];
rz(-2.0025496) q[1];
sx q[1];
rz(1.4686718) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.136841) q[0];
sx q[0];
rz(-0.84293619) q[0];
sx q[0];
rz(2.1086295) q[0];
x q[1];
rz(-1.0650915) q[2];
sx q[2];
rz(-1.347216) q[2];
sx q[2];
rz(-0.93709842) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1717559) q[1];
sx q[1];
rz(-1.3984658) q[1];
sx q[1];
rz(1.2660162) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0861868) q[3];
sx q[3];
rz(-2.1089777) q[3];
sx q[3];
rz(2.2937867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4520182) q[2];
sx q[2];
rz(-1.8147899) q[2];
sx q[2];
rz(-0.16239521) q[2];
rz(-1.1299805) q[3];
sx q[3];
rz(-0.88310784) q[3];
sx q[3];
rz(1.1348772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75055403) q[0];
sx q[0];
rz(-0.51893187) q[0];
sx q[0];
rz(-0.050405141) q[0];
rz(-2.9734036) q[1];
sx q[1];
rz(-1.5805809) q[1];
sx q[1];
rz(-0.87356299) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7058327) q[0];
sx q[0];
rz(-1.8144023) q[0];
sx q[0];
rz(0.63023366) q[0];
rz(-pi) q[1];
rz(1.7123772) q[2];
sx q[2];
rz(-1.537286) q[2];
sx q[2];
rz(-0.17591116) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.5205914) q[1];
sx q[1];
rz(-2.0157923) q[1];
sx q[1];
rz(0.59153647) q[1];
rz(-2.843552) q[3];
sx q[3];
rz(-1.3750769) q[3];
sx q[3];
rz(-0.51575553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1217338) q[2];
sx q[2];
rz(-0.21616082) q[2];
sx q[2];
rz(0.38062322) q[2];
rz(-0.56626433) q[3];
sx q[3];
rz(-1.7411313) q[3];
sx q[3];
rz(-2.3316021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64067632) q[0];
sx q[0];
rz(-2.930142) q[0];
sx q[0];
rz(-2.7389615) q[0];
rz(1.3817878) q[1];
sx q[1];
rz(-2.333162) q[1];
sx q[1];
rz(-0.056338739) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23898174) q[0];
sx q[0];
rz(-1.1711811) q[0];
sx q[0];
rz(2.7863471) q[0];
rz(0.041002657) q[2];
sx q[2];
rz(-1.1004538) q[2];
sx q[2];
rz(3.0029763) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.496324) q[1];
sx q[1];
rz(-2.5773415) q[1];
sx q[1];
rz(-1.5936046) q[1];
x q[2];
rz(1.0545441) q[3];
sx q[3];
rz(-2.1813381) q[3];
sx q[3];
rz(2.1757389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9075883) q[2];
sx q[2];
rz(-1.7981217) q[2];
sx q[2];
rz(-0.50055707) q[2];
rz(3.0808926) q[3];
sx q[3];
rz(-1.7874291) q[3];
sx q[3];
rz(0.48262706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57182264) q[0];
sx q[0];
rz(-1.3615384) q[0];
sx q[0];
rz(1.7806336) q[0];
rz(-2.3979777) q[1];
sx q[1];
rz(-1.4185602) q[1];
sx q[1];
rz(-3.0027711) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25538975) q[0];
sx q[0];
rz(-2.3950044) q[0];
sx q[0];
rz(-2.3476178) q[0];
x q[1];
rz(2.2108364) q[2];
sx q[2];
rz(-0.70022196) q[2];
sx q[2];
rz(-0.28565684) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.946859) q[1];
sx q[1];
rz(-1.7725962) q[1];
sx q[1];
rz(-1.9029593) q[1];
x q[2];
rz(-2.9337594) q[3];
sx q[3];
rz(-1.1869988) q[3];
sx q[3];
rz(1.919534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.6737785) q[2];
sx q[2];
rz(-1.0445341) q[2];
sx q[2];
rz(-2.4274801) q[2];
rz(2.1821187) q[3];
sx q[3];
rz(-1.6474612) q[3];
sx q[3];
rz(2.1625471) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4176843) q[0];
sx q[0];
rz(-2.4651616) q[0];
sx q[0];
rz(0.12829256) q[0];
rz(-2.7543606) q[1];
sx q[1];
rz(-0.59806824) q[1];
sx q[1];
rz(-1.8181575) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.58777) q[0];
sx q[0];
rz(-0.21311305) q[0];
sx q[0];
rz(1.0928667) q[0];
x q[1];
rz(0.88487423) q[2];
sx q[2];
rz(-1.0909683) q[2];
sx q[2];
rz(0.72140933) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1897498) q[1];
sx q[1];
rz(-1.4268759) q[1];
sx q[1];
rz(0.014149498) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3714482) q[3];
sx q[3];
rz(-2.2888765) q[3];
sx q[3];
rz(-1.9950641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4667751) q[2];
sx q[2];
rz(-1.6206436) q[2];
sx q[2];
rz(-2.9368994) q[2];
rz(-1.8681059) q[3];
sx q[3];
rz(-2.4114362) q[3];
sx q[3];
rz(-1.5761121) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650763) q[0];
sx q[0];
rz(-2.5193494) q[0];
sx q[0];
rz(-0.21328558) q[0];
rz(0.22008303) q[1];
sx q[1];
rz(-1.8890231) q[1];
sx q[1];
rz(-1.782104) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8187788) q[0];
sx q[0];
rz(-1.5867763) q[0];
sx q[0];
rz(3.131991) q[0];
rz(-pi) q[1];
rz(0.47021265) q[2];
sx q[2];
rz(-0.5222975) q[2];
sx q[2];
rz(3.0015415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2343341) q[1];
sx q[1];
rz(-1.2266907) q[1];
sx q[1];
rz(-2.1109778) q[1];
x q[2];
rz(2.1207934) q[3];
sx q[3];
rz(-0.74899835) q[3];
sx q[3];
rz(-1.4326381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.517211) q[2];
sx q[2];
rz(-1.5959847) q[2];
sx q[2];
rz(-0.33779302) q[2];
rz(-1.4126011) q[3];
sx q[3];
rz(-1.1510886) q[3];
sx q[3];
rz(-0.82144773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4008041) q[0];
sx q[0];
rz(-0.82397006) q[0];
sx q[0];
rz(-1.9116221) q[0];
rz(-0.38372718) q[1];
sx q[1];
rz(-1.6576672) q[1];
sx q[1];
rz(-2.3794649) q[1];
rz(2.4356213) q[2];
sx q[2];
rz(-2.6298475) q[2];
sx q[2];
rz(2.0012326) q[2];
rz(1.0596342) q[3];
sx q[3];
rz(-1.7851024) q[3];
sx q[3];
rz(2.057597) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
