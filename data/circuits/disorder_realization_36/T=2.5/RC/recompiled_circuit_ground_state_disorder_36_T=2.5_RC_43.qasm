OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.8776956) q[0];
sx q[0];
rz(-2.1030302) q[0];
sx q[0];
rz(-2.7466018) q[0];
rz(0.78139961) q[1];
sx q[1];
rz(-1.7698987) q[1];
sx q[1];
rz(-0.79010195) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2443197) q[0];
sx q[0];
rz(-1.6023912) q[0];
sx q[0];
rz(2.5665166) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1320791) q[2];
sx q[2];
rz(-2.2971675) q[2];
sx q[2];
rz(1.4216636) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.61617) q[1];
sx q[1];
rz(-1.425263) q[1];
sx q[1];
rz(-1.0305745) q[1];
rz(-pi) q[2];
rz(-1.9856307) q[3];
sx q[3];
rz(-2.0409763) q[3];
sx q[3];
rz(-0.0068723504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6140952) q[2];
sx q[2];
rz(-0.082753269) q[2];
sx q[2];
rz(0.56438524) q[2];
rz(-0.91974059) q[3];
sx q[3];
rz(-0.65110937) q[3];
sx q[3];
rz(1.9981599) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1564002) q[0];
sx q[0];
rz(-2.1512478) q[0];
sx q[0];
rz(1.1457957) q[0];
rz(1.0076373) q[1];
sx q[1];
rz(-2.2538908) q[1];
sx q[1];
rz(0.064858286) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6211557) q[0];
sx q[0];
rz(-1.1722992) q[0];
sx q[0];
rz(0.99834672) q[0];
rz(-0.62297946) q[2];
sx q[2];
rz(-2.2588135) q[2];
sx q[2];
rz(2.1114897) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.3895208) q[1];
sx q[1];
rz(-0.18587843) q[1];
sx q[1];
rz(-0.13711547) q[1];
x q[2];
rz(2.3228538) q[3];
sx q[3];
rz(-2.4609339) q[3];
sx q[3];
rz(-2.7622032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8811983) q[2];
sx q[2];
rz(-0.10513267) q[2];
sx q[2];
rz(-0.51327389) q[2];
rz(-1.5848292) q[3];
sx q[3];
rz(-1.1791752) q[3];
sx q[3];
rz(1.56196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.8455115) q[0];
sx q[0];
rz(-3.0873612) q[0];
sx q[0];
rz(-0.84501141) q[0];
rz(-2.4301279) q[1];
sx q[1];
rz(-0.80159801) q[1];
sx q[1];
rz(2.5149288) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1820647) q[0];
sx q[0];
rz(-1.3747038) q[0];
sx q[0];
rz(-0.0053868731) q[0];
rz(-2.0925631) q[2];
sx q[2];
rz(-0.43704068) q[2];
sx q[2];
rz(-2.4573987) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6271403) q[1];
sx q[1];
rz(-2.4042292) q[1];
sx q[1];
rz(-1.8335349) q[1];
rz(0.0073284433) q[3];
sx q[3];
rz(-1.0798608) q[3];
sx q[3];
rz(-2.1027264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6381548) q[2];
sx q[2];
rz(-1.902709) q[2];
sx q[2];
rz(-0.16866355) q[2];
rz(-0.1564129) q[3];
sx q[3];
rz(-0.63786879) q[3];
sx q[3];
rz(-2.8580581) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4446568) q[0];
sx q[0];
rz(-2.038027) q[0];
sx q[0];
rz(-2.7705833) q[0];
rz(1.8095398) q[1];
sx q[1];
rz(-1.5468372) q[1];
sx q[1];
rz(0.36088774) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9661117) q[0];
sx q[0];
rz(-1.4752441) q[0];
sx q[0];
rz(3.0555834) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0118594) q[2];
sx q[2];
rz(-0.12309821) q[2];
sx q[2];
rz(2.77746) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8239221) q[1];
sx q[1];
rz(-0.63613632) q[1];
sx q[1];
rz(-2.0348861) q[1];
x q[2];
rz(2.4655409) q[3];
sx q[3];
rz(-2.0094487) q[3];
sx q[3];
rz(-1.0880409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0648301) q[2];
sx q[2];
rz(-0.37134376) q[2];
sx q[2];
rz(-1.4999464) q[2];
rz(-2.0867945) q[3];
sx q[3];
rz(-1.0902371) q[3];
sx q[3];
rz(-1.9487901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.763279) q[0];
sx q[0];
rz(-1.9739456) q[0];
sx q[0];
rz(1.2953229) q[0];
rz(-3.0150705) q[1];
sx q[1];
rz(-2.449072) q[1];
sx q[1];
rz(-1.893938) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0678789) q[0];
sx q[0];
rz(-2.824221) q[0];
sx q[0];
rz(0.48010357) q[0];
x q[1];
rz(-0.40389668) q[2];
sx q[2];
rz(-0.81012175) q[2];
sx q[2];
rz(-0.4357117) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3646255) q[1];
sx q[1];
rz(-2.7086621) q[1];
sx q[1];
rz(-2.6294699) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9548847) q[3];
sx q[3];
rz(-1.8481585) q[3];
sx q[3];
rz(1.385078) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.7724472) q[2];
sx q[2];
rz(-2.7771066) q[2];
sx q[2];
rz(-2.6920953) q[2];
rz(-2.7033324) q[3];
sx q[3];
rz(-2.0354383) q[3];
sx q[3];
rz(-1.0615758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448369) q[0];
sx q[0];
rz(-2.0773092) q[0];
sx q[0];
rz(2.7979895) q[0];
rz(-2.4876439) q[1];
sx q[1];
rz(-0.75075722) q[1];
sx q[1];
rz(-0.59069815) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1794248) q[0];
sx q[0];
rz(-0.54467595) q[0];
sx q[0];
rz(2.0228128) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4015437) q[2];
sx q[2];
rz(-0.97566665) q[2];
sx q[2];
rz(0.33444946) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.5370767) q[1];
sx q[1];
rz(-0.13128102) q[1];
sx q[1];
rz(1.4907565) q[1];
rz(-0.77883522) q[3];
sx q[3];
rz(-2.2766354) q[3];
sx q[3];
rz(2.6293628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5609453) q[2];
sx q[2];
rz(-2.2717768) q[2];
sx q[2];
rz(0.84442863) q[2];
rz(2.1010418) q[3];
sx q[3];
rz(-1.8620164) q[3];
sx q[3];
rz(-1.1738663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(1.7548673) q[0];
sx q[0];
rz(-2.6193021) q[0];
sx q[0];
rz(-2.0544384) q[0];
rz(-0.48671752) q[1];
sx q[1];
rz(-1.4954647) q[1];
sx q[1];
rz(-1.4901935) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5958811) q[0];
sx q[0];
rz(-2.427248) q[0];
sx q[0];
rz(-1.7313135) q[0];
x q[1];
rz(-0.59972024) q[2];
sx q[2];
rz(-2.2968946) q[2];
sx q[2];
rz(-3.0350189) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.5512276) q[1];
sx q[1];
rz(-1.9842973) q[1];
sx q[1];
rz(1.9848003) q[1];
rz(-pi) q[2];
rz(2.9341615) q[3];
sx q[3];
rz(-1.279502) q[3];
sx q[3];
rz(2.8136611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3460441) q[2];
sx q[2];
rz(-0.86980021) q[2];
sx q[2];
rz(-0.33934936) q[2];
rz(-0.36335534) q[3];
sx q[3];
rz(-2.8715869) q[3];
sx q[3];
rz(-1.5615162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.617008) q[0];
sx q[0];
rz(-1.120485) q[0];
sx q[0];
rz(-2.8560915) q[0];
rz(-0.78954804) q[1];
sx q[1];
rz(-1.6098846) q[1];
sx q[1];
rz(-1.550536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94211468) q[0];
sx q[0];
rz(-1.3923595) q[0];
sx q[0];
rz(-2.8392544) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9934814) q[2];
sx q[2];
rz(-2.092871) q[2];
sx q[2];
rz(-2.8294308) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.62666368) q[1];
sx q[1];
rz(-1.0486992) q[1];
sx q[1];
rz(2.057328) q[1];
rz(-pi) q[2];
rz(1.490216) q[3];
sx q[3];
rz(-2.2364607) q[3];
sx q[3];
rz(-2.5192266) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8153814) q[2];
sx q[2];
rz(-2.9823494) q[2];
sx q[2];
rz(2.2747269) q[2];
rz(2.4032812) q[3];
sx q[3];
rz(-1.5840931) q[3];
sx q[3];
rz(-0.90832925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
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
rz(-0.065780491) q[0];
sx q[0];
rz(-2.3920238) q[0];
sx q[0];
rz(-0.69920364) q[0];
rz(0.54222822) q[1];
sx q[1];
rz(-1.3424073) q[1];
sx q[1];
rz(0.29621616) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2767773) q[0];
sx q[0];
rz(-2.3010198) q[0];
sx q[0];
rz(3.0481799) q[0];
rz(-pi) q[1];
x q[1];
rz(0.64546236) q[2];
sx q[2];
rz(-0.92593595) q[2];
sx q[2];
rz(2.1301816) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.49668) q[1];
sx q[1];
rz(-2.1261922) q[1];
sx q[1];
rz(0.72873656) q[1];
x q[2];
rz(-1.9486197) q[3];
sx q[3];
rz(-1.6512031) q[3];
sx q[3];
rz(0.37624826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.88973033) q[2];
sx q[2];
rz(-2.6984062) q[2];
sx q[2];
rz(-0.52081338) q[2];
rz(-3.1098747) q[3];
sx q[3];
rz(-2.7364276) q[3];
sx q[3];
rz(-0.051137663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6601335) q[0];
sx q[0];
rz(-2.0429459) q[0];
sx q[0];
rz(-0.81137401) q[0];
rz(2.7157281) q[1];
sx q[1];
rz(-0.90161294) q[1];
sx q[1];
rz(-1.3070377) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9872914) q[0];
sx q[0];
rz(-1.5679006) q[0];
sx q[0];
rz(-1.5737246) q[0];
x q[1];
rz(1.1063031) q[2];
sx q[2];
rz(-1.5896322) q[2];
sx q[2];
rz(2.1003758) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2073686) q[1];
sx q[1];
rz(-1.1282053) q[1];
sx q[1];
rz(1.2138483) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8167231) q[3];
sx q[3];
rz(-1.4063971) q[3];
sx q[3];
rz(-1.4332958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1618774) q[2];
sx q[2];
rz(-2.3780509) q[2];
sx q[2];
rz(-0.36583501) q[2];
rz(-2.2063935) q[3];
sx q[3];
rz(-1.5465982) q[3];
sx q[3];
rz(-0.38177761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81015051) q[0];
sx q[0];
rz(-1.7773542) q[0];
sx q[0];
rz(-1.6553028) q[0];
rz(-1.4797795) q[1];
sx q[1];
rz(-1.9098837) q[1];
sx q[1];
rz(2.2178537) q[1];
rz(-1.6681832) q[2];
sx q[2];
rz(-2.591572) q[2];
sx q[2];
rz(-0.77629065) q[2];
rz(0.52025683) q[3];
sx q[3];
rz(-2.8388966) q[3];
sx q[3];
rz(2.7404529) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
