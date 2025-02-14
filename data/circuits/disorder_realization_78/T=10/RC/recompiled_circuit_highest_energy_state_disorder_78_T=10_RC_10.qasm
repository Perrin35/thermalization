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
rz(-1.7595093) q[0];
sx q[0];
rz(-2.3926662) q[0];
sx q[0];
rz(1.393526) q[0];
rz(-2.3387609) q[1];
sx q[1];
rz(3.9376942) q[1];
sx q[1];
rz(12.035523) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8061764) q[0];
sx q[0];
rz(-2.5160976) q[0];
sx q[0];
rz(-0.92910398) q[0];
rz(-0.10239281) q[2];
sx q[2];
rz(-2.9739485) q[2];
sx q[2];
rz(1.2861811) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.4208274) q[1];
sx q[1];
rz(-1.6457575) q[1];
sx q[1];
rz(1.0383181) q[1];
rz(-pi) q[2];
rz(3.0396456) q[3];
sx q[3];
rz(-1.0539712) q[3];
sx q[3];
rz(-1.0727796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.406245) q[2];
sx q[2];
rz(-2.2082059) q[2];
sx q[2];
rz(2.1437342) q[2];
rz(-0.88456279) q[3];
sx q[3];
rz(-1.9622784) q[3];
sx q[3];
rz(-2.4295889) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.74580055) q[0];
sx q[0];
rz(-0.53676787) q[0];
sx q[0];
rz(2.0641548) q[0];
rz(-2.5224345) q[1];
sx q[1];
rz(-1.2666603) q[1];
sx q[1];
rz(-1.9177297) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3676476) q[0];
sx q[0];
rz(-1.2804693) q[0];
sx q[0];
rz(0.18240697) q[0];
rz(-pi) q[1];
rz(3.139758) q[2];
sx q[2];
rz(-1.6009838) q[2];
sx q[2];
rz(2.052553) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6757766) q[1];
sx q[1];
rz(-0.73938939) q[1];
sx q[1];
rz(2.5104816) q[1];
rz(-pi) q[2];
rz(-2.4183351) q[3];
sx q[3];
rz(-0.92953909) q[3];
sx q[3];
rz(1.0783522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.18664843) q[2];
sx q[2];
rz(-0.84543219) q[2];
sx q[2];
rz(-2.1050982) q[2];
rz(1.1896108) q[3];
sx q[3];
rz(-1.2396953) q[3];
sx q[3];
rz(-0.069124393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67967296) q[0];
sx q[0];
rz(-0.24149495) q[0];
sx q[0];
rz(-1.6756206) q[0];
rz(-1.2566603) q[1];
sx q[1];
rz(-2.4938221) q[1];
sx q[1];
rz(-0.57674903) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.758848) q[0];
sx q[0];
rz(-0.85475105) q[0];
sx q[0];
rz(2.3780253) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6079712) q[2];
sx q[2];
rz(-2.0293183) q[2];
sx q[2];
rz(-0.22500817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.27611342) q[1];
sx q[1];
rz(-2.7250368) q[1];
sx q[1];
rz(1.4426484) q[1];
rz(-1.4068094) q[3];
sx q[3];
rz(-1.6583558) q[3];
sx q[3];
rz(-1.5045741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74167788) q[2];
sx q[2];
rz(-2.8027813) q[2];
sx q[2];
rz(-2.3728288) q[2];
rz(-2.9703043) q[3];
sx q[3];
rz(-1.7007217) q[3];
sx q[3];
rz(0.35458529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6866368) q[0];
sx q[0];
rz(-0.88382116) q[0];
sx q[0];
rz(-2.6052642) q[0];
rz(-0.80279154) q[1];
sx q[1];
rz(-1.2840459) q[1];
sx q[1];
rz(0.098043052) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2659543) q[0];
sx q[0];
rz(-1.2438626) q[0];
sx q[0];
rz(-3.0578259) q[0];
rz(-0.98020245) q[2];
sx q[2];
rz(-0.88104311) q[2];
sx q[2];
rz(2.142327) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.72934276) q[1];
sx q[1];
rz(-2.5855484) q[1];
sx q[1];
rz(1.1792609) q[1];
rz(-pi) q[2];
x q[2];
rz(0.40880568) q[3];
sx q[3];
rz(-1.9681566) q[3];
sx q[3];
rz(1.8026343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1596277) q[2];
sx q[2];
rz(-2.1106796) q[2];
sx q[2];
rz(0.31943303) q[2];
rz(0.56644136) q[3];
sx q[3];
rz(-2.3947377) q[3];
sx q[3];
rz(-1.9802861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7403858) q[0];
sx q[0];
rz(-1.8219319) q[0];
sx q[0];
rz(2.5886986) q[0];
rz(2.3623908) q[1];
sx q[1];
rz(-0.82256493) q[1];
sx q[1];
rz(2.353277) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6471651) q[0];
sx q[0];
rz(-2.8719423) q[0];
sx q[0];
rz(-2.8799876) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.20878883) q[2];
sx q[2];
rz(-1.5475071) q[2];
sx q[2];
rz(1.67729) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.3723046) q[1];
sx q[1];
rz(-1.2944229) q[1];
sx q[1];
rz(0.45120542) q[1];
rz(0.68627091) q[3];
sx q[3];
rz(-2.7618558) q[3];
sx q[3];
rz(0.59829955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.39764443) q[2];
sx q[2];
rz(-1.1279736) q[2];
sx q[2];
rz(-2.9359342) q[2];
rz(1.7913943) q[3];
sx q[3];
rz(-0.74068991) q[3];
sx q[3];
rz(3.1309483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.81269294) q[0];
sx q[0];
rz(-0.44726547) q[0];
sx q[0];
rz(1.6269667) q[0];
rz(-0.80500066) q[1];
sx q[1];
rz(-1.6945508) q[1];
sx q[1];
rz(2.8256493) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63271967) q[0];
sx q[0];
rz(-2.609786) q[0];
sx q[0];
rz(2.4200632) q[0];
rz(-pi) q[1];
rz(-2.6348389) q[2];
sx q[2];
rz(-2.4656418) q[2];
sx q[2];
rz(1.6981704) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.915357) q[1];
sx q[1];
rz(-1.2958044) q[1];
sx q[1];
rz(0.26696856) q[1];
rz(-pi) q[2];
x q[2];
rz(0.31031761) q[3];
sx q[3];
rz(-1.4253221) q[3];
sx q[3];
rz(0.014257243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.67779764) q[2];
sx q[2];
rz(-0.69910502) q[2];
sx q[2];
rz(0.15764906) q[2];
rz(-2.6080103) q[3];
sx q[3];
rz(-1.491051) q[3];
sx q[3];
rz(1.6238448) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59876281) q[0];
sx q[0];
rz(-0.49648008) q[0];
sx q[0];
rz(-2.1606523) q[0];
rz(-2.6988103) q[1];
sx q[1];
rz(-1.1863703) q[1];
sx q[1];
rz(-0.47169366) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6169051) q[0];
sx q[0];
rz(-1.1393095) q[0];
sx q[0];
rz(0.3226852) q[0];
rz(-0.27916698) q[2];
sx q[2];
rz(-1.4067603) q[2];
sx q[2];
rz(3.0085473) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0761169) q[1];
sx q[1];
rz(-0.13607506) q[1];
sx q[1];
rz(0.10175609) q[1];
x q[2];
rz(1.7552492) q[3];
sx q[3];
rz(-1.9058936) q[3];
sx q[3];
rz(1.9012251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.20785759) q[2];
sx q[2];
rz(-0.8693049) q[2];
sx q[2];
rz(0.11824879) q[2];
rz(2.5737428) q[3];
sx q[3];
rz(-1.3563145) q[3];
sx q[3];
rz(-2.3387199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0589013) q[0];
sx q[0];
rz(-1.4014129) q[0];
sx q[0];
rz(-2.4543104) q[0];
rz(-1.8742689) q[1];
sx q[1];
rz(-1.9272389) q[1];
sx q[1];
rz(-2.5158688) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5289297) q[0];
sx q[0];
rz(-1.5497396) q[0];
sx q[0];
rz(0.0046757129) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5753046) q[2];
sx q[2];
rz(-2.1293961) q[2];
sx q[2];
rz(-1.9273014) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3753452) q[1];
sx q[1];
rz(-1.5832003) q[1];
sx q[1];
rz(0.9961025) q[1];
rz(2.8318303) q[3];
sx q[3];
rz(-1.0132257) q[3];
sx q[3];
rz(2.6738809) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.9132793) q[2];
sx q[2];
rz(-1.8931171) q[2];
sx q[2];
rz(2.0466364) q[2];
rz(2.6878808) q[3];
sx q[3];
rz(-2.3504421) q[3];
sx q[3];
rz(2.9554844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5943282) q[0];
sx q[0];
rz(-2.6481977) q[0];
sx q[0];
rz(-3.0739947) q[0];
rz(1.0293845) q[1];
sx q[1];
rz(-1.8600347) q[1];
sx q[1];
rz(3.0060815) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1757723) q[0];
sx q[0];
rz(-2.1012573) q[0];
sx q[0];
rz(0.6483174) q[0];
rz(2.4004647) q[2];
sx q[2];
rz(-1.4504045) q[2];
sx q[2];
rz(0.84726221) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.9635591) q[1];
sx q[1];
rz(-1.2894703) q[1];
sx q[1];
rz(2.540349) q[1];
rz(1.3247847) q[3];
sx q[3];
rz(-1.8455077) q[3];
sx q[3];
rz(1.5266071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4542666) q[2];
sx q[2];
rz(-2.8886075) q[2];
sx q[2];
rz(-0.48736408) q[2];
rz(-2.7519915) q[3];
sx q[3];
rz(-0.95249683) q[3];
sx q[3];
rz(-3.0756557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.6479263) q[0];
sx q[0];
rz(-2.0559897) q[0];
sx q[0];
rz(-1.7568463) q[0];
rz(1.0549649) q[1];
sx q[1];
rz(-2.2672548) q[1];
sx q[1];
rz(-2.8816282) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1459634) q[0];
sx q[0];
rz(-0.71486799) q[0];
sx q[0];
rz(-1.9674106) q[0];
rz(1.7021263) q[2];
sx q[2];
rz(-1.471278) q[2];
sx q[2];
rz(3.1414349) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1321226) q[1];
sx q[1];
rz(-2.1984221) q[1];
sx q[1];
rz(-3.0982137) q[1];
x q[2];
rz(0.72885363) q[3];
sx q[3];
rz(-2.8480004) q[3];
sx q[3];
rz(2.7135682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7197623) q[2];
sx q[2];
rz(-2.7037342) q[2];
sx q[2];
rz(-1.7913294) q[2];
rz(-1.0849902) q[3];
sx q[3];
rz(-1.2144054) q[3];
sx q[3];
rz(2.0124281) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0706901) q[0];
sx q[0];
rz(-2.4546843) q[0];
sx q[0];
rz(-0.51921459) q[0];
rz(-1.6593973) q[1];
sx q[1];
rz(-1.2980325) q[1];
sx q[1];
rz(1.4549805) q[1];
rz(-2.5605911) q[2];
sx q[2];
rz(-2.0145363) q[2];
sx q[2];
rz(1.2546652) q[2];
rz(1.0709892) q[3];
sx q[3];
rz(-1.8355814) q[3];
sx q[3];
rz(1.7818835) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
