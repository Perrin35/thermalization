OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7933529) q[0];
sx q[0];
rz(-1.5433595) q[0];
sx q[0];
rz(1.5016851) q[0];
rz(-1.3950672) q[1];
sx q[1];
rz(2.7434064) q[1];
sx q[1];
rz(12.413496) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22963755) q[0];
sx q[0];
rz(-2.9792333) q[0];
sx q[0];
rz(2.1549757) q[0];
rz(-pi) q[1];
rz(0.42793226) q[2];
sx q[2];
rz(-0.38201354) q[2];
sx q[2];
rz(-2.0150507) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.0331206) q[1];
sx q[1];
rz(-1.7594928) q[1];
sx q[1];
rz(-0.061339247) q[1];
rz(-2.128162) q[3];
sx q[3];
rz(-2.3881222) q[3];
sx q[3];
rz(2.4581152) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.07831002) q[2];
sx q[2];
rz(-1.4274884) q[2];
sx q[2];
rz(2.5417292) q[2];
rz(-2.6317224) q[3];
sx q[3];
rz(-0.32250914) q[3];
sx q[3];
rz(-2.9738284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.584499) q[0];
sx q[0];
rz(-0.85894132) q[0];
sx q[0];
rz(2.9993045) q[0];
rz(-1.2917057) q[1];
sx q[1];
rz(-2.6085491) q[1];
sx q[1];
rz(2.5407015) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93705015) q[0];
sx q[0];
rz(-1.1536404) q[0];
sx q[0];
rz(-2.8931505) q[0];
x q[1];
rz(-0.63670701) q[2];
sx q[2];
rz(-1.8786624) q[2];
sx q[2];
rz(-0.40154469) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8411202) q[1];
sx q[1];
rz(-2.4266647) q[1];
sx q[1];
rz(0.27473533) q[1];
x q[2];
rz(0.51436794) q[3];
sx q[3];
rz(-1.4050438) q[3];
sx q[3];
rz(2.7837798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.94704023) q[2];
sx q[2];
rz(-1.9779454) q[2];
sx q[2];
rz(0.8645424) q[2];
rz(0.38327992) q[3];
sx q[3];
rz(-2.7195103) q[3];
sx q[3];
rz(0.88467902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8466723) q[0];
sx q[0];
rz(-0.28691322) q[0];
sx q[0];
rz(0.82557803) q[0];
rz(-2.3041252) q[1];
sx q[1];
rz(-2.1223964) q[1];
sx q[1];
rz(0.79140633) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.74293711) q[0];
sx q[0];
rz(-2.3283132) q[0];
sx q[0];
rz(0.11095993) q[0];
rz(-1.7010274) q[2];
sx q[2];
rz(-1.4716867) q[2];
sx q[2];
rz(-2.2203022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7792862) q[1];
sx q[1];
rz(-2.0842881) q[1];
sx q[1];
rz(-1.1720285) q[1];
rz(-pi) q[2];
rz(-0.43413991) q[3];
sx q[3];
rz(-2.0159147) q[3];
sx q[3];
rz(-0.49830084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.5241663) q[2];
sx q[2];
rz(-2.7238621) q[2];
sx q[2];
rz(1.8903271) q[2];
rz(-1.9620126) q[3];
sx q[3];
rz(-1.1359295) q[3];
sx q[3];
rz(-2.484926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(0.11611045) q[0];
sx q[0];
rz(-1.2283093) q[0];
sx q[0];
rz(-2.3226698) q[0];
rz(0.081309155) q[1];
sx q[1];
rz(-2.5550877) q[1];
sx q[1];
rz(2.3611045) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1063447) q[0];
sx q[0];
rz(-2.7625933) q[0];
sx q[0];
rz(2.9826791) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.444818) q[2];
sx q[2];
rz(-0.65444817) q[2];
sx q[2];
rz(2.12754) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38950237) q[1];
sx q[1];
rz(-0.47359514) q[1];
sx q[1];
rz(0.69382967) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1673766) q[3];
sx q[3];
rz(-1.8634923) q[3];
sx q[3];
rz(-0.044805275) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9235733) q[2];
sx q[2];
rz(-1.6417445) q[2];
sx q[2];
rz(1.2884595) q[2];
rz(-1.4675379) q[3];
sx q[3];
rz(-0.54687423) q[3];
sx q[3];
rz(0.43681496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9913637) q[0];
sx q[0];
rz(-2.4910091) q[0];
sx q[0];
rz(2.9184166) q[0];
rz(-3.0617833) q[1];
sx q[1];
rz(-0.97750074) q[1];
sx q[1];
rz(-2.3056183) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8337934) q[0];
sx q[0];
rz(-2.6386441) q[0];
sx q[0];
rz(-1.0566684) q[0];
rz(-pi) q[1];
rz(-2.4673632) q[2];
sx q[2];
rz(-0.59665702) q[2];
sx q[2];
rz(2.3643189) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.32415307) q[1];
sx q[1];
rz(-2.0803703) q[1];
sx q[1];
rz(-0.47199366) q[1];
rz(-pi) q[2];
rz(-0.23578819) q[3];
sx q[3];
rz(-0.97321586) q[3];
sx q[3];
rz(-1.8862629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23433805) q[2];
sx q[2];
rz(-0.85504389) q[2];
sx q[2];
rz(-1.7503395) q[2];
rz(1.2587345) q[3];
sx q[3];
rz(-1.7573059) q[3];
sx q[3];
rz(-2.7425227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(1.2788972) q[0];
sx q[0];
rz(-2.6628222) q[0];
sx q[0];
rz(-0.71267772) q[0];
rz(2.0172987) q[1];
sx q[1];
rz(-1.5713888) q[1];
sx q[1];
rz(-1.0363151) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46902114) q[0];
sx q[0];
rz(-0.9378792) q[0];
sx q[0];
rz(1.4255037) q[0];
rz(-pi) q[1];
rz(-0.89547248) q[2];
sx q[2];
rz(-2.7717441) q[2];
sx q[2];
rz(2.0167882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.24806286) q[1];
sx q[1];
rz(-1.7772733) q[1];
sx q[1];
rz(2.1724387) q[1];
rz(-pi) q[2];
x q[2];
rz(0.18671496) q[3];
sx q[3];
rz(-2.2590504) q[3];
sx q[3];
rz(2.6783239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0816281) q[2];
sx q[2];
rz(-0.57264239) q[2];
sx q[2];
rz(0.51679483) q[2];
rz(1.146727) q[3];
sx q[3];
rz(-2.1490993) q[3];
sx q[3];
rz(-3.0920046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0741172) q[0];
sx q[0];
rz(-0.68454409) q[0];
sx q[0];
rz(1.462498) q[0];
rz(1.0985724) q[1];
sx q[1];
rz(-1.699828) q[1];
sx q[1];
rz(0.77902737) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1027101) q[0];
sx q[0];
rz(-1.5691367) q[0];
sx q[0];
rz(0.05430211) q[0];
rz(0.99188478) q[2];
sx q[2];
rz(-1.9702458) q[2];
sx q[2];
rz(-2.7766984) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.9866619) q[1];
sx q[1];
rz(-2.8720461) q[1];
sx q[1];
rz(0.69828548) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38052151) q[3];
sx q[3];
rz(-2.2448934) q[3];
sx q[3];
rz(3.1269011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.74759185) q[2];
sx q[2];
rz(-0.70351768) q[2];
sx q[2];
rz(-0.35959378) q[2];
rz(-0.30073419) q[3];
sx q[3];
rz(-2.3263003) q[3];
sx q[3];
rz(-0.99641478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4312298) q[0];
sx q[0];
rz(-1.9332941) q[0];
sx q[0];
rz(2.8408458) q[0];
rz(0.53384471) q[1];
sx q[1];
rz(-2.0678935) q[1];
sx q[1];
rz(-3.0278382) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9917002) q[0];
sx q[0];
rz(-1.0731369) q[0];
sx q[0];
rz(2.0927621) q[0];
rz(-3.0871307) q[2];
sx q[2];
rz(-0.60945933) q[2];
sx q[2];
rz(-0.0527339) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.1492474) q[1];
sx q[1];
rz(-2.2015044) q[1];
sx q[1];
rz(2.135732) q[1];
x q[2];
rz(-1.9451009) q[3];
sx q[3];
rz(-1.7262207) q[3];
sx q[3];
rz(-0.96159305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6834695) q[2];
sx q[2];
rz(-0.14294954) q[2];
sx q[2];
rz(-2.6911823) q[2];
rz(-2.2266375) q[3];
sx q[3];
rz(-0.92792845) q[3];
sx q[3];
rz(-1.8042867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27110711) q[0];
sx q[0];
rz(-2.4021689) q[0];
sx q[0];
rz(-2.7258605) q[0];
rz(0.43139002) q[1];
sx q[1];
rz(-0.58512551) q[1];
sx q[1];
rz(0.77692568) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3142648) q[0];
sx q[0];
rz(-1.5175455) q[0];
sx q[0];
rz(-3.1287249) q[0];
x q[1];
rz(-0.81468366) q[2];
sx q[2];
rz(-2.2534342) q[2];
sx q[2];
rz(-1.4614568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5239656) q[1];
sx q[1];
rz(-1.0203779) q[1];
sx q[1];
rz(-1.074076) q[1];
x q[2];
rz(-2.611972) q[3];
sx q[3];
rz(-2.6130015) q[3];
sx q[3];
rz(1.1666672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0100157) q[2];
sx q[2];
rz(-0.78993979) q[2];
sx q[2];
rz(-1.5124403) q[2];
rz(-0.43859628) q[3];
sx q[3];
rz(-0.83493835) q[3];
sx q[3];
rz(-0.51630539) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5877391) q[0];
sx q[0];
rz(-2.6022311) q[0];
sx q[0];
rz(1.3158276) q[0];
rz(3.059803) q[1];
sx q[1];
rz(-1.9606699) q[1];
sx q[1];
rz(1.2710424) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91453856) q[0];
sx q[0];
rz(-2.7167121) q[0];
sx q[0];
rz(-1.4999313) q[0];
rz(-pi) q[1];
rz(-3.0767303) q[2];
sx q[2];
rz(-2.2998875) q[2];
sx q[2];
rz(-0.87769485) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.20746729) q[1];
sx q[1];
rz(-2.3015056) q[1];
sx q[1];
rz(-0.39189696) q[1];
x q[2];
rz(1.5457235) q[3];
sx q[3];
rz(-1.6524602) q[3];
sx q[3];
rz(-1.9586399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3424993) q[2];
sx q[2];
rz(-1.4405595) q[2];
sx q[2];
rz(-1.052617) q[2];
rz(-1.8140225) q[3];
sx q[3];
rz(-1.9460287) q[3];
sx q[3];
rz(-2.6322406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7256182) q[0];
sx q[0];
rz(-1.7617891) q[0];
sx q[0];
rz(-1.2431385) q[0];
rz(1.3432518) q[1];
sx q[1];
rz(-2.6814798) q[1];
sx q[1];
rz(0.36569256) q[1];
rz(-1.2058291) q[2];
sx q[2];
rz(-0.50508026) q[2];
sx q[2];
rz(2.6618367) q[2];
rz(2.7300077) q[3];
sx q[3];
rz(-1.5829908) q[3];
sx q[3];
rz(-0.79394059) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
