OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.055846) q[0];
sx q[0];
rz(-3.0598109) q[0];
sx q[0];
rz(2.6401289) q[0];
rz(-1.6429098) q[1];
sx q[1];
rz(-0.39615762) q[1];
sx q[1];
rz(-2.8191541) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0816406) q[0];
sx q[0];
rz(-1.8917221) q[0];
sx q[0];
rz(3.0517464) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85876043) q[2];
sx q[2];
rz(-0.81958629) q[2];
sx q[2];
rz(-2.8263) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4336366) q[1];
sx q[1];
rz(-0.70322137) q[1];
sx q[1];
rz(2.1232848) q[1];
x q[2];
rz(2.1665855) q[3];
sx q[3];
rz(-0.72837225) q[3];
sx q[3];
rz(-1.4260074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.50513187) q[2];
sx q[2];
rz(-2.5487066) q[2];
sx q[2];
rz(-0.55603975) q[2];
rz(2.3089144) q[3];
sx q[3];
rz(-1.6502389) q[3];
sx q[3];
rz(2.1957943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44822025) q[0];
sx q[0];
rz(-1.6813261) q[0];
sx q[0];
rz(-2.9843176) q[0];
rz(2.8804624) q[1];
sx q[1];
rz(-1.7938679) q[1];
sx q[1];
rz(-3.0325586) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5599247) q[0];
sx q[0];
rz(-1.7413057) q[0];
sx q[0];
rz(1.3397564) q[0];
rz(-1.3729587) q[2];
sx q[2];
rz(-1.3126144) q[2];
sx q[2];
rz(-0.82222647) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1034531) q[1];
sx q[1];
rz(-0.98587576) q[1];
sx q[1];
rz(-1.000688) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3974959) q[3];
sx q[3];
rz(-1.0565851) q[3];
sx q[3];
rz(1.7931995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9033501) q[2];
sx q[2];
rz(-1.976333) q[2];
sx q[2];
rz(-1.8781352) q[2];
rz(-0.3271099) q[3];
sx q[3];
rz(-1.5644904) q[3];
sx q[3];
rz(-1.9272778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.0771714) q[0];
sx q[0];
rz(-3.0922958) q[0];
sx q[0];
rz(-1.3431312) q[0];
rz(0.24761565) q[1];
sx q[1];
rz(-0.74664465) q[1];
sx q[1];
rz(-0.48167357) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7307229) q[0];
sx q[0];
rz(-2.0559089) q[0];
sx q[0];
rz(-2.7184125) q[0];
rz(-2.8086497) q[2];
sx q[2];
rz(-2.1579086) q[2];
sx q[2];
rz(-1.3981896) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.18322769) q[1];
sx q[1];
rz(-1.0291568) q[1];
sx q[1];
rz(0.72802131) q[1];
rz(-pi) q[2];
rz(-2.5211469) q[3];
sx q[3];
rz(-0.9471604) q[3];
sx q[3];
rz(-2.7321531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.3383011) q[2];
sx q[2];
rz(-2.3241966) q[2];
sx q[2];
rz(-2.6417007) q[2];
rz(2.5806184) q[3];
sx q[3];
rz(-1.8818972) q[3];
sx q[3];
rz(1.4612173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9445779) q[0];
sx q[0];
rz(-2.975583) q[0];
sx q[0];
rz(-2.5894077) q[0];
rz(-1.5532956) q[1];
sx q[1];
rz(-0.89893666) q[1];
sx q[1];
rz(1.8968556) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9075605) q[0];
sx q[0];
rz(-1.2313156) q[0];
sx q[0];
rz(1.5725122) q[0];
rz(-pi) q[1];
x q[1];
rz(0.99889465) q[2];
sx q[2];
rz(-0.23795393) q[2];
sx q[2];
rz(1.9902802) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6914122) q[1];
sx q[1];
rz(-1.0346864) q[1];
sx q[1];
rz(2.8944573) q[1];
rz(-pi) q[2];
rz(2.5707385) q[3];
sx q[3];
rz(-1.7613162) q[3];
sx q[3];
rz(-1.966147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84919471) q[2];
sx q[2];
rz(-1.2595824) q[2];
sx q[2];
rz(-1.1506895) q[2];
rz(1.6644647) q[3];
sx q[3];
rz(-1.5090347) q[3];
sx q[3];
rz(2.6586444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0634336) q[0];
sx q[0];
rz(-0.76197356) q[0];
sx q[0];
rz(-3.0601236) q[0];
rz(-3.07913) q[1];
sx q[1];
rz(-2.000258) q[1];
sx q[1];
rz(1.5030456) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8836356) q[0];
sx q[0];
rz(-0.99146087) q[0];
sx q[0];
rz(0.15276076) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91707768) q[2];
sx q[2];
rz(-3.0055025) q[2];
sx q[2];
rz(1.2166785) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.690783) q[1];
sx q[1];
rz(-0.329031) q[1];
sx q[1];
rz(-1.8915218) q[1];
rz(1.8893858) q[3];
sx q[3];
rz(-1.0831523) q[3];
sx q[3];
rz(2.3209751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.5082671) q[2];
sx q[2];
rz(-2.1990364) q[2];
sx q[2];
rz(1.903669) q[2];
rz(-2.0189019) q[3];
sx q[3];
rz(-2.4653547) q[3];
sx q[3];
rz(2.6200263) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6102819) q[0];
sx q[0];
rz(-2.1370482) q[0];
sx q[0];
rz(-2.8748728) q[0];
rz(-2.5807014) q[1];
sx q[1];
rz(-1.8436878) q[1];
sx q[1];
rz(0.7985324) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5312936) q[0];
sx q[0];
rz(-1.4302505) q[0];
sx q[0];
rz(1.848624) q[0];
rz(-pi) q[1];
rz(1.2217667) q[2];
sx q[2];
rz(-2.0418228) q[2];
sx q[2];
rz(0.40086056) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.16520271) q[1];
sx q[1];
rz(-2.5678647) q[1];
sx q[1];
rz(1.4742673) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1875601) q[3];
sx q[3];
rz(-0.87152374) q[3];
sx q[3];
rz(-1.9642533) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9810527) q[2];
sx q[2];
rz(-1.8926228) q[2];
sx q[2];
rz(-0.36671656) q[2];
rz(-1.8803053) q[3];
sx q[3];
rz(-1.6882608) q[3];
sx q[3];
rz(0.096207531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40903184) q[0];
sx q[0];
rz(-2.2213187) q[0];
sx q[0];
rz(-0.60638705) q[0];
rz(-0.19730332) q[1];
sx q[1];
rz(-2.0154672) q[1];
sx q[1];
rz(-0.46404776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89550971) q[0];
sx q[0];
rz(-1.9946788) q[0];
sx q[0];
rz(-2.258582) q[0];
rz(-pi) q[1];
rz(3.1214141) q[2];
sx q[2];
rz(-1.3300465) q[2];
sx q[2];
rz(-1.9764331) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0276427) q[1];
sx q[1];
rz(-2.0538035) q[1];
sx q[1];
rz(-2.7999858) q[1];
rz(-pi) q[2];
rz(1.4887605) q[3];
sx q[3];
rz(-1.335252) q[3];
sx q[3];
rz(-1.1788648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.93418926) q[2];
sx q[2];
rz(-1.0031909) q[2];
sx q[2];
rz(2.8835473) q[2];
rz(-1.9559654) q[3];
sx q[3];
rz(-1.5153171) q[3];
sx q[3];
rz(-0.0035088249) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.946452) q[0];
sx q[0];
rz(-1.8608681) q[0];
sx q[0];
rz(0.38129693) q[0];
rz(0.095245846) q[1];
sx q[1];
rz(-0.97243273) q[1];
sx q[1];
rz(1.4415178) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9593175) q[0];
sx q[0];
rz(-1.5556941) q[0];
sx q[0];
rz(0.063477593) q[0];
rz(-pi) q[1];
rz(-1.6091789) q[2];
sx q[2];
rz(-0.95853171) q[2];
sx q[2];
rz(1.6910764) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4691094) q[1];
sx q[1];
rz(-1.4454578) q[1];
sx q[1];
rz(-0.13087665) q[1];
rz(1.5278682) q[3];
sx q[3];
rz(-1.2371484) q[3];
sx q[3];
rz(1.4294525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1429446) q[2];
sx q[2];
rz(-2.7286077) q[2];
sx q[2];
rz(-2.9150035) q[2];
rz(0.44858027) q[3];
sx q[3];
rz(-1.5357163) q[3];
sx q[3];
rz(2.3118238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-1.3806234) q[0];
sx q[0];
rz(-0.7614823) q[0];
sx q[0];
rz(-1.7425849) q[0];
rz(0.31708583) q[1];
sx q[1];
rz(-1.6665019) q[1];
sx q[1];
rz(-2.1549966) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87514585) q[0];
sx q[0];
rz(-2.6212924) q[0];
sx q[0];
rz(0.32169028) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.87395845) q[2];
sx q[2];
rz(-1.9947589) q[2];
sx q[2];
rz(-1.6800113) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7762737) q[1];
sx q[1];
rz(-2.1734108) q[1];
sx q[1];
rz(-1.3990632) q[1];
rz(-pi) q[2];
rz(-1.1499693) q[3];
sx q[3];
rz(-1.7289366) q[3];
sx q[3];
rz(1.4827673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.2150779) q[2];
sx q[2];
rz(-0.72700095) q[2];
sx q[2];
rz(-2.731936) q[2];
rz(-2.8783197) q[3];
sx q[3];
rz(-1.3116838) q[3];
sx q[3];
rz(-0.51945654) q[3];
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
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39919329) q[0];
sx q[0];
rz(-0.078646794) q[0];
sx q[0];
rz(-1.7364527) q[0];
rz(2.3204904) q[1];
sx q[1];
rz(-0.91870538) q[1];
sx q[1];
rz(1.4155037) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4840568) q[0];
sx q[0];
rz(-2.7298379) q[0];
sx q[0];
rz(2.5176237) q[0];
x q[1];
rz(1.9344994) q[2];
sx q[2];
rz(-1.8038245) q[2];
sx q[2];
rz(-1.9050913) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0701329) q[1];
sx q[1];
rz(-1.3576641) q[1];
sx q[1];
rz(0.15503426) q[1];
x q[2];
rz(1.0584162) q[3];
sx q[3];
rz(-1.4231764) q[3];
sx q[3];
rz(-0.83422134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4225509) q[2];
sx q[2];
rz(-2.8406403) q[2];
sx q[2];
rz(0.12410513) q[2];
rz(2.1758046) q[3];
sx q[3];
rz(-1.4744759) q[3];
sx q[3];
rz(0.66108274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.538095) q[0];
sx q[0];
rz(-2.8932543) q[0];
sx q[0];
rz(2.2809991) q[0];
rz(2.8339236) q[1];
sx q[1];
rz(-1.2528906) q[1];
sx q[1];
rz(1.2045592) q[1];
rz(-1.3901426) q[2];
sx q[2];
rz(-1.8271108) q[2];
sx q[2];
rz(1.9432632) q[2];
rz(-2.9028805) q[3];
sx q[3];
rz(-2.5720027) q[3];
sx q[3];
rz(-3.0726074) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
