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
rz(0.48756227) q[0];
sx q[0];
rz(-0.3806448) q[0];
sx q[0];
rz(-3.0523377) q[0];
rz(-0.46269497) q[1];
sx q[1];
rz(-2.8609639) q[1];
sx q[1];
rz(1.7791003) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2089522) q[0];
sx q[0];
rz(-2.7679043) q[0];
sx q[0];
rz(-0.057500827) q[0];
rz(-pi) q[1];
rz(1.6073741) q[2];
sx q[2];
rz(-1.4223411) q[2];
sx q[2];
rz(-2.0853993) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.67586865) q[1];
sx q[1];
rz(-1.3480519) q[1];
sx q[1];
rz(-0.59091076) q[1];
rz(-pi) q[2];
rz(-2.9611601) q[3];
sx q[3];
rz(-0.46262925) q[3];
sx q[3];
rz(-0.53354665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.51916277) q[2];
sx q[2];
rz(-2.3795655) q[2];
sx q[2];
rz(2.2117174) q[2];
rz(0.30396384) q[3];
sx q[3];
rz(-1.8094939) q[3];
sx q[3];
rz(-2.9728594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.25794432) q[0];
sx q[0];
rz(-1.5644263) q[0];
sx q[0];
rz(-1.2516578) q[0];
rz(-2.9257863) q[1];
sx q[1];
rz(-2.1025175) q[1];
sx q[1];
rz(-0.13914093) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5128327) q[0];
sx q[0];
rz(-0.91673512) q[0];
sx q[0];
rz(-0.20166986) q[0];
x q[1];
rz(1.2982803) q[2];
sx q[2];
rz(-2.0534229) q[2];
sx q[2];
rz(-2.9543205) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.5522668) q[1];
sx q[1];
rz(-1.3331474) q[1];
sx q[1];
rz(2.4219805) q[1];
x q[2];
rz(0.71897935) q[3];
sx q[3];
rz(-0.96299473) q[3];
sx q[3];
rz(1.3083576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.4591878) q[2];
sx q[2];
rz(-0.46899691) q[2];
sx q[2];
rz(-2.0501308) q[2];
rz(-1.546321) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(-0.44927621) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0794373) q[0];
sx q[0];
rz(-1.9154444) q[0];
sx q[0];
rz(-0.011818258) q[0];
rz(-2.2305409) q[1];
sx q[1];
rz(-1.751519) q[1];
sx q[1];
rz(-1.3284838) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9244316) q[0];
sx q[0];
rz(-1.7329441) q[0];
sx q[0];
rz(0.52442201) q[0];
rz(0.53880589) q[2];
sx q[2];
rz(-0.22141117) q[2];
sx q[2];
rz(0.087269727) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6360259) q[1];
sx q[1];
rz(-0.88249245) q[1];
sx q[1];
rz(2.0783552) q[1];
rz(-pi) q[2];
rz(0.073224501) q[3];
sx q[3];
rz(-1.4595881) q[3];
sx q[3];
rz(-1.2014845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.38982424) q[2];
sx q[2];
rz(-1.5247034) q[2];
sx q[2];
rz(0.41979182) q[2];
rz(-1.7828364) q[3];
sx q[3];
rz(-0.32521453) q[3];
sx q[3];
rz(1.9512008) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.457394) q[0];
sx q[0];
rz(-1.1379108) q[0];
sx q[0];
rz(0.5140636) q[0];
rz(2.7534516) q[1];
sx q[1];
rz(-0.61182794) q[1];
sx q[1];
rz(-1.2302037) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3183751) q[0];
sx q[0];
rz(-1.8288695) q[0];
sx q[0];
rz(0.31320546) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.72429652) q[2];
sx q[2];
rz(-1.8593374) q[2];
sx q[2];
rz(-2.9674825) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.99707264) q[1];
sx q[1];
rz(-1.9341661) q[1];
sx q[1];
rz(2.4097869) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3421273) q[3];
sx q[3];
rz(-2.2284751) q[3];
sx q[3];
rz(-0.87008892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6396883) q[2];
sx q[2];
rz(-2.1268667) q[2];
sx q[2];
rz(0.91628966) q[2];
rz(-2.6472951) q[3];
sx q[3];
rz(-1.9435792) q[3];
sx q[3];
rz(0.77427197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2275527) q[0];
sx q[0];
rz(-2.1844449) q[0];
sx q[0];
rz(2.6555632) q[0];
rz(-2.5388429) q[1];
sx q[1];
rz(-1.9662247) q[1];
sx q[1];
rz(-1.1154307) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1971216) q[0];
sx q[0];
rz(-0.96707464) q[0];
sx q[0];
rz(-2.05462) q[0];
rz(-0.40256315) q[2];
sx q[2];
rz(-2.3057115) q[2];
sx q[2];
rz(2.1083567) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.045687606) q[1];
sx q[1];
rz(-2.3153911) q[1];
sx q[1];
rz(-0.5788486) q[1];
rz(-1.0904457) q[3];
sx q[3];
rz(-1.8060038) q[3];
sx q[3];
rz(3.1344828) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.1439765) q[2];
sx q[2];
rz(-2.82085) q[2];
sx q[2];
rz(1.0616659) q[2];
rz(1.5205421) q[3];
sx q[3];
rz(-1.942626) q[3];
sx q[3];
rz(-2.2475713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0735556) q[0];
sx q[0];
rz(-3.1017922) q[0];
sx q[0];
rz(-1.2233446) q[0];
rz(-2.6262737) q[1];
sx q[1];
rz(-0.66910187) q[1];
sx q[1];
rz(1.7759148) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13785744) q[0];
sx q[0];
rz(-0.23085871) q[0];
sx q[0];
rz(-2.0876711) q[0];
rz(-pi) q[1];
x q[1];
rz(0.093528259) q[2];
sx q[2];
rz(-1.3892738) q[2];
sx q[2];
rz(-0.43147555) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.9881043) q[1];
sx q[1];
rz(-1.0968913) q[1];
sx q[1];
rz(1.4010282) q[1];
x q[2];
rz(-0.11752085) q[3];
sx q[3];
rz(-0.29942808) q[3];
sx q[3];
rz(-2.260331) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.078865163) q[2];
sx q[2];
rz(-1.3389503) q[2];
sx q[2];
rz(-1.8143181) q[2];
rz(1.7257388) q[3];
sx q[3];
rz(-3.0684107) q[3];
sx q[3];
rz(2.2590526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0075204) q[0];
sx q[0];
rz(-0.53711689) q[0];
sx q[0];
rz(-0.19675955) q[0];
rz(-0.21601954) q[1];
sx q[1];
rz(-2.1327503) q[1];
sx q[1];
rz(0.78741995) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8926933) q[0];
sx q[0];
rz(-1.3266047) q[0];
sx q[0];
rz(-2.9466183) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.51047275) q[2];
sx q[2];
rz(-1.831448) q[2];
sx q[2];
rz(0.54768054) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.52975028) q[1];
sx q[1];
rz(-2.5482436) q[1];
sx q[1];
rz(0.081562138) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.18936746) q[3];
sx q[3];
rz(-1.8872617) q[3];
sx q[3];
rz(0.53839436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.6009723) q[2];
sx q[2];
rz(-1.5121907) q[2];
sx q[2];
rz(-0.38360325) q[2];
rz(0.90841928) q[3];
sx q[3];
rz(-2.3288265) q[3];
sx q[3];
rz(0.15373357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.28021321) q[0];
sx q[0];
rz(-2.6160243) q[0];
sx q[0];
rz(2.5049765) q[0];
rz(2.5909297) q[1];
sx q[1];
rz(-1.6267585) q[1];
sx q[1];
rz(-0.47264636) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7028684) q[0];
sx q[0];
rz(-1.4819549) q[0];
sx q[0];
rz(2.4636561) q[0];
rz(-1.6238975) q[2];
sx q[2];
rz(-1.4776973) q[2];
sx q[2];
rz(-1.129221) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3797261) q[1];
sx q[1];
rz(-1.7032924) q[1];
sx q[1];
rz(-1.5195058) q[1];
rz(0.63153679) q[3];
sx q[3];
rz(-1.3270006) q[3];
sx q[3];
rz(-2.7840028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.4281281) q[2];
sx q[2];
rz(-1.2423923) q[2];
sx q[2];
rz(-0.098527519) q[2];
rz(-3.1108917) q[3];
sx q[3];
rz(-1.5701598) q[3];
sx q[3];
rz(-0.15392412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-2.5310265) q[0];
sx q[0];
rz(-2.1838146) q[0];
sx q[0];
rz(-2.4089693) q[0];
rz(0.0064119617) q[1];
sx q[1];
rz(-1.0259722) q[1];
sx q[1];
rz(1.4220672) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32860562) q[0];
sx q[0];
rz(-1.3805132) q[0];
sx q[0];
rz(-1.3178131) q[0];
x q[1];
rz(-1.9344011) q[2];
sx q[2];
rz(-1.0997314) q[2];
sx q[2];
rz(0.94719145) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8796685) q[1];
sx q[1];
rz(-1.4170483) q[1];
sx q[1];
rz(1.6052482) q[1];
rz(-pi) q[2];
rz(-1.645817) q[3];
sx q[3];
rz(-2.0891694) q[3];
sx q[3];
rz(3.1278267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0942568) q[2];
sx q[2];
rz(-1.7862659) q[2];
sx q[2];
rz(2.6195841) q[2];
rz(-0.020307288) q[3];
sx q[3];
rz(-0.69965196) q[3];
sx q[3];
rz(0.29683963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8727528) q[0];
sx q[0];
rz(-1.5248542) q[0];
sx q[0];
rz(-1.1312477) q[0];
rz(-2.3616683) q[1];
sx q[1];
rz(-1.9386539) q[1];
sx q[1];
rz(-2.6663229) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49233046) q[0];
sx q[0];
rz(-1.2075675) q[0];
sx q[0];
rz(0.80357768) q[0];
rz(-pi) q[1];
rz(-0.23481253) q[2];
sx q[2];
rz(-2.2249892) q[2];
sx q[2];
rz(-1.7476817) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.1058029) q[1];
sx q[1];
rz(-2.8052108) q[1];
sx q[1];
rz(-0.85806429) q[1];
rz(-0.17367878) q[3];
sx q[3];
rz(-2.4376737) q[3];
sx q[3];
rz(-0.46797637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.37107006) q[2];
sx q[2];
rz(-1.5658242) q[2];
sx q[2];
rz(-1.2062262) q[2];
rz(1.7462339) q[3];
sx q[3];
rz(-2.5985056) q[3];
sx q[3];
rz(-0.079308184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4974925) q[0];
sx q[0];
rz(-0.19484367) q[0];
sx q[0];
rz(2.3051443) q[0];
rz(-2.1438228) q[1];
sx q[1];
rz(-1.0918959) q[1];
sx q[1];
rz(0.087654884) q[1];
rz(-3.072425) q[2];
sx q[2];
rz(-1.7887049) q[2];
sx q[2];
rz(-0.6520581) q[2];
rz(-3.0042778) q[3];
sx q[3];
rz(-0.68344322) q[3];
sx q[3];
rz(-2.4333011) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
