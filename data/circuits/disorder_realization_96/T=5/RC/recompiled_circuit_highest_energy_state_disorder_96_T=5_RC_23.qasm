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
rz(0.30492914) q[0];
sx q[0];
rz(-0.066216901) q[0];
sx q[0];
rz(0.15596341) q[0];
rz(-1.9744385) q[1];
sx q[1];
rz(-0.65216291) q[1];
sx q[1];
rz(-2.6551533) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32549511) q[0];
sx q[0];
rz(-0.98113546) q[0];
sx q[0];
rz(-2.4611121) q[0];
rz(-0.42066853) q[2];
sx q[2];
rz(-2.0528057) q[2];
sx q[2];
rz(-1.9276305) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9984879) q[1];
sx q[1];
rz(-1.2799026) q[1];
sx q[1];
rz(2.8105633) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3399038) q[3];
sx q[3];
rz(-2.3844815) q[3];
sx q[3];
rz(0.40386236) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8638986) q[2];
sx q[2];
rz(-0.70964491) q[2];
sx q[2];
rz(-2.7126183) q[2];
rz(-3.0947558) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(-2.528791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2752537) q[0];
sx q[0];
rz(-2.1493122) q[0];
sx q[0];
rz(2.476165) q[0];
rz(1.1412507) q[1];
sx q[1];
rz(-1.7612061) q[1];
sx q[1];
rz(0.64249396) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8905971) q[0];
sx q[0];
rz(-0.83775508) q[0];
sx q[0];
rz(-1.6392073) q[0];
rz(2.7879509) q[2];
sx q[2];
rz(-2.414915) q[2];
sx q[2];
rz(-0.35219372) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4027109) q[1];
sx q[1];
rz(-2.4428574) q[1];
sx q[1];
rz(0.043234392) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1118494) q[3];
sx q[3];
rz(-1.7223889) q[3];
sx q[3];
rz(0.21546396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.82681727) q[2];
sx q[2];
rz(-0.78709698) q[2];
sx q[2];
rz(0.55244201) q[2];
rz(-1.6636482) q[3];
sx q[3];
rz(-1.843957) q[3];
sx q[3];
rz(-0.67599952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33573547) q[0];
sx q[0];
rz(-0.72000802) q[0];
sx q[0];
rz(-0.82255256) q[0];
rz(0.81126732) q[1];
sx q[1];
rz(-2.8370116) q[1];
sx q[1];
rz(-1.4623581) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84483355) q[0];
sx q[0];
rz(-0.61723304) q[0];
sx q[0];
rz(-2.4485641) q[0];
x q[1];
rz(0.41303582) q[2];
sx q[2];
rz(-1.5025284) q[2];
sx q[2];
rz(-2.3562246) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.18204) q[1];
sx q[1];
rz(-0.55110747) q[1];
sx q[1];
rz(2.1564547) q[1];
x q[2];
rz(-0.7861761) q[3];
sx q[3];
rz(-2.7179931) q[3];
sx q[3];
rz(-0.15555233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.86639577) q[2];
sx q[2];
rz(-2.3123645) q[2];
sx q[2];
rz(0.59494507) q[2];
rz(-2.7368937) q[3];
sx q[3];
rz(-2.0345104) q[3];
sx q[3];
rz(1.2987202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73612708) q[0];
sx q[0];
rz(-1.8647702) q[0];
sx q[0];
rz(0.24293105) q[0];
rz(2.2566707) q[1];
sx q[1];
rz(-0.72598571) q[1];
sx q[1];
rz(-0.0028217908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76317518) q[0];
sx q[0];
rz(-2.2522914) q[0];
sx q[0];
rz(0.50294089) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6153938) q[2];
sx q[2];
rz(-1.0164653) q[2];
sx q[2];
rz(2.2458011) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2828212) q[1];
sx q[1];
rz(-2.438218) q[1];
sx q[1];
rz(1.4618631) q[1];
rz(-2.3798928) q[3];
sx q[3];
rz(-2.3683067) q[3];
sx q[3];
rz(-2.9536216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.60488492) q[2];
sx q[2];
rz(-1.1515836) q[2];
sx q[2];
rz(-2.4082157) q[2];
rz(-0.65507656) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(-0.56183279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4578399) q[0];
sx q[0];
rz(-0.34128749) q[0];
sx q[0];
rz(0.14232464) q[0];
rz(1.7798452) q[1];
sx q[1];
rz(-0.35265499) q[1];
sx q[1];
rz(0.49837643) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67923421) q[0];
sx q[0];
rz(-0.7158044) q[0];
sx q[0];
rz(1.4767893) q[0];
x q[1];
rz(-0.75888486) q[2];
sx q[2];
rz(-1.5336131) q[2];
sx q[2];
rz(2.5868724) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.26311647) q[1];
sx q[1];
rz(-2.2848087) q[1];
sx q[1];
rz(-1.962838) q[1];
x q[2];
rz(2.8041995) q[3];
sx q[3];
rz(-0.80885799) q[3];
sx q[3];
rz(2.0615426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0087697) q[2];
sx q[2];
rz(-1.885773) q[2];
sx q[2];
rz(-2.447017) q[2];
rz(2.3092367) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(-2.8913403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4271127) q[0];
sx q[0];
rz(-0.50943333) q[0];
sx q[0];
rz(2.4795649) q[0];
rz(1.9561249) q[1];
sx q[1];
rz(-1.0292425) q[1];
sx q[1];
rz(1.1659291) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87905721) q[0];
sx q[0];
rz(-1.3692807) q[0];
sx q[0];
rz(2.461947) q[0];
rz(-1.7244965) q[2];
sx q[2];
rz(-0.96818334) q[2];
sx q[2];
rz(3.1069047) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7699189) q[1];
sx q[1];
rz(-0.087962463) q[1];
sx q[1];
rz(0.96925737) q[1];
rz(2.5218368) q[3];
sx q[3];
rz(-2.0888076) q[3];
sx q[3];
rz(-0.56700804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.379443) q[2];
sx q[2];
rz(-1.2078441) q[2];
sx q[2];
rz(-0.70972788) q[2];
rz(0.48192561) q[3];
sx q[3];
rz(-0.47526264) q[3];
sx q[3];
rz(-0.019439241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4703281) q[0];
sx q[0];
rz(-0.98099357) q[0];
sx q[0];
rz(-1.574466) q[0];
rz(-1.4944685) q[1];
sx q[1];
rz(-2.6946805) q[1];
sx q[1];
rz(-2.2692197) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0013638) q[0];
sx q[0];
rz(-1.2622381) q[0];
sx q[0];
rz(-0.17805992) q[0];
rz(-pi) q[1];
rz(2.9296604) q[2];
sx q[2];
rz(-0.733558) q[2];
sx q[2];
rz(-0.19679815) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.7561314) q[1];
sx q[1];
rz(-1.4646936) q[1];
sx q[1];
rz(3.1110292) q[1];
rz(-pi) q[2];
rz(-2.2254785) q[3];
sx q[3];
rz(-1.5903842) q[3];
sx q[3];
rz(0.71924984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.81277728) q[2];
sx q[2];
rz(-0.8605364) q[2];
sx q[2];
rz(-2.013773) q[2];
rz(-2.7382216) q[3];
sx q[3];
rz(-1.4453459) q[3];
sx q[3];
rz(-2.4283714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9047852) q[0];
sx q[0];
rz(-1.1431563) q[0];
sx q[0];
rz(1.8444201) q[0];
rz(-0.5212658) q[1];
sx q[1];
rz(-1.8788985) q[1];
sx q[1];
rz(-3.0788132) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32048827) q[0];
sx q[0];
rz(-2.0135499) q[0];
sx q[0];
rz(-0.073320079) q[0];
rz(-pi) q[1];
rz(1.3387483) q[2];
sx q[2];
rz(-1.9221483) q[2];
sx q[2];
rz(0.50594375) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.52555841) q[1];
sx q[1];
rz(-1.8174606) q[1];
sx q[1];
rz(-0.14825578) q[1];
x q[2];
rz(-2.6046329) q[3];
sx q[3];
rz(-2.1459422) q[3];
sx q[3];
rz(-1.6292269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.89821833) q[2];
sx q[2];
rz(-2.2369907) q[2];
sx q[2];
rz(-0.96837366) q[2];
rz(-0.51356703) q[3];
sx q[3];
rz(-1.3526724) q[3];
sx q[3];
rz(-2.8106522) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6043855) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(-2.3576417) q[0];
rz(-1.9369269) q[1];
sx q[1];
rz(-0.37310633) q[1];
sx q[1];
rz(1.7663667) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88343745) q[0];
sx q[0];
rz(-2.7984507) q[0];
sx q[0];
rz(2.457042) q[0];
rz(1.398574) q[2];
sx q[2];
rz(-0.15744124) q[2];
sx q[2];
rz(-2.1724608) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8158477) q[1];
sx q[1];
rz(-2.1195852) q[1];
sx q[1];
rz(0.71760787) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0635942) q[3];
sx q[3];
rz(-0.35420277) q[3];
sx q[3];
rz(-2.2627047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.84460008) q[2];
sx q[2];
rz(-0.30868369) q[2];
sx q[2];
rz(-2.6958579) q[2];
rz(2.5388057) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(-2.2980237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
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
rz(-0.26850253) q[0];
sx q[0];
rz(-0.52274811) q[0];
sx q[0];
rz(0.52421808) q[0];
rz(-1.2007319) q[1];
sx q[1];
rz(-1.2501161) q[1];
sx q[1];
rz(1.1842309) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4864593) q[0];
sx q[0];
rz(-0.39748991) q[0];
sx q[0];
rz(0.50661202) q[0];
rz(-pi) q[1];
rz(0.068151926) q[2];
sx q[2];
rz(-1.8859204) q[2];
sx q[2];
rz(-1.5308876) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6683015) q[1];
sx q[1];
rz(-1.4379825) q[1];
sx q[1];
rz(-1.5584438) q[1];
rz(-pi) q[2];
rz(2.5126825) q[3];
sx q[3];
rz(-1.5987608) q[3];
sx q[3];
rz(0.6599676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.9253917) q[2];
sx q[2];
rz(-0.4500469) q[2];
sx q[2];
rz(1.0518543) q[2];
rz(0.22107407) q[3];
sx q[3];
rz(-1.3007921) q[3];
sx q[3];
rz(-2.2033447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
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
rz(-1.5939519) q[0];
sx q[0];
rz(-1.5008391) q[0];
sx q[0];
rz(0.048901625) q[0];
rz(-2.7652057) q[1];
sx q[1];
rz(-0.86224894) q[1];
sx q[1];
rz(-1.1663762) q[1];
rz(-0.57009956) q[2];
sx q[2];
rz(-1.5975633) q[2];
sx q[2];
rz(-0.057889197) q[2];
rz(-0.71394271) q[3];
sx q[3];
rz(-1.9458013) q[3];
sx q[3];
rz(-1.8251606) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
