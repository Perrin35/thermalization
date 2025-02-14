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
rz(3.0753758) q[0];
sx q[0];
rz(9.2688146) q[0];
rz(-1.9744385) q[1];
sx q[1];
rz(-0.65216291) q[1];
sx q[1];
rz(-2.6551533) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8160975) q[0];
sx q[0];
rz(-0.98113546) q[0];
sx q[0];
rz(-0.6804806) q[0];
x q[1];
rz(-2.0912284) q[2];
sx q[2];
rz(-1.2005521) q[2];
sx q[2];
rz(2.5802719) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6156494) q[1];
sx q[1];
rz(-1.8874223) q[1];
sx q[1];
rz(-1.2642045) q[1];
rz(2.9286195) q[3];
sx q[3];
rz(-0.83847441) q[3];
sx q[3];
rz(2.4249083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.277694) q[2];
sx q[2];
rz(-2.4319477) q[2];
sx q[2];
rz(-0.42897439) q[2];
rz(-0.046836827) q[3];
sx q[3];
rz(-2.7497079) q[3];
sx q[3];
rz(-0.61280167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2752537) q[0];
sx q[0];
rz(-0.99228042) q[0];
sx q[0];
rz(2.476165) q[0];
rz(-2.0003419) q[1];
sx q[1];
rz(-1.7612061) q[1];
sx q[1];
rz(0.64249396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8905971) q[0];
sx q[0];
rz(-2.3038376) q[0];
sx q[0];
rz(1.6392073) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69509721) q[2];
sx q[2];
rz(-1.3386269) q[2];
sx q[2];
rz(0.94939453) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.86502581) q[1];
sx q[1];
rz(-1.5986018) q[1];
sx q[1];
rz(-0.69827484) q[1];
x q[2];
rz(-2.965224) q[3];
sx q[3];
rz(-2.1049728) q[3];
sx q[3];
rz(-1.4458223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3147754) q[2];
sx q[2];
rz(-0.78709698) q[2];
sx q[2];
rz(0.55244201) q[2];
rz(1.4779444) q[3];
sx q[3];
rz(-1.843957) q[3];
sx q[3];
rz(-0.67599952) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
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
rz(0.33573547) q[0];
sx q[0];
rz(-0.72000802) q[0];
sx q[0];
rz(2.3190401) q[0];
rz(0.81126732) q[1];
sx q[1];
rz(-0.30458105) q[1];
sx q[1];
rz(1.4623581) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3212292) q[0];
sx q[0];
rz(-1.1920394) q[0];
sx q[0];
rz(2.6418153) q[0];
rz(-pi) q[1];
rz(0.16872318) q[2];
sx q[2];
rz(-0.41831917) q[2];
sx q[2];
rz(2.5105798) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87440675) q[1];
sx q[1];
rz(-1.864434) q[1];
sx q[1];
rz(-2.044157) q[1];
rz(2.8331746) q[3];
sx q[3];
rz(-1.2756516) q[3];
sx q[3];
rz(-0.67527387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2751969) q[2];
sx q[2];
rz(-2.3123645) q[2];
sx q[2];
rz(-2.5466476) q[2];
rz(-0.40469894) q[3];
sx q[3];
rz(-2.0345104) q[3];
sx q[3];
rz(-1.2987202) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4054656) q[0];
sx q[0];
rz(-1.8647702) q[0];
sx q[0];
rz(2.8986616) q[0];
rz(-2.2566707) q[1];
sx q[1];
rz(-2.4156069) q[1];
sx q[1];
rz(-0.0028217908) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76317518) q[0];
sx q[0];
rz(-0.88930128) q[0];
sx q[0];
rz(-2.6386518) q[0];
x q[1];
rz(2.6153938) q[2];
sx q[2];
rz(-2.1251273) q[2];
sx q[2];
rz(-0.89579158) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.62880781) q[1];
sx q[1];
rz(-1.6411726) q[1];
sx q[1];
rz(-2.2712399) q[1];
rz(-pi) q[2];
rz(-2.5266227) q[3];
sx q[3];
rz(-1.0677862) q[3];
sx q[3];
rz(0.78395432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.60488492) q[2];
sx q[2];
rz(-1.1515836) q[2];
sx q[2];
rz(2.4082157) q[2];
rz(0.65507656) q[3];
sx q[3];
rz(-2.8337182) q[3];
sx q[3];
rz(0.56183279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(2.4578399) q[0];
sx q[0];
rz(-2.8003052) q[0];
sx q[0];
rz(-2.999268) q[0];
rz(1.7798452) q[1];
sx q[1];
rz(-0.35265499) q[1];
sx q[1];
rz(-2.6432162) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55492586) q[0];
sx q[0];
rz(-2.2827671) q[0];
sx q[0];
rz(3.0601383) q[0];
x q[1];
rz(0.75888486) q[2];
sx q[2];
rz(-1.6079796) q[2];
sx q[2];
rz(-0.55472022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8784762) q[1];
sx q[1];
rz(-2.2848087) q[1];
sx q[1];
rz(-1.962838) q[1];
rz(-pi) q[2];
rz(1.2368535) q[3];
sx q[3];
rz(-2.32226) q[3];
sx q[3];
rz(-0.60988319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.132823) q[2];
sx q[2];
rz(-1.885773) q[2];
sx q[2];
rz(2.447017) q[2];
rz(-0.83235598) q[3];
sx q[3];
rz(-1.1135626) q[3];
sx q[3];
rz(-2.8913403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4271127) q[0];
sx q[0];
rz(-0.50943333) q[0];
sx q[0];
rz(0.66202778) q[0];
rz(-1.1854677) q[1];
sx q[1];
rz(-1.0292425) q[1];
sx q[1];
rz(1.1659291) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93454516) q[0];
sx q[0];
rz(-2.4372792) q[0];
sx q[0];
rz(2.8273316) q[0];
x q[1];
rz(2.922631) q[2];
sx q[2];
rz(-0.6195407) q[2];
sx q[2];
rz(-0.23212405) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.37167376) q[1];
sx q[1];
rz(-3.0536302) q[1];
sx q[1];
rz(-0.96925737) q[1];
rz(0.77591578) q[3];
sx q[3];
rz(-0.78523038) q[3];
sx q[3];
rz(-1.531158) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.379443) q[2];
sx q[2];
rz(-1.2078441) q[2];
sx q[2];
rz(0.70972788) q[2];
rz(-2.659667) q[3];
sx q[3];
rz(-0.47526264) q[3];
sx q[3];
rz(3.1221534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67126453) q[0];
sx q[0];
rz(-0.98099357) q[0];
sx q[0];
rz(-1.574466) q[0];
rz(1.4944685) q[1];
sx q[1];
rz(-2.6946805) q[1];
sx q[1];
rz(-0.87237298) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14022889) q[0];
sx q[0];
rz(-1.2622381) q[0];
sx q[0];
rz(-2.9635327) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4192737) q[2];
sx q[2];
rz(-1.4294942) q[2];
sx q[2];
rz(-1.6091122) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.0371767) q[1];
sx q[1];
rz(-3.0311916) q[1];
sx q[1];
rz(1.2913741) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.116901) q[3];
sx q[3];
rz(-0.91626142) q[3];
sx q[3];
rz(0.86658044) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3288154) q[2];
sx q[2];
rz(-0.8605364) q[2];
sx q[2];
rz(-2.013773) q[2];
rz(-0.40337107) q[3];
sx q[3];
rz(-1.6962467) q[3];
sx q[3];
rz(0.71322125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.990898) q[0];
sx q[0];
rz(-2.6932062) q[0];
sx q[0];
rz(-1.4175116) q[0];
x q[1];
rz(1.3387483) q[2];
sx q[2];
rz(-1.9221483) q[2];
sx q[2];
rz(0.50594375) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.52555841) q[1];
sx q[1];
rz(-1.8174606) q[1];
sx q[1];
rz(0.14825578) q[1];
rz(-pi) q[2];
rz(0.90274685) q[3];
sx q[3];
rz(-0.76551907) q[3];
sx q[3];
rz(-2.3428041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.89821833) q[2];
sx q[2];
rz(-2.2369907) q[2];
sx q[2];
rz(-0.96837366) q[2];
rz(0.51356703) q[3];
sx q[3];
rz(-1.7889203) q[3];
sx q[3];
rz(0.33094049) q[3];
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
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53720713) q[0];
sx q[0];
rz(-2.4810915) q[0];
sx q[0];
rz(0.78395098) q[0];
rz(-1.2046658) q[1];
sx q[1];
rz(-0.37310633) q[1];
sx q[1];
rz(1.3752259) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.169302) q[0];
sx q[0];
rz(-1.307104) q[0];
sx q[0];
rz(1.792981) q[0];
rz(-pi) q[1];
rz(1.4156467) q[2];
sx q[2];
rz(-1.5976693) q[2];
sx q[2];
rz(0.43153119) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.8177796) q[1];
sx q[1];
rz(-2.166417) q[1];
sx q[1];
rz(0.88902529) q[1];
rz(-pi) q[2];
rz(1.8834388) q[3];
sx q[3];
rz(-1.4015159) q[3];
sx q[3];
rz(1.9692242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.84460008) q[2];
sx q[2];
rz(-0.30868369) q[2];
sx q[2];
rz(-2.6958579) q[2];
rz(2.5388057) q[3];
sx q[3];
rz(-1.2647537) q[3];
sx q[3];
rz(0.84356892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8730901) q[0];
sx q[0];
rz(-0.52274811) q[0];
sx q[0];
rz(-2.6173746) q[0];
rz(1.9408608) q[1];
sx q[1];
rz(-1.2501161) q[1];
sx q[1];
rz(-1.9573617) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94468695) q[0];
sx q[0];
rz(-1.2254929) q[0];
sx q[0];
rz(1.3698335) q[0];
rz(-pi) q[1];
rz(-3.0734407) q[2];
sx q[2];
rz(-1.2556723) q[2];
sx q[2];
rz(1.5308876) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0457234) q[1];
sx q[1];
rz(-1.5585526) q[1];
sx q[1];
rz(0.13282383) q[1];
rz(-pi) q[2];
rz(-3.0940787) q[3];
sx q[3];
rz(-2.5121452) q[3];
sx q[3];
rz(-2.269182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9253917) q[2];
sx q[2];
rz(-0.4500469) q[2];
sx q[2];
rz(1.0518543) q[2];
rz(-0.22107407) q[3];
sx q[3];
rz(-1.3007921) q[3];
sx q[3];
rz(2.2033447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5476407) q[0];
sx q[0];
rz(-1.5008391) q[0];
sx q[0];
rz(0.048901625) q[0];
rz(2.7652057) q[1];
sx q[1];
rz(-2.2793437) q[1];
sx q[1];
rz(1.9752165) q[1];
rz(-0.57009956) q[2];
sx q[2];
rz(-1.5975633) q[2];
sx q[2];
rz(-0.057889197) q[2];
rz(2.0509649) q[3];
sx q[3];
rz(-2.2259897) q[3];
sx q[3];
rz(2.5799354) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
