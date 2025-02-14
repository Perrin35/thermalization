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
rz(0.25207818) q[0];
sx q[0];
rz(3.7171465) q[0];
sx q[0];
rz(9.3024749) q[0];
rz(-2.0343434) q[1];
sx q[1];
rz(-2.1802433) q[1];
sx q[1];
rz(-3.1032739) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3275546) q[0];
sx q[0];
rz(-0.71056847) q[0];
sx q[0];
rz(-2.1585805) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0228268) q[2];
sx q[2];
rz(-1.6890172) q[2];
sx q[2];
rz(1.880065) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.9895175) q[1];
sx q[1];
rz(-0.65836009) q[1];
sx q[1];
rz(2.1170298) q[1];
rz(0.24738048) q[3];
sx q[3];
rz(-1.5508964) q[3];
sx q[3];
rz(0.54075356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.4965839) q[2];
sx q[2];
rz(-1.4270447) q[2];
sx q[2];
rz(1.6928147) q[2];
rz(-2.951176) q[3];
sx q[3];
rz(-2.0864291) q[3];
sx q[3];
rz(-0.14934389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.2038302) q[0];
sx q[0];
rz(-1.8730524) q[0];
sx q[0];
rz(-2.8835836) q[0];
rz(1.5915271) q[1];
sx q[1];
rz(-1.9015046) q[1];
sx q[1];
rz(-0.16664997) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44010559) q[0];
sx q[0];
rz(-0.99033725) q[0];
sx q[0];
rz(-1.6850182) q[0];
x q[1];
rz(2.7562253) q[2];
sx q[2];
rz(-1.8358942) q[2];
sx q[2];
rz(0.4601269) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.8434324) q[1];
sx q[1];
rz(-0.86151228) q[1];
sx q[1];
rz(-1.405836) q[1];
rz(-pi) q[2];
x q[2];
rz(0.54661481) q[3];
sx q[3];
rz(-2.0677462) q[3];
sx q[3];
rz(0.95595804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0717281) q[2];
sx q[2];
rz(-2.0191777) q[2];
sx q[2];
rz(1.9094763) q[2];
rz(-0.92464906) q[3];
sx q[3];
rz(-0.93622127) q[3];
sx q[3];
rz(-1.1054976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.6919747) q[0];
sx q[0];
rz(-1.469935) q[0];
sx q[0];
rz(-0.0019419226) q[0];
rz(3.107403) q[1];
sx q[1];
rz(-1.9296153) q[1];
sx q[1];
rz(1.5431822) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6337834) q[0];
sx q[0];
rz(-0.52716053) q[0];
sx q[0];
rz(-2.6109004) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4056881) q[2];
sx q[2];
rz(-1.2703478) q[2];
sx q[2];
rz(1.4020593) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1605121) q[1];
sx q[1];
rz(-0.35582572) q[1];
sx q[1];
rz(-0.11248223) q[1];
rz(2.8974124) q[3];
sx q[3];
rz(-0.93726087) q[3];
sx q[3];
rz(1.4118495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.2446642) q[2];
sx q[2];
rz(-0.43411532) q[2];
sx q[2];
rz(1.0373235) q[2];
rz(-1.1050998) q[3];
sx q[3];
rz(-1.5367855) q[3];
sx q[3];
rz(-2.0028116) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8741375) q[0];
sx q[0];
rz(-0.48478165) q[0];
sx q[0];
rz(-1.4111891) q[0];
rz(-0.63938582) q[1];
sx q[1];
rz(-1.6849898) q[1];
sx q[1];
rz(3.0944518) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61091238) q[0];
sx q[0];
rz(-0.99448863) q[0];
sx q[0];
rz(1.9082101) q[0];
x q[1];
rz(0.36836715) q[2];
sx q[2];
rz(-0.50058156) q[2];
sx q[2];
rz(0.49886242) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2915708) q[1];
sx q[1];
rz(-2.154989) q[1];
sx q[1];
rz(2.7136346) q[1];
rz(-1.239226) q[3];
sx q[3];
rz(-2.3904388) q[3];
sx q[3];
rz(-2.0560665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.822927) q[2];
sx q[2];
rz(-1.1456127) q[2];
sx q[2];
rz(1.9226496) q[2];
rz(-2.8259891) q[3];
sx q[3];
rz(-3.0867519) q[3];
sx q[3];
rz(-2.992673) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8197935) q[0];
sx q[0];
rz(-2.5892374) q[0];
sx q[0];
rz(0.34859443) q[0];
rz(-1.2527342) q[1];
sx q[1];
rz(-1.1628954) q[1];
sx q[1];
rz(-1.3016275) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7093703) q[0];
sx q[0];
rz(-2.3314868) q[0];
sx q[0];
rz(3.0322984) q[0];
rz(-pi) q[1];
rz(-1.9867861) q[2];
sx q[2];
rz(-2.092053) q[2];
sx q[2];
rz(-0.74300569) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.0685454) q[1];
sx q[1];
rz(-1.7918192) q[1];
sx q[1];
rz(-0.18494341) q[1];
x q[2];
rz(1.0318448) q[3];
sx q[3];
rz(-2.5205118) q[3];
sx q[3];
rz(-0.10824848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18315135) q[2];
sx q[2];
rz(-2.8330467) q[2];
sx q[2];
rz(-1.4754254) q[2];
rz(-2.4557377) q[3];
sx q[3];
rz(-0.97969046) q[3];
sx q[3];
rz(0.53387749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0297861) q[0];
sx q[0];
rz(-0.15203467) q[0];
sx q[0];
rz(1.4923805) q[0];
rz(0.97914186) q[1];
sx q[1];
rz(-1.5359595) q[1];
sx q[1];
rz(-1.917256) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5021073) q[0];
sx q[0];
rz(-1.3728791) q[0];
sx q[0];
rz(1.60198) q[0];
x q[1];
rz(-2.5923205) q[2];
sx q[2];
rz(-1.3830997) q[2];
sx q[2];
rz(-2.930738) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.2883207) q[1];
sx q[1];
rz(-1.5832807) q[1];
sx q[1];
rz(2.1846143) q[1];
rz(0.34747261) q[3];
sx q[3];
rz(-2.1712448) q[3];
sx q[3];
rz(-0.97755177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.41902038) q[2];
sx q[2];
rz(-1.6442862) q[2];
sx q[2];
rz(-0.9160308) q[2];
rz(1.3845059) q[3];
sx q[3];
rz(-1.8839096) q[3];
sx q[3];
rz(-2.4345583) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1322121) q[0];
sx q[0];
rz(-3.0896602) q[0];
sx q[0];
rz(-0.95426553) q[0];
rz(0.43680278) q[1];
sx q[1];
rz(-1.5635468) q[1];
sx q[1];
rz(3.0314235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5075275) q[0];
sx q[0];
rz(-1.332009) q[0];
sx q[0];
rz(-2.8077543) q[0];
rz(-pi) q[1];
rz(-1.0344347) q[2];
sx q[2];
rz(-2.1948994) q[2];
sx q[2];
rz(1.0076866) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.71437826) q[1];
sx q[1];
rz(-0.34571474) q[1];
sx q[1];
rz(1.5835254) q[1];
x q[2];
rz(-2.8875868) q[3];
sx q[3];
rz(-2.0585103) q[3];
sx q[3];
rz(-0.28764492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7097912) q[2];
sx q[2];
rz(-0.0063889901) q[2];
sx q[2];
rz(-0.94376454) q[2];
rz(2.5778095) q[3];
sx q[3];
rz(-1.0958593) q[3];
sx q[3];
rz(-1.110466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3540045) q[0];
sx q[0];
rz(-2.8546951) q[0];
sx q[0];
rz(-0.34550825) q[0];
rz(1.0954789) q[1];
sx q[1];
rz(-1.4778719) q[1];
sx q[1];
rz(1.5617721) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45392198) q[0];
sx q[0];
rz(-0.78886388) q[0];
sx q[0];
rz(0.044290941) q[0];
rz(0.91362915) q[2];
sx q[2];
rz(-1.5217178) q[2];
sx q[2];
rz(0.14684453) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9672444) q[1];
sx q[1];
rz(-2.1617608) q[1];
sx q[1];
rz(2.6978542) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6437279) q[3];
sx q[3];
rz(-1.5988873) q[3];
sx q[3];
rz(1.3974578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.3788508) q[2];
sx q[2];
rz(-0.26143917) q[2];
sx q[2];
rz(1.7108819) q[2];
rz(-0.53449574) q[3];
sx q[3];
rz(-1.675324) q[3];
sx q[3];
rz(-0.9526332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5296103) q[0];
sx q[0];
rz(-3.072325) q[0];
sx q[0];
rz(-1.3105422) q[0];
rz(-2.2546841) q[1];
sx q[1];
rz(-2.3014018) q[1];
sx q[1];
rz(1.2695674) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3536772) q[0];
sx q[0];
rz(-1.6991109) q[0];
sx q[0];
rz(-1.9246035) q[0];
rz(-pi) q[1];
rz(0.15523703) q[2];
sx q[2];
rz(-1.198215) q[2];
sx q[2];
rz(-0.22107228) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3514401) q[1];
sx q[1];
rz(-1.8428486) q[1];
sx q[1];
rz(2.9012009) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.19143243) q[3];
sx q[3];
rz(-2.9342954) q[3];
sx q[3];
rz(1.3216922) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.612192) q[2];
sx q[2];
rz(-1.8058913) q[2];
sx q[2];
rz(0.23294918) q[2];
rz(1.285078) q[3];
sx q[3];
rz(-2.464747) q[3];
sx q[3];
rz(2.7555833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1583629) q[0];
sx q[0];
rz(-2.5322999) q[0];
sx q[0];
rz(-3.0443211) q[0];
rz(-2.3376047) q[1];
sx q[1];
rz(-2.4318047) q[1];
sx q[1];
rz(-0.21673094) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0618503) q[0];
sx q[0];
rz(-1.5155795) q[0];
sx q[0];
rz(1.400007) q[0];
rz(-pi) q[1];
rz(-2.8267536) q[2];
sx q[2];
rz(-2.9926569) q[2];
sx q[2];
rz(-0.51477942) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8179124) q[1];
sx q[1];
rz(-0.92880946) q[1];
sx q[1];
rz(0.85425185) q[1];
rz(-2.4794934) q[3];
sx q[3];
rz(-1.2939014) q[3];
sx q[3];
rz(-2.3014042) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.29297605) q[2];
sx q[2];
rz(-2.8605707) q[2];
sx q[2];
rz(2.6463553) q[2];
rz(-1.5826591) q[3];
sx q[3];
rz(-1.0844743) q[3];
sx q[3];
rz(-2.9787279) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0804629) q[0];
sx q[0];
rz(-1.5123788) q[0];
sx q[0];
rz(2.6176591) q[0];
rz(2.0551266) q[1];
sx q[1];
rz(-2.6230984) q[1];
sx q[1];
rz(-2.7883504) q[1];
rz(-1.3475781) q[2];
sx q[2];
rz(-0.032608727) q[2];
sx q[2];
rz(-0.41667117) q[2];
rz(-0.090251017) q[3];
sx q[3];
rz(-0.8606438) q[3];
sx q[3];
rz(0.46910486) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
