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
rz(1.264313) q[0];
sx q[0];
rz(-2.8180583) q[0];
sx q[0];
rz(-2.6927595) q[0];
rz(-1.1558865) q[1];
sx q[1];
rz(1.356025) q[1];
sx q[1];
rz(7.4577509) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6960775) q[0];
sx q[0];
rz(-2.4598274) q[0];
sx q[0];
rz(2.5755139) q[0];
rz(-2.817906) q[2];
sx q[2];
rz(-1.5059221) q[2];
sx q[2];
rz(-1.6271568) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5331363) q[1];
sx q[1];
rz(-1.605833) q[1];
sx q[1];
rz(-0.94989454) q[1];
rz(-pi) q[2];
rz(-1.6714736) q[3];
sx q[3];
rz(-2.0824471) q[3];
sx q[3];
rz(-0.6938405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.590098) q[2];
sx q[2];
rz(-1.3000725) q[2];
sx q[2];
rz(-2.2287915) q[2];
rz(1.270795) q[3];
sx q[3];
rz(-1.0400892) q[3];
sx q[3];
rz(-2.2916268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.0047282334) q[0];
sx q[0];
rz(-2.5623463) q[0];
sx q[0];
rz(2.7194523) q[0];
rz(-0.36960754) q[1];
sx q[1];
rz(-1.0969176) q[1];
sx q[1];
rz(0.94013989) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5135315) q[0];
sx q[0];
rz(-1.2413176) q[0];
sx q[0];
rz(0.0026436289) q[0];
x q[1];
rz(1.0375848) q[2];
sx q[2];
rz(-1.4681446) q[2];
sx q[2];
rz(-0.29666049) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.861189) q[1];
sx q[1];
rz(-0.94377667) q[1];
sx q[1];
rz(-2.1908733) q[1];
x q[2];
rz(2.8690954) q[3];
sx q[3];
rz(-0.91236189) q[3];
sx q[3];
rz(-2.8373425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.9422354) q[2];
sx q[2];
rz(-2.4713559) q[2];
sx q[2];
rz(0.90636903) q[2];
rz(-0.97936112) q[3];
sx q[3];
rz(-2.7370079) q[3];
sx q[3];
rz(1.2649068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1449428) q[0];
sx q[0];
rz(-3.0746089) q[0];
sx q[0];
rz(1.8970733) q[0];
rz(2.7527346) q[1];
sx q[1];
rz(-1.8617947) q[1];
sx q[1];
rz(2.8188474) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30814442) q[0];
sx q[0];
rz(-0.9128154) q[0];
sx q[0];
rz(1.4288155) q[0];
x q[1];
rz(1.1006196) q[2];
sx q[2];
rz(-0.61264804) q[2];
sx q[2];
rz(2.8620811) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1908993) q[1];
sx q[1];
rz(-2.3436693) q[1];
sx q[1];
rz(1.8969593) q[1];
x q[2];
rz(-0.1348008) q[3];
sx q[3];
rz(-0.35405891) q[3];
sx q[3];
rz(-1.9446789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9129703) q[2];
sx q[2];
rz(-1.3801489) q[2];
sx q[2];
rz(1.3529533) q[2];
rz(2.8690423) q[3];
sx q[3];
rz(-1.291178) q[3];
sx q[3];
rz(-1.6779617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54101855) q[0];
sx q[0];
rz(-2.8193642) q[0];
sx q[0];
rz(-1.898265) q[0];
rz(0.86878949) q[1];
sx q[1];
rz(-2.181535) q[1];
sx q[1];
rz(2.0713846) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7420827) q[0];
sx q[0];
rz(-2.1993981) q[0];
sx q[0];
rz(0.6357155) q[0];
x q[1];
rz(0.79480537) q[2];
sx q[2];
rz(-1.897395) q[2];
sx q[2];
rz(0.87233938) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.29954545) q[1];
sx q[1];
rz(-1.4562325) q[1];
sx q[1];
rz(-1.8946527) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7620874) q[3];
sx q[3];
rz(-1.9076348) q[3];
sx q[3];
rz(-1.9302521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.59226817) q[2];
sx q[2];
rz(-0.78846875) q[2];
sx q[2];
rz(-2.2198524) q[2];
rz(0.24374572) q[3];
sx q[3];
rz(-1.7071525) q[3];
sx q[3];
rz(-2.8488979) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9859966) q[0];
sx q[0];
rz(-1.8159741) q[0];
sx q[0];
rz(-0.50088125) q[0];
rz(-1.1174551) q[1];
sx q[1];
rz(-1.3331579) q[1];
sx q[1];
rz(-2.8270328) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68334333) q[0];
sx q[0];
rz(-1.2309845) q[0];
sx q[0];
rz(2.763624) q[0];
x q[1];
rz(-1.5273845) q[2];
sx q[2];
rz(-1.5752666) q[2];
sx q[2];
rz(-2.4624766) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52175888) q[1];
sx q[1];
rz(-1.9653826) q[1];
sx q[1];
rz(-2.16741) q[1];
rz(1.3077478) q[3];
sx q[3];
rz(-2.7125053) q[3];
sx q[3];
rz(1.7744296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.306281) q[2];
sx q[2];
rz(-1.1113144) q[2];
sx q[2];
rz(1.8446946) q[2];
rz(-0.48526192) q[3];
sx q[3];
rz(-1.5819712) q[3];
sx q[3];
rz(-1.1135134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6473963) q[0];
sx q[0];
rz(-1.2911456) q[0];
sx q[0];
rz(0.099763481) q[0];
rz(-1.8708723) q[1];
sx q[1];
rz(-2.2427509) q[1];
sx q[1];
rz(1.3894003) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64369624) q[0];
sx q[0];
rz(-0.31539105) q[0];
sx q[0];
rz(-2.5853378) q[0];
rz(1.6758789) q[2];
sx q[2];
rz(-1.7438403) q[2];
sx q[2];
rz(-1.6478754) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9854765) q[1];
sx q[1];
rz(-0.96806541) q[1];
sx q[1];
rz(0.74405693) q[1];
rz(-0.12677315) q[3];
sx q[3];
rz(-0.78713464) q[3];
sx q[3];
rz(-1.2490369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.49326593) q[2];
sx q[2];
rz(-0.47097012) q[2];
sx q[2];
rz(-0.89034447) q[2];
rz(1.1912311) q[3];
sx q[3];
rz(-1.3072562) q[3];
sx q[3];
rz(0.90203917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4297727) q[0];
sx q[0];
rz(-0.044782488) q[0];
sx q[0];
rz(-0.032715948) q[0];
rz(0.50321594) q[1];
sx q[1];
rz(-1.0992173) q[1];
sx q[1];
rz(0.12282664) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5635934) q[0];
sx q[0];
rz(-1.2482572) q[0];
sx q[0];
rz(0.75046993) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5678504) q[2];
sx q[2];
rz(-2.5649568) q[2];
sx q[2];
rz(1.6078311) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.87007755) q[1];
sx q[1];
rz(-2.6598499) q[1];
sx q[1];
rz(-1.355257) q[1];
rz(2.8005573) q[3];
sx q[3];
rz(-0.37217316) q[3];
sx q[3];
rz(-2.2256874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7922625) q[2];
sx q[2];
rz(-0.96668875) q[2];
sx q[2];
rz(1.9897423) q[2];
rz(2.2564015) q[3];
sx q[3];
rz(-1.3693634) q[3];
sx q[3];
rz(-1.7159897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.952878) q[0];
sx q[0];
rz(-2.2080244) q[0];
sx q[0];
rz(2.2123912) q[0];
rz(-0.62905351) q[1];
sx q[1];
rz(-1.355143) q[1];
sx q[1];
rz(0.74310511) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5261054) q[0];
sx q[0];
rz(-2.0172523) q[0];
sx q[0];
rz(-0.087421405) q[0];
rz(1.6665206) q[2];
sx q[2];
rz(-2.0572955) q[2];
sx q[2];
rz(-0.41991389) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.4854606) q[1];
sx q[1];
rz(-0.4194057) q[1];
sx q[1];
rz(2.9784883) q[1];
rz(0.93502126) q[3];
sx q[3];
rz(-1.1639769) q[3];
sx q[3];
rz(2.0182899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.5653845) q[2];
sx q[2];
rz(-2.3614063) q[2];
sx q[2];
rz(-0.72594491) q[2];
rz(-0.37096008) q[3];
sx q[3];
rz(-1.5434664) q[3];
sx q[3];
rz(-2.7788739) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.87042701) q[0];
sx q[0];
rz(-2.3112264) q[0];
sx q[0];
rz(-0.56394947) q[0];
rz(-1.14934) q[1];
sx q[1];
rz(-0.96200689) q[1];
sx q[1];
rz(-0.74251485) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71096984) q[0];
sx q[0];
rz(-1.1046243) q[0];
sx q[0];
rz(-0.52180334) q[0];
x q[1];
rz(2.7269823) q[2];
sx q[2];
rz(-1.8413723) q[2];
sx q[2];
rz(-2.0913578) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5919784) q[1];
sx q[1];
rz(-1.5266298) q[1];
sx q[1];
rz(-2.7034034) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5768859) q[3];
sx q[3];
rz(-2.2032524) q[3];
sx q[3];
rz(-2.6868827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.62770647) q[2];
sx q[2];
rz(-1.0739001) q[2];
sx q[2];
rz(-2.7970496) q[2];
rz(-0.44376093) q[3];
sx q[3];
rz(-2.0659476) q[3];
sx q[3];
rz(1.3226604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(-0.76170707) q[0];
sx q[0];
rz(-2.8590617) q[0];
sx q[0];
rz(-0.66201061) q[0];
rz(-0.33540353) q[1];
sx q[1];
rz(-1.6796203) q[1];
sx q[1];
rz(-0.9368771) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7279908) q[0];
sx q[0];
rz(-1.9303891) q[0];
sx q[0];
rz(-0.92443621) q[0];
rz(1.30055) q[2];
sx q[2];
rz(-1.593538) q[2];
sx q[2];
rz(1.8128527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.7171927) q[1];
sx q[1];
rz(-0.48312995) q[1];
sx q[1];
rz(-2.4880954) q[1];
rz(-pi) q[2];
x q[2];
rz(0.85031894) q[3];
sx q[3];
rz(-2.8788044) q[3];
sx q[3];
rz(-2.9966054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.4235437) q[2];
sx q[2];
rz(-0.23423883) q[2];
sx q[2];
rz(2.6424109) q[2];
rz(1.6792123) q[3];
sx q[3];
rz(-1.1464109) q[3];
sx q[3];
rz(-1.4596938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725631) q[0];
sx q[0];
rz(-1.220663) q[0];
sx q[0];
rz(2.2882373) q[0];
rz(-0.60494963) q[1];
sx q[1];
rz(-1.1440944) q[1];
sx q[1];
rz(0.49857421) q[1];
rz(2.5685979) q[2];
sx q[2];
rz(-1.6934494) q[2];
sx q[2];
rz(2.1394503) q[2];
rz(-0.50157401) q[3];
sx q[3];
rz(-1.4091103) q[3];
sx q[3];
rz(1.4179358) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
