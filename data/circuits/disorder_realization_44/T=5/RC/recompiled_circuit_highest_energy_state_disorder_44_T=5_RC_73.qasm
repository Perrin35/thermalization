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
rz(-0.15260881) q[0];
sx q[0];
rz(-1.1996562) q[0];
sx q[0];
rz(-0.3056404) q[0];
rz(0.29023728) q[1];
sx q[1];
rz(-1.6936392) q[1];
sx q[1];
rz(1.0040959) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476513) q[0];
sx q[0];
rz(-3.0783639) q[0];
sx q[0];
rz(2.7020523) q[0];
x q[1];
rz(-1.8199241) q[2];
sx q[2];
rz(-1.7580571) q[2];
sx q[2];
rz(1.934777) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6184147) q[1];
sx q[1];
rz(-2.6175099) q[1];
sx q[1];
rz(-2.2487703) q[1];
x q[2];
rz(-3.0707487) q[3];
sx q[3];
rz(-0.53991441) q[3];
sx q[3];
rz(1.2666463) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4940138) q[2];
sx q[2];
rz(-0.61654377) q[2];
sx q[2];
rz(2.5971863) q[2];
rz(0.38862774) q[3];
sx q[3];
rz(-1.2231239) q[3];
sx q[3];
rz(2.2012034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87043864) q[0];
sx q[0];
rz(-2.1156023) q[0];
sx q[0];
rz(2.0197268) q[0];
rz(2.0815966) q[1];
sx q[1];
rz(-0.50299877) q[1];
sx q[1];
rz(-2.2559135) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1141027) q[0];
sx q[0];
rz(-1.8654279) q[0];
sx q[0];
rz(-2.824371) q[0];
rz(-pi) q[1];
rz(-0.76876872) q[2];
sx q[2];
rz(-2.2115006) q[2];
sx q[2];
rz(-2.6219683) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.52351511) q[1];
sx q[1];
rz(-1.4898058) q[1];
sx q[1];
rz(-0.063132719) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4526414) q[3];
sx q[3];
rz(-1.6756578) q[3];
sx q[3];
rz(-2.855043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1765959) q[2];
sx q[2];
rz(-1.6792038) q[2];
sx q[2];
rz(-3.1122567) q[2];
rz(-0.36597478) q[3];
sx q[3];
rz(-2.6656606) q[3];
sx q[3];
rz(-0.16987814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[3];
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
rz(-2.733424) q[0];
sx q[0];
rz(-1.3329196) q[0];
sx q[0];
rz(-2.9881706) q[0];
rz(0.19639213) q[1];
sx q[1];
rz(-1.6729665) q[1];
sx q[1];
rz(-2.3207655) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9388537) q[0];
sx q[0];
rz(-1.0717055) q[0];
sx q[0];
rz(-0.75627723) q[0];
x q[1];
rz(2.1961741) q[2];
sx q[2];
rz(-2.4074005) q[2];
sx q[2];
rz(-1.1933094) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.094999639) q[1];
sx q[1];
rz(-0.41727704) q[1];
sx q[1];
rz(-0.098578171) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5215932) q[3];
sx q[3];
rz(-1.9370034) q[3];
sx q[3];
rz(-2.5871426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.11744943) q[2];
sx q[2];
rz(-1.3261745) q[2];
sx q[2];
rz(-2.2085025) q[2];
rz(2.7945331) q[3];
sx q[3];
rz(-2.0324028) q[3];
sx q[3];
rz(-0.6412653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1953122) q[0];
sx q[0];
rz(-1.1271789) q[0];
sx q[0];
rz(2.0830925) q[0];
rz(3.0463386) q[1];
sx q[1];
rz(-1.1493827) q[1];
sx q[1];
rz(0.20406318) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0515186) q[0];
sx q[0];
rz(-1.7885889) q[0];
sx q[0];
rz(-0.15845944) q[0];
rz(-pi) q[1];
rz(-0.45943164) q[2];
sx q[2];
rz(-2.0028015) q[2];
sx q[2];
rz(2.9615796) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.82641027) q[1];
sx q[1];
rz(-0.99821222) q[1];
sx q[1];
rz(-1.8099643) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7768562) q[3];
sx q[3];
rz(-2.2584174) q[3];
sx q[3];
rz(-1.7182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.05693398) q[2];
sx q[2];
rz(-1.441842) q[2];
sx q[2];
rz(1.8208246) q[2];
rz(-1.002958) q[3];
sx q[3];
rz(-1.2489677) q[3];
sx q[3];
rz(-2.5761719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4249307) q[0];
sx q[0];
rz(-1.2053763) q[0];
sx q[0];
rz(0.64111125) q[0];
rz(0.59116108) q[1];
sx q[1];
rz(-1.4407651) q[1];
sx q[1];
rz(1.444918) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4137694) q[0];
sx q[0];
rz(-1.1846847) q[0];
sx q[0];
rz(-1.5408433) q[0];
x q[1];
rz(-2.9282753) q[2];
sx q[2];
rz(-0.89175311) q[2];
sx q[2];
rz(0.17903331) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.843335) q[1];
sx q[1];
rz(-1.3239256) q[1];
sx q[1];
rz(2.6554537) q[1];
rz(1.0459445) q[3];
sx q[3];
rz(-1.305745) q[3];
sx q[3];
rz(1.8432528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.558202) q[2];
sx q[2];
rz(-2.148197) q[2];
sx q[2];
rz(2.2885585) q[2];
rz(-2.9663626) q[3];
sx q[3];
rz(-1.3220738) q[3];
sx q[3];
rz(2.5186553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.919642) q[0];
sx q[0];
rz(-0.80606824) q[0];
sx q[0];
rz(0.44575086) q[0];
rz(0.99459612) q[1];
sx q[1];
rz(-1.0044731) q[1];
sx q[1];
rz(-1.6533143) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22791187) q[0];
sx q[0];
rz(-0.74800038) q[0];
sx q[0];
rz(2.8196003) q[0];
rz(-pi) q[1];
rz(-1.4966665) q[2];
sx q[2];
rz(-1.526381) q[2];
sx q[2];
rz(1.6629459) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.52998972) q[1];
sx q[1];
rz(-1.7771136) q[1];
sx q[1];
rz(1.2633268) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7986695) q[3];
sx q[3];
rz(-2.8066435) q[3];
sx q[3];
rz(0.11996525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.17066869) q[2];
sx q[2];
rz(-1.2430151) q[2];
sx q[2];
rz(-1.0260014) q[2];
rz(0.79801997) q[3];
sx q[3];
rz(-2.3012216) q[3];
sx q[3];
rz(0.26871267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84665027) q[0];
sx q[0];
rz(-1.7468528) q[0];
sx q[0];
rz(-0.41644874) q[0];
rz(1.7448447) q[1];
sx q[1];
rz(-1.6114707) q[1];
sx q[1];
rz(-1.4769311) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0357405) q[0];
sx q[0];
rz(-2.1243411) q[0];
sx q[0];
rz(-0.043364924) q[0];
x q[1];
rz(-0.91628051) q[2];
sx q[2];
rz(-2.6323292) q[2];
sx q[2];
rz(-1.3167574) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1394964) q[1];
sx q[1];
rz(-2.0130035) q[1];
sx q[1];
rz(-0.79331974) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6536413) q[3];
sx q[3];
rz(-1.9085064) q[3];
sx q[3];
rz(-3.0520671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.2485409) q[2];
sx q[2];
rz(-0.38022843) q[2];
sx q[2];
rz(-2.5377972) q[2];
rz(0.90886146) q[3];
sx q[3];
rz(-1.4330319) q[3];
sx q[3];
rz(-2.042167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73868442) q[0];
sx q[0];
rz(-1.7576907) q[0];
sx q[0];
rz(0.72147328) q[0];
rz(-1.8062704) q[1];
sx q[1];
rz(-1.1386917) q[1];
sx q[1];
rz(-1.8849751) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4271256) q[0];
sx q[0];
rz(-0.79475105) q[0];
sx q[0];
rz(2.7156272) q[0];
rz(0.5783028) q[2];
sx q[2];
rz(-1.1546557) q[2];
sx q[2];
rz(1.5612891) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.2366888) q[1];
sx q[1];
rz(-1.8042548) q[1];
sx q[1];
rz(1.9953362) q[1];
rz(-1.6001892) q[3];
sx q[3];
rz(-1.5252068) q[3];
sx q[3];
rz(-0.33908333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.98715544) q[2];
sx q[2];
rz(-2.4989231) q[2];
sx q[2];
rz(-1.2745693) q[2];
rz(-0.68862033) q[3];
sx q[3];
rz(-1.4911115) q[3];
sx q[3];
rz(3.137818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(-2.1275198) q[0];
sx q[0];
rz(-1.8919683) q[0];
sx q[0];
rz(0.67181146) q[0];
rz(-1.6067243) q[1];
sx q[1];
rz(-2.913919) q[1];
sx q[1];
rz(-2.6188376) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3776457) q[0];
sx q[0];
rz(-1.3113759) q[0];
sx q[0];
rz(2.460145) q[0];
rz(0.050809697) q[2];
sx q[2];
rz(-1.9138971) q[2];
sx q[2];
rz(-1.1468514) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6784489) q[1];
sx q[1];
rz(-1.1511973) q[1];
sx q[1];
rz(1.9897919) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32355752) q[3];
sx q[3];
rz(-1.3731628) q[3];
sx q[3];
rz(1.4909286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25040024) q[2];
sx q[2];
rz(-0.6468536) q[2];
sx q[2];
rz(2.7817173) q[2];
rz(1.8187836) q[3];
sx q[3];
rz(-1.6836932) q[3];
sx q[3];
rz(0.047164269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5322402) q[0];
sx q[0];
rz(-2.1873964) q[0];
sx q[0];
rz(-0.94616079) q[0];
rz(0.040945176) q[1];
sx q[1];
rz(-1.8649201) q[1];
sx q[1];
rz(-1.8765705) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0909983) q[0];
sx q[0];
rz(-0.48830214) q[0];
sx q[0];
rz(0.22240345) q[0];
rz(-pi) q[1];
rz(0.47750116) q[2];
sx q[2];
rz(-2.6386119) q[2];
sx q[2];
rz(1.8790085) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.6835306) q[1];
sx q[1];
rz(-0.89408656) q[1];
sx q[1];
rz(-2.2735734) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6947097) q[3];
sx q[3];
rz(-0.37624761) q[3];
sx q[3];
rz(-0.28722426) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.9762743) q[2];
sx q[2];
rz(-1.9860622) q[2];
sx q[2];
rz(-2.9602642) q[2];
rz(2.2881962) q[3];
sx q[3];
rz(-2.3386164) q[3];
sx q[3];
rz(1.3388504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4637909) q[0];
sx q[0];
rz(-1.811857) q[0];
sx q[0];
rz(-1.3743251) q[0];
rz(1.6865267) q[1];
sx q[1];
rz(-1.68119) q[1];
sx q[1];
rz(2.841058) q[1];
rz(0.64114943) q[2];
sx q[2];
rz(-1.1417626) q[2];
sx q[2];
rz(-0.7456197) q[2];
rz(2.8541722) q[3];
sx q[3];
rz(-1.8183553) q[3];
sx q[3];
rz(2.9498065) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
