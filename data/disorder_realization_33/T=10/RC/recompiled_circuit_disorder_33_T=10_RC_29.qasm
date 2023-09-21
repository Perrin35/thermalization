OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(-1.3077868) q[0];
rz(-2.0055327) q[1];
sx q[1];
rz(-2.2059031) q[1];
sx q[1];
rz(1.5712665) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64496541) q[0];
sx q[0];
rz(-0.49155203) q[0];
sx q[0];
rz(0.77792032) q[0];
rz(0.76531305) q[2];
sx q[2];
rz(-2.1047449) q[2];
sx q[2];
rz(-3.090976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1882602) q[1];
sx q[1];
rz(-0.42020513) q[1];
sx q[1];
rz(2.5145867) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9340431) q[3];
sx q[3];
rz(-1.3875291) q[3];
sx q[3];
rz(2.7024384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.87542614) q[2];
sx q[2];
rz(-2.8484919) q[2];
sx q[2];
rz(-1.1323294) q[2];
rz(1.6752361) q[3];
sx q[3];
rz(-1.3365859) q[3];
sx q[3];
rz(-2.1291389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(2.9557513) q[0];
rz(2.5813685) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(0.21683189) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8502055) q[0];
sx q[0];
rz(-0.7190401) q[0];
sx q[0];
rz(-2.015381) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47302834) q[2];
sx q[2];
rz(-2.0748667) q[2];
sx q[2];
rz(-2.8216528) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.195897) q[1];
sx q[1];
rz(-0.90598124) q[1];
sx q[1];
rz(-0.2701387) q[1];
rz(0.84083765) q[3];
sx q[3];
rz(-0.84078046) q[3];
sx q[3];
rz(-2.9527612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.8314787) q[2];
sx q[2];
rz(-0.82565132) q[2];
sx q[2];
rz(1.8537834) q[2];
rz(-2.3790322) q[3];
sx q[3];
rz(-1.972714) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4644311) q[0];
sx q[0];
rz(-0.34496775) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(-1.3263946) q[1];
sx q[1];
rz(-1.7809968) q[1];
sx q[1];
rz(-2.2089829) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1754477) q[0];
sx q[0];
rz(-1.6245337) q[0];
sx q[0];
rz(1.8885814) q[0];
rz(-0.89838018) q[2];
sx q[2];
rz(-2.8542238) q[2];
sx q[2];
rz(0.76295602) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0886791) q[1];
sx q[1];
rz(-1.0148078) q[1];
sx q[1];
rz(-1.710379) q[1];
rz(-pi) q[2];
rz(-2.2266085) q[3];
sx q[3];
rz(-2.4768156) q[3];
sx q[3];
rz(-1.2802326) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.9937667) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-2.0489342) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-2.0565624) q[3];
sx q[3];
rz(0.96737635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7820691) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(2.6413667) q[0];
rz(2.3362828) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(1.4979699) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1456137) q[0];
sx q[0];
rz(-1.0756452) q[0];
sx q[0];
rz(-2.8061295) q[0];
rz(-pi) q[1];
rz(3.0849886) q[2];
sx q[2];
rz(-1.3805693) q[2];
sx q[2];
rz(-2.9375926) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8151617) q[1];
sx q[1];
rz(-1.8306499) q[1];
sx q[1];
rz(-1.3304779) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3718932) q[3];
sx q[3];
rz(-0.61004988) q[3];
sx q[3];
rz(-1.9609914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.74636373) q[2];
sx q[2];
rz(-2.5791898) q[2];
sx q[2];
rz(-0.70181075) q[2];
rz(2.3102405) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(-0.62197661) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2410626) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(-0.81533122) q[0];
rz(1.5218081) q[1];
sx q[1];
rz(-0.83414572) q[1];
sx q[1];
rz(-2.0933847) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5883023) q[0];
sx q[0];
rz(-2.7748845) q[0];
sx q[0];
rz(0.81272965) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.65152119) q[2];
sx q[2];
rz(-1.5805575) q[2];
sx q[2];
rz(-0.62192813) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9427467) q[1];
sx q[1];
rz(-1.2578576) q[1];
sx q[1];
rz(-1.320977) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7311677) q[3];
sx q[3];
rz(-0.99567185) q[3];
sx q[3];
rz(-0.23469532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(1.099951) q[2];
rz(2.3163017) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(-0.88551372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0734171) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(2.2391879) q[0];
rz(-2.1249318) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(-3.0117603) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79486217) q[0];
sx q[0];
rz(-2.4299893) q[0];
sx q[0];
rz(2.5559588) q[0];
rz(-pi) q[1];
rz(2.1461357) q[2];
sx q[2];
rz(-1.8838922) q[2];
sx q[2];
rz(0.42195937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1631158) q[1];
sx q[1];
rz(-0.69677959) q[1];
sx q[1];
rz(-2.5549868) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3744266) q[3];
sx q[3];
rz(-2.11103) q[3];
sx q[3];
rz(1.2353209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(-0.20425805) q[2];
rz(-1.9355109) q[3];
sx q[3];
rz(-1.5217425) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(1.4468505) q[0];
rz(-1.2591259) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(-0.68626219) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9335564) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(-2.3963388) q[0];
rz(2.1967728) q[2];
sx q[2];
rz(-2.2959024) q[2];
sx q[2];
rz(3.0996029) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8441019) q[1];
sx q[1];
rz(-0.19980783) q[1];
sx q[1];
rz(-1.954345) q[1];
rz(-pi) q[2];
rz(-1.2495743) q[3];
sx q[3];
rz(-1.9739082) q[3];
sx q[3];
rz(1.5954799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.69616047) q[2];
sx q[2];
rz(-1.3779209) q[2];
sx q[2];
rz(-0.0017722842) q[2];
rz(-0.56162515) q[3];
sx q[3];
rz(-0.91149819) q[3];
sx q[3];
rz(1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6034265) q[0];
sx q[0];
rz(-0.68646938) q[0];
sx q[0];
rz(-1.6954533) q[0];
rz(-0.7810477) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.5015645) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43270375) q[0];
sx q[0];
rz(-2.6484657) q[0];
sx q[0];
rz(1.5296442) q[0];
rz(-pi) q[1];
rz(3.115032) q[2];
sx q[2];
rz(-1.5090669) q[2];
sx q[2];
rz(-0.82690566) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6540263) q[1];
sx q[1];
rz(-1.8974202) q[1];
sx q[1];
rz(-1.6649654) q[1];
rz(2.1498508) q[3];
sx q[3];
rz(-1.0240882) q[3];
sx q[3];
rz(-2.9255097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.35187307) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(-1.8224576) q[2];
rz(1.9296648) q[3];
sx q[3];
rz(-1.2865678) q[3];
sx q[3];
rz(-0.31931988) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8050352) q[0];
sx q[0];
rz(-2.5890077) q[0];
sx q[0];
rz(-1.9375027) q[0];
rz(2.7583292) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(-2.7899172) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8512745) q[0];
sx q[0];
rz(-1.9418678) q[0];
sx q[0];
rz(-2.0593658) q[0];
x q[1];
rz(0.36231626) q[2];
sx q[2];
rz(-3*pi/13) q[2];
sx q[2];
rz(2.97646) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.6955399) q[1];
sx q[1];
rz(-0.76247588) q[1];
sx q[1];
rz(1.6474849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44378186) q[3];
sx q[3];
rz(-3.0203331) q[3];
sx q[3];
rz(1.4951984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.3433156) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(-1.2822255) q[2];
rz(1.6451689) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4984109) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(-2.9472651) q[0];
rz(2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(-2.1077572) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68034222) q[0];
sx q[0];
rz(-0.67078062) q[0];
sx q[0];
rz(-1.781342) q[0];
rz(-0.50844426) q[2];
sx q[2];
rz(-0.93548453) q[2];
sx q[2];
rz(-2.9490162) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2227877) q[1];
sx q[1];
rz(-2.5606887) q[1];
sx q[1];
rz(-1.3851628) q[1];
x q[2];
rz(-0.25063534) q[3];
sx q[3];
rz(-2.4126629) q[3];
sx q[3];
rz(-0.38754101) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0620492) q[2];
sx q[2];
rz(-0.94576183) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(-0.27030269) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.6132145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939659) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(-1.7059965) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(0.031899115) q[2];
sx q[2];
rz(-0.96822856) q[2];
sx q[2];
rz(-0.4005489) q[2];
rz(0.070449645) q[3];
sx q[3];
rz(-1.9000713) q[3];
sx q[3];
rz(-2.6303359) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];