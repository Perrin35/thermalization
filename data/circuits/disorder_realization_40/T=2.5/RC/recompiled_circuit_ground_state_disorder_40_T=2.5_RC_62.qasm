OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3482398) q[0];
sx q[0];
rz(-1.5982331) q[0];
sx q[0];
rz(-1.5016851) q[0];
rz(-1.3950672) q[1];
sx q[1];
rz(-0.39818624) q[1];
sx q[1];
rz(0.15287486) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76305394) q[0];
sx q[0];
rz(-1.4815271) q[0];
sx q[0];
rz(-1.7065926) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7136604) q[2];
sx q[2];
rz(-0.38201354) q[2];
sx q[2];
rz(-1.126542) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0331206) q[1];
sx q[1];
rz(-1.7594928) q[1];
sx q[1];
rz(0.061339247) q[1];
x q[2];
rz(2.6809815) q[3];
sx q[3];
rz(-2.1902927) q[3];
sx q[3];
rz(-1.7509489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.07831002) q[2];
sx q[2];
rz(-1.4274884) q[2];
sx q[2];
rz(0.59986344) q[2];
rz(-0.50987023) q[3];
sx q[3];
rz(-2.8190835) q[3];
sx q[3];
rz(0.16776423) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
sx q[2];
rz(-pi) q[2];
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
rz(-0.53304356) q[1];
sx q[1];
rz(0.60089111) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4054234) q[0];
sx q[0];
rz(-1.3440508) q[0];
sx q[0];
rz(1.9996243) q[0];
x q[1];
rz(2.6504806) q[2];
sx q[2];
rz(-2.4437856) q[2];
sx q[2];
rz(1.5581824) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0809853) q[1];
sx q[1];
rz(-1.391996) q[1];
sx q[1];
rz(-2.4455322) q[1];
rz(-2.8138312) q[3];
sx q[3];
rz(-0.53812384) q[3];
sx q[3];
rz(-1.6445352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.1945524) q[2];
sx q[2];
rz(-1.1636473) q[2];
sx q[2];
rz(-2.2770503) q[2];
rz(0.38327992) q[3];
sx q[3];
rz(-0.42208233) q[3];
sx q[3];
rz(2.2569136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29492038) q[0];
sx q[0];
rz(-2.8546794) q[0];
sx q[0];
rz(-0.82557803) q[0];
rz(0.83746743) q[1];
sx q[1];
rz(-1.0191963) q[1];
sx q[1];
rz(-0.79140633) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.237903) q[0];
sx q[0];
rz(-2.3775953) q[0];
sx q[0];
rz(-1.6873515) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7010274) q[2];
sx q[2];
rz(-1.4716867) q[2];
sx q[2];
rz(-2.2203022) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1371932) q[1];
sx q[1];
rz(-1.9158084) q[1];
sx q[1];
rz(2.5924204) q[1];
rz(-pi) q[2];
rz(2.0548925) q[3];
sx q[3];
rz(-1.181385) q[3];
sx q[3];
rz(2.2661346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.5241663) q[2];
sx q[2];
rz(-2.7238621) q[2];
sx q[2];
rz(1.2512655) q[2];
rz(1.9620126) q[3];
sx q[3];
rz(-1.1359295) q[3];
sx q[3];
rz(2.484926) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11611045) q[0];
sx q[0];
rz(-1.9132834) q[0];
sx q[0];
rz(0.81892282) q[0];
rz(0.081309155) q[1];
sx q[1];
rz(-0.58650494) q[1];
sx q[1];
rz(-2.3611045) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7538502) q[0];
sx q[0];
rz(-1.6293793) q[0];
sx q[0];
rz(-2.7669319) q[0];
x q[1];
rz(1.6967746) q[2];
sx q[2];
rz(-0.65444817) q[2];
sx q[2];
rz(-1.0140527) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7520903) q[1];
sx q[1];
rz(-2.6679975) q[1];
sx q[1];
rz(-2.447763) q[1];
rz(-1.9742161) q[3];
sx q[3];
rz(-1.8634923) q[3];
sx q[3];
rz(-3.0967874) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9235733) q[2];
sx q[2];
rz(-1.4998481) q[2];
sx q[2];
rz(1.8531331) q[2];
rz(1.4675379) q[3];
sx q[3];
rz(-0.54687423) q[3];
sx q[3];
rz(-0.43681496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9913637) q[0];
sx q[0];
rz(-0.65058351) q[0];
sx q[0];
rz(-0.22317602) q[0];
rz(3.0617833) q[1];
sx q[1];
rz(-0.97750074) q[1];
sx q[1];
rz(2.3056183) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88038091) q[0];
sx q[0];
rz(-1.1376842) q[0];
sx q[0];
rz(-0.26421896) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4673632) q[2];
sx q[2];
rz(-0.59665702) q[2];
sx q[2];
rz(0.7772738) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.0093352) q[1];
sx q[1];
rz(-2.461488) q[1];
sx q[1];
rz(-2.2538005) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96009843) q[3];
sx q[3];
rz(-1.7651422) q[3];
sx q[3];
rz(2.9604908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.23433805) q[2];
sx q[2];
rz(-2.2865488) q[2];
sx q[2];
rz(-1.3912531) q[2];
rz(1.2587345) q[3];
sx q[3];
rz(-1.3842868) q[3];
sx q[3];
rz(-0.39906991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.8626955) q[0];
sx q[0];
rz(-2.6628222) q[0];
sx q[0];
rz(-0.71267772) q[0];
rz(-1.124294) q[1];
sx q[1];
rz(-1.5702039) q[1];
sx q[1];
rz(-2.1052776) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9150805) q[0];
sx q[0];
rz(-2.4944502) q[0];
sx q[0];
rz(0.1948561) q[0];
rz(-pi) q[1];
rz(-2.9038139) q[2];
sx q[2];
rz(-1.8568123) q[2];
sx q[2];
rz(0.41508383) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.24806286) q[1];
sx q[1];
rz(-1.3643193) q[1];
sx q[1];
rz(-2.1724387) q[1];
rz(-pi) q[2];
x q[2];
rz(0.87392636) q[3];
sx q[3];
rz(-1.4269265) q[3];
sx q[3];
rz(1.914639) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0599646) q[2];
sx q[2];
rz(-2.5689503) q[2];
sx q[2];
rz(2.6247978) q[2];
rz(-1.146727) q[3];
sx q[3];
rz(-0.99249339) q[3];
sx q[3];
rz(0.049588047) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
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
rz(2.0674755) q[0];
sx q[0];
rz(-0.68454409) q[0];
sx q[0];
rz(1.6790947) q[0];
rz(-2.0430203) q[1];
sx q[1];
rz(-1.4417646) q[1];
sx q[1];
rz(2.3625653) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.532004) q[0];
sx q[0];
rz(-1.6250984) q[0];
sx q[0];
rz(-1.5691343) q[0];
rz(-pi) q[1];
rz(-2.674496) q[2];
sx q[2];
rz(-2.0990685) q[2];
sx q[2];
rz(-1.4548276) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.56151831) q[1];
sx q[1];
rz(-1.3653879) q[1];
sx q[1];
rz(-1.7465786) q[1];
rz(-2.0059465) q[3];
sx q[3];
rz(-2.3823934) q[3];
sx q[3];
rz(0.58461207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.74759185) q[2];
sx q[2];
rz(-0.70351768) q[2];
sx q[2];
rz(-0.35959378) q[2];
rz(-2.8408585) q[3];
sx q[3];
rz(-0.81529236) q[3];
sx q[3];
rz(2.1451779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71036285) q[0];
sx q[0];
rz(-1.2082986) q[0];
sx q[0];
rz(0.30074686) q[0];
rz(-2.6077479) q[1];
sx q[1];
rz(-1.0736991) q[1];
sx q[1];
rz(3.0278382) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1136459) q[0];
sx q[0];
rz(-0.7048713) q[0];
sx q[0];
rz(2.3990554) q[0];
rz(-pi) q[1];
rz(0.054461976) q[2];
sx q[2];
rz(-2.5321333) q[2];
sx q[2];
rz(-3.0888588) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99234527) q[1];
sx q[1];
rz(-0.94008821) q[1];
sx q[1];
rz(-2.135732) q[1];
rz(-2.9748125) q[3];
sx q[3];
rz(-1.2012225) q[3];
sx q[3];
rz(2.5931234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.45812312) q[2];
sx q[2];
rz(-2.9986431) q[2];
sx q[2];
rz(-2.6911823) q[2];
rz(-0.9149552) q[3];
sx q[3];
rz(-2.2136642) q[3];
sx q[3];
rz(1.3373059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8704855) q[0];
sx q[0];
rz(-2.4021689) q[0];
sx q[0];
rz(-0.41573218) q[0];
rz(0.43139002) q[1];
sx q[1];
rz(-0.58512551) q[1];
sx q[1];
rz(-2.364667) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3142648) q[0];
sx q[0];
rz(-1.6240472) q[0];
sx q[0];
rz(-3.1287249) q[0];
x q[1];
rz(0.84085744) q[2];
sx q[2];
rz(-2.1321725) q[2];
sx q[2];
rz(-0.64602588) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7207234) q[1];
sx q[1];
rz(-2.4178079) q[1];
sx q[1];
rz(-0.66026824) q[1];
rz(-2.6747392) q[3];
sx q[3];
rz(-1.8284208) q[3];
sx q[3];
rz(0.063916884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1315769) q[2];
sx q[2];
rz(-2.3516529) q[2];
sx q[2];
rz(-1.6291523) q[2];
rz(-0.43859628) q[3];
sx q[3];
rz(-2.3066543) q[3];
sx q[3];
rz(0.51630539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55385357) q[0];
sx q[0];
rz(-2.6022311) q[0];
sx q[0];
rz(-1.3158276) q[0];
rz(-3.059803) q[1];
sx q[1];
rz(-1.1809228) q[1];
sx q[1];
rz(1.2710424) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5499175) q[0];
sx q[0];
rz(-1.5999874) q[0];
sx q[0];
rz(1.1468588) q[0];
rz(-3.0767303) q[2];
sx q[2];
rz(-0.84170512) q[2];
sx q[2];
rz(0.87769485) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0473905) q[1];
sx q[1];
rz(-1.8592111) q[1];
sx q[1];
rz(0.80077313) q[1];
rz(-pi) q[2];
rz(0.2972352) q[3];
sx q[3];
rz(-0.085418022) q[3];
sx q[3];
rz(-0.88469317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.79909331) q[2];
sx q[2];
rz(-1.4405595) q[2];
sx q[2];
rz(-2.0889757) q[2];
rz(-1.3275702) q[3];
sx q[3];
rz(-1.9460287) q[3];
sx q[3];
rz(-0.50935203) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4159745) q[0];
sx q[0];
rz(-1.3798036) q[0];
sx q[0];
rz(1.8984541) q[0];
rz(1.3432518) q[1];
sx q[1];
rz(-2.6814798) q[1];
sx q[1];
rz(0.36569256) q[1];
rz(2.0475564) q[2];
sx q[2];
rz(-1.3972211) q[2];
sx q[2];
rz(0.76835189) q[2];
rz(0.030473564) q[3];
sx q[3];
rz(-0.41175523) q[3];
sx q[3];
rz(-2.336809) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
