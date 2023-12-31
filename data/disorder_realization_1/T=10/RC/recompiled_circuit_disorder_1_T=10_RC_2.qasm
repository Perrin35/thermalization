OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3172265) q[0];
sx q[0];
rz(-2.0269725) q[0];
sx q[0];
rz(0.00014076509) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(-1.1934086) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0187459) q[0];
sx q[0];
rz(-1.393264) q[0];
sx q[0];
rz(1.2794504) q[0];
rz(-pi) q[1];
rz(-2.675406) q[2];
sx q[2];
rz(-2.5417915) q[2];
sx q[2];
rz(-2.8592062) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.57230091) q[1];
sx q[1];
rz(-0.83218677) q[1];
sx q[1];
rz(-0.65170793) q[1];
rz(-1.7883349) q[3];
sx q[3];
rz(-1.6780403) q[3];
sx q[3];
rz(-1.6579962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.45941916) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(-1.2288644) q[2];
rz(-1.7284283) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(-1.4878954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5380149) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(-2.1287825) q[0];
rz(3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(1.123463) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9119308) q[0];
sx q[0];
rz(-2.6768885) q[0];
sx q[0];
rz(1.6994516) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.77983071) q[2];
sx q[2];
rz(-1.5260328) q[2];
sx q[2];
rz(-1.1080527) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.0151129) q[1];
sx q[1];
rz(-1.0547332) q[1];
sx q[1];
rz(0.59241398) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1334322) q[3];
sx q[3];
rz(-2.1481272) q[3];
sx q[3];
rz(1.2853704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.79364395) q[2];
sx q[2];
rz(-2.0517893) q[2];
sx q[2];
rz(0.91903764) q[2];
rz(-0.67409003) q[3];
sx q[3];
rz(-0.6522817) q[3];
sx q[3];
rz(-1.6154217) q[3];
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
rz(-pi) q[0];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27750257) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(1.8664237) q[0];
rz(0.69349849) q[1];
sx q[1];
rz(-1.2561412) q[1];
sx q[1];
rz(1.1330053) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4988574) q[0];
sx q[0];
rz(-0.14232902) q[0];
sx q[0];
rz(-1.5475153) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6373367) q[2];
sx q[2];
rz(-2.2849053) q[2];
sx q[2];
rz(2.2222663) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7743467) q[1];
sx q[1];
rz(-2.5173752) q[1];
sx q[1];
rz(-1.063785) q[1];
rz(-pi) q[2];
rz(2.5127701) q[3];
sx q[3];
rz(-1.0477133) q[3];
sx q[3];
rz(2.3716225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8901849) q[2];
sx q[2];
rz(-0.79139411) q[2];
sx q[2];
rz(1.8481002) q[2];
rz(-0.039316468) q[3];
sx q[3];
rz(-1.2189564) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.9807293) q[0];
rz(2.2456031) q[1];
sx q[1];
rz(-1.4410102) q[1];
sx q[1];
rz(-0.13555759) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7318152) q[0];
sx q[0];
rz(-1.8827794) q[0];
sx q[0];
rz(2.413495) q[0];
rz(-pi) q[1];
rz(0.074684871) q[2];
sx q[2];
rz(-1.9539781) q[2];
sx q[2];
rz(-0.89786868) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.63771026) q[1];
sx q[1];
rz(-1.4117068) q[1];
sx q[1];
rz(1.7590982) q[1];
x q[2];
rz(2.0849864) q[3];
sx q[3];
rz(-2.8286472) q[3];
sx q[3];
rz(1.3514951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.9049412) q[2];
sx q[2];
rz(-0.94649482) q[2];
sx q[2];
rz(2.2616852) q[2];
rz(3.0974292) q[3];
sx q[3];
rz(-1.5019838) q[3];
sx q[3];
rz(-0.28863171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1039466) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(-2.1283545) q[0];
rz(0.049731072) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(2.0577046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.097682) q[0];
sx q[0];
rz(-1.8143166) q[0];
sx q[0];
rz(2.9527412) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82917825) q[2];
sx q[2];
rz(-0.36311705) q[2];
sx q[2];
rz(1.5311637) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6127527) q[1];
sx q[1];
rz(-1.4626979) q[1];
sx q[1];
rz(-1.0428863) q[1];
rz(-pi) q[2];
rz(-2.6277222) q[3];
sx q[3];
rz(-1.6146982) q[3];
sx q[3];
rz(1.1843137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.23285398) q[2];
sx q[2];
rz(-0.32662699) q[2];
sx q[2];
rz(-2.8971635) q[2];
rz(-0.43236732) q[3];
sx q[3];
rz(-1.7418539) q[3];
sx q[3];
rz(2.6385245) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.4215707) q[0];
sx q[0];
rz(3.0474512) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-2.005902) q[1];
sx q[1];
rz(-2.24618) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46195128) q[0];
sx q[0];
rz(-1.2286751) q[0];
sx q[0];
rz(-1.530594) q[0];
x q[1];
rz(-2.5080639) q[2];
sx q[2];
rz(-1.3760566) q[2];
sx q[2];
rz(-2.5031236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0370334) q[1];
sx q[1];
rz(-1.4022786) q[1];
sx q[1];
rz(-3.0969572) q[1];
rz(-pi) q[2];
x q[2];
rz(1.496109) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(-2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.133693) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(2.3383979) q[2];
rz(-1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(2.7289594) q[3];
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
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-2.9769653) q[0];
sx q[0];
rz(0.51914006) q[0];
rz(0.58147645) q[1];
sx q[1];
rz(-1.1053718) q[1];
sx q[1];
rz(-1.2566459) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4706659) q[0];
sx q[0];
rz(-1.5304655) q[0];
sx q[0];
rz(-1.3192024) q[0];
rz(-pi) q[1];
rz(1.16876) q[2];
sx q[2];
rz(-2.6425344) q[2];
sx q[2];
rz(-1.4025584) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.7278442) q[1];
sx q[1];
rz(-1.2894948) q[1];
sx q[1];
rz(-1.6378228) q[1];
rz(-pi) q[2];
rz(2.7729633) q[3];
sx q[3];
rz(-2.3978524) q[3];
sx q[3];
rz(1.3524692) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1533623) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(-1.3640277) q[2];
rz(2.2310232) q[3];
sx q[3];
rz(-1.986859) q[3];
sx q[3];
rz(-1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0664739) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(2.8334154) q[0];
rz(-0.072487436) q[1];
sx q[1];
rz(-1.0132353) q[1];
sx q[1];
rz(0.38696188) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9403801) q[0];
sx q[0];
rz(-1.7965172) q[0];
sx q[0];
rz(2.1696027) q[0];
rz(-0.96078028) q[2];
sx q[2];
rz(-0.79723251) q[2];
sx q[2];
rz(-2.0955992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.19883979) q[1];
sx q[1];
rz(-1.6457335) q[1];
sx q[1];
rz(2.4692175) q[1];
rz(-pi) q[2];
rz(-2.7510795) q[3];
sx q[3];
rz(-0.53005855) q[3];
sx q[3];
rz(0.9583677) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.62362921) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(0.44000885) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0004262) q[0];
sx q[0];
rz(-0.74581242) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(0.72775841) q[1];
sx q[1];
rz(-2.7658503) q[1];
sx q[1];
rz(0.049302014) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4816149) q[0];
sx q[0];
rz(-1.0064631) q[0];
sx q[0];
rz(-0.83156395) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.49528893) q[2];
sx q[2];
rz(-2.1914748) q[2];
sx q[2];
rz(-2.4923313) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.1972563) q[1];
sx q[1];
rz(-2.0347056) q[1];
sx q[1];
rz(1.0524366) q[1];
rz(-pi) q[2];
rz(-2.6155869) q[3];
sx q[3];
rz(-0.81323871) q[3];
sx q[3];
rz(0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0631642) q[2];
sx q[2];
rz(-2.5421263) q[2];
sx q[2];
rz(0.72193974) q[2];
rz(-2.1980964) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(-0.25434428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9938875) q[0];
sx q[0];
rz(-1.1575971) q[0];
sx q[0];
rz(2.06185) q[0];
rz(-2.0823157) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(-1.7396897) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7125177) q[0];
sx q[0];
rz(-1.4445462) q[0];
sx q[0];
rz(1.7849126) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9548168) q[2];
sx q[2];
rz(-1.597607) q[2];
sx q[2];
rz(1.9865799) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.55587308) q[1];
sx q[1];
rz(-1.266529) q[1];
sx q[1];
rz(2.2833706) q[1];
rz(-2.510342) q[3];
sx q[3];
rz(-1.7389986) q[3];
sx q[3];
rz(-1.2390618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.4828651) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(1.6213017) q[2];
rz(0.55082095) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-2.4441392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.993492) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(-0.91611721) q[1];
sx q[1];
rz(-2.5506908) q[1];
sx q[1];
rz(2.5509902) q[1];
rz(-0.82952164) q[2];
sx q[2];
rz(-0.98486949) q[2];
sx q[2];
rz(-2.9688901) q[2];
rz(-0.017833088) q[3];
sx q[3];
rz(-2.087839) q[3];
sx q[3];
rz(2.9057333) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
