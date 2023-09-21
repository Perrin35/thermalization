OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(1.8811037) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(-1.7927875) q[1];
sx q[1];
rz(-0.92372149) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7605654) q[0];
sx q[0];
rz(-2.0187223) q[0];
sx q[0];
rz(0.3737803) q[0];
x q[1];
rz(1.107723) q[2];
sx q[2];
rz(-2.827364) q[2];
sx q[2];
rz(-2.3067834) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.87286283) q[1];
sx q[1];
rz(-1.9951207) q[1];
sx q[1];
rz(-2.5904168) q[1];
rz(-pi) q[2];
x q[2];
rz(1.9507017) q[3];
sx q[3];
rz(-1.1067179) q[3];
sx q[3];
rz(-0.78117785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(2.9795734) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0682003) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(-0.67990047) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.686036) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50549492) q[0];
sx q[0];
rz(-2.539145) q[0];
sx q[0];
rz(1.2732182) q[0];
rz(2.7669719) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(-2.3945216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1728954) q[1];
sx q[1];
rz(-2.6186133) q[1];
sx q[1];
rz(-1.5437267) q[1];
rz(1.5442113) q[3];
sx q[3];
rz(-1.9823325) q[3];
sx q[3];
rz(-0.4707903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.42852795) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(-0.18243608) q[3];
sx q[3];
rz(-0.97674102) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(2.341111) q[0];
rz(-3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.172539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7045672) q[0];
sx q[0];
rz(-0.66172681) q[0];
sx q[0];
rz(-2.6252803) q[0];
rz(1.8981947) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(2.3936405) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0769656) q[1];
sx q[1];
rz(-2.5588227) q[1];
sx q[1];
rz(1.1175734) q[1];
rz(-pi) q[2];
rz(-1.3200687) q[3];
sx q[3];
rz(-1.6276976) q[3];
sx q[3];
rz(2.1123321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(0.95101142) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(-0.12810853) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(-2.6180843) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2142221) q[0];
sx q[0];
rz(-2.4282051) q[0];
sx q[0];
rz(0.58332304) q[0];
rz(-pi) q[1];
x q[1];
rz(0.33004327) q[2];
sx q[2];
rz(-1.651262) q[2];
sx q[2];
rz(0.45804322) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6701811) q[1];
sx q[1];
rz(-1.8243316) q[1];
sx q[1];
rz(-2.2815435) q[1];
rz(-pi) q[2];
rz(-1.773049) q[3];
sx q[3];
rz(-2.5407255) q[3];
sx q[3];
rz(-2.3893389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(0.28856746) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(0.55571663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(2.2619757) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-2.1496444) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0846227) q[0];
sx q[0];
rz(-1.7449433) q[0];
sx q[0];
rz(1.4625545) q[0];
x q[1];
rz(-1.5993824) q[2];
sx q[2];
rz(-0.30297849) q[2];
sx q[2];
rz(0.44705331) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2111152) q[1];
sx q[1];
rz(-0.31146295) q[1];
sx q[1];
rz(-3.009797) q[1];
rz(1.4335853) q[3];
sx q[3];
rz(-2.3793594) q[3];
sx q[3];
rz(1.2587794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(0.27080718) q[2];
rz(2.9233542) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72702423) q[0];
sx q[0];
rz(-0.71838656) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.7165002) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093229175) q[0];
sx q[0];
rz(-1.7943802) q[0];
sx q[0];
rz(2.6327052) q[0];
x q[1];
rz(-1.233333) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(2.595682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.9565935) q[1];
sx q[1];
rz(-1.6374267) q[1];
sx q[1];
rz(-2.0208298) q[1];
rz(-pi) q[2];
x q[2];
rz(0.96958843) q[3];
sx q[3];
rz(-2.0242656) q[3];
sx q[3];
rz(-0.78905247) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(2.9610736) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(-0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6329704) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1497027) q[0];
sx q[0];
rz(-2.2737962) q[0];
sx q[0];
rz(-2.4292612) q[0];
x q[1];
rz(-0.49634883) q[2];
sx q[2];
rz(-0.90663547) q[2];
sx q[2];
rz(-1.6830483) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.51470876) q[1];
sx q[1];
rz(-0.85299546) q[1];
sx q[1];
rz(1.2674598) q[1];
x q[2];
rz(-0.52845593) q[3];
sx q[3];
rz(-2.4848357) q[3];
sx q[3];
rz(-2.9627851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(1.5130419) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.725175) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-2.3186671) q[0];
rz(-2.8322463) q[1];
sx q[1];
rz(-1.7495218) q[1];
sx q[1];
rz(-1.3051422) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6179498) q[0];
sx q[0];
rz(-2.4542913) q[0];
sx q[0];
rz(-2.6108517) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.42758503) q[2];
sx q[2];
rz(-0.75140778) q[2];
sx q[2];
rz(2.4497355) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1888684) q[1];
sx q[1];
rz(-0.75718588) q[1];
sx q[1];
rz(-2.5785179) q[1];
x q[2];
rz(1.5042138) q[3];
sx q[3];
rz(-0.3399907) q[3];
sx q[3];
rz(0.41223994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(-1.4902327) q[2];
rz(2.0643318) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(-0.13154496) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3867144) q[0];
sx q[0];
rz(-1.8443549) q[0];
sx q[0];
rz(0.3219147) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-2.4386491) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2962869) q[0];
sx q[0];
rz(-1.1830813) q[0];
sx q[0];
rz(1.182343) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6417575) q[2];
sx q[2];
rz(-0.5103726) q[2];
sx q[2];
rz(-0.77529782) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.96303899) q[1];
sx q[1];
rz(-2.5759765) q[1];
sx q[1];
rz(1.8815243) q[1];
rz(-1.1144936) q[3];
sx q[3];
rz(-0.66702402) q[3];
sx q[3];
rz(-1.7175355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-2.8618331) q[2];
sx q[2];
rz(-1.3396324) q[2];
rz(-2.83589) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(-1.8113177) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16383485) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(-1.9158069) q[0];
rz(-2.2380791) q[1];
sx q[1];
rz(-2.5279896) q[1];
sx q[1];
rz(-0.46863619) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7077431) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(-1.7953403) q[0];
rz(-pi) q[1];
rz(-0.87601985) q[2];
sx q[2];
rz(-0.48601905) q[2];
sx q[2];
rz(-1.3181869) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.2965282) q[1];
sx q[1];
rz(-1.5025286) q[1];
sx q[1];
rz(1.3045842) q[1];
rz(-pi) q[2];
rz(0.43379421) q[3];
sx q[3];
rz(-1.9926096) q[3];
sx q[3];
rz(-1.0292366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.48352155) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(-1.998385) q[2];
rz(0.11463595) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(-1.5293998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6464012) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(0.62190965) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(-0.86482277) q[2];
sx q[2];
rz(-2.13158) q[2];
sx q[2];
rz(-2.0958015) q[2];
rz(0.91924304) q[3];
sx q[3];
rz(-1.6332492) q[3];
sx q[3];
rz(1.9132683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
