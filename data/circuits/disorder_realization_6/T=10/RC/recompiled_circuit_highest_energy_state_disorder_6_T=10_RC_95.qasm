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
rz(-2.1498635) q[0];
sx q[0];
rz(-2.5340134) q[0];
sx q[0];
rz(2.1460331) q[0];
rz(-2.6140656) q[1];
sx q[1];
rz(-1.0310643) q[1];
sx q[1];
rz(-0.48479015) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51895307) q[0];
sx q[0];
rz(-0.52871176) q[0];
sx q[0];
rz(-1.3081119) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.44102168) q[2];
sx q[2];
rz(-0.82953757) q[2];
sx q[2];
rz(-1.5953568) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.97580662) q[1];
sx q[1];
rz(-0.96250223) q[1];
sx q[1];
rz(3.0700141) q[1];
rz(-2.8504478) q[3];
sx q[3];
rz(-1.640794) q[3];
sx q[3];
rz(-1.1394046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.506044) q[2];
sx q[2];
rz(-0.94749331) q[2];
sx q[2];
rz(-0.49723899) q[2];
rz(-2.5908568) q[3];
sx q[3];
rz(-2.4516055) q[3];
sx q[3];
rz(-1.1373038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.9841442) q[0];
sx q[0];
rz(-1.4009615) q[0];
sx q[0];
rz(-2.3781811) q[0];
rz(2.5570671) q[1];
sx q[1];
rz(-1.3817363) q[1];
sx q[1];
rz(2.1302628) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4268734) q[0];
sx q[0];
rz(-1.701864) q[0];
sx q[0];
rz(0.1161954) q[0];
x q[1];
rz(-2.9052034) q[2];
sx q[2];
rz(-1.4589785) q[2];
sx q[2];
rz(-2.7674302) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4262169) q[1];
sx q[1];
rz(-1.9131887) q[1];
sx q[1];
rz(2.9222548) q[1];
x q[2];
rz(0.39803466) q[3];
sx q[3];
rz(-1.0938489) q[3];
sx q[3];
rz(-1.8534891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6611019) q[2];
sx q[2];
rz(-1.2031518) q[2];
sx q[2];
rz(2.1413546) q[2];
rz(-2.5601322) q[3];
sx q[3];
rz(-2.4088819) q[3];
sx q[3];
rz(1.2015517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8345399) q[0];
sx q[0];
rz(-0.59891278) q[0];
sx q[0];
rz(-0.55759984) q[0];
rz(-2.8343976) q[1];
sx q[1];
rz(-2.5418044) q[1];
sx q[1];
rz(-2.6934326) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0791866) q[0];
sx q[0];
rz(-2.6074479) q[0];
sx q[0];
rz(2.8791761) q[0];
x q[1];
rz(-1.0501625) q[2];
sx q[2];
rz(-1.7562281) q[2];
sx q[2];
rz(-0.38927573) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3002172) q[1];
sx q[1];
rz(-2.5102477) q[1];
sx q[1];
rz(-1.5247099) q[1];
x q[2];
rz(0.056164785) q[3];
sx q[3];
rz(-0.43496385) q[3];
sx q[3];
rz(0.54612225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8223411) q[2];
sx q[2];
rz(-0.71270573) q[2];
sx q[2];
rz(2.5437497) q[2];
rz(-0.44761014) q[3];
sx q[3];
rz(-1.8748583) q[3];
sx q[3];
rz(0.058535695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7287801) q[0];
sx q[0];
rz(-3.1106115) q[0];
sx q[0];
rz(-0.70984167) q[0];
rz(-2.1624883) q[1];
sx q[1];
rz(-0.85936463) q[1];
sx q[1];
rz(-1.3213347) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13793547) q[0];
sx q[0];
rz(-0.70323479) q[0];
sx q[0];
rz(-0.08589311) q[0];
x q[1];
rz(-2.3611907) q[2];
sx q[2];
rz(-1.681904) q[2];
sx q[2];
rz(1.358913) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.85985124) q[1];
sx q[1];
rz(-2.4415589) q[1];
sx q[1];
rz(1.9938049) q[1];
rz(-pi) q[2];
rz(2.0632486) q[3];
sx q[3];
rz(-0.60090099) q[3];
sx q[3];
rz(2.7067824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.1506302) q[2];
sx q[2];
rz(-1.4272828) q[2];
sx q[2];
rz(-2.4347351) q[2];
rz(1.7957211) q[3];
sx q[3];
rz(-1.6790877) q[3];
sx q[3];
rz(2.4367387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3020346) q[0];
sx q[0];
rz(-0.78080451) q[0];
sx q[0];
rz(-1.5020405) q[0];
rz(1.5976277) q[1];
sx q[1];
rz(-0.85999703) q[1];
sx q[1];
rz(1.4942716) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2199545) q[0];
sx q[0];
rz(-1.0003547) q[0];
sx q[0];
rz(-3.1396237) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5780894) q[2];
sx q[2];
rz(-2.1661758) q[2];
sx q[2];
rz(0.0059758107) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0922904) q[1];
sx q[1];
rz(-0.43152025) q[1];
sx q[1];
rz(-0.79705654) q[1];
rz(-pi) q[2];
rz(-0.95559883) q[3];
sx q[3];
rz(-1.5892059) q[3];
sx q[3];
rz(1.2332724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.8813701) q[2];
sx q[2];
rz(-1.2726731) q[2];
sx q[2];
rz(-0.55848813) q[2];
rz(2.6327366) q[3];
sx q[3];
rz(-1.6130092) q[3];
sx q[3];
rz(-1.2671027) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59638554) q[0];
sx q[0];
rz(-1.243243) q[0];
sx q[0];
rz(-0.16045706) q[0];
rz(-0.67333737) q[1];
sx q[1];
rz(-1.9247749) q[1];
sx q[1];
rz(2.014726) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50235275) q[0];
sx q[0];
rz(-1.9316088) q[0];
sx q[0];
rz(2.8176184) q[0];
rz(1.8931863) q[2];
sx q[2];
rz(-0.83288542) q[2];
sx q[2];
rz(-0.59925848) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4194132) q[1];
sx q[1];
rz(-0.99687885) q[1];
sx q[1];
rz(2.9001321) q[1];
x q[2];
rz(1.2282727) q[3];
sx q[3];
rz(-2.5728101) q[3];
sx q[3];
rz(0.82297135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0141853) q[2];
sx q[2];
rz(-1.5421474) q[2];
sx q[2];
rz(-2.3678153) q[2];
rz(-2.1524147) q[3];
sx q[3];
rz(-0.4072322) q[3];
sx q[3];
rz(-1.0729084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9147375) q[0];
sx q[0];
rz(-1.6586774) q[0];
sx q[0];
rz(0.7731272) q[0];
rz(1.5559366) q[1];
sx q[1];
rz(-1.2890041) q[1];
sx q[1];
rz(1.9645436) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0814514) q[0];
sx q[0];
rz(-2.3075275) q[0];
sx q[0];
rz(-3.0990776) q[0];
rz(-pi) q[1];
rz(-2.8297205) q[2];
sx q[2];
rz(-2.1215028) q[2];
sx q[2];
rz(0.68974173) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8428264) q[1];
sx q[1];
rz(-0.93888679) q[1];
sx q[1];
rz(-2.7999987) q[1];
rz(-pi) q[2];
rz(-2.9639951) q[3];
sx q[3];
rz(-1.7081722) q[3];
sx q[3];
rz(-3.0876189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.3561463) q[2];
sx q[2];
rz(-2.7768504) q[2];
sx q[2];
rz(-2.6981603) q[2];
rz(-2.7555079) q[3];
sx q[3];
rz(-1.723879) q[3];
sx q[3];
rz(-0.20080565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3767553) q[0];
sx q[0];
rz(-1.5926462) q[0];
sx q[0];
rz(3.1078597) q[0];
rz(2.4434166) q[1];
sx q[1];
rz(-2.5444784) q[1];
sx q[1];
rz(-0.66506344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9586104) q[0];
sx q[0];
rz(-2.0183051) q[0];
sx q[0];
rz(-0.29194269) q[0];
rz(-1.0475137) q[2];
sx q[2];
rz(-2.3436597) q[2];
sx q[2];
rz(2.8454859) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.60282502) q[1];
sx q[1];
rz(-2.1884349) q[1];
sx q[1];
rz(0.62299872) q[1];
rz(-pi) q[2];
rz(-0.31954664) q[3];
sx q[3];
rz(-1.5906153) q[3];
sx q[3];
rz(1.457959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.72208059) q[2];
sx q[2];
rz(-0.91700143) q[2];
sx q[2];
rz(1.3631932) q[2];
rz(2.6204387) q[3];
sx q[3];
rz(-1.8173953) q[3];
sx q[3];
rz(0.49191973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9687965) q[0];
sx q[0];
rz(-1.7954614) q[0];
sx q[0];
rz(0.81934339) q[0];
rz(2.4567545) q[1];
sx q[1];
rz(-1.2872773) q[1];
sx q[1];
rz(-0.12231621) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31560311) q[0];
sx q[0];
rz(-1.9902744) q[0];
sx q[0];
rz(2.2537838) q[0];
rz(-pi) q[1];
x q[1];
rz(2.165602) q[2];
sx q[2];
rz(-1.8907818) q[2];
sx q[2];
rz(-0.94727635) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6717447) q[1];
sx q[1];
rz(-1.2472595) q[1];
sx q[1];
rz(-2.5128415) q[1];
rz(-pi) q[2];
rz(-0.58845406) q[3];
sx q[3];
rz(-1.9663484) q[3];
sx q[3];
rz(1.5423397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.43716064) q[2];
sx q[2];
rz(-1.8081343) q[2];
sx q[2];
rz(-0.86755794) q[2];
rz(-1.7433172) q[3];
sx q[3];
rz(-0.38845348) q[3];
sx q[3];
rz(2.5625663) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45359465) q[0];
sx q[0];
rz(-1.9168251) q[0];
sx q[0];
rz(-1.0892185) q[0];
rz(-0.047529686) q[1];
sx q[1];
rz(-2.0952974) q[1];
sx q[1];
rz(0.68738031) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6733765) q[0];
sx q[0];
rz(-0.49096732) q[0];
sx q[0];
rz(-1.6223905) q[0];
x q[1];
rz(-2.0014761) q[2];
sx q[2];
rz(-1.3146558) q[2];
sx q[2];
rz(1.3821951) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.43722758) q[1];
sx q[1];
rz(-0.95672551) q[1];
sx q[1];
rz(2.1350056) q[1];
rz(3.0484588) q[3];
sx q[3];
rz(-2.2054005) q[3];
sx q[3];
rz(-0.30090162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1462732) q[2];
sx q[2];
rz(-0.61890382) q[2];
sx q[2];
rz(-1.6590365) q[2];
rz(-0.65685529) q[3];
sx q[3];
rz(-1.1494136) q[3];
sx q[3];
rz(2.759554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.252608) q[0];
sx q[0];
rz(-1.1532619) q[0];
sx q[0];
rz(1.1070195) q[0];
rz(-0.57869115) q[1];
sx q[1];
rz(-0.83732579) q[1];
sx q[1];
rz(2.8566828) q[1];
rz(0.35619932) q[2];
sx q[2];
rz(-1.648633) q[2];
sx q[2];
rz(1.6666694) q[2];
rz(-2.9978776) q[3];
sx q[3];
rz(-2.6201709) q[3];
sx q[3];
rz(-1.3377671) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
