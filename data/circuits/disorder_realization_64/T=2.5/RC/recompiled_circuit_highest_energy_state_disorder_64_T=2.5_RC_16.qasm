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
rz(3.1177899) q[0];
sx q[0];
rz(-1.1060214) q[0];
sx q[0];
rz(-0.7769146) q[0];
rz(-1.7493526) q[1];
sx q[1];
rz(-1.8267781) q[1];
sx q[1];
rz(-2.1652752) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2433246) q[0];
sx q[0];
rz(-1.4372361) q[0];
sx q[0];
rz(1.1569174) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0357321) q[2];
sx q[2];
rz(-0.56625329) q[2];
sx q[2];
rz(0.98041269) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.0529671) q[1];
sx q[1];
rz(-1.8421116) q[1];
sx q[1];
rz(-2.0120828) q[1];
rz(-3.1093978) q[3];
sx q[3];
rz(-2.1436084) q[3];
sx q[3];
rz(-2.9762852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.50600791) q[2];
sx q[2];
rz(-3.0878461) q[2];
sx q[2];
rz(0.40679833) q[2];
rz(2.9721416) q[3];
sx q[3];
rz(-2.6120766) q[3];
sx q[3];
rz(1.0725526) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2605543) q[0];
sx q[0];
rz(-2.9091703) q[0];
sx q[0];
rz(-0.0037923092) q[0];
rz(3.0637528) q[1];
sx q[1];
rz(-2.4816315) q[1];
sx q[1];
rz(2.8357764) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90214848) q[0];
sx q[0];
rz(-1.3864707) q[0];
sx q[0];
rz(1.0499766) q[0];
x q[1];
rz(-1.4364834) q[2];
sx q[2];
rz(-0.82342734) q[2];
sx q[2];
rz(0.01625492) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.393049) q[1];
sx q[1];
rz(-2.3094588) q[1];
sx q[1];
rz(-2.119675) q[1];
rz(-pi) q[2];
rz(1.5068077) q[3];
sx q[3];
rz(-1.631569) q[3];
sx q[3];
rz(-1.9072541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.99016142) q[2];
sx q[2];
rz(-1.3913245) q[2];
sx q[2];
rz(-2.4988417) q[2];
rz(-3.0691872) q[3];
sx q[3];
rz(-2.0824771) q[3];
sx q[3];
rz(-1.36093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0533503) q[0];
sx q[0];
rz(-0.14277661) q[0];
sx q[0];
rz(-0.15637583) q[0];
rz(-3.1001672) q[1];
sx q[1];
rz(-2.5138469) q[1];
sx q[1];
rz(-1.5511537) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0254733) q[0];
sx q[0];
rz(-1.0295233) q[0];
sx q[0];
rz(-1.812029) q[0];
rz(2.4855108) q[2];
sx q[2];
rz(-1.5170043) q[2];
sx q[2];
rz(-1.0840067) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.52122766) q[1];
sx q[1];
rz(-1.4583734) q[1];
sx q[1];
rz(0.53037723) q[1];
rz(-0.82156397) q[3];
sx q[3];
rz(-2.2872777) q[3];
sx q[3];
rz(-0.91148538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7330043) q[2];
sx q[2];
rz(-0.53885794) q[2];
sx q[2];
rz(3.1206701) q[2];
rz(2.9544592) q[3];
sx q[3];
rz(-0.20550607) q[3];
sx q[3];
rz(0.11370295) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.1784096) q[0];
sx q[0];
rz(-0.33247501) q[0];
sx q[0];
rz(-0.43854976) q[0];
rz(1.5248388) q[1];
sx q[1];
rz(-0.33477819) q[1];
sx q[1];
rz(2.8964892) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49798508) q[0];
sx q[0];
rz(-1.6542572) q[0];
sx q[0];
rz(0.85026922) q[0];
rz(-pi) q[1];
rz(2.7510277) q[2];
sx q[2];
rz(-2.7465944) q[2];
sx q[2];
rz(-2.4633138) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41293834) q[1];
sx q[1];
rz(-2.1082343) q[1];
sx q[1];
rz(-0.16220233) q[1];
rz(-0.33632261) q[3];
sx q[3];
rz(-0.78236474) q[3];
sx q[3];
rz(-2.4049644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.601292) q[2];
sx q[2];
rz(-0.43593323) q[2];
sx q[2];
rz(2.8098246) q[2];
rz(0.48745421) q[3];
sx q[3];
rz(-1.0468227) q[3];
sx q[3];
rz(0.99307466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8496534) q[0];
sx q[0];
rz(-1.6878457) q[0];
sx q[0];
rz(-0.77350235) q[0];
rz(1.1812814) q[1];
sx q[1];
rz(-2.9995194) q[1];
sx q[1];
rz(1.389651) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5050801) q[0];
sx q[0];
rz(-0.55547041) q[0];
sx q[0];
rz(-2.4019) q[0];
rz(-pi) q[1];
rz(0.90221407) q[2];
sx q[2];
rz(-2.3093975) q[2];
sx q[2];
rz(0.26022831) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4331095) q[1];
sx q[1];
rz(-0.481284) q[1];
sx q[1];
rz(-0.17848707) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.23791194) q[3];
sx q[3];
rz(-1.7470659) q[3];
sx q[3];
rz(-0.59559866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.2191849) q[2];
sx q[2];
rz(-1.9318523) q[2];
sx q[2];
rz(-0.47214559) q[2];
rz(1.8426497) q[3];
sx q[3];
rz(-1.2932212) q[3];
sx q[3];
rz(0.81056547) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9206813) q[0];
sx q[0];
rz(-2.8784316) q[0];
sx q[0];
rz(-2.8780908) q[0];
rz(1.1031411) q[1];
sx q[1];
rz(-1.8270854) q[1];
sx q[1];
rz(-0.37364328) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220396) q[0];
sx q[0];
rz(-0.094498903) q[0];
sx q[0];
rz(2.2631133) q[0];
rz(-pi) q[1];
rz(-1.4336072) q[2];
sx q[2];
rz(-2.409916) q[2];
sx q[2];
rz(2.8525994) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.3195575) q[1];
sx q[1];
rz(-0.6614092) q[1];
sx q[1];
rz(-2.103785) q[1];
rz(-pi) q[2];
rz(-0.5881891) q[3];
sx q[3];
rz(-1.08687) q[3];
sx q[3];
rz(1.5671135) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0282447) q[2];
sx q[2];
rz(-2.9714606) q[2];
sx q[2];
rz(-0.55383468) q[2];
rz(1.7438186) q[3];
sx q[3];
rz(-0.60269409) q[3];
sx q[3];
rz(0.3127313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58436191) q[0];
sx q[0];
rz(-1.0065684) q[0];
sx q[0];
rz(-1.2297909) q[0];
rz(-2.9025485) q[1];
sx q[1];
rz(-1.6300647) q[1];
sx q[1];
rz(2.8412433) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2984788) q[0];
sx q[0];
rz(-1.4550147) q[0];
sx q[0];
rz(-2.9779469) q[0];
rz(0.70234583) q[2];
sx q[2];
rz(-1.1804198) q[2];
sx q[2];
rz(2.3927488) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9253941) q[1];
sx q[1];
rz(-1.555519) q[1];
sx q[1];
rz(3.141204) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2836779) q[3];
sx q[3];
rz(-1.5893717) q[3];
sx q[3];
rz(-1.1934848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.131669) q[2];
sx q[2];
rz(-1.4574304) q[2];
sx q[2];
rz(-0.24492502) q[2];
rz(0.51982546) q[3];
sx q[3];
rz(-0.86331415) q[3];
sx q[3];
rz(0.68827099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-0.35172611) q[0];
sx q[0];
rz(-2.7930197) q[0];
sx q[0];
rz(-1.1630195) q[0];
rz(3.0746958) q[1];
sx q[1];
rz(-1.4935378) q[1];
sx q[1];
rz(-1.012872) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5236008) q[0];
sx q[0];
rz(-0.80192425) q[0];
sx q[0];
rz(-2.2973934) q[0];
rz(0.17872058) q[2];
sx q[2];
rz(-2.131049) q[2];
sx q[2];
rz(2.5541039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2248963) q[1];
sx q[1];
rz(-2.6768497) q[1];
sx q[1];
rz(-0.76239379) q[1];
rz(-pi) q[2];
rz(0.49533923) q[3];
sx q[3];
rz(-2.3916349) q[3];
sx q[3];
rz(2.0928252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.11671994) q[2];
sx q[2];
rz(-0.97815424) q[2];
sx q[2];
rz(2.8187974) q[2];
rz(2.5358477) q[3];
sx q[3];
rz(-2.3492458) q[3];
sx q[3];
rz(0.34887031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.054319687) q[0];
sx q[0];
rz(-0.1467341) q[0];
sx q[0];
rz(-0.012454575) q[0];
rz(-2.3948578) q[1];
sx q[1];
rz(-0.92339271) q[1];
sx q[1];
rz(0.27997231) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0736948) q[0];
sx q[0];
rz(-1.520074) q[0];
sx q[0];
rz(1.6832666) q[0];
rz(2.8146663) q[2];
sx q[2];
rz(-1.9844311) q[2];
sx q[2];
rz(1.7978316) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.88290434) q[1];
sx q[1];
rz(-1.8778442) q[1];
sx q[1];
rz(-0.76489246) q[1];
rz(-pi) q[2];
rz(2.8948726) q[3];
sx q[3];
rz(-1.9233875) q[3];
sx q[3];
rz(0.22724928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3296457) q[2];
sx q[2];
rz(-0.75355607) q[2];
sx q[2];
rz(0.21198708) q[2];
rz(-0.82344615) q[3];
sx q[3];
rz(-1.4444838) q[3];
sx q[3];
rz(-2.8823891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.14281808) q[0];
sx q[0];
rz(-3.0840254) q[0];
sx q[0];
rz(0.69277358) q[0];
rz(-0.57299262) q[1];
sx q[1];
rz(-1.7968105) q[1];
sx q[1];
rz(-0.43100345) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13700542) q[0];
sx q[0];
rz(-1.0197029) q[0];
sx q[0];
rz(-1.5020834) q[0];
rz(-2.0153322) q[2];
sx q[2];
rz(-2.1241786) q[2];
sx q[2];
rz(1.044342) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.17934701) q[1];
sx q[1];
rz(-1.822346) q[1];
sx q[1];
rz(-0.72466447) q[1];
x q[2];
rz(-2.1206843) q[3];
sx q[3];
rz(-1.8326933) q[3];
sx q[3];
rz(-2.6481215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7196322) q[2];
sx q[2];
rz(-2.8675291) q[2];
sx q[2];
rz(-2.5813622) q[2];
rz(-0.46960056) q[3];
sx q[3];
rz(-0.40265366) q[3];
sx q[3];
rz(0.71389055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2847168) q[0];
sx q[0];
rz(-1.7354043) q[0];
sx q[0];
rz(2.0866557) q[0];
rz(-0.84125413) q[1];
sx q[1];
rz(-2.0396736) q[1];
sx q[1];
rz(3.094818) q[1];
rz(-0.31568676) q[2];
sx q[2];
rz(-1.9236947) q[2];
sx q[2];
rz(-2.535939) q[2];
rz(0.026997707) q[3];
sx q[3];
rz(-2.229573) q[3];
sx q[3];
rz(0.98462478) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
