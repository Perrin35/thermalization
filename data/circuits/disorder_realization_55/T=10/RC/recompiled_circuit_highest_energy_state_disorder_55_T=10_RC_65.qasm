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
rz(-2.6540304) q[0];
sx q[0];
rz(-2.7609479) q[0];
sx q[0];
rz(-0.08925499) q[0];
rz(2.6788977) q[1];
sx q[1];
rz(-0.2806288) q[1];
sx q[1];
rz(1.3624924) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2089522) q[0];
sx q[0];
rz(-0.37368837) q[0];
sx q[0];
rz(0.057500827) q[0];
rz(-pi) q[1];
rz(-0.2398165) q[2];
sx q[2];
rz(-0.15286286) q[2];
sx q[2];
rz(-2.3279362) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.67586865) q[1];
sx q[1];
rz(-1.7935408) q[1];
sx q[1];
rz(2.5506819) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4815349) q[3];
sx q[3];
rz(-2.0253455) q[3];
sx q[3];
rz(0.7346357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51916277) q[2];
sx q[2];
rz(-0.76202718) q[2];
sx q[2];
rz(0.92987522) q[2];
rz(2.8376288) q[3];
sx q[3];
rz(-1.8094939) q[3];
sx q[3];
rz(2.9728594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(2.8836483) q[0];
sx q[0];
rz(-1.5771663) q[0];
sx q[0];
rz(-1.8899348) q[0];
rz(2.9257863) q[1];
sx q[1];
rz(-2.1025175) q[1];
sx q[1];
rz(0.13914093) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.62876) q[0];
sx q[0];
rz(-0.91673512) q[0];
sx q[0];
rz(2.9399228) q[0];
rz(-pi) q[1];
rz(-0.49824841) q[2];
sx q[2];
rz(-1.3300657) q[2];
sx q[2];
rz(-1.5125076) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.58932585) q[1];
sx q[1];
rz(-1.3331474) q[1];
sx q[1];
rz(-0.71961211) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.81273791) q[3];
sx q[3];
rz(-2.2366282) q[3];
sx q[3];
rz(-2.8259371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.4591878) q[2];
sx q[2];
rz(-0.46899691) q[2];
sx q[2];
rz(-2.0501308) q[2];
rz(1.5952716) q[3];
sx q[3];
rz(-1.073444) q[3];
sx q[3];
rz(-0.44927621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.062155398) q[0];
sx q[0];
rz(-1.2261483) q[0];
sx q[0];
rz(0.011818258) q[0];
rz(0.91105175) q[1];
sx q[1];
rz(-1.751519) q[1];
sx q[1];
rz(1.8131088) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8810711) q[0];
sx q[0];
rz(-2.0876472) q[0];
sx q[0];
rz(-1.3840186) q[0];
rz(-2.6027868) q[2];
sx q[2];
rz(-2.9201815) q[2];
sx q[2];
rz(3.0543229) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6360259) q[1];
sx q[1];
rz(-0.88249245) q[1];
sx q[1];
rz(-2.0783552) q[1];
x q[2];
rz(-0.99081466) q[3];
sx q[3];
rz(-3.0085251) q[3];
sx q[3];
rz(-1.7855438) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.38982424) q[2];
sx q[2];
rz(-1.6168892) q[2];
sx q[2];
rz(-0.41979182) q[2];
rz(1.3587562) q[3];
sx q[3];
rz(-0.32521453) q[3];
sx q[3];
rz(1.9512008) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6841986) q[0];
sx q[0];
rz(-2.0036819) q[0];
sx q[0];
rz(2.6275291) q[0];
rz(-2.7534516) q[1];
sx q[1];
rz(-2.5297647) q[1];
sx q[1];
rz(1.911389) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3183751) q[0];
sx q[0];
rz(-1.8288695) q[0];
sx q[0];
rz(-0.31320546) q[0];
rz(-pi) q[1];
rz(2.7204334) q[2];
sx q[2];
rz(-0.76986662) q[2];
sx q[2];
rz(-2.056096) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.19691218) q[1];
sx q[1];
rz(-2.3397603) q[1];
sx q[1];
rz(-0.51735984) q[1];
rz(1.3421273) q[3];
sx q[3];
rz(-0.91311753) q[3];
sx q[3];
rz(2.2715037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.5019044) q[2];
sx q[2];
rz(-2.1268667) q[2];
sx q[2];
rz(0.91628966) q[2];
rz(0.4942975) q[3];
sx q[3];
rz(-1.9435792) q[3];
sx q[3];
rz(0.77427197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(0.91404) q[0];
sx q[0];
rz(-0.95714772) q[0];
sx q[0];
rz(-2.6555632) q[0];
rz(-2.5388429) q[1];
sx q[1];
rz(-1.9662247) q[1];
sx q[1];
rz(2.026162) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94447105) q[0];
sx q[0];
rz(-0.96707464) q[0];
sx q[0];
rz(2.05462) q[0];
rz(-2.7390295) q[2];
sx q[2];
rz(-0.83588119) q[2];
sx q[2];
rz(-1.0332359) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9420227) q[1];
sx q[1];
rz(-1.1567819) q[1];
sx q[1];
rz(-2.4041818) q[1];
x q[2];
rz(2.0511469) q[3];
sx q[3];
rz(-1.8060038) q[3];
sx q[3];
rz(-0.0071098162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9976161) q[2];
sx q[2];
rz(-0.3207427) q[2];
sx q[2];
rz(-2.0799267) q[2];
rz(-1.6210506) q[3];
sx q[3];
rz(-1.1989667) q[3];
sx q[3];
rz(-0.89402136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0735556) q[0];
sx q[0];
rz(-3.1017922) q[0];
sx q[0];
rz(1.2233446) q[0];
rz(2.6262737) q[1];
sx q[1];
rz(-2.4724908) q[1];
sx q[1];
rz(-1.3656778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66636234) q[0];
sx q[0];
rz(-1.7710553) q[0];
sx q[0];
rz(0.11563511) q[0];
x q[1];
rz(-1.1001105) q[2];
sx q[2];
rz(-2.9376279) q[2];
sx q[2];
rz(2.230913) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.4953782) q[1];
sx q[1];
rz(-1.7217024) q[1];
sx q[1];
rz(2.6617828) q[1];
rz(-2.8441098) q[3];
sx q[3];
rz(-1.5362036) q[3];
sx q[3];
rz(0.57719798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.078865163) q[2];
sx q[2];
rz(-1.8026423) q[2];
sx q[2];
rz(-1.3272746) q[2];
rz(1.4158538) q[3];
sx q[3];
rz(-3.0684107) q[3];
sx q[3];
rz(0.88254005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0075204) q[0];
sx q[0];
rz(-2.6044758) q[0];
sx q[0];
rz(2.9448331) q[0];
rz(-0.21601954) q[1];
sx q[1];
rz(-1.0088423) q[1];
sx q[1];
rz(2.3541727) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2077825) q[0];
sx q[0];
rz(-0.31126198) q[0];
sx q[0];
rz(-0.90988501) q[0];
rz(-pi) q[1];
rz(-0.49968736) q[2];
sx q[2];
rz(-0.56788194) q[2];
sx q[2];
rz(0.5918006) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1087139) q[1];
sx q[1];
rz(-1.6163663) q[1];
sx q[1];
rz(-0.5918064) q[1];
rz(-pi) q[2];
rz(-0.18936746) q[3];
sx q[3];
rz(-1.8872617) q[3];
sx q[3];
rz(0.53839436) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.54062033) q[2];
sx q[2];
rz(-1.629402) q[2];
sx q[2];
rz(-2.7579894) q[2];
rz(-2.2331734) q[3];
sx q[3];
rz(-0.81276613) q[3];
sx q[3];
rz(-0.15373357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28021321) q[0];
sx q[0];
rz(-0.52556831) q[0];
sx q[0];
rz(-0.63661611) q[0];
rz(0.55066291) q[1];
sx q[1];
rz(-1.6267585) q[1];
sx q[1];
rz(0.47264636) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.11926) q[0];
sx q[0];
rz(-0.68281931) q[0];
sx q[0];
rz(0.14108087) q[0];
x q[1];
rz(-1.6238975) q[2];
sx q[2];
rz(-1.6638954) q[2];
sx q[2];
rz(1.129221) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.9437406) q[1];
sx q[1];
rz(-1.5199558) q[1];
sx q[1];
rz(3.0089241) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63153679) q[3];
sx q[3];
rz(-1.8145921) q[3];
sx q[3];
rz(-2.7840028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7134646) q[2];
sx q[2];
rz(-1.8992004) q[2];
sx q[2];
rz(0.098527519) q[2];
rz(0.030700961) q[3];
sx q[3];
rz(-1.5701598) q[3];
sx q[3];
rz(2.9876685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6105662) q[0];
sx q[0];
rz(-0.95777804) q[0];
sx q[0];
rz(2.4089693) q[0];
rz(3.1351807) q[1];
sx q[1];
rz(-2.1156204) q[1];
sx q[1];
rz(1.4220672) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32860562) q[0];
sx q[0];
rz(-1.7610794) q[0];
sx q[0];
rz(-1.3178131) q[0];
rz(-pi) q[1];
rz(0.60958013) q[2];
sx q[2];
rz(-0.58660075) q[2];
sx q[2];
rz(-1.4966436) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8274424) q[1];
sx q[1];
rz(-1.6048417) q[1];
sx q[1];
rz(2.9877547) q[1];
x q[2];
rz(1.4957756) q[3];
sx q[3];
rz(-1.0524233) q[3];
sx q[3];
rz(-3.1278267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0942568) q[2];
sx q[2];
rz(-1.3553268) q[2];
sx q[2];
rz(-2.6195841) q[2];
rz(-3.1212854) q[3];
sx q[3];
rz(-2.4419407) q[3];
sx q[3];
rz(0.29683963) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2688399) q[0];
sx q[0];
rz(-1.5248542) q[0];
sx q[0];
rz(1.1312477) q[0];
rz(-0.77992431) q[1];
sx q[1];
rz(-1.9386539) q[1];
sx q[1];
rz(-0.47526971) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49233046) q[0];
sx q[0];
rz(-1.9340252) q[0];
sx q[0];
rz(2.338015) q[0];
rz(2.9067801) q[2];
sx q[2];
rz(-0.91660344) q[2];
sx q[2];
rz(1.7476817) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.14932528) q[1];
sx q[1];
rz(-1.7883446) q[1];
sx q[1];
rz(1.8294249) q[1];
rz(-pi) q[2];
rz(2.9679139) q[3];
sx q[3];
rz(-2.4376737) q[3];
sx q[3];
rz(2.6736163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.7705226) q[2];
sx q[2];
rz(-1.5757685) q[2];
sx q[2];
rz(1.9353665) q[2];
rz(1.7462339) q[3];
sx q[3];
rz(-2.5985056) q[3];
sx q[3];
rz(3.0622845) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64410011) q[0];
sx q[0];
rz(-0.19484367) q[0];
sx q[0];
rz(2.3051443) q[0];
rz(-2.1438228) q[1];
sx q[1];
rz(-1.0918959) q[1];
sx q[1];
rz(0.087654884) q[1];
rz(1.2682512) q[2];
sx q[2];
rz(-2.9131363) q[2];
sx q[2];
rz(2.7996488) q[2];
rz(-0.67881696) q[3];
sx q[3];
rz(-1.6573418) q[3];
sx q[3];
rz(-0.96924622) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
