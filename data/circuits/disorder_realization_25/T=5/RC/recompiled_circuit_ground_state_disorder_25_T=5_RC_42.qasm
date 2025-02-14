OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.3397665) q[0];
sx q[0];
rz(-0.56639329) q[0];
sx q[0];
rz(-1.9422148) q[0];
rz(-0.83130032) q[1];
sx q[1];
rz(4.755862) q[1];
sx q[1];
rz(8.9250467) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3008376) q[0];
sx q[0];
rz(-0.84512701) q[0];
sx q[0];
rz(-0.27780224) q[0];
rz(-1.209916) q[2];
sx q[2];
rz(-1.7583876) q[2];
sx q[2];
rz(-0.23439342) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9341683) q[1];
sx q[1];
rz(-0.8189924) q[1];
sx q[1];
rz(2.7374205) q[1];
rz(1.5152895) q[3];
sx q[3];
rz(-2.422547) q[3];
sx q[3];
rz(0.87874246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1685593) q[2];
sx q[2];
rz(-1.3848105) q[2];
sx q[2];
rz(-2.2200269) q[2];
rz(1.5714931) q[3];
sx q[3];
rz(-2.470033) q[3];
sx q[3];
rz(-2.4422372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(1.8916931) q[0];
sx q[0];
rz(-1.0279011) q[0];
sx q[0];
rz(1.6722884) q[0];
rz(-1.6647313) q[1];
sx q[1];
rz(-1.7988484) q[1];
sx q[1];
rz(0.5161759) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3625534) q[0];
sx q[0];
rz(-1.6677471) q[0];
sx q[0];
rz(-1.1753814) q[0];
rz(2.8968229) q[2];
sx q[2];
rz(-1.2483856) q[2];
sx q[2];
rz(2.4800194) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4975289) q[1];
sx q[1];
rz(-1.6402771) q[1];
sx q[1];
rz(-0.25150464) q[1];
x q[2];
rz(-2.1352876) q[3];
sx q[3];
rz(-0.7745452) q[3];
sx q[3];
rz(1.2513127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3790562) q[2];
sx q[2];
rz(-1.2006589) q[2];
sx q[2];
rz(0.29736796) q[2];
rz(2.6250046) q[3];
sx q[3];
rz(-1.9818431) q[3];
sx q[3];
rz(0.62492257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6487938) q[0];
sx q[0];
rz(-0.53557098) q[0];
sx q[0];
rz(-0.18185644) q[0];
rz(-0.44405538) q[1];
sx q[1];
rz(-2.2477138) q[1];
sx q[1];
rz(0.081238834) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67760175) q[0];
sx q[0];
rz(-0.39636546) q[0];
sx q[0];
rz(3.0667846) q[0];
rz(2.2611558) q[2];
sx q[2];
rz(-1.7882344) q[2];
sx q[2];
rz(-1.4354998) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4887152) q[1];
sx q[1];
rz(-1.3849568) q[1];
sx q[1];
rz(-2.04791) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0100097) q[3];
sx q[3];
rz(-1.6688235) q[3];
sx q[3];
rz(1.3658226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.9265201) q[2];
sx q[2];
rz(-0.36538616) q[2];
sx q[2];
rz(-0.60058769) q[2];
rz(-1.8961204) q[3];
sx q[3];
rz(-2.5966817) q[3];
sx q[3];
rz(3.1020402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0092225) q[0];
sx q[0];
rz(-2.4339269) q[0];
sx q[0];
rz(-0.02455499) q[0];
rz(0.37697667) q[1];
sx q[1];
rz(-1.8507277) q[1];
sx q[1];
rz(2.7797508) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.637476) q[0];
sx q[0];
rz(-1.2607122) q[0];
sx q[0];
rz(1.4092567) q[0];
x q[1];
rz(-1.3738427) q[2];
sx q[2];
rz(-0.3817238) q[2];
sx q[2];
rz(2.2854855) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6201219) q[1];
sx q[1];
rz(-1.748834) q[1];
sx q[1];
rz(-2.5012052) q[1];
rz(1.0933769) q[3];
sx q[3];
rz(-1.9481812) q[3];
sx q[3];
rz(1.6536825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.61818608) q[2];
sx q[2];
rz(-1.2647102) q[2];
sx q[2];
rz(-2.3183909) q[2];
rz(-0.19436714) q[3];
sx q[3];
rz(-0.40209642) q[3];
sx q[3];
rz(2.226734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0219367) q[0];
sx q[0];
rz(-1.5384262) q[0];
sx q[0];
rz(-2.5198779) q[0];
rz(-0.74553982) q[1];
sx q[1];
rz(-2.3108683) q[1];
sx q[1];
rz(-3.0527557) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1232077) q[0];
sx q[0];
rz(-1.6980217) q[0];
sx q[0];
rz(-1.5036262) q[0];
rz(0.88314573) q[2];
sx q[2];
rz(-2.7252203) q[2];
sx q[2];
rz(-1.0036381) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2076513) q[1];
sx q[1];
rz(-0.7867032) q[1];
sx q[1];
rz(-0.3982597) q[1];
rz(-pi) q[2];
rz(-1.9078498) q[3];
sx q[3];
rz(-1.4160857) q[3];
sx q[3];
rz(-2.0706994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.9125774) q[2];
sx q[2];
rz(-1.503016) q[2];
sx q[2];
rz(1.6055321) q[2];
rz(0.80794263) q[3];
sx q[3];
rz(-0.84351051) q[3];
sx q[3];
rz(-1.9690751) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.262893) q[0];
sx q[0];
rz(-1.7182588) q[0];
sx q[0];
rz(2.189157) q[0];
rz(1.7772504) q[1];
sx q[1];
rz(-1.9352501) q[1];
sx q[1];
rz(-0.84699026) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5085707) q[0];
sx q[0];
rz(-1.428033) q[0];
sx q[0];
rz(-1.4748013) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79361992) q[2];
sx q[2];
rz(-1.8066896) q[2];
sx q[2];
rz(-1.2737887) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.0274732) q[1];
sx q[1];
rz(-0.67016685) q[1];
sx q[1];
rz(0.7188188) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0383864) q[3];
sx q[3];
rz(-1.5852734) q[3];
sx q[3];
rz(-0.38179427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.18818894) q[2];
sx q[2];
rz(-2.0764515) q[2];
sx q[2];
rz(0.38189253) q[2];
rz(2.0415107) q[3];
sx q[3];
rz(-2.3298658) q[3];
sx q[3];
rz(0.41128099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(0.57322684) q[0];
sx q[0];
rz(-1.8170284) q[0];
sx q[0];
rz(-0.49760094) q[0];
rz(2.7865903) q[1];
sx q[1];
rz(-1.6604796) q[1];
sx q[1];
rz(-2.3752046) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8615211) q[0];
sx q[0];
rz(-2.6277271) q[0];
sx q[0];
rz(0.070529672) q[0];
x q[1];
rz(0.42986912) q[2];
sx q[2];
rz(-2.1235222) q[2];
sx q[2];
rz(-1.6195219) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.16535266) q[1];
sx q[1];
rz(-1.8596949) q[1];
sx q[1];
rz(1.6331178) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.4637632) q[3];
sx q[3];
rz(-1.7683709) q[3];
sx q[3];
rz(3.1238926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7327205) q[2];
sx q[2];
rz(-1.6072105) q[2];
sx q[2];
rz(-1.4917779) q[2];
rz(1.07897) q[3];
sx q[3];
rz(-1.3078728) q[3];
sx q[3];
rz(0.61416793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33477467) q[0];
sx q[0];
rz(-0.65009999) q[0];
sx q[0];
rz(0.41710576) q[0];
rz(-0.053622309) q[1];
sx q[1];
rz(-1.6157179) q[1];
sx q[1];
rz(-2.0448304) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8452783) q[0];
sx q[0];
rz(-2.1597549) q[0];
sx q[0];
rz(-2.9968569) q[0];
rz(-pi) q[1];
rz(-1.4545733) q[2];
sx q[2];
rz(-0.46346482) q[2];
sx q[2];
rz(1.6348686) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.027710304) q[1];
sx q[1];
rz(-1.7358297) q[1];
sx q[1];
rz(-1.5636958) q[1];
x q[2];
rz(-1.4194059) q[3];
sx q[3];
rz(-0.89821363) q[3];
sx q[3];
rz(-2.0167493) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2199478) q[2];
sx q[2];
rz(-1.3741477) q[2];
sx q[2];
rz(-2.5448223) q[2];
rz(2.2611332) q[3];
sx q[3];
rz(-1.7170649) q[3];
sx q[3];
rz(-0.60976353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6841458) q[0];
sx q[0];
rz(-0.54231751) q[0];
sx q[0];
rz(-0.31998262) q[0];
rz(0.34842247) q[1];
sx q[1];
rz(-0.44487822) q[1];
sx q[1];
rz(2.6382823) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.704858) q[0];
sx q[0];
rz(-2.9665369) q[0];
sx q[0];
rz(2.1268093) q[0];
rz(-pi) q[1];
rz(-2.636049) q[2];
sx q[2];
rz(-2.2525666) q[2];
sx q[2];
rz(0.56150061) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.383068) q[1];
sx q[1];
rz(-0.43979859) q[1];
sx q[1];
rz(-1.3851993) q[1];
x q[2];
rz(0.52127922) q[3];
sx q[3];
rz(-1.67799) q[3];
sx q[3];
rz(0.85338795) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2549501) q[2];
sx q[2];
rz(-2.6889668) q[2];
sx q[2];
rz(2.5060999) q[2];
rz(1.5358216) q[3];
sx q[3];
rz(-1.7255892) q[3];
sx q[3];
rz(-3.0543069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6814293) q[0];
sx q[0];
rz(-0.59578139) q[0];
sx q[0];
rz(-1.7728565) q[0];
rz(-2.6112556) q[1];
sx q[1];
rz(-1.6531205) q[1];
sx q[1];
rz(-2.368685) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3221098) q[0];
sx q[0];
rz(-0.81429309) q[0];
sx q[0];
rz(0.54634516) q[0];
x q[1];
rz(0.10112986) q[2];
sx q[2];
rz(-0.79532184) q[2];
sx q[2];
rz(0.9507319) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3858531) q[1];
sx q[1];
rz(-0.47702152) q[1];
sx q[1];
rz(-1.9128837) q[1];
rz(-2.021455) q[3];
sx q[3];
rz(-1.3067596) q[3];
sx q[3];
rz(2.9024189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3717926) q[2];
sx q[2];
rz(-2.1830406) q[2];
sx q[2];
rz(-2.8690763) q[2];
rz(0.46485999) q[3];
sx q[3];
rz(-1.9404989) q[3];
sx q[3];
rz(2.4242937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7335325) q[0];
sx q[0];
rz(-1.536674) q[0];
sx q[0];
rz(1.6996171) q[0];
rz(-1.5245262) q[1];
sx q[1];
rz(-2.1433612) q[1];
sx q[1];
rz(-2.8467766) q[1];
rz(1.3140368) q[2];
sx q[2];
rz(-1.2645559) q[2];
sx q[2];
rz(3.1380646) q[2];
rz(2.0130264) q[3];
sx q[3];
rz(-1.3696522) q[3];
sx q[3];
rz(-0.29396653) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
