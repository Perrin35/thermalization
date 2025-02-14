OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.1281066) q[0];
sx q[0];
rz(3.9689316) q[0];
sx q[0];
rz(7.6710424) q[0];
rz(0.85343051) q[1];
sx q[1];
rz(5.8813385) q[1];
sx q[1];
rz(11.452236) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46060503) q[0];
sx q[0];
rz(-1.3189684) q[0];
sx q[0];
rz(2.2973209) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9164785) q[2];
sx q[2];
rz(-0.71768453) q[2];
sx q[2];
rz(2.8717741) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6932363) q[1];
sx q[1];
rz(-0.020388842) q[1];
sx q[1];
rz(-1.8113813) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2306289) q[3];
sx q[3];
rz(-1.0940043) q[3];
sx q[3];
rz(1.7309675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9806597) q[2];
sx q[2];
rz(-0.42155835) q[2];
sx q[2];
rz(-1.4483615) q[2];
rz(-2.8990922) q[3];
sx q[3];
rz(-0.91553965) q[3];
sx q[3];
rz(-1.8814794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5777957) q[0];
sx q[0];
rz(-1.3160492) q[0];
sx q[0];
rz(1.0003566) q[0];
rz(-0.73161221) q[1];
sx q[1];
rz(-1.4294521) q[1];
sx q[1];
rz(1.9614722) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2041572) q[0];
sx q[0];
rz(-1.674386) q[0];
sx q[0];
rz(-1.4803463) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6890516) q[2];
sx q[2];
rz(-0.27980556) q[2];
sx q[2];
rz(0.56820936) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0958671) q[1];
sx q[1];
rz(-1.5689227) q[1];
sx q[1];
rz(1.4366524) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3445156) q[3];
sx q[3];
rz(-2.7555572) q[3];
sx q[3];
rz(-0.71938709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.09329) q[2];
sx q[2];
rz(-0.90128171) q[2];
sx q[2];
rz(1.0789336) q[2];
rz(-0.087873936) q[3];
sx q[3];
rz(-1.7885957) q[3];
sx q[3];
rz(2.8972076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3946149) q[0];
sx q[0];
rz(-1.9786388) q[0];
sx q[0];
rz(-1.8633457) q[0];
rz(-2.6784015) q[1];
sx q[1];
rz(-1.3563503) q[1];
sx q[1];
rz(2.3203704) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0671135) q[0];
sx q[0];
rz(-2.4601106) q[0];
sx q[0];
rz(-1.2740178) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6326007) q[2];
sx q[2];
rz(-1.9414177) q[2];
sx q[2];
rz(2.6573617) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1248904) q[1];
sx q[1];
rz(-1.115106) q[1];
sx q[1];
rz(-1.7348921) q[1];
rz(-pi) q[2];
rz(-2.3825112) q[3];
sx q[3];
rz(-2.0584752) q[3];
sx q[3];
rz(1.9239192) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.897573) q[2];
sx q[2];
rz(-1.579318) q[2];
sx q[2];
rz(-2.139034) q[2];
rz(-1.2282486) q[3];
sx q[3];
rz(-1.7398261) q[3];
sx q[3];
rz(1.7206515) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9035852) q[0];
sx q[0];
rz(-1.0500195) q[0];
sx q[0];
rz(0.28582698) q[0];
rz(0.26126513) q[1];
sx q[1];
rz(-1.8087872) q[1];
sx q[1];
rz(-3.1108943) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0945902) q[0];
sx q[0];
rz(-1.6436716) q[0];
sx q[0];
rz(-2.8776309) q[0];
x q[1];
rz(-1.4723563) q[2];
sx q[2];
rz(-2.4623496) q[2];
sx q[2];
rz(-1.6427276) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.6469137) q[1];
sx q[1];
rz(-2.621197) q[1];
sx q[1];
rz(-2.097165) q[1];
rz(-0.081153374) q[3];
sx q[3];
rz(-1.8834682) q[3];
sx q[3];
rz(2.816538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3280481) q[2];
sx q[2];
rz(-1.8009461) q[2];
sx q[2];
rz(0.0030585232) q[2];
rz(2.3004153) q[3];
sx q[3];
rz(-0.58924651) q[3];
sx q[3];
rz(2.7868748) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67741126) q[0];
sx q[0];
rz(-2.2941636) q[0];
sx q[0];
rz(1.6743976) q[0];
rz(-1.5122308) q[1];
sx q[1];
rz(-2.1758175) q[1];
sx q[1];
rz(-0.95219749) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5291572) q[0];
sx q[0];
rz(-1.5780492) q[0];
sx q[0];
rz(-3.0713505) q[0];
rz(1.1811851) q[2];
sx q[2];
rz(-0.91105538) q[2];
sx q[2];
rz(-0.3241186) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.72022382) q[1];
sx q[1];
rz(-1.6744057) q[1];
sx q[1];
rz(3.1006579) q[1];
rz(2.5467993) q[3];
sx q[3];
rz(-2.1457971) q[3];
sx q[3];
rz(-1.2269542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3516922) q[2];
sx q[2];
rz(-1.5396427) q[2];
sx q[2];
rz(2.6969625) q[2];
rz(2.6897258) q[3];
sx q[3];
rz(-0.63932747) q[3];
sx q[3];
rz(-0.062084559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43279466) q[0];
sx q[0];
rz(-2.4169156) q[0];
sx q[0];
rz(-2.2213347) q[0];
rz(2.1814749) q[1];
sx q[1];
rz(-1.3939539) q[1];
sx q[1];
rz(-2.6422909) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0408142) q[0];
sx q[0];
rz(-1.4199323) q[0];
sx q[0];
rz(2.8013381) q[0];
x q[1];
rz(-0.66125691) q[2];
sx q[2];
rz(-1.7098797) q[2];
sx q[2];
rz(-2.5360803) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.074452049) q[1];
sx q[1];
rz(-2.4765091) q[1];
sx q[1];
rz(0.52560271) q[1];
rz(-pi) q[2];
rz(-0.74120993) q[3];
sx q[3];
rz(-1.831372) q[3];
sx q[3];
rz(1.8678566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3906735) q[2];
sx q[2];
rz(-1.6030739) q[2];
sx q[2];
rz(1.0323367) q[2];
rz(-2.27683) q[3];
sx q[3];
rz(-1.7059749) q[3];
sx q[3];
rz(-0.56149703) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.432935) q[0];
sx q[0];
rz(-1.189804) q[0];
sx q[0];
rz(-1.3849965) q[0];
rz(2.2191091) q[1];
sx q[1];
rz(-1.9360767) q[1];
sx q[1];
rz(-0.14499697) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2534197) q[0];
sx q[0];
rz(-1.7494666) q[0];
sx q[0];
rz(2.9736142) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.85741557) q[2];
sx q[2];
rz(-1.5930297) q[2];
sx q[2];
rz(-2.2043101) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0516982) q[1];
sx q[1];
rz(-1.6015918) q[1];
sx q[1];
rz(-1.4975862) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0078405) q[3];
sx q[3];
rz(-0.87742311) q[3];
sx q[3];
rz(-2.106732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.9292235) q[2];
sx q[2];
rz(-0.23304686) q[2];
sx q[2];
rz(-1.2082427) q[2];
rz(0.82331795) q[3];
sx q[3];
rz(-1.1813141) q[3];
sx q[3];
rz(1.7821144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.061148297) q[0];
sx q[0];
rz(-2.3539982) q[0];
sx q[0];
rz(-2.0620692) q[0];
rz(2.1615255) q[1];
sx q[1];
rz(-1.020224) q[1];
sx q[1];
rz(-0.10990873) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82775527) q[0];
sx q[0];
rz(-1.2775363) q[0];
sx q[0];
rz(-0.37773962) q[0];
rz(-pi) q[1];
rz(-2.8632322) q[2];
sx q[2];
rz(-1.5794069) q[2];
sx q[2];
rz(1.4573163) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.76131436) q[1];
sx q[1];
rz(-1.2732693) q[1];
sx q[1];
rz(0.44400906) q[1];
x q[2];
rz(-1.9426467) q[3];
sx q[3];
rz(-2.5805485) q[3];
sx q[3];
rz(-2.4435465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93426934) q[2];
sx q[2];
rz(-1.0730275) q[2];
sx q[2];
rz(0.70460021) q[2];
rz(-0.28471026) q[3];
sx q[3];
rz(-2.4094818) q[3];
sx q[3];
rz(0.43340161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1339742) q[0];
sx q[0];
rz(-1.3173137) q[0];
sx q[0];
rz(0.91868573) q[0];
rz(-2.6559415) q[1];
sx q[1];
rz(-0.44356569) q[1];
sx q[1];
rz(2.0408911) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31195413) q[0];
sx q[0];
rz(-1.432278) q[0];
sx q[0];
rz(-1.7111219) q[0];
rz(-pi) q[1];
rz(-1.8728016) q[2];
sx q[2];
rz(-1.7922316) q[2];
sx q[2];
rz(1.0332274) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.13362836) q[1];
sx q[1];
rz(-1.4660264) q[1];
sx q[1];
rz(1.973335) q[1];
rz(-pi) q[2];
rz(1.507255) q[3];
sx q[3];
rz(-2.2054234) q[3];
sx q[3];
rz(-1.5295446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6653768) q[2];
sx q[2];
rz(-2.2874338) q[2];
sx q[2];
rz(-0.7594792) q[2];
rz(-2.1935943) q[3];
sx q[3];
rz(-1.4576603) q[3];
sx q[3];
rz(-1.6215526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.043902472) q[0];
sx q[0];
rz(-2.6444785) q[0];
sx q[0];
rz(-0.33101606) q[0];
rz(-2.379592) q[1];
sx q[1];
rz(-2.8490729) q[1];
sx q[1];
rz(-0.62823546) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2822489) q[0];
sx q[0];
rz(-0.64447908) q[0];
sx q[0];
rz(2.6656239) q[0];
x q[1];
rz(-1.3359137) q[2];
sx q[2];
rz(-1.8705436) q[2];
sx q[2];
rz(-0.91938726) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3207631) q[1];
sx q[1];
rz(-2.4270128) q[1];
sx q[1];
rz(2.5933867) q[1];
rz(-pi) q[2];
rz(-1.2248433) q[3];
sx q[3];
rz(-1.1802434) q[3];
sx q[3];
rz(-1.8486763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1198173) q[2];
sx q[2];
rz(-2.6675318) q[2];
sx q[2];
rz(-2.9277756) q[2];
rz(2.1553433) q[3];
sx q[3];
rz(-1.2912368) q[3];
sx q[3];
rz(-3.0333062) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4719791) q[0];
sx q[0];
rz(-1.6983953) q[0];
sx q[0];
rz(-1.4629913) q[0];
rz(-1.1769453) q[1];
sx q[1];
rz(-1.4851478) q[1];
sx q[1];
rz(-3.1335395) q[1];
rz(2.0811002) q[2];
sx q[2];
rz(-2.7955187) q[2];
sx q[2];
rz(-2.1375755) q[2];
rz(3.0358992) q[3];
sx q[3];
rz(-0.65423818) q[3];
sx q[3];
rz(-0.26442179) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
