OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.174515) q[0];
sx q[0];
rz(-1.6835901) q[0];
sx q[0];
rz(2.8414677) q[0];
rz(1.1974273) q[1];
sx q[1];
rz(4.6680968) q[1];
sx q[1];
rz(11.065281) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5434721) q[0];
sx q[0];
rz(-1.5556922) q[0];
sx q[0];
rz(-1.5746818) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9674703) q[2];
sx q[2];
rz(-0.80840092) q[2];
sx q[2];
rz(1.8044745) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.38084322) q[1];
sx q[1];
rz(-2.2135206) q[1];
sx q[1];
rz(1.2973143) q[1];
x q[2];
rz(-0.80880161) q[3];
sx q[3];
rz(-1.5743269) q[3];
sx q[3];
rz(0.88909066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.7370721) q[2];
sx q[2];
rz(-1.9998735) q[2];
sx q[2];
rz(2.4386151) q[2];
rz(-2.5718555) q[3];
sx q[3];
rz(-2.2629786) q[3];
sx q[3];
rz(1.5841293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0074145929) q[0];
sx q[0];
rz(-2.5226722) q[0];
sx q[0];
rz(-1.2331569) q[0];
rz(2.5470219) q[1];
sx q[1];
rz(-1.0392799) q[1];
sx q[1];
rz(-1.0979244) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1933408) q[0];
sx q[0];
rz(-1.7714714) q[0];
sx q[0];
rz(-3.0737123) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6498943) q[2];
sx q[2];
rz(-1.227081) q[2];
sx q[2];
rz(-0.51174405) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.2757551) q[1];
sx q[1];
rz(-3.0601353) q[1];
sx q[1];
rz(-1.8862731) q[1];
rz(0.091836496) q[3];
sx q[3];
rz(-0.79376924) q[3];
sx q[3];
rz(-2.5240999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7132831) q[2];
sx q[2];
rz(-1.2998591) q[2];
sx q[2];
rz(-1.606288) q[2];
rz(-0.71896583) q[3];
sx q[3];
rz(-2.1956367) q[3];
sx q[3];
rz(-2.9276221) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83800256) q[0];
sx q[0];
rz(-0.8135697) q[0];
sx q[0];
rz(-2.549262) q[0];
rz(1.8968286) q[1];
sx q[1];
rz(-2.0015621) q[1];
sx q[1];
rz(0.14561428) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1942822) q[0];
sx q[0];
rz(-2.6366451) q[0];
sx q[0];
rz(0.28316747) q[0];
x q[1];
rz(-1.4136366) q[2];
sx q[2];
rz(-2.4470235) q[2];
sx q[2];
rz(-0.63169152) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1864724) q[1];
sx q[1];
rz(-1.8545517) q[1];
sx q[1];
rz(-0.99631817) q[1];
x q[2];
rz(1.0823194) q[3];
sx q[3];
rz(-1.8104189) q[3];
sx q[3];
rz(2.3828196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.40272063) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(1.2325475) q[2];
rz(0.23972073) q[3];
sx q[3];
rz(-2.0870049) q[3];
sx q[3];
rz(2.6565552) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0204912) q[0];
sx q[0];
rz(-0.70206577) q[0];
sx q[0];
rz(-0.030601587) q[0];
rz(0.55514151) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(2.0487002) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7858148) q[0];
sx q[0];
rz(-0.91453881) q[0];
sx q[0];
rz(2.2758765) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6299155) q[2];
sx q[2];
rz(-2.335304) q[2];
sx q[2];
rz(1.3313683) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38072571) q[1];
sx q[1];
rz(-0.66752258) q[1];
sx q[1];
rz(-1.236626) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.67375848) q[3];
sx q[3];
rz(-2.1762037) q[3];
sx q[3];
rz(-2.377569) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0838202) q[2];
sx q[2];
rz(-1.3670992) q[2];
sx q[2];
rz(2.8687381) q[2];
rz(-1.3911635) q[3];
sx q[3];
rz(-2.302156) q[3];
sx q[3];
rz(-2.1896037) q[3];
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
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6735753) q[0];
sx q[0];
rz(-1.3382358) q[0];
sx q[0];
rz(-1.8433628) q[0];
rz(-0.6908373) q[1];
sx q[1];
rz(-1.1996484) q[1];
sx q[1];
rz(1.615049) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95193255) q[0];
sx q[0];
rz(-2.2554923) q[0];
sx q[0];
rz(-0.71345274) q[0];
rz(-pi) q[1];
rz(-2.7886018) q[2];
sx q[2];
rz(-1.304875) q[2];
sx q[2];
rz(-0.28133767) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4936798) q[1];
sx q[1];
rz(-0.70131749) q[1];
sx q[1];
rz(1.5307242) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.2722716) q[3];
sx q[3];
rz(-1.2236475) q[3];
sx q[3];
rz(2.783028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.4096628) q[2];
sx q[2];
rz(-0.80594984) q[2];
sx q[2];
rz(-0.16732495) q[2];
rz(1.3903728) q[3];
sx q[3];
rz(-2.6087587) q[3];
sx q[3];
rz(2.4840241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.14199) q[0];
sx q[0];
rz(-0.36407343) q[0];
sx q[0];
rz(-1.2784736) q[0];
rz(-1.4356042) q[1];
sx q[1];
rz(-1.4878788) q[1];
sx q[1];
rz(-2.7659168) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3669691) q[0];
sx q[0];
rz(-1.7967829) q[0];
sx q[0];
rz(1.8480728) q[0];
rz(-0.85083811) q[2];
sx q[2];
rz(-0.22432835) q[2];
sx q[2];
rz(1.4590603) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.7573812) q[1];
sx q[1];
rz(-2.7486292) q[1];
sx q[1];
rz(2.1860366) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.2578854) q[3];
sx q[3];
rz(-0.77270672) q[3];
sx q[3];
rz(-1.4290641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0756691) q[2];
sx q[2];
rz(-2.0772987) q[2];
sx q[2];
rz(-2.2806878) q[2];
rz(-1.143645) q[3];
sx q[3];
rz(-0.30503169) q[3];
sx q[3];
rz(-2.8633269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.025573108) q[0];
sx q[0];
rz(-1.8208068) q[0];
sx q[0];
rz(-2.386911) q[0];
rz(2.2017551) q[1];
sx q[1];
rz(-0.87516963) q[1];
sx q[1];
rz(1.2831203) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1719753) q[0];
sx q[0];
rz(-1.2640868) q[0];
sx q[0];
rz(0.0071010751) q[0];
rz(-pi) q[1];
rz(-2.2290984) q[2];
sx q[2];
rz(-2.5844838) q[2];
sx q[2];
rz(2.0832286) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2312113) q[1];
sx q[1];
rz(-1.2802135) q[1];
sx q[1];
rz(2.742275) q[1];
rz(-pi) q[2];
rz(-0.51605172) q[3];
sx q[3];
rz(-2.1911804) q[3];
sx q[3];
rz(-0.32143264) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.54520404) q[2];
sx q[2];
rz(-2.0241006) q[2];
sx q[2];
rz(-2.0580573) q[2];
rz(-0.74294535) q[3];
sx q[3];
rz(-1.5321621) q[3];
sx q[3];
rz(-1.7090181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1525986) q[0];
sx q[0];
rz(-2.3567663) q[0];
sx q[0];
rz(0.75491828) q[0];
rz(0.49916357) q[1];
sx q[1];
rz(-0.9639591) q[1];
sx q[1];
rz(-0.69854936) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38235086) q[0];
sx q[0];
rz(-2.1672591) q[0];
sx q[0];
rz(-1.9475219) q[0];
rz(-pi) q[1];
rz(-1.4261943) q[2];
sx q[2];
rz(-2.0315665) q[2];
sx q[2];
rz(-0.20587155) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9883635) q[1];
sx q[1];
rz(-0.53255845) q[1];
sx q[1];
rz(-1.531324) q[1];
rz(-pi) q[2];
rz(-1.9898765) q[3];
sx q[3];
rz(-1.5718096) q[3];
sx q[3];
rz(0.59194293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2373206) q[2];
sx q[2];
rz(-0.99232435) q[2];
sx q[2];
rz(-0.80238706) q[2];
rz(-0.55195105) q[3];
sx q[3];
rz(-2.3410083) q[3];
sx q[3];
rz(-2.5210181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17860831) q[0];
sx q[0];
rz(-1.7010138) q[0];
sx q[0];
rz(-1.1075903) q[0];
rz(1.7604609) q[1];
sx q[1];
rz(-1.7778722) q[1];
sx q[1];
rz(1.6345056) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7716146) q[0];
sx q[0];
rz(-0.42015892) q[0];
sx q[0];
rz(-0.3492979) q[0];
rz(0.92169833) q[2];
sx q[2];
rz(-0.49921152) q[2];
sx q[2];
rz(-1.8550903) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6927106) q[1];
sx q[1];
rz(-1.4254942) q[1];
sx q[1];
rz(-0.36621817) q[1];
rz(-pi) q[2];
rz(-1.4816557) q[3];
sx q[3];
rz(-2.3957068) q[3];
sx q[3];
rz(2.8583877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0845118) q[2];
sx q[2];
rz(-0.19597404) q[2];
sx q[2];
rz(-1.2723119) q[2];
rz(-1.6614527) q[3];
sx q[3];
rz(-1.1353761) q[3];
sx q[3];
rz(-0.60012668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24499527) q[0];
sx q[0];
rz(-0.56140459) q[0];
sx q[0];
rz(-1.2835314) q[0];
rz(-0.01783477) q[1];
sx q[1];
rz(-2.5051038) q[1];
sx q[1];
rz(-0.12558118) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.784362) q[0];
sx q[0];
rz(-2.1249558) q[0];
sx q[0];
rz(-0.016655075) q[0];
rz(-pi) q[1];
rz(-0.67809503) q[2];
sx q[2];
rz(-1.9481648) q[2];
sx q[2];
rz(-1.8630149) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6148147) q[1];
sx q[1];
rz(-1.1452951) q[1];
sx q[1];
rz(-2.6761901) q[1];
rz(2.2108881) q[3];
sx q[3];
rz(-0.32880201) q[3];
sx q[3];
rz(-2.9931675) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3632726) q[2];
sx q[2];
rz(-1.8724226) q[2];
sx q[2];
rz(2.3385091) q[2];
rz(-0.17656365) q[3];
sx q[3];
rz(-1.3253788) q[3];
sx q[3];
rz(-2.363502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17978996) q[0];
sx q[0];
rz(-2.640124) q[0];
sx q[0];
rz(-1.5928706) q[0];
rz(-1.8190307) q[1];
sx q[1];
rz(-1.3643199) q[1];
sx q[1];
rz(1.551569) q[1];
rz(-1.6336468) q[2];
sx q[2];
rz(-1.812915) q[2];
sx q[2];
rz(-1.6356638) q[2];
rz(0.17038067) q[3];
sx q[3];
rz(-0.34210404) q[3];
sx q[3];
rz(0.31755527) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
