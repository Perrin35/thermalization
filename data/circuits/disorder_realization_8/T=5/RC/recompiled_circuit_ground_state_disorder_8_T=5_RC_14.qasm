OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.3437164) q[0];
sx q[0];
rz(-2.1581082) q[0];
sx q[0];
rz(-2.9498192) q[0];
rz(-2.6311488) q[1];
sx q[1];
rz(-1.3671083) q[1];
sx q[1];
rz(1.3956611) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5001427) q[0];
sx q[0];
rz(-0.035273835) q[0];
sx q[0];
rz(-2.7462237) q[0];
rz(0.67740719) q[2];
sx q[2];
rz(-2.8500994) q[2];
sx q[2];
rz(3.0335333) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8366791) q[1];
sx q[1];
rz(-2.3354451) q[1];
sx q[1];
rz(-2.8668001) q[1];
rz(-pi) q[2];
rz(1.9590501) q[3];
sx q[3];
rz(-2.7388529) q[3];
sx q[3];
rz(-2.5975063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6338966) q[2];
sx q[2];
rz(-0.32749978) q[2];
sx q[2];
rz(-2.2155649) q[2];
rz(-1.6294468) q[3];
sx q[3];
rz(-0.67073268) q[3];
sx q[3];
rz(1.9984455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.73285) q[0];
sx q[0];
rz(-2.4091305) q[0];
sx q[0];
rz(2.9090885) q[0];
rz(-2.0022424) q[1];
sx q[1];
rz(-1.5683697) q[1];
sx q[1];
rz(0.92612902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8208198) q[0];
sx q[0];
rz(-2.8348587) q[0];
sx q[0];
rz(-2.308918) q[0];
rz(-pi) q[1];
x q[1];
rz(2.861169) q[2];
sx q[2];
rz(-0.33992514) q[2];
sx q[2];
rz(-0.47001335) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2567443) q[1];
sx q[1];
rz(-1.7769281) q[1];
sx q[1];
rz(0.54289674) q[1];
rz(2.0235463) q[3];
sx q[3];
rz(-2.0914075) q[3];
sx q[3];
rz(-1.3829766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6050379) q[2];
sx q[2];
rz(-1.6190517) q[2];
sx q[2];
rz(-0.79263318) q[2];
rz(-2.7301181) q[3];
sx q[3];
rz(-2.1992407) q[3];
sx q[3];
rz(-2.3365432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3181535) q[0];
sx q[0];
rz(-1.0018438) q[0];
sx q[0];
rz(0.26082984) q[0];
rz(-0.28314319) q[1];
sx q[1];
rz(-1.4500376) q[1];
sx q[1];
rz(-2.0859437) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2973141) q[0];
sx q[0];
rz(-1.3998654) q[0];
sx q[0];
rz(-2.2173082) q[0];
rz(-1.5670384) q[2];
sx q[2];
rz(-0.93766391) q[2];
sx q[2];
rz(2.1830677) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.187589) q[1];
sx q[1];
rz(-2.2936645) q[1];
sx q[1];
rz(1.4299117) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6840213) q[3];
sx q[3];
rz(-0.49834305) q[3];
sx q[3];
rz(-1.6087974) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.61567125) q[2];
sx q[2];
rz(-2.0522223) q[2];
sx q[2];
rz(1.2582568) q[2];
rz(-1.0387756) q[3];
sx q[3];
rz(-2.153331) q[3];
sx q[3];
rz(1.725215) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1309758) q[0];
sx q[0];
rz(-1.5357786) q[0];
sx q[0];
rz(-0.11949874) q[0];
rz(2.9833228) q[1];
sx q[1];
rz(-2.6893078) q[1];
sx q[1];
rz(-1.4403042) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7410204) q[0];
sx q[0];
rz(-1.9306679) q[0];
sx q[0];
rz(-2.3163356) q[0];
x q[1];
rz(-1.794852) q[2];
sx q[2];
rz(-0.34287004) q[2];
sx q[2];
rz(2.1488239) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.54771891) q[1];
sx q[1];
rz(-0.55127326) q[1];
sx q[1];
rz(2.6216595) q[1];
rz(-pi) q[2];
rz(-2.1705008) q[3];
sx q[3];
rz(-2.5673742) q[3];
sx q[3];
rz(-0.80050877) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.4629472) q[2];
sx q[2];
rz(-3.0061649) q[2];
sx q[2];
rz(2.7079008) q[2];
rz(2.4667451) q[3];
sx q[3];
rz(-1.0899455) q[3];
sx q[3];
rz(1.9851782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(2.7334412) q[0];
sx q[0];
rz(-2.0858522) q[0];
sx q[0];
rz(2.5422886) q[0];
rz(0.76404461) q[1];
sx q[1];
rz(-2.1115477) q[1];
sx q[1];
rz(0.61242551) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9117994) q[0];
sx q[0];
rz(-0.62055885) q[0];
sx q[0];
rz(1.8259274) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5715309) q[2];
sx q[2];
rz(-2.2876566) q[2];
sx q[2];
rz(1.5518722) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.99701559) q[1];
sx q[1];
rz(-1.738228) q[1];
sx q[1];
rz(2.5232878) q[1];
x q[2];
rz(-1.0507842) q[3];
sx q[3];
rz(-2.0377167) q[3];
sx q[3];
rz(2.0869521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4446438) q[2];
sx q[2];
rz(-0.79978839) q[2];
sx q[2];
rz(-0.4701699) q[2];
rz(2.7759975) q[3];
sx q[3];
rz(-1.9496893) q[3];
sx q[3];
rz(2.5308334) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81569165) q[0];
sx q[0];
rz(-2.3453562) q[0];
sx q[0];
rz(-0.66194397) q[0];
rz(-3.0603307) q[1];
sx q[1];
rz(-2.6632402) q[1];
sx q[1];
rz(-3.001396) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7707342) q[0];
sx q[0];
rz(-2.6261289) q[0];
sx q[0];
rz(-0.45592842) q[0];
x q[1];
rz(-2.0070932) q[2];
sx q[2];
rz(-0.66777523) q[2];
sx q[2];
rz(2.2420924) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6684456) q[1];
sx q[1];
rz(-2.6700017) q[1];
sx q[1];
rz(1.5890597) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.46886985) q[3];
sx q[3];
rz(-0.38772407) q[3];
sx q[3];
rz(-2.455251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7774272) q[2];
sx q[2];
rz(-1.825793) q[2];
sx q[2];
rz(-0.3248471) q[2];
rz(-0.15803629) q[3];
sx q[3];
rz(-0.41301781) q[3];
sx q[3];
rz(-1.0154999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3196816) q[0];
sx q[0];
rz(-0.076229036) q[0];
sx q[0];
rz(-0.88777375) q[0];
rz(0.38331389) q[1];
sx q[1];
rz(-0.8822459) q[1];
sx q[1];
rz(-0.15957889) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68113806) q[0];
sx q[0];
rz(-2.7260906) q[0];
sx q[0];
rz(0.88232188) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80990661) q[2];
sx q[2];
rz(-0.72945539) q[2];
sx q[2];
rz(-0.73939656) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.23160203) q[1];
sx q[1];
rz(-1.195915) q[1];
sx q[1];
rz(-1.7085307) q[1];
x q[2];
rz(1.4042312) q[3];
sx q[3];
rz(-1.8565531) q[3];
sx q[3];
rz(0.33657246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5923656) q[2];
sx q[2];
rz(-1.1575907) q[2];
sx q[2];
rz(-1.7249736) q[2];
rz(1.2688515) q[3];
sx q[3];
rz(-0.78059355) q[3];
sx q[3];
rz(-2.1835073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6855327) q[0];
sx q[0];
rz(-2.0262418) q[0];
sx q[0];
rz(0.90150315) q[0];
rz(0.69425663) q[1];
sx q[1];
rz(-1.6777104) q[1];
sx q[1];
rz(-1.136397) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0944871) q[0];
sx q[0];
rz(-1.4102954) q[0];
sx q[0];
rz(-3.1094527) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.55113422) q[2];
sx q[2];
rz(-1.8922085) q[2];
sx q[2];
rz(1.2642191) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8759114) q[1];
sx q[1];
rz(-2.2655575) q[1];
sx q[1];
rz(3.0819703) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0497007) q[3];
sx q[3];
rz(-0.9113833) q[3];
sx q[3];
rz(-1.5707317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14094341) q[2];
sx q[2];
rz(-1.8724172) q[2];
sx q[2];
rz(2.3285274) q[2];
rz(-0.55417577) q[3];
sx q[3];
rz(-1.7289836) q[3];
sx q[3];
rz(-0.36177844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
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
rz(-2.5799705) q[0];
sx q[0];
rz(-2.0273835) q[0];
sx q[0];
rz(0.17350523) q[0];
rz(-0.42501998) q[1];
sx q[1];
rz(-1.605502) q[1];
sx q[1];
rz(-0.82957155) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7873142) q[0];
sx q[0];
rz(-1.4256546) q[0];
sx q[0];
rz(-2.4046242) q[0];
rz(-3.095171) q[2];
sx q[2];
rz(-1.5299462) q[2];
sx q[2];
rz(-2.61039) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8977938) q[1];
sx q[1];
rz(-1.3659371) q[1];
sx q[1];
rz(0.61989354) q[1];
rz(-pi) q[2];
rz(2.9584269) q[3];
sx q[3];
rz(-2.7999827) q[3];
sx q[3];
rz(-2.6721725) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.35036626) q[2];
sx q[2];
rz(-1.841505) q[2];
sx q[2];
rz(-1.680797) q[2];
rz(-2.1469877) q[3];
sx q[3];
rz(-2.1456783) q[3];
sx q[3];
rz(0.88609707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9643672) q[0];
sx q[0];
rz(-2.565964) q[0];
sx q[0];
rz(3.0395569) q[0];
rz(-0.61965865) q[1];
sx q[1];
rz(-1.4980059) q[1];
sx q[1];
rz(0.59725753) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6031979) q[0];
sx q[0];
rz(-1.4369597) q[0];
sx q[0];
rz(-0.60026818) q[0];
rz(2.7623457) q[2];
sx q[2];
rz(-0.37988099) q[2];
sx q[2];
rz(-0.82373649) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.958739) q[1];
sx q[1];
rz(-2.5748061) q[1];
sx q[1];
rz(-1.6394407) q[1];
x q[2];
rz(-3.0773874) q[3];
sx q[3];
rz(-2.6693025) q[3];
sx q[3];
rz(2.9407168) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9336046) q[2];
sx q[2];
rz(-2.7808069) q[2];
sx q[2];
rz(-0.92850816) q[2];
rz(-1.9814631) q[3];
sx q[3];
rz(-1.7103633) q[3];
sx q[3];
rz(2.3742356) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3930897) q[0];
sx q[0];
rz(-0.7001644) q[0];
sx q[0];
rz(-2.9099332) q[0];
rz(1.1275935) q[1];
sx q[1];
rz(-2.6078106) q[1];
sx q[1];
rz(-0.003905205) q[1];
rz(-1.4814346) q[2];
sx q[2];
rz(-2.4632774) q[2];
sx q[2];
rz(2.6954271) q[2];
rz(-0.87455187) q[3];
sx q[3];
rz(-1.3037494) q[3];
sx q[3];
rz(-1.3791374) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
