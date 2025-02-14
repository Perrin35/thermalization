OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8018262) q[0];
sx q[0];
rz(-2.5751994) q[0];
sx q[0];
rz(1.9422148) q[0];
rz(2.3102923) q[1];
sx q[1];
rz(-1.6142694) q[1];
sx q[1];
rz(0.49973127) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89497772) q[0];
sx q[0];
rz(-2.3736989) q[0];
sx q[0];
rz(1.270986) q[0];
rz(1.9316767) q[2];
sx q[2];
rz(-1.3832051) q[2];
sx q[2];
rz(-2.9071992) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.079165212) q[1];
sx q[1];
rz(-1.8621596) q[1];
sx q[1];
rz(-2.3645762) q[1];
x q[2];
rz(-0.048526564) q[3];
sx q[3];
rz(-2.2884946) q[3];
sx q[3];
rz(-0.80503073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.1685593) q[2];
sx q[2];
rz(-1.7567822) q[2];
sx q[2];
rz(0.92156571) q[2];
rz(-1.5714931) q[3];
sx q[3];
rz(-0.67155963) q[3];
sx q[3];
rz(0.69935548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8916931) q[0];
sx q[0];
rz(-1.0279011) q[0];
sx q[0];
rz(-1.6722884) q[0];
rz(-1.4768614) q[1];
sx q[1];
rz(-1.3427443) q[1];
sx q[1];
rz(-2.6254168) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0196386) q[0];
sx q[0];
rz(-0.40651822) q[0];
sx q[0];
rz(1.818114) q[0];
rz(-2.8968229) q[2];
sx q[2];
rz(-1.2483856) q[2];
sx q[2];
rz(-2.4800194) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8044334) q[1];
sx q[1];
rz(-0.26072956) q[1];
sx q[1];
rz(2.8689137) q[1];
rz(2.2616181) q[3];
sx q[3];
rz(-1.954292) q[3];
sx q[3];
rz(-3.0360706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.76253647) q[2];
sx q[2];
rz(-1.9409337) q[2];
sx q[2];
rz(2.8442247) q[2];
rz(-0.51658806) q[3];
sx q[3];
rz(-1.9818431) q[3];
sx q[3];
rz(0.62492257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6487938) q[0];
sx q[0];
rz(-2.6060217) q[0];
sx q[0];
rz(-2.9597362) q[0];
rz(0.44405538) q[1];
sx q[1];
rz(-2.2477138) q[1];
sx q[1];
rz(-0.081238834) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4639909) q[0];
sx q[0];
rz(-0.39636546) q[0];
sx q[0];
rz(-0.074808077) q[0];
x q[1];
rz(2.2611558) q[2];
sx q[2];
rz(-1.3533583) q[2];
sx q[2];
rz(-1.7060929) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4887152) q[1];
sx q[1];
rz(-1.7566359) q[1];
sx q[1];
rz(1.0936827) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3879561) q[3];
sx q[3];
rz(-0.5683848) q[3];
sx q[3];
rz(-3.0912085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.2150725) q[2];
sx q[2];
rz(-0.36538616) q[2];
sx q[2];
rz(-0.60058769) q[2];
rz(1.8961204) q[3];
sx q[3];
rz(-2.5966817) q[3];
sx q[3];
rz(-3.1020402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-1.1323701) q[0];
sx q[0];
rz(-2.4339269) q[0];
sx q[0];
rz(3.1170377) q[0];
rz(-0.37697667) q[1];
sx q[1];
rz(-1.2908649) q[1];
sx q[1];
rz(-0.36184186) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1279917) q[0];
sx q[0];
rz(-0.34843081) q[0];
sx q[0];
rz(-2.6763787) q[0];
rz(-pi) q[1];
rz(0.078388647) q[2];
sx q[2];
rz(-1.9447717) q[2];
sx q[2];
rz(-1.0678991) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1805318) q[1];
sx q[1];
rz(-2.1994563) q[1];
sx q[1];
rz(-1.3500477) q[1];
rz(-pi) q[2];
rz(-0.85899701) q[3];
sx q[3];
rz(-0.59938369) q[3];
sx q[3];
rz(0.70181018) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5234066) q[2];
sx q[2];
rz(-1.8768825) q[2];
sx q[2];
rz(-0.82320172) q[2];
rz(-0.19436714) q[3];
sx q[3];
rz(-0.40209642) q[3];
sx q[3];
rz(-0.91485867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.119656) q[0];
sx q[0];
rz(-1.6031665) q[0];
sx q[0];
rz(-2.5198779) q[0];
rz(2.3960528) q[1];
sx q[1];
rz(-2.3108683) q[1];
sx q[1];
rz(0.088836975) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6357036) q[0];
sx q[0];
rz(-2.997809) q[0];
sx q[0];
rz(-0.48322387) q[0];
rz(-2.2584469) q[2];
sx q[2];
rz(-0.41637233) q[2];
sx q[2];
rz(1.0036381) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6703265) q[1];
sx q[1];
rz(-0.85976344) q[1];
sx q[1];
rz(1.1999576) q[1];
rz(-pi) q[2];
rz(-2.0114548) q[3];
sx q[3];
rz(-0.36962907) q[3];
sx q[3];
rz(-0.91418788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9125774) q[2];
sx q[2];
rz(-1.6385767) q[2];
sx q[2];
rz(-1.5360606) q[2];
rz(-2.33365) q[3];
sx q[3];
rz(-0.84351051) q[3];
sx q[3];
rz(-1.9690751) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87869969) q[0];
sx q[0];
rz(-1.4233339) q[0];
sx q[0];
rz(2.189157) q[0];
rz(1.3643422) q[1];
sx q[1];
rz(-1.9352501) q[1];
sx q[1];
rz(-2.2946024) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1901189) q[0];
sx q[0];
rz(-1.4757809) q[0];
sx q[0];
rz(0.14341469) q[0];
rz(-0.79361992) q[2];
sx q[2];
rz(-1.334903) q[2];
sx q[2];
rz(-1.867804) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.8678363) q[1];
sx q[1];
rz(-1.0844007) q[1];
sx q[1];
rz(2.0517906) q[1];
rz(-pi) q[2];
rz(-1.5562418) q[3];
sx q[3];
rz(-1.6739917) q[3];
sx q[3];
rz(-1.1905014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.18818894) q[2];
sx q[2];
rz(-2.0764515) q[2];
sx q[2];
rz(-0.38189253) q[2];
rz(-2.0415107) q[3];
sx q[3];
rz(-0.81172687) q[3];
sx q[3];
rz(-2.7303117) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5683658) q[0];
sx q[0];
rz(-1.8170284) q[0];
sx q[0];
rz(0.49760094) q[0];
rz(-0.35500232) q[1];
sx q[1];
rz(-1.6604796) q[1];
sx q[1];
rz(-2.3752046) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9123133) q[0];
sx q[0];
rz(-1.5361495) q[0];
sx q[0];
rz(-2.628792) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7117235) q[2];
sx q[2];
rz(-2.1235222) q[2];
sx q[2];
rz(-1.5220707) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.38098225) q[1];
sx q[1];
rz(-2.8462324) q[1];
sx q[1];
rz(-2.9350314) q[1];
x q[2];
rz(0.67782948) q[3];
sx q[3];
rz(-1.7683709) q[3];
sx q[3];
rz(-0.017700087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4088722) q[2];
sx q[2];
rz(-1.5343821) q[2];
sx q[2];
rz(-1.4917779) q[2];
rz(-2.0626227) q[3];
sx q[3];
rz(-1.8337199) q[3];
sx q[3];
rz(-0.61416793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.806818) q[0];
sx q[0];
rz(-2.4914927) q[0];
sx q[0];
rz(-2.7244869) q[0];
rz(3.0879703) q[1];
sx q[1];
rz(-1.5258748) q[1];
sx q[1];
rz(2.0448304) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7863207) q[0];
sx q[0];
rz(-1.4505761) q[0];
sx q[0];
rz(0.97697301) q[0];
rz(-pi) q[1];
rz(0.057889537) q[2];
sx q[2];
rz(-1.1107003) q[2];
sx q[2];
rz(1.636508) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1138823) q[1];
sx q[1];
rz(-1.405763) q[1];
sx q[1];
rz(1.5636958) q[1];
x q[2];
rz(-0.67819579) q[3];
sx q[3];
rz(-1.4525529) q[3];
sx q[3];
rz(2.6008796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2199478) q[2];
sx q[2];
rz(-1.3741477) q[2];
sx q[2];
rz(2.5448223) q[2];
rz(-0.88045949) q[3];
sx q[3];
rz(-1.7170649) q[3];
sx q[3];
rz(-0.60976353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-2.6841458) q[0];
sx q[0];
rz(-0.54231751) q[0];
sx q[0];
rz(-0.31998262) q[0];
rz(-0.34842247) q[1];
sx q[1];
rz(-0.44487822) q[1];
sx q[1];
rz(-2.6382823) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1419174) q[0];
sx q[0];
rz(-1.7192695) q[0];
sx q[0];
rz(0.093080892) q[0];
rz(-pi) q[1];
rz(1.032802) q[2];
sx q[2];
rz(-0.82394407) q[2];
sx q[2];
rz(1.8593553) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7585246) q[1];
sx q[1];
rz(-0.43979859) q[1];
sx q[1];
rz(-1.3851993) q[1];
rz(1.4473404) q[3];
sx q[3];
rz(-2.0887825) q[3];
sx q[3];
rz(-0.65604612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8866426) q[2];
sx q[2];
rz(-2.6889668) q[2];
sx q[2];
rz(2.5060999) q[2];
rz(1.5358216) q[3];
sx q[3];
rz(-1.7255892) q[3];
sx q[3];
rz(0.087285727) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6814293) q[0];
sx q[0];
rz(-2.5458113) q[0];
sx q[0];
rz(-1.3687362) q[0];
rz(2.6112556) q[1];
sx q[1];
rz(-1.6531205) q[1];
sx q[1];
rz(2.368685) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0470974) q[0];
sx q[0];
rz(-0.90029923) q[0];
sx q[0];
rz(-2.0740202) q[0];
x q[1];
rz(2.3488316) q[2];
sx q[2];
rz(-1.498641) q[2];
sx q[2];
rz(-2.4506086) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.0046282) q[1];
sx q[1];
rz(-2.0180948) q[1];
sx q[1];
rz(-0.17166673) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8498296) q[3];
sx q[3];
rz(-2.004753) q[3];
sx q[3];
rz(-1.6843554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.7698001) q[2];
sx q[2];
rz(-2.1830406) q[2];
sx q[2];
rz(-2.8690763) q[2];
rz(0.46485999) q[3];
sx q[3];
rz(-1.9404989) q[3];
sx q[3];
rz(-0.71729898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40806017) q[0];
sx q[0];
rz(-1.6049186) q[0];
sx q[0];
rz(-1.4419755) q[0];
rz(-1.5245262) q[1];
sx q[1];
rz(-2.1433612) q[1];
sx q[1];
rz(-2.8467766) q[1];
rz(2.8256399) q[2];
sx q[2];
rz(-1.8153594) q[2];
sx q[2];
rz(1.6462576) q[2];
rz(1.1285662) q[3];
sx q[3];
rz(-1.7719405) q[3];
sx q[3];
rz(2.8476261) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
