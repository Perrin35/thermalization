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
rz(0.022194447) q[0];
sx q[0];
rz(3.7137346) q[0];
sx q[0];
rz(9.619286) q[0];
rz(2.1394849) q[1];
sx q[1];
rz(3.9245457) q[1];
sx q[1];
rz(9.4406162) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5436264) q[0];
sx q[0];
rz(-0.84613325) q[0];
sx q[0];
rz(1.0084125) q[0];
rz(2.7439666) q[2];
sx q[2];
rz(-0.83774746) q[2];
sx q[2];
rz(-2.308297) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.3091649) q[1];
sx q[1];
rz(-2.2057071) q[1];
sx q[1];
rz(-3.0135703) q[1];
rz(-pi) q[2];
rz(1.3531308) q[3];
sx q[3];
rz(-2.7870745) q[3];
sx q[3];
rz(-2.5237172) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.1281434) q[2];
sx q[2];
rz(-1.6703419) q[2];
sx q[2];
rz(-0.61689287) q[2];
rz(-1.8938176) q[3];
sx q[3];
rz(-0.82472491) q[3];
sx q[3];
rz(0.29243803) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1013252) q[0];
sx q[0];
rz(-1.1722925) q[0];
sx q[0];
rz(0.030666703) q[0];
rz(1.4847633) q[1];
sx q[1];
rz(-2.8430884) q[1];
sx q[1];
rz(-2.8224077) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46542612) q[0];
sx q[0];
rz(-2.5990708) q[0];
sx q[0];
rz(1.2529621) q[0];
rz(-pi) q[1];
rz(1.1775374) q[2];
sx q[2];
rz(-1.3280091) q[2];
sx q[2];
rz(-0.78906203) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.044925) q[1];
sx q[1];
rz(-2.6003692) q[1];
sx q[1];
rz(-0.9344395) q[1];
rz(-pi) q[2];
rz(-1.033554) q[3];
sx q[3];
rz(-1.1059009) q[3];
sx q[3];
rz(-3.1102095) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9926051) q[2];
sx q[2];
rz(-0.4036029) q[2];
sx q[2];
rz(-2.4888424) q[2];
rz(0.8695237) q[3];
sx q[3];
rz(-1.0966938) q[3];
sx q[3];
rz(2.7147527) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2000548) q[0];
sx q[0];
rz(-0.24987276) q[0];
sx q[0];
rz(3.0043434) q[0];
rz(-0.50357729) q[1];
sx q[1];
rz(-2.5027687) q[1];
sx q[1];
rz(-0.30127475) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5001016) q[0];
sx q[0];
rz(-1.8805686) q[0];
sx q[0];
rz(-0.70112438) q[0];
rz(-pi) q[1];
rz(-3.0820222) q[2];
sx q[2];
rz(-1.2822617) q[2];
sx q[2];
rz(2.3602006) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2259146) q[1];
sx q[1];
rz(-1.7748194) q[1];
sx q[1];
rz(1.7197439) q[1];
rz(-pi) q[2];
rz(-1.7295763) q[3];
sx q[3];
rz(-2.232312) q[3];
sx q[3];
rz(0.56453913) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0990024) q[2];
sx q[2];
rz(-2.4799535) q[2];
sx q[2];
rz(-2.674687) q[2];
rz(-0.0070988797) q[3];
sx q[3];
rz(-3.0341798) q[3];
sx q[3];
rz(-0.98889178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.057509363) q[0];
sx q[0];
rz(-3.1355661) q[0];
sx q[0];
rz(0.65473336) q[0];
rz(1.5714802) q[1];
sx q[1];
rz(-1.4288158) q[1];
sx q[1];
rz(2.7007801) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8424555) q[0];
sx q[0];
rz(-1.5464325) q[0];
sx q[0];
rz(0.024473622) q[0];
rz(-pi) q[1];
rz(-0.47870584) q[2];
sx q[2];
rz(-1.5619593) q[2];
sx q[2];
rz(-2.1610726) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.310439) q[1];
sx q[1];
rz(-0.37292994) q[1];
sx q[1];
rz(2.4937477) q[1];
rz(-pi) q[2];
rz(2.6648947) q[3];
sx q[3];
rz(-1.293382) q[3];
sx q[3];
rz(-1.6200844) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.64678994) q[2];
sx q[2];
rz(-1.9014634) q[2];
sx q[2];
rz(2.6498762) q[2];
rz(-1.803319) q[3];
sx q[3];
rz(-1.0889784) q[3];
sx q[3];
rz(1.4664388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76984513) q[0];
sx q[0];
rz(-2.0553698) q[0];
sx q[0];
rz(2.0628498) q[0];
rz(-0.046401333) q[1];
sx q[1];
rz(-1.6175783) q[1];
sx q[1];
rz(-2.0414415) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.28957) q[0];
sx q[0];
rz(-1.4542113) q[0];
sx q[0];
rz(-1.5800716) q[0];
x q[1];
rz(-2.4920667) q[2];
sx q[2];
rz(-0.99156556) q[2];
sx q[2];
rz(1.9686101) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.8535159) q[1];
sx q[1];
rz(-1.3259246) q[1];
sx q[1];
rz(2.4254786) q[1];
x q[2];
rz(-0.38918503) q[3];
sx q[3];
rz(-0.2901623) q[3];
sx q[3];
rz(-0.25012514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.84755737) q[2];
sx q[2];
rz(-0.39176771) q[2];
sx q[2];
rz(-0.34226391) q[2];
rz(1.3058454) q[3];
sx q[3];
rz(-1.9394453) q[3];
sx q[3];
rz(2.0391298) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24868988) q[0];
sx q[0];
rz(-1.9547639) q[0];
sx q[0];
rz(-2.7500395) q[0];
rz(2.4941817) q[1];
sx q[1];
rz(-2.7029111) q[1];
sx q[1];
rz(2.2364) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12654141) q[0];
sx q[0];
rz(-0.76989664) q[0];
sx q[0];
rz(-2.4036744) q[0];
rz(-0.91692704) q[2];
sx q[2];
rz(-0.77007896) q[2];
sx q[2];
rz(-0.21756141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.74108142) q[1];
sx q[1];
rz(-1.834474) q[1];
sx q[1];
rz(2.5225894) q[1];
x q[2];
rz(1.5213826) q[3];
sx q[3];
rz(-1.2826398) q[3];
sx q[3];
rz(-2.2505983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.68760005) q[2];
sx q[2];
rz(-0.31012055) q[2];
sx q[2];
rz(-1.0529168) q[2];
rz(0.29371253) q[3];
sx q[3];
rz(-2.3070344) q[3];
sx q[3];
rz(-0.53148758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(2.3783136) q[0];
sx q[0];
rz(-1.2792307) q[0];
sx q[0];
rz(3.1148425) q[0];
rz(1.9552975) q[1];
sx q[1];
rz(-0.18272884) q[1];
sx q[1];
rz(-2.9881086) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59236103) q[0];
sx q[0];
rz(-0.57825297) q[0];
sx q[0];
rz(1.0628257) q[0];
rz(1.5181958) q[2];
sx q[2];
rz(-2.1204877) q[2];
sx q[2];
rz(-2.652183) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.19928923) q[1];
sx q[1];
rz(-2.2484712) q[1];
sx q[1];
rz(0.79395644) q[1];
rz(-pi) q[2];
rz(-1.1304705) q[3];
sx q[3];
rz(-3.0743361) q[3];
sx q[3];
rz(-2.7471144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.89580506) q[2];
sx q[2];
rz(-0.20496002) q[2];
sx q[2];
rz(1.8765571) q[2];
rz(-2.9686019) q[3];
sx q[3];
rz(-1.750662) q[3];
sx q[3];
rz(2.0632108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.10385253) q[0];
sx q[0];
rz(-2.9321892) q[0];
sx q[0];
rz(-3.0331392) q[0];
rz(1.3718038) q[1];
sx q[1];
rz(-1.3550974) q[1];
sx q[1];
rz(3.1350737) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9193503) q[0];
sx q[0];
rz(-2.0323638) q[0];
sx q[0];
rz(2.190175) q[0];
rz(-1.4638607) q[2];
sx q[2];
rz(-1.7902014) q[2];
sx q[2];
rz(2.1768513) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7951736) q[1];
sx q[1];
rz(-1.0246207) q[1];
sx q[1];
rz(2.4743028) q[1];
x q[2];
rz(2.3586999) q[3];
sx q[3];
rz(-1.8777579) q[3];
sx q[3];
rz(1.3166733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45695496) q[2];
sx q[2];
rz(-2.4943887) q[2];
sx q[2];
rz(-1.4400488) q[2];
rz(0.35259926) q[3];
sx q[3];
rz(-2.2779901) q[3];
sx q[3];
rz(2.6889804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.088293485) q[0];
sx q[0];
rz(-2.0501917) q[0];
sx q[0];
rz(-1.1540867) q[0];
rz(-1.6905009) q[1];
sx q[1];
rz(-0.68203753) q[1];
sx q[1];
rz(-2.3114204) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.80487554) q[0];
sx q[0];
rz(-1.5600388) q[0];
sx q[0];
rz(-1.5675926) q[0];
x q[1];
rz(1.4373892) q[2];
sx q[2];
rz(-2.0683943) q[2];
sx q[2];
rz(-1.069151) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8389924) q[1];
sx q[1];
rz(-2.5328089) q[1];
sx q[1];
rz(2.3997612) q[1];
x q[2];
rz(-1.3187879) q[3];
sx q[3];
rz(-2.4652836) q[3];
sx q[3];
rz(2.0700434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.95149583) q[2];
sx q[2];
rz(-2.2647965) q[2];
sx q[2];
rz(-2.4693176) q[2];
rz(-0.43664524) q[3];
sx q[3];
rz(-1.838622) q[3];
sx q[3];
rz(2.8521027) q[3];
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
x q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3960246) q[0];
sx q[0];
rz(-2.4085299) q[0];
sx q[0];
rz(2.720604) q[0];
rz(1.1517395) q[1];
sx q[1];
rz(-0.20944171) q[1];
sx q[1];
rz(-0.71796012) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1678667) q[0];
sx q[0];
rz(-1.5930678) q[0];
sx q[0];
rz(1.2931965) q[0];
rz(-pi) q[1];
x q[1];
rz(0.46109445) q[2];
sx q[2];
rz(-2.3873219) q[2];
sx q[2];
rz(-0.2939156) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.1970994) q[1];
sx q[1];
rz(-1.2596895) q[1];
sx q[1];
rz(-0.13469805) q[1];
rz(-pi) q[2];
rz(1.1208543) q[3];
sx q[3];
rz(-1.8186343) q[3];
sx q[3];
rz(-0.42603394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1154321) q[2];
sx q[2];
rz(-2.5371976) q[2];
sx q[2];
rz(1.7660512) q[2];
rz(0.17624217) q[3];
sx q[3];
rz(-2.3459489) q[3];
sx q[3];
rz(-3.0680883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7701223) q[0];
sx q[0];
rz(-1.4513411) q[0];
sx q[0];
rz(2.0392192) q[0];
rz(-0.86028987) q[1];
sx q[1];
rz(-1.6382891) q[1];
sx q[1];
rz(2.0232497) q[1];
rz(-2.3762351) q[2];
sx q[2];
rz(-1.5258963) q[2];
sx q[2];
rz(-0.85276251) q[2];
rz(1.1841487) q[3];
sx q[3];
rz(-1.2332698) q[3];
sx q[3];
rz(-1.5982066) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
