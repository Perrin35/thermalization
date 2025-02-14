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
rz(0.30000609) q[0];
sx q[0];
rz(-1.0862779) q[0];
sx q[0];
rz(2.3334184) q[0];
rz(-1.1440682) q[1];
sx q[1];
rz(-0.77163982) q[1];
sx q[1];
rz(-1.5789403) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8573541) q[0];
sx q[0];
rz(-2.1942217) q[0];
sx q[0];
rz(-0.68238544) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0879395) q[2];
sx q[2];
rz(-2.8364193) q[2];
sx q[2];
rz(2.9834753) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4889604) q[1];
sx q[1];
rz(-0.51981407) q[1];
sx q[1];
rz(1.7955154) q[1];
rz(-1.5366301) q[3];
sx q[3];
rz(-1.6164268) q[3];
sx q[3];
rz(0.51866097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.19848862) q[2];
sx q[2];
rz(-2.9297332) q[2];
sx q[2];
rz(2.1073821) q[2];
rz(-0.69283038) q[3];
sx q[3];
rz(-1.0537909) q[3];
sx q[3];
rz(2.0809765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3028054) q[0];
sx q[0];
rz(-0.028258709) q[0];
sx q[0];
rz(-2.5651108) q[0];
rz(0.020596404) q[1];
sx q[1];
rz(-2.7437904) q[1];
sx q[1];
rz(-1.0284665) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6860424) q[0];
sx q[0];
rz(-2.8032113) q[0];
sx q[0];
rz(2.4421921) q[0];
rz(1.9443545) q[2];
sx q[2];
rz(-1.8578863) q[2];
sx q[2];
rz(-0.32880983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2814863) q[1];
sx q[1];
rz(-0.53412837) q[1];
sx q[1];
rz(-0.87546443) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64739689) q[3];
sx q[3];
rz(-1.8893554) q[3];
sx q[3];
rz(2.7872374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.2241406) q[2];
sx q[2];
rz(-2.1805306) q[2];
sx q[2];
rz(-1.8769598) q[2];
rz(-2.1680016) q[3];
sx q[3];
rz(-1.6907938) q[3];
sx q[3];
rz(2.8784331) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50315404) q[0];
sx q[0];
rz(-1.3042903) q[0];
sx q[0];
rz(2.2880182) q[0];
rz(-2.0897934) q[1];
sx q[1];
rz(-0.97426668) q[1];
sx q[1];
rz(3.1265756) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6081469) q[0];
sx q[0];
rz(-0.86020532) q[0];
sx q[0];
rz(0.011562499) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.096108561) q[2];
sx q[2];
rz(-1.5495566) q[2];
sx q[2];
rz(2.6526895) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.7370809) q[1];
sx q[1];
rz(-0.85857262) q[1];
sx q[1];
rz(-1.668307) q[1];
x q[2];
rz(0.28983966) q[3];
sx q[3];
rz(-1.6320758) q[3];
sx q[3];
rz(1.0672399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.0217648) q[2];
sx q[2];
rz(-1.652635) q[2];
sx q[2];
rz(-0.28142288) q[2];
rz(-1.4540539) q[3];
sx q[3];
rz(-1.931793) q[3];
sx q[3];
rz(-0.76930261) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5948831) q[0];
sx q[0];
rz(-1.2458206) q[0];
sx q[0];
rz(2.967714) q[0];
rz(1.8100544) q[1];
sx q[1];
rz(-1.4215819) q[1];
sx q[1];
rz(1.7132267) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4177307) q[0];
sx q[0];
rz(-1.763146) q[0];
sx q[0];
rz(3.1270315) q[0];
rz(-pi) q[1];
rz(2.7816781) q[2];
sx q[2];
rz(-1.6891306) q[2];
sx q[2];
rz(-1.182076) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.54879872) q[1];
sx q[1];
rz(-2.4830635) q[1];
sx q[1];
rz(3.0786773) q[1];
rz(-2.3263116) q[3];
sx q[3];
rz(-0.51183701) q[3];
sx q[3];
rz(2.3533604) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.7017158) q[2];
sx q[2];
rz(-2.3150847) q[2];
sx q[2];
rz(-2.7244549) q[2];
rz(-2.9719628) q[3];
sx q[3];
rz(-0.093234213) q[3];
sx q[3];
rz(2.4197742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7222662) q[0];
sx q[0];
rz(-0.37501431) q[0];
sx q[0];
rz(0.30935031) q[0];
rz(2.3720692) q[1];
sx q[1];
rz(-2.0517495) q[1];
sx q[1];
rz(-1.5187029) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5681388) q[0];
sx q[0];
rz(-2.4208957) q[0];
sx q[0];
rz(1.1842404) q[0];
rz(-pi) q[1];
rz(2.4028087) q[2];
sx q[2];
rz(-1.0157983) q[2];
sx q[2];
rz(-2.8679297) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.9936693) q[1];
sx q[1];
rz(-1.6003612) q[1];
sx q[1];
rz(1.3815341) q[1];
rz(-pi) q[2];
rz(-0.098747323) q[3];
sx q[3];
rz(-1.4374466) q[3];
sx q[3];
rz(2.6790706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3994483) q[2];
sx q[2];
rz(-0.98230201) q[2];
sx q[2];
rz(2.1461416) q[2];
rz(0.71850145) q[3];
sx q[3];
rz(-1.3281497) q[3];
sx q[3];
rz(-2.3464581) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8566078) q[0];
sx q[0];
rz(-1.6866848) q[0];
sx q[0];
rz(1.6108151) q[0];
rz(-1.4467622) q[1];
sx q[1];
rz(-1.6981354) q[1];
sx q[1];
rz(-1.4379427) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7890678) q[0];
sx q[0];
rz(-1.97856) q[0];
sx q[0];
rz(-3.1192794) q[0];
rz(-1.9447359) q[2];
sx q[2];
rz(-1.0122006) q[2];
sx q[2];
rz(-2.2064759) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.6561778) q[1];
sx q[1];
rz(-0.5159157) q[1];
sx q[1];
rz(2.9713216) q[1];
rz(-1.5678649) q[3];
sx q[3];
rz(-2.3303707) q[3];
sx q[3];
rz(2.7377759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1630254) q[2];
sx q[2];
rz(-1.7820396) q[2];
sx q[2];
rz(-1.2320409) q[2];
rz(-1.9949404) q[3];
sx q[3];
rz(-1.8012643) q[3];
sx q[3];
rz(-0.31057096) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6335886) q[0];
sx q[0];
rz(-0.81729832) q[0];
sx q[0];
rz(-1.1401796) q[0];
rz(-2.2445402) q[1];
sx q[1];
rz(-1.8042118) q[1];
sx q[1];
rz(-3.1111029) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6935007) q[0];
sx q[0];
rz(-2.8909978) q[0];
sx q[0];
rz(1.084024) q[0];
rz(-1.2435799) q[2];
sx q[2];
rz(-2.1874965) q[2];
sx q[2];
rz(-1.9685352) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97059879) q[1];
sx q[1];
rz(-2.355651) q[1];
sx q[1];
rz(-1.1517609) q[1];
x q[2];
rz(-2.6875911) q[3];
sx q[3];
rz(-0.89278614) q[3];
sx q[3];
rz(2.0583722) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9474779) q[2];
sx q[2];
rz(-1.0741445) q[2];
sx q[2];
rz(-1.2987632) q[2];
rz(1.4878368) q[3];
sx q[3];
rz(-2.1066809) q[3];
sx q[3];
rz(-1.3207159) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3461935) q[0];
sx q[0];
rz(-0.29618707) q[0];
sx q[0];
rz(1.3854223) q[0];
rz(-1.2162544) q[1];
sx q[1];
rz(-1.4477891) q[1];
sx q[1];
rz(2.6522955) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84626694) q[0];
sx q[0];
rz(-1.5510484) q[0];
sx q[0];
rz(-0.67157816) q[0];
rz(-pi) q[1];
rz(1.6750181) q[2];
sx q[2];
rz(-2.8036593) q[2];
sx q[2];
rz(-0.80160415) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.7304436) q[1];
sx q[1];
rz(-2.659982) q[1];
sx q[1];
rz(-1.3661307) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.10189773) q[3];
sx q[3];
rz(-1.8188261) q[3];
sx q[3];
rz(2.0534648) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.74799246) q[2];
sx q[2];
rz(-1.9990498) q[2];
sx q[2];
rz(-2.6753814) q[2];
rz(1.5318058) q[3];
sx q[3];
rz(-1.6377056) q[3];
sx q[3];
rz(2.8209749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4894678) q[0];
sx q[0];
rz(-0.89972275) q[0];
sx q[0];
rz(-0.049064431) q[0];
rz(-0.24066726) q[1];
sx q[1];
rz(-1.0180232) q[1];
sx q[1];
rz(3.1191471) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72180333) q[0];
sx q[0];
rz(-1.5468441) q[0];
sx q[0];
rz(-2.1502058) q[0];
x q[1];
rz(0.38752611) q[2];
sx q[2];
rz(-1.8948613) q[2];
sx q[2];
rz(1.9112196) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.55521527) q[1];
sx q[1];
rz(-2.0421861) q[1];
sx q[1];
rz(-2.4280598) q[1];
x q[2];
rz(-2.4565036) q[3];
sx q[3];
rz(-1.9307319) q[3];
sx q[3];
rz(2.0956831) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-3.0040794) q[2];
sx q[2];
rz(-1.235032) q[2];
sx q[2];
rz(2.1732886) q[2];
rz(0.83640313) q[3];
sx q[3];
rz(-2.8439549) q[3];
sx q[3];
rz(-1.77232) q[3];
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
rz(-pi/2) q[0];
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
rz(-1.7264929) q[0];
sx q[0];
rz(-1.8510171) q[0];
sx q[0];
rz(2.7628164) q[0];
rz(3.1355766) q[1];
sx q[1];
rz(-0.52979398) q[1];
sx q[1];
rz(-2.5078497) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0519052) q[0];
sx q[0];
rz(-1.2335334) q[0];
sx q[0];
rz(-2.428672) q[0];
rz(-pi) q[1];
rz(2.9815282) q[2];
sx q[2];
rz(-2.1383921) q[2];
sx q[2];
rz(-2.4307323) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.56713518) q[1];
sx q[1];
rz(-2.0525041) q[1];
sx q[1];
rz(-2.5888799) q[1];
x q[2];
rz(2.2489088) q[3];
sx q[3];
rz(-2.0555946) q[3];
sx q[3];
rz(-2.3546653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.70574957) q[2];
sx q[2];
rz(-1.5219995) q[2];
sx q[2];
rz(-0.61326927) q[2];
rz(2.4663726) q[3];
sx q[3];
rz(-2.7628511) q[3];
sx q[3];
rz(0.81775445) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
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
rz(-2.2262065) q[0];
sx q[0];
rz(-1.7788667) q[0];
sx q[0];
rz(-1.9052196) q[0];
rz(2.9551103) q[1];
sx q[1];
rz(-1.7541371) q[1];
sx q[1];
rz(-1.6288155) q[1];
rz(-0.45966799) q[2];
sx q[2];
rz(-1.5528233) q[2];
sx q[2];
rz(-0.96832392) q[2];
rz(-0.9486089) q[3];
sx q[3];
rz(-1.4962248) q[3];
sx q[3];
rz(0.20635508) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
