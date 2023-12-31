OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.4424326) q[0];
sx q[0];
rz(-1.3843098) q[0];
sx q[0];
rz(-1.260489) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(-2.2178712) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0246575) q[0];
sx q[0];
rz(-0.57514656) q[0];
sx q[0];
rz(0.92098178) q[0];
rz(-1.853763) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(-1.9622918) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.10780653) q[1];
sx q[1];
rz(-2.4596655) q[1];
sx q[1];
rz(2.4297907) q[1];
x q[2];
rz(-0.63763036) q[3];
sx q[3];
rz(-0.59083592) q[3];
sx q[3];
rz(-1.5096111) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(0.16201924) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-2.155442) q[3];
sx q[3];
rz(-0.71301618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0682003) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.1967999) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(1.4555567) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6360977) q[0];
sx q[0];
rz(-2.539145) q[0];
sx q[0];
rz(1.8683744) q[0];
x q[1];
rz(1.4886841) q[2];
sx q[2];
rz(-1.9442691) q[2];
sx q[2];
rz(2.3478846) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1728954) q[1];
sx q[1];
rz(-2.6186133) q[1];
sx q[1];
rz(-1.5978659) q[1];
rz(-pi) q[2];
rz(-0.41166572) q[3];
sx q[3];
rz(-1.5951612) q[3];
sx q[3];
rz(2.0522232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.42852795) q[2];
sx q[2];
rz(-1.6899127) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(-2.9591566) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3882554) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(0.80048168) q[0];
rz(-3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(-1.172539) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3282933) q[0];
sx q[0];
rz(-1.0070224) q[0];
sx q[0];
rz(-1.203712) q[0];
x q[1];
rz(-1.8981947) q[2];
sx q[2];
rz(-0.84257579) q[2];
sx q[2];
rz(0.74795216) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.4634358) q[1];
sx q[1];
rz(-2.0883745) q[1];
sx q[1];
rz(-2.8606158) q[1];
x q[2];
rz(1.7964732) q[3];
sx q[3];
rz(-0.25697069) q[3];
sx q[3];
rz(-0.32303177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(-0.91119901) q[2];
rz(2.1905812) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-2.2495911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38055414) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(0.12810853) q[0];
rz(0.076106636) q[1];
sx q[1];
rz(-1.2144621) q[1];
sx q[1];
rz(-2.6180843) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2142221) q[0];
sx q[0];
rz(-2.4282051) q[0];
sx q[0];
rz(-0.58332304) q[0];
rz(1.4857616) q[2];
sx q[2];
rz(-1.2418613) q[2];
sx q[2];
rz(-2.0563682) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7586786) q[1];
sx q[1];
rz(-2.3944693) q[1];
sx q[1];
rz(1.1927356) q[1];
rz(-pi) q[2];
rz(-2.162096) q[3];
sx q[3];
rz(-1.4569836) q[3];
sx q[3];
rz(-2.1554961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(-0.564044) q[2];
rz(-0.28856746) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(-2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(0.87961698) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(0.99194828) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6447727) q[0];
sx q[0];
rz(-2.9368375) q[0];
sx q[0];
rz(0.55069189) q[0];
rz(-pi) q[1];
x q[1];
rz(0.0089346272) q[2];
sx q[2];
rz(-1.8736471) q[2];
sx q[2];
rz(-0.47700275) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79210287) q[1];
sx q[1];
rz(-1.8794685) q[1];
sx q[1];
rz(-1.6130788) q[1];
rz(-pi) q[2];
rz(3.0117412) q[3];
sx q[3];
rz(-0.81749812) q[3];
sx q[3];
rz(-2.0714456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.1296967) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-0.27080718) q[2];
rz(2.9233542) q[3];
sx q[3];
rz(-1.821358) q[3];
sx q[3];
rz(-0.22578421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72702423) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.7927992) q[0];
rz(-0.38189608) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(-1.4250925) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2857367) q[0];
sx q[0];
rz(-0.55186134) q[0];
sx q[0];
rz(-0.43666552) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9082597) q[2];
sx q[2];
rz(-2.87185) q[2];
sx q[2];
rz(-2.595682) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9565935) q[1];
sx q[1];
rz(-1.6374267) q[1];
sx q[1];
rz(2.0208298) q[1];
rz(-2.1720042) q[3];
sx q[3];
rz(-1.117327) q[3];
sx q[3];
rz(-2.3525402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0075334) q[2];
sx q[2];
rz(-0.39742658) q[2];
sx q[2];
rz(2.5777204) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.6231977) q[3];
sx q[3];
rz(-2.738651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5086223) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(1.58889) q[1];
sx q[1];
rz(-1.8808552) q[1];
sx q[1];
rz(-2.3197876) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2215295) q[0];
sx q[0];
rz(-0.95525817) q[0];
sx q[0];
rz(0.91381844) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2981811) q[2];
sx q[2];
rz(-1.1864098) q[2];
sx q[2];
rz(-0.20993983) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2591178) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(0.74101733) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.58737289) q[3];
sx q[3];
rz(-1.2578739) q[3];
sx q[3];
rz(-1.3164933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8043148) q[2];
sx q[2];
rz(-0.75411212) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(-0.129536) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.725175) q[0];
sx q[0];
rz(-3.1224407) q[0];
sx q[0];
rz(-2.3186671) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(-1.3051422) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96852036) q[0];
sx q[0];
rz(-2.1497796) q[0];
sx q[0];
rz(-1.9645683) q[0];
rz(-2.4370286) q[2];
sx q[2];
rz(-1.2837871) q[2];
sx q[2];
rz(-1.2003843) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.23755632) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(-1.1035641) q[1];
x q[2];
rz(0.023530258) q[3];
sx q[3];
rz(-1.9100034) q[3];
sx q[3];
rz(-0.48285218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.3858162) q[2];
sx q[2];
rz(1.4902327) q[2];
rz(1.0772609) q[3];
sx q[3];
rz(-0.96499413) q[3];
sx q[3];
rz(-3.0100477) q[3];
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
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-2.819678) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-0.70294356) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84530572) q[0];
sx q[0];
rz(-1.1830813) q[0];
sx q[0];
rz(-1.182343) q[0];
rz(-0.039673294) q[2];
sx q[2];
rz(-1.0618321) q[2];
sx q[2];
rz(2.2850125) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.798588) q[1];
sx q[1];
rz(-1.4061905) q[1];
sx q[1];
rz(-2.1144457) q[1];
rz(-2.0270991) q[3];
sx q[3];
rz(-0.66702402) q[3];
sx q[3];
rz(1.7175355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.90074173) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.3396324) q[2];
rz(-0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(-1.3302749) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16383485) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(1.9158069) q[0];
rz(-2.2380791) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(-2.6729565) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0271284) q[0];
sx q[0];
rz(-2.7002324) q[0];
sx q[0];
rz(-2.6370185) q[0];
rz(-pi) q[1];
rz(0.87601985) q[2];
sx q[2];
rz(-2.6555736) q[2];
sx q[2];
rz(1.8234058) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.84506449) q[1];
sx q[1];
rz(-1.5025286) q[1];
sx q[1];
rz(-1.3045842) q[1];
x q[2];
rz(-2.3234899) q[3];
sx q[3];
rz(-0.59554282) q[3];
sx q[3];
rz(-2.9593352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6580711) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(1.6121929) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4951915) q[0];
sx q[0];
rz(-1.9468745) q[0];
sx q[0];
rz(-0.68328802) q[0];
rz(-2.519683) q[1];
sx q[1];
rz(-1.4629296) q[1];
sx q[1];
rz(-0.32348979) q[1];
rz(-0.8016349) q[2];
sx q[2];
rz(-2.2710706) q[2];
sx q[2];
rz(2.0588277) q[2];
rz(-0.078483742) q[3];
sx q[3];
rz(-0.92072903) q[3];
sx q[3];
rz(-2.846684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
