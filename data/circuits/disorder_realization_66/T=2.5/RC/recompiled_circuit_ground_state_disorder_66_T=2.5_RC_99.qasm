OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8074789) q[0];
sx q[0];
rz(2.4456094) q[0];
sx q[0];
rz(13.453475) q[0];
rz(2.4298985) q[1];
sx q[1];
rz(-2.0790172) q[1];
sx q[1];
rz(0.33023155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45107156) q[0];
sx q[0];
rz(-0.94342782) q[0];
sx q[0];
rz(-0.55212195) q[0];
rz(-pi) q[1];
rz(-2.5862972) q[2];
sx q[2];
rz(-1.2568297) q[2];
sx q[2];
rz(-1.2295251) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.2621999) q[1];
sx q[1];
rz(-2.6385251) q[1];
sx q[1];
rz(-2.9285953) q[1];
rz(-pi) q[2];
rz(-3.0631271) q[3];
sx q[3];
rz(-2.0002504) q[3];
sx q[3];
rz(2.0450908) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.3868788) q[2];
sx q[2];
rz(-0.42749307) q[2];
sx q[2];
rz(-1.0270366) q[2];
rz(2.2586281) q[3];
sx q[3];
rz(-1.797978) q[3];
sx q[3];
rz(-0.55356717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3565732) q[0];
sx q[0];
rz(-0.078332575) q[0];
sx q[0];
rz(-2.8739492) q[0];
rz(1.74125) q[1];
sx q[1];
rz(-0.61336556) q[1];
sx q[1];
rz(-2.5852481) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2927383) q[0];
sx q[0];
rz(-1.8195099) q[0];
sx q[0];
rz(2.0746887) q[0];
rz(-pi) q[1];
x q[1];
rz(0.39926932) q[2];
sx q[2];
rz(-2.396691) q[2];
sx q[2];
rz(2.5585554) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.90168629) q[1];
sx q[1];
rz(-1.0178029) q[1];
sx q[1];
rz(0.4497952) q[1];
rz(-pi) q[2];
rz(2.2351867) q[3];
sx q[3];
rz(-1.0944546) q[3];
sx q[3];
rz(1.2211459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.9050425) q[2];
sx q[2];
rz(-1.4507797) q[2];
sx q[2];
rz(-1.8897918) q[2];
rz(1.8768138) q[3];
sx q[3];
rz(-0.71139657) q[3];
sx q[3];
rz(-1.0714162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4816137) q[0];
sx q[0];
rz(-0.55584207) q[0];
sx q[0];
rz(-2.3023093) q[0];
rz(-1.6370157) q[1];
sx q[1];
rz(-1.7103651) q[1];
sx q[1];
rz(-1.7617216) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8913075) q[0];
sx q[0];
rz(-1.4472741) q[0];
sx q[0];
rz(0.082300622) q[0];
x q[1];
rz(-1.9154869) q[2];
sx q[2];
rz(-0.57248964) q[2];
sx q[2];
rz(-0.16333157) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.4288001) q[1];
sx q[1];
rz(-1.2512491) q[1];
sx q[1];
rz(-1.0733114) q[1];
rz(1.0187418) q[3];
sx q[3];
rz(-2.0716487) q[3];
sx q[3];
rz(0.93397442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.442349) q[2];
sx q[2];
rz(-1.2171429) q[2];
sx q[2];
rz(0.71845636) q[2];
rz(2.5213304) q[3];
sx q[3];
rz(-2.2060427) q[3];
sx q[3];
rz(-0.81234318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.45850596) q[0];
sx q[0];
rz(-1.9106671) q[0];
sx q[0];
rz(1.7339535) q[0];
rz(-0.55157026) q[1];
sx q[1];
rz(-3.0307814) q[1];
sx q[1];
rz(1.3161906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16770076) q[0];
sx q[0];
rz(-0.97588723) q[0];
sx q[0];
rz(0.3428726) q[0];
rz(3.1291601) q[2];
sx q[2];
rz(-1.2913449) q[2];
sx q[2];
rz(-2.3258727) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.11034695) q[1];
sx q[1];
rz(-1.7025885) q[1];
sx q[1];
rz(0.40459569) q[1];
rz(-pi) q[2];
rz(-0.82877036) q[3];
sx q[3];
rz(-2.2743524) q[3];
sx q[3];
rz(0.32876536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.59336415) q[2];
sx q[2];
rz(-2.7061988) q[2];
sx q[2];
rz(-1.0080053) q[2];
rz(-2.8625782) q[3];
sx q[3];
rz(-0.81441003) q[3];
sx q[3];
rz(-0.45849714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27610436) q[0];
sx q[0];
rz(-1.8115598) q[0];
sx q[0];
rz(-2.8651067) q[0];
rz(-0.28929389) q[1];
sx q[1];
rz(-1.1951059) q[1];
sx q[1];
rz(-0.2624661) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2967178) q[0];
sx q[0];
rz(-0.27063677) q[0];
sx q[0];
rz(2.4225745) q[0];
x q[1];
rz(-1.5868411) q[2];
sx q[2];
rz(-2.5501921) q[2];
sx q[2];
rz(0.41447869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.9346859) q[1];
sx q[1];
rz(-1.0240203) q[1];
sx q[1];
rz(-1.9711167) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5632202) q[3];
sx q[3];
rz(-1.4422073) q[3];
sx q[3];
rz(1.2728387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3929954) q[2];
sx q[2];
rz(-2.5454919) q[2];
sx q[2];
rz(1.7945012) q[2];
rz(3.0387481) q[3];
sx q[3];
rz(-1.9072073) q[3];
sx q[3];
rz(1.6245406) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1384001) q[0];
sx q[0];
rz(-1.5035368) q[0];
sx q[0];
rz(-0.060977161) q[0];
rz(-1.2334476) q[1];
sx q[1];
rz(-1.3950709) q[1];
sx q[1];
rz(1.6159509) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71919854) q[0];
sx q[0];
rz(-0.53421796) q[0];
sx q[0];
rz(-3.0052276) q[0];
rz(-pi) q[1];
x q[1];
rz(0.032268957) q[2];
sx q[2];
rz(-2.319558) q[2];
sx q[2];
rz(-0.2781333) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.7394289) q[1];
sx q[1];
rz(-1.3203748) q[1];
sx q[1];
rz(3.0274505) q[1];
rz(-0.16651972) q[3];
sx q[3];
rz(-2.2158479) q[3];
sx q[3];
rz(1.695961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6700651) q[2];
sx q[2];
rz(-1.7933041) q[2];
sx q[2];
rz(1.3304375) q[2];
rz(-2.0165675) q[3];
sx q[3];
rz(-2.2671813) q[3];
sx q[3];
rz(3.103638) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1838945) q[0];
sx q[0];
rz(-1.4426458) q[0];
sx q[0];
rz(-2.884602) q[0];
rz(-2.7067302) q[1];
sx q[1];
rz(-0.75006524) q[1];
sx q[1];
rz(-0.90352568) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6174451) q[0];
sx q[0];
rz(-3.0492231) q[0];
sx q[0];
rz(1.6156107) q[0];
rz(-2.9952094) q[2];
sx q[2];
rz(-0.75136649) q[2];
sx q[2];
rz(-2.9379972) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3147887) q[1];
sx q[1];
rz(-1.5233524) q[1];
sx q[1];
rz(-0.06158365) q[1];
rz(-pi) q[2];
x q[2];
rz(0.89374505) q[3];
sx q[3];
rz(-0.66679685) q[3];
sx q[3];
rz(-2.9725985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.32288185) q[2];
sx q[2];
rz(-1.633753) q[2];
sx q[2];
rz(1.4212849) q[2];
rz(-2.4317702) q[3];
sx q[3];
rz(-1.1387419) q[3];
sx q[3];
rz(2.6074954) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1555971) q[0];
sx q[0];
rz(-0.43102145) q[0];
sx q[0];
rz(2.0565597) q[0];
rz(0.64708465) q[1];
sx q[1];
rz(-1.1754464) q[1];
sx q[1];
rz(3.0771902) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8540809) q[0];
sx q[0];
rz(-2.4123544) q[0];
sx q[0];
rz(-0.50042787) q[0];
rz(-2.4738381) q[2];
sx q[2];
rz(-2.1760094) q[2];
sx q[2];
rz(2.3507694) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6820087) q[1];
sx q[1];
rz(-2.0279851) q[1];
sx q[1];
rz(0.80322916) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2256718) q[3];
sx q[3];
rz(-0.43817156) q[3];
sx q[3];
rz(-2.8000591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.881968) q[2];
sx q[2];
rz(-2.7590064) q[2];
sx q[2];
rz(3.0982223) q[2];
rz(0.67443332) q[3];
sx q[3];
rz(-1.3573656) q[3];
sx q[3];
rz(1.9969214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6944273) q[0];
sx q[0];
rz(-0.58037102) q[0];
sx q[0];
rz(-2.8149783) q[0];
rz(1.0940374) q[1];
sx q[1];
rz(-2.2973165) q[1];
sx q[1];
rz(-0.48386595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.992386) q[0];
sx q[0];
rz(-1.0629553) q[0];
sx q[0];
rz(0.81224982) q[0];
rz(0.47205527) q[2];
sx q[2];
rz(-2.6723571) q[2];
sx q[2];
rz(-3.0735441) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.30712515) q[1];
sx q[1];
rz(-0.79453429) q[1];
sx q[1];
rz(-1.8737206) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8362379) q[3];
sx q[3];
rz(-0.66514665) q[3];
sx q[3];
rz(-1.0843474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.6924374) q[2];
sx q[2];
rz(-1.038895) q[2];
sx q[2];
rz(0.078744002) q[2];
rz(2.3732843) q[3];
sx q[3];
rz(-1.8470304) q[3];
sx q[3];
rz(-0.19101846) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67097265) q[0];
sx q[0];
rz(-2.5363531) q[0];
sx q[0];
rz(-3.0249) q[0];
rz(2.281588) q[1];
sx q[1];
rz(-1.3415033) q[1];
sx q[1];
rz(-2.2681627) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78055843) q[0];
sx q[0];
rz(-1.3300899) q[0];
sx q[0];
rz(-0.81674256) q[0];
rz(1.8928746) q[2];
sx q[2];
rz(-2.0587741) q[2];
sx q[2];
rz(-2.8805594) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.22120093) q[1];
sx q[1];
rz(-1.6947692) q[1];
sx q[1];
rz(-0.50235475) q[1];
rz(-pi) q[2];
rz(1.6386581) q[3];
sx q[3];
rz(-2.181119) q[3];
sx q[3];
rz(2.7009956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.42150911) q[2];
sx q[2];
rz(-2.8953711) q[2];
sx q[2];
rz(0.949238) q[2];
rz(-1.2468437) q[3];
sx q[3];
rz(-1.4693762) q[3];
sx q[3];
rz(2.5505572) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.128189) q[0];
sx q[0];
rz(-1.0008151) q[0];
sx q[0];
rz(-1.6304954) q[0];
rz(0.36698256) q[1];
sx q[1];
rz(-1.3584247) q[1];
sx q[1];
rz(-0.57986837) q[1];
rz(0.30324528) q[2];
sx q[2];
rz(-2.0187536) q[2];
sx q[2];
rz(3.0222859) q[2];
rz(0.49574481) q[3];
sx q[3];
rz(-2.3875925) q[3];
sx q[3];
rz(2.7238764) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
