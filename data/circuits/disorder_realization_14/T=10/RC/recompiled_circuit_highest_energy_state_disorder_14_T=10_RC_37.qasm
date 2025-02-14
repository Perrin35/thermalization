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
rz(-2.2351216) q[0];
sx q[0];
rz(-1.6989166) q[0];
sx q[0];
rz(-0.28191167) q[0];
rz(0.52892041) q[1];
sx q[1];
rz(4.79098) q[1];
sx q[1];
rz(11.00287) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.603005) q[0];
sx q[0];
rz(-1.9360124) q[0];
sx q[0];
rz(2.9486548) q[0];
rz(-1.5761539) q[2];
sx q[2];
rz(-1.845635) q[2];
sx q[2];
rz(1.9071867) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.2927478) q[1];
sx q[1];
rz(-2.1969921) q[1];
sx q[1];
rz(1.4153773) q[1];
x q[2];
rz(-2.4535577) q[3];
sx q[3];
rz(-1.8092938) q[3];
sx q[3];
rz(-1.9515338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6866744) q[2];
sx q[2];
rz(-2.8439971) q[2];
sx q[2];
rz(2.5285524) q[2];
rz(0.4736627) q[3];
sx q[3];
rz(-1.1981755) q[3];
sx q[3];
rz(-1.4325498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4050201) q[0];
sx q[0];
rz(-0.18017811) q[0];
sx q[0];
rz(2.3905684) q[0];
rz(0.48149064) q[1];
sx q[1];
rz(-1.0844237) q[1];
sx q[1];
rz(-2.1717333) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4879726) q[0];
sx q[0];
rz(-0.75838415) q[0];
sx q[0];
rz(0.55413306) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.18377797) q[2];
sx q[2];
rz(-2.922195) q[2];
sx q[2];
rz(0.32599005) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4837326) q[1];
sx q[1];
rz(-0.93440234) q[1];
sx q[1];
rz(0.32245335) q[1];
rz(-pi) q[2];
rz(2.0131575) q[3];
sx q[3];
rz(-2.1214888) q[3];
sx q[3];
rz(-0.1503508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.91886175) q[2];
sx q[2];
rz(-1.4758045) q[2];
sx q[2];
rz(0.53885031) q[2];
rz(3.1239964) q[3];
sx q[3];
rz(-2.9774057) q[3];
sx q[3];
rz(0.9084475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.4410412) q[0];
sx q[0];
rz(-2.7563162) q[0];
sx q[0];
rz(0.061263099) q[0];
rz(1.6368658) q[1];
sx q[1];
rz(-2.442339) q[1];
sx q[1];
rz(-1.9452852) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29522592) q[0];
sx q[0];
rz(-2.2714473) q[0];
sx q[0];
rz(-2.8267751) q[0];
rz(0.26232403) q[2];
sx q[2];
rz(-1.6379426) q[2];
sx q[2];
rz(3.0284428) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6711025) q[1];
sx q[1];
rz(-1.2467978) q[1];
sx q[1];
rz(-2.1598706) q[1];
rz(0.34516224) q[3];
sx q[3];
rz(-0.90800873) q[3];
sx q[3];
rz(1.2488641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4521744) q[2];
sx q[2];
rz(-1.8345366) q[2];
sx q[2];
rz(2.7090731) q[2];
rz(-2.050926) q[3];
sx q[3];
rz(-2.5011823) q[3];
sx q[3];
rz(-2.045491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5525621) q[0];
sx q[0];
rz(-2.2046389) q[0];
sx q[0];
rz(2.9402148) q[0];
rz(0.51678139) q[1];
sx q[1];
rz(-0.38025451) q[1];
sx q[1];
rz(0.94863272) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7761013) q[0];
sx q[0];
rz(-0.94697111) q[0];
sx q[0];
rz(2.7073949) q[0];
rz(-1.2327551) q[2];
sx q[2];
rz(-2.4043407) q[2];
sx q[2];
rz(-2.8056932) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.81056726) q[1];
sx q[1];
rz(-2.0750891) q[1];
sx q[1];
rz(1.585998) q[1];
rz(-1.6679156) q[3];
sx q[3];
rz(-2.0896974) q[3];
sx q[3];
rz(3.1176381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8676694) q[2];
sx q[2];
rz(-0.40931585) q[2];
sx q[2];
rz(2.5841827) q[2];
rz(-3.0607767) q[3];
sx q[3];
rz(-1.2405453) q[3];
sx q[3];
rz(1.2062937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7252561) q[0];
sx q[0];
rz(-1.6660322) q[0];
sx q[0];
rz(-1.0303372) q[0];
rz(2.985785) q[1];
sx q[1];
rz(-2.0741597) q[1];
sx q[1];
rz(-2.9373998) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0263173) q[0];
sx q[0];
rz(-1.5020619) q[0];
sx q[0];
rz(-0.24253775) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.5459909) q[2];
sx q[2];
rz(-0.24790774) q[2];
sx q[2];
rz(-1.966983) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.1487919) q[1];
sx q[1];
rz(-0.53314236) q[1];
sx q[1];
rz(-0.41457446) q[1];
x q[2];
rz(0.77301003) q[3];
sx q[3];
rz(-1.7436308) q[3];
sx q[3];
rz(-0.59316778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.92945176) q[2];
sx q[2];
rz(-1.6484304) q[2];
sx q[2];
rz(-2.7280651) q[2];
rz(2.2996969) q[3];
sx q[3];
rz(-1.3827518) q[3];
sx q[3];
rz(2.1690185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
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
rz(-3.0464762) q[0];
sx q[0];
rz(-1.0914509) q[0];
sx q[0];
rz(-3.1360151) q[0];
rz(1.3439517) q[1];
sx q[1];
rz(-2.9426212) q[1];
sx q[1];
rz(2.4406348) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13667915) q[0];
sx q[0];
rz(-2.4304996) q[0];
sx q[0];
rz(-0.85791608) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43918886) q[2];
sx q[2];
rz(-0.97676313) q[2];
sx q[2];
rz(-0.33719815) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8474104) q[1];
sx q[1];
rz(-2.8632616) q[1];
sx q[1];
rz(-0.14974447) q[1];
x q[2];
rz(2.4180321) q[3];
sx q[3];
rz(-2.4306707) q[3];
sx q[3];
rz(1.5462745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.31412101) q[2];
sx q[2];
rz(-2.8077329) q[2];
sx q[2];
rz(-1.1673048) q[2];
rz(-2.2589034) q[3];
sx q[3];
rz(-2.3736931) q[3];
sx q[3];
rz(-2.816443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2348787) q[0];
sx q[0];
rz(-2.0218847) q[0];
sx q[0];
rz(0.085302189) q[0];
rz(-0.75434297) q[1];
sx q[1];
rz(-0.5032379) q[1];
sx q[1];
rz(2.1139961) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7639973) q[0];
sx q[0];
rz(-0.21276424) q[0];
sx q[0];
rz(-0.048325267) q[0];
x q[1];
rz(1.5761887) q[2];
sx q[2];
rz(-0.9741592) q[2];
sx q[2];
rz(-2.3062458) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4484549) q[1];
sx q[1];
rz(-1.692728) q[1];
sx q[1];
rz(-1.2454525) q[1];
x q[2];
rz(0.87966921) q[3];
sx q[3];
rz(-2.4604049) q[3];
sx q[3];
rz(-0.78024125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.1520285) q[2];
sx q[2];
rz(-0.98735183) q[2];
sx q[2];
rz(-0.41934553) q[2];
rz(-0.63465214) q[3];
sx q[3];
rz(-0.80497634) q[3];
sx q[3];
rz(1.1254719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80970508) q[0];
sx q[0];
rz(-2.2655847) q[0];
sx q[0];
rz(2.4494655) q[0];
rz(0.57811111) q[1];
sx q[1];
rz(-0.41936857) q[1];
sx q[1];
rz(-1.7587761) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.064474) q[0];
sx q[0];
rz(-3.0713284) q[0];
sx q[0];
rz(1.8514567) q[0];
rz(-pi) q[1];
rz(0.67202576) q[2];
sx q[2];
rz(-1.0254854) q[2];
sx q[2];
rz(1.3915075) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.18305138) q[1];
sx q[1];
rz(-1.019283) q[1];
sx q[1];
rz(2.8666878) q[1];
rz(-2.0973849) q[3];
sx q[3];
rz(-2.3447737) q[3];
sx q[3];
rz(0.45180333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.77070037) q[2];
sx q[2];
rz(-2.6110677) q[2];
sx q[2];
rz(1.4250866) q[2];
rz(-2.6254081) q[3];
sx q[3];
rz(-0.97847146) q[3];
sx q[3];
rz(1.2396575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6462964) q[0];
sx q[0];
rz(-1.4343867) q[0];
sx q[0];
rz(-2.7700951) q[0];
rz(2.3067572) q[1];
sx q[1];
rz(-0.8131665) q[1];
sx q[1];
rz(-3.1239948) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4478701) q[0];
sx q[0];
rz(-1.2945064) q[0];
sx q[0];
rz(1.0280861) q[0];
x q[1];
rz(-1.8210635) q[2];
sx q[2];
rz(-2.2907718) q[2];
sx q[2];
rz(1.7257476) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0068215) q[1];
sx q[1];
rz(-2.4398514) q[1];
sx q[1];
rz(-0.59533466) q[1];
x q[2];
rz(0.60815717) q[3];
sx q[3];
rz(-1.3038692) q[3];
sx q[3];
rz(-1.9723231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0355494) q[2];
sx q[2];
rz(-0.86842662) q[2];
sx q[2];
rz(3.0412728) q[2];
rz(0.22942461) q[3];
sx q[3];
rz(-1.0474297) q[3];
sx q[3];
rz(0.44982287) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.68596524) q[0];
sx q[0];
rz(-2.8534511) q[0];
sx q[0];
rz(-2.567754) q[0];
rz(-1.0539184) q[1];
sx q[1];
rz(-1.046448) q[1];
sx q[1];
rz(0.67646772) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8343) q[0];
sx q[0];
rz(-0.3737475) q[0];
sx q[0];
rz(0.71889241) q[0];
rz(-pi) q[1];
rz(0.25679882) q[2];
sx q[2];
rz(-1.1780103) q[2];
sx q[2];
rz(1.2744255) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.4852462) q[1];
sx q[1];
rz(-2.606483) q[1];
sx q[1];
rz(2.9675304) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0384239) q[3];
sx q[3];
rz(-1.7422973) q[3];
sx q[3];
rz(-0.72491437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4544868) q[2];
sx q[2];
rz(-0.7578308) q[2];
sx q[2];
rz(-1.494361) q[2];
rz(2.2729661) q[3];
sx q[3];
rz(-0.70610154) q[3];
sx q[3];
rz(0.47006616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.070521991) q[0];
sx q[0];
rz(-1.1916397) q[0];
sx q[0];
rz(1.2217039) q[0];
rz(-1.8078049) q[1];
sx q[1];
rz(-1.318327) q[1];
sx q[1];
rz(-0.64073906) q[1];
rz(2.376198) q[2];
sx q[2];
rz(-1.6366048) q[2];
sx q[2];
rz(-0.66726782) q[2];
rz(1.815769) q[3];
sx q[3];
rz(-0.180937) q[3];
sx q[3];
rz(-2.141249) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
