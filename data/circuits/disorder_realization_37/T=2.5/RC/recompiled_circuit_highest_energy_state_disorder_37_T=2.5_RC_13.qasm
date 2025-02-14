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
rz(1.2020741) q[0];
sx q[0];
rz(-2.095686) q[0];
sx q[0];
rz(2.3910971) q[0];
rz(6.5039492) q[1];
sx q[1];
rz(1.060744) q[1];
sx q[1];
rz(4.3252601) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.409371) q[0];
sx q[0];
rz(-2.3187162) q[0];
sx q[0];
rz(2.5636682) q[0];
rz(-pi) q[1];
rz(-1.7375101) q[2];
sx q[2];
rz(-1.6958456) q[2];
sx q[2];
rz(2.2956583) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.88088502) q[1];
sx q[1];
rz(-1.5591677) q[1];
sx q[1];
rz(-2.4616918) q[1];
rz(-pi) q[2];
rz(-0.20211192) q[3];
sx q[3];
rz(-1.8732866) q[3];
sx q[3];
rz(-1.2260395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.815328) q[2];
sx q[2];
rz(-1.999819) q[2];
sx q[2];
rz(-0.092078837) q[2];
rz(2.0103256) q[3];
sx q[3];
rz(-2.0720033) q[3];
sx q[3];
rz(2.2856975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49694127) q[0];
sx q[0];
rz(-0.48321378) q[0];
sx q[0];
rz(-2.605873) q[0];
rz(0.016853111) q[1];
sx q[1];
rz(-0.87937513) q[1];
sx q[1];
rz(-0.0171612) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6494228) q[0];
sx q[0];
rz(-2.0606406) q[0];
sx q[0];
rz(-0.36618284) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13249915) q[2];
sx q[2];
rz(-2.6713058) q[2];
sx q[2];
rz(2.7072226) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8507599) q[1];
sx q[1];
rz(-2.8597745) q[1];
sx q[1];
rz(2.5435104) q[1];
rz(-0.593556) q[3];
sx q[3];
rz(-2.9800804) q[3];
sx q[3];
rz(-2.7167729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0993593) q[2];
sx q[2];
rz(-1.6681654) q[2];
sx q[2];
rz(0.24822203) q[2];
rz(0.92784268) q[3];
sx q[3];
rz(-0.78301269) q[3];
sx q[3];
rz(-0.46970126) q[3];
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
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27702734) q[0];
sx q[0];
rz(-1.1772573) q[0];
sx q[0];
rz(0.16265854) q[0];
rz(0.76205572) q[1];
sx q[1];
rz(-0.22480741) q[1];
sx q[1];
rz(0.83388296) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.57788173) q[0];
sx q[0];
rz(-2.1378008) q[0];
sx q[0];
rz(-1.7162697) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.32270785) q[2];
sx q[2];
rz(-1.6889179) q[2];
sx q[2];
rz(2.8594227) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(3.0030886) q[1];
sx q[1];
rz(-2.2416267) q[1];
sx q[1];
rz(-0.53536949) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0224946) q[3];
sx q[3];
rz(-2.6980711) q[3];
sx q[3];
rz(1.8325072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.062408058) q[2];
sx q[2];
rz(-1.1391613) q[2];
sx q[2];
rz(0.2641693) q[2];
rz(-2.0832113) q[3];
sx q[3];
rz(-2.1975785) q[3];
sx q[3];
rz(-0.03037608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77085483) q[0];
sx q[0];
rz(-0.83468947) q[0];
sx q[0];
rz(-1.7133065) q[0];
rz(-2.3772073) q[1];
sx q[1];
rz(-1.1420219) q[1];
sx q[1];
rz(1.8199325) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8115615) q[0];
sx q[0];
rz(-1.9178277) q[0];
sx q[0];
rz(0.42591947) q[0];
rz(-1.9025397) q[2];
sx q[2];
rz(-1.52196) q[2];
sx q[2];
rz(-1.8845173) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.22239599) q[1];
sx q[1];
rz(-0.42342227) q[1];
sx q[1];
rz(-0.75091305) q[1];
x q[2];
rz(-1.0760078) q[3];
sx q[3];
rz(-0.77829276) q[3];
sx q[3];
rz(-2.8552804) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.6617714) q[2];
sx q[2];
rz(-1.7281374) q[2];
sx q[2];
rz(-0.28523764) q[2];
rz(-1.2626922) q[3];
sx q[3];
rz(-1.9243536) q[3];
sx q[3];
rz(1.7181905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8406463) q[0];
sx q[0];
rz(-2.3941289) q[0];
sx q[0];
rz(-0.88229156) q[0];
rz(2.4709985) q[1];
sx q[1];
rz(-2.0457485) q[1];
sx q[1];
rz(-0.98463279) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6961215) q[0];
sx q[0];
rz(-1.6021671) q[0];
sx q[0];
rz(-0.062792129) q[0];
rz(-0.042858275) q[2];
sx q[2];
rz(-1.0295261) q[2];
sx q[2];
rz(-2.0861911) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.4337816) q[1];
sx q[1];
rz(-1.0278152) q[1];
sx q[1];
rz(0.94306268) q[1];
rz(2.4299942) q[3];
sx q[3];
rz(-2.4715373) q[3];
sx q[3];
rz(-1.7473451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.275445) q[2];
sx q[2];
rz(-2.6068942) q[2];
sx q[2];
rz(2.5246485) q[2];
rz(-2.3894737) q[3];
sx q[3];
rz(-1.813846) q[3];
sx q[3];
rz(2.6822958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
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
rz(0.98944703) q[0];
sx q[0];
rz(-2.4247657) q[0];
sx q[0];
rz(0.77051198) q[0];
rz(0.72692263) q[1];
sx q[1];
rz(-0.22218552) q[1];
sx q[1];
rz(-2.4368584) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73913888) q[0];
sx q[0];
rz(-2.1069035) q[0];
sx q[0];
rz(-2.2741063) q[0];
rz(-pi) q[1];
x q[1];
rz(0.36811916) q[2];
sx q[2];
rz(-1.0794753) q[2];
sx q[2];
rz(1.4368601) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.26624666) q[1];
sx q[1];
rz(-0.35142144) q[1];
sx q[1];
rz(0.09697341) q[1];
rz(-pi) q[2];
rz(-2.8877657) q[3];
sx q[3];
rz(-2.1269264) q[3];
sx q[3];
rz(1.1185916) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.7695693) q[2];
sx q[2];
rz(-1.4906887) q[2];
sx q[2];
rz(-2.7083569) q[2];
rz(2.9406934) q[3];
sx q[3];
rz(-1.9379987) q[3];
sx q[3];
rz(2.6065629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
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
rz(0.83605003) q[0];
sx q[0];
rz(-1.9925646) q[0];
sx q[0];
rz(-2.7591163) q[0];
rz(1.0410694) q[1];
sx q[1];
rz(-1.6761227) q[1];
sx q[1];
rz(-2.6952851) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8167717) q[0];
sx q[0];
rz(-1.8292301) q[0];
sx q[0];
rz(1.6245317) q[0];
x q[1];
rz(-0.52526229) q[2];
sx q[2];
rz(-2.7064311) q[2];
sx q[2];
rz(-2.9829142) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.32532) q[1];
sx q[1];
rz(-2.153646) q[1];
sx q[1];
rz(2.7144949) q[1];
x q[2];
rz(-0.84844037) q[3];
sx q[3];
rz(-0.99096837) q[3];
sx q[3];
rz(-2.8715796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.81991601) q[2];
sx q[2];
rz(-2.0240462) q[2];
sx q[2];
rz(2.4134911) q[2];
rz(2.8150832) q[3];
sx q[3];
rz(-1.6312586) q[3];
sx q[3];
rz(0.85730332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9051055) q[0];
sx q[0];
rz(-2.7323854) q[0];
sx q[0];
rz(-1.2976728) q[0];
rz(2.1406651) q[1];
sx q[1];
rz(-0.74556723) q[1];
sx q[1];
rz(-2.7330858) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9169642) q[0];
sx q[0];
rz(-0.014212) q[0];
sx q[0];
rz(1.5111708) q[0];
rz(-pi) q[1];
rz(0.82168545) q[2];
sx q[2];
rz(-0.61640451) q[2];
sx q[2];
rz(2.6954755) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0169819) q[1];
sx q[1];
rz(-1.7045741) q[1];
sx q[1];
rz(0.11511187) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85902416) q[3];
sx q[3];
rz(-2.6005496) q[3];
sx q[3];
rz(2.4699901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-3.0268176) q[2];
sx q[2];
rz(-1.5207745) q[2];
sx q[2];
rz(0.82359037) q[2];
rz(-2.1987727) q[3];
sx q[3];
rz(-1.4225682) q[3];
sx q[3];
rz(-1.6387117) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5470062) q[0];
sx q[0];
rz(-2.6972045) q[0];
sx q[0];
rz(-1.2521) q[0];
rz(1.3638672) q[1];
sx q[1];
rz(-1.6412647) q[1];
sx q[1];
rz(-1.3346671) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20864381) q[0];
sx q[0];
rz(-1.2675537) q[0];
sx q[0];
rz(-1.6088845) q[0];
rz(-pi) q[1];
rz(1.7450856) q[2];
sx q[2];
rz(-2.2054923) q[2];
sx q[2];
rz(-1.2668934) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.081055911) q[1];
sx q[1];
rz(-1.4976275) q[1];
sx q[1];
rz(2.6899978) q[1];
x q[2];
rz(1.077721) q[3];
sx q[3];
rz(-2.6209313) q[3];
sx q[3];
rz(2.9021104) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-3.1198662) q[2];
sx q[2];
rz(-0.29325565) q[2];
sx q[2];
rz(-0.29652706) q[2];
rz(-1.6622539) q[3];
sx q[3];
rz(-1.5115503) q[3];
sx q[3];
rz(0.1571981) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5088365) q[0];
sx q[0];
rz(-0.65082508) q[0];
sx q[0];
rz(0.79488361) q[0];
rz(0.55566135) q[1];
sx q[1];
rz(-1.2763005) q[1];
sx q[1];
rz(1.4478252) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5710766) q[0];
sx q[0];
rz(-0.85502842) q[0];
sx q[0];
rz(-1.060673) q[0];
rz(-1.6332288) q[2];
sx q[2];
rz(-1.0737906) q[2];
sx q[2];
rz(1.7071498) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9196133) q[1];
sx q[1];
rz(-1.7697502) q[1];
sx q[1];
rz(2.0556021) q[1];
rz(-2.1673739) q[3];
sx q[3];
rz(-2.0250626) q[3];
sx q[3];
rz(0.59859401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.29622233) q[2];
sx q[2];
rz(-1.6310548) q[2];
sx q[2];
rz(1.4307865) q[2];
rz(-0.75302643) q[3];
sx q[3];
rz(-0.81792653) q[3];
sx q[3];
rz(0.16451612) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0059148) q[0];
sx q[0];
rz(-0.68886859) q[0];
sx q[0];
rz(-1.9266358) q[0];
rz(1.6750492) q[1];
sx q[1];
rz(-2.2129682) q[1];
sx q[1];
rz(0.60563544) q[1];
rz(2.3382414) q[2];
sx q[2];
rz(-2.4016082) q[2];
sx q[2];
rz(-0.27863816) q[2];
rz(-0.10726992) q[3];
sx q[3];
rz(-1.4013877) q[3];
sx q[3];
rz(0.79584484) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
