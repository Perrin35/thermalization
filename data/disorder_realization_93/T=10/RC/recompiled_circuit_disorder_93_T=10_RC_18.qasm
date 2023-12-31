OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6991601) q[0];
sx q[0];
rz(-1.7572829) q[0];
sx q[0];
rz(1.260489) q[0];
rz(-1.0386382) q[1];
sx q[1];
rz(4.4903978) q[1];
sx q[1];
rz(8.5010565) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7605654) q[0];
sx q[0];
rz(-2.0187223) q[0];
sx q[0];
rz(2.7678124) q[0];
rz(-1.853763) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(-1.9622918) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.10780653) q[1];
sx q[1];
rz(-0.68192712) q[1];
sx q[1];
rz(-0.71180196) q[1];
rz(0.49433319) q[3];
sx q[3];
rz(-1.2327854) q[3];
sx q[3];
rz(2.5288343) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.8929409) q[2];
sx q[2];
rz(-2.9795734) q[2];
rz(-0.93531936) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(2.4285765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.686036) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3128132) q[0];
sx q[0];
rz(-1.7377186) q[0];
sx q[0];
rz(-2.1524327) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7669719) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(0.74707109) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.96869722) q[1];
sx q[1];
rz(-0.52297938) q[1];
sx q[1];
rz(1.5978659) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5442113) q[3];
sx q[3];
rz(-1.1592602) q[3];
sx q[3];
rz(-2.6708024) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(-1.7896174) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75333726) q[0];
sx q[0];
rz(-0.68080807) q[0];
sx q[0];
rz(-0.80048168) q[0];
rz(-3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(1.9690537) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5867509) q[0];
sx q[0];
rz(-1.2625853) q[0];
sx q[0];
rz(0.59535938) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34611361) q[2];
sx q[2];
rz(-2.3556404) q[2];
sx q[2];
rz(-1.2197989) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8923924) q[1];
sx q[1];
rz(-1.8141659) q[1];
sx q[1];
rz(2.1057486) q[1];
rz(-pi) q[2];
rz(-1.8215239) q[3];
sx q[3];
rz(-1.6276976) q[3];
sx q[3];
rz(1.0292605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0744434) q[2];
sx q[2];
rz(-1.477244) q[2];
sx q[2];
rz(-0.91119901) q[2];
rz(-0.95101142) q[3];
sx q[3];
rz(-2.337303) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610385) q[0];
sx q[0];
rz(-0.13042139) q[0];
sx q[0];
rz(-3.0134841) q[0];
rz(-3.065486) q[1];
sx q[1];
rz(-1.9271306) q[1];
sx q[1];
rz(2.6180843) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9273705) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(-0.58332304) q[0];
rz(1.6558311) q[2];
sx q[2];
rz(-1.2418613) q[2];
sx q[2];
rz(-1.0852244) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.88672968) q[1];
sx q[1];
rz(-0.88725315) q[1];
sx q[1];
rz(-0.32943326) q[1];
rz(1.773049) q[3];
sx q[3];
rz(-0.60086717) q[3];
sx q[3];
rz(0.7522538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.6161502) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(2.5775487) q[2];
rz(0.28856746) q[3];
sx q[3];
rz(-2.7189062) q[3];
sx q[3];
rz(-2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6600835) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(-1.4915285) q[0];
rz(-2.2619757) q[1];
sx q[1];
rz(-1.8938226) q[1];
sx q[1];
rz(0.99194828) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0569699) q[0];
sx q[0];
rz(-1.3966494) q[0];
sx q[0];
rz(1.6790381) q[0];
rz(-pi) q[1];
rz(1.8736585) q[2];
sx q[2];
rz(-1.5793243) q[2];
sx q[2];
rz(-1.0964583) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3757513) q[1];
sx q[1];
rz(-1.5305133) q[1];
sx q[1];
rz(-0.30893107) q[1];
rz(-pi) q[2];
rz(-0.12985142) q[3];
sx q[3];
rz(-0.81749812) q[3];
sx q[3];
rz(1.070147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-2.8707855) q[2];
rz(2.9233542) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(-2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72702423) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.7927992) q[0];
rz(2.7596966) q[1];
sx q[1];
rz(-0.31612879) q[1];
sx q[1];
rz(1.7165002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3544918) q[0];
sx q[0];
rz(-1.0757425) q[0];
sx q[0];
rz(1.8255193) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.233333) q[2];
sx q[2];
rz(-0.26974264) q[2];
sx q[2];
rz(-2.595682) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1849991) q[1];
sx q[1];
rz(-1.6374267) q[1];
sx q[1];
rz(2.0208298) q[1];
x q[2];
rz(0.53374966) q[3];
sx q[3];
rz(-1.0373877) q[3];
sx q[3];
rz(-0.48983869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.0075334) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(2.5777204) q[2];
rz(-2.9610736) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-0.40294161) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5086223) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(-1.58889) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99188995) q[0];
sx q[0];
rz(-0.8677965) q[0];
sx q[0];
rz(0.71233149) q[0];
rz(-pi) q[1];
rz(-0.84341151) q[2];
sx q[2];
rz(-1.9551829) q[2];
sx q[2];
rz(2.9316528) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0710443) q[1];
sx q[1];
rz(-0.76862915) q[1];
sx q[1];
rz(-2.8119836) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6131367) q[3];
sx q[3];
rz(-0.65675694) q[3];
sx q[3];
rz(0.17880759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(-2.896893) q[2];
rz(0.129536) q[3];
sx q[3];
rz(-1.9774388) q[3];
sx q[3];
rz(-1.6285508) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
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
rz(1.8364505) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37874052) q[0];
sx q[0];
rz(-1.2438602) q[0];
sx q[0];
rz(-0.61600323) q[0];
x q[1];
rz(-0.70456409) q[2];
sx q[2];
rz(-1.2837871) q[2];
sx q[2];
rz(-1.9412083) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9040363) q[1];
sx q[1];
rz(-2.190553) q[1];
sx q[1];
rz(2.0380286) q[1];
x q[2];
rz(-1.9100902) q[3];
sx q[3];
rz(-1.5486071) q[3];
sx q[3];
rz(2.0458178) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.0017073) q[2];
sx q[2];
rz(-1.7557764) q[2];
sx q[2];
rz(-1.4902327) q[2];
rz(-1.0772609) q[3];
sx q[3];
rz(-2.1765985) q[3];
sx q[3];
rz(-3.0100477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3867144) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(0.3219147) q[0];
rz(-1.6053258) q[1];
sx q[1];
rz(-1.221311) q[1];
sx q[1];
rz(-2.4386491) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.020333175) q[0];
sx q[0];
rz(-2.5998305) q[0];
sx q[0];
rz(-0.74777491) q[0];
rz(0.039673294) q[2];
sx q[2];
rz(-2.0797605) q[2];
sx q[2];
rz(2.2850125) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.326509) q[1];
sx q[1];
rz(-2.106296) q[1];
sx q[1];
rz(-2.9498847) q[1];
rz(-pi) q[2];
rz(-0.33396696) q[3];
sx q[3];
rz(-2.1595862) q[3];
sx q[3];
rz(-0.86563084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(1.8019603) q[2];
rz(-0.30570269) q[3];
sx q[3];
rz(-1.327508) q[3];
sx q[3];
rz(1.8113177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9777578) q[0];
sx q[0];
rz(-0.75755388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(-0.90351358) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(-0.46863619) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9194473) q[0];
sx q[0];
rz(-1.3627909) q[0];
sx q[0];
rz(2.749445) q[0];
rz(-1.956316) q[2];
sx q[2];
rz(-1.8744933) q[2];
sx q[2];
rz(-0.38244837) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.2965282) q[1];
sx q[1];
rz(-1.6390641) q[1];
sx q[1];
rz(1.8370085) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7077984) q[3];
sx q[3];
rz(-1.9926096) q[3];
sx q[3];
rz(2.1123561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.6580711) q[2];
sx q[2];
rz(-1.2926241) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(-3.0269567) q[3];
sx q[3];
rz(-0.95364037) q[3];
sx q[3];
rz(-1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(-2.3399578) q[2];
sx q[2];
rz(-0.87052204) q[2];
sx q[2];
rz(-1.082765) q[2];
rz(-0.91924304) q[3];
sx q[3];
rz(-1.5083434) q[3];
sx q[3];
rz(-1.2283243) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
