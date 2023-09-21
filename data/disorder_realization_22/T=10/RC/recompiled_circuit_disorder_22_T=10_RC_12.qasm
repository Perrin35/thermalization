OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7093631) q[0];
sx q[0];
rz(-2.1837283) q[0];
sx q[0];
rz(-0.14444484) q[0];
rz(-2.5748409) q[1];
sx q[1];
rz(-2.6161939) q[1];
sx q[1];
rz(2.1638343) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.190783) q[0];
sx q[0];
rz(-1.3124183) q[0];
sx q[0];
rz(1.6590614) q[0];
rz(-pi) q[1];
rz(1.058504) q[2];
sx q[2];
rz(-0.57527486) q[2];
sx q[2];
rz(-0.011205999) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6482918) q[1];
sx q[1];
rz(-2.9420605) q[1];
sx q[1];
rz(2.3471911) q[1];
rz(1.7498738) q[3];
sx q[3];
rz(-0.79154166) q[3];
sx q[3];
rz(2.4812428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4686761) q[2];
sx q[2];
rz(-1.1489931) q[2];
sx q[2];
rz(-2.2093175) q[2];
rz(-2.9428234) q[3];
sx q[3];
rz(-1.1105744) q[3];
sx q[3];
rz(-0.96536243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(-0.7070049) q[0];
sx q[0];
rz(-2.2362237) q[0];
sx q[0];
rz(2.7804651) q[0];
rz(-1.4350285) q[1];
sx q[1];
rz(-1.3577434) q[1];
sx q[1];
rz(-2.3235869) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.60048238) q[0];
sx q[0];
rz(-2.6594909) q[0];
sx q[0];
rz(0.62645285) q[0];
x q[1];
rz(-1.2641764) q[2];
sx q[2];
rz(-2.7536256) q[2];
sx q[2];
rz(2.1167133) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9243014) q[1];
sx q[1];
rz(-1.1326619) q[1];
sx q[1];
rz(2.2469254) q[1];
x q[2];
rz(-0.057095842) q[3];
sx q[3];
rz(-1.207587) q[3];
sx q[3];
rz(2.925194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.25257418) q[2];
sx q[2];
rz(-0.44712862) q[2];
sx q[2];
rz(1.7209631) q[2];
rz(1.3160926) q[3];
sx q[3];
rz(-0.75794739) q[3];
sx q[3];
rz(2.7533598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4017568) q[0];
sx q[0];
rz(-0.49210423) q[0];
sx q[0];
rz(-0.87093583) q[0];
rz(-0.31618205) q[1];
sx q[1];
rz(-0.28156391) q[1];
sx q[1];
rz(-0.2972163) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0270099) q[0];
sx q[0];
rz(-1.6342589) q[0];
sx q[0];
rz(-2.7364536) q[0];
x q[1];
rz(2.259953) q[2];
sx q[2];
rz(-2.214553) q[2];
sx q[2];
rz(-1.8635441) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3618468) q[1];
sx q[1];
rz(-1.6184813) q[1];
sx q[1];
rz(1.7576799) q[1];
rz(-pi) q[2];
rz(-0.8244332) q[3];
sx q[3];
rz(-1.0472877) q[3];
sx q[3];
rz(-1.6214961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3729942) q[2];
sx q[2];
rz(-2.3186389) q[2];
sx q[2];
rz(-2.7139943) q[2];
rz(-1.9528495) q[3];
sx q[3];
rz(-2.5116428) q[3];
sx q[3];
rz(-2.6141613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7052085) q[0];
sx q[0];
rz(-1.5773062) q[0];
sx q[0];
rz(0.77392459) q[0];
rz(-2.4286843) q[1];
sx q[1];
rz(-1.0293101) q[1];
sx q[1];
rz(-2.4598222) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49127689) q[0];
sx q[0];
rz(-0.33704764) q[0];
sx q[0];
rz(-2.2469673) q[0];
rz(-pi) q[1];
rz(-1.8800456) q[2];
sx q[2];
rz(-2.5144858) q[2];
sx q[2];
rz(-2.4583465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.75979739) q[1];
sx q[1];
rz(-2.0054381) q[1];
sx q[1];
rz(0.97310193) q[1];
rz(1.2069615) q[3];
sx q[3];
rz(-1.2483276) q[3];
sx q[3];
rz(2.913167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1327847) q[2];
sx q[2];
rz(-2.4286353) q[2];
sx q[2];
rz(-2.183389) q[2];
rz(-1.0664252) q[3];
sx q[3];
rz(-1.8582148) q[3];
sx q[3];
rz(-2.7619894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.585007) q[0];
sx q[0];
rz(-1.7946694) q[0];
sx q[0];
rz(0.56030309) q[0];
rz(2.141748) q[1];
sx q[1];
rz(-0.20345774) q[1];
sx q[1];
rz(-1.6220185) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1513838) q[0];
sx q[0];
rz(-0.68400331) q[0];
sx q[0];
rz(0.74599501) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.6117758) q[2];
sx q[2];
rz(-2.0423186) q[2];
sx q[2];
rz(-0.39786354) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.49415627) q[1];
sx q[1];
rz(-1.7687706) q[1];
sx q[1];
rz(-0.76311771) q[1];
x q[2];
rz(-2.1836957) q[3];
sx q[3];
rz(-1.6727722) q[3];
sx q[3];
rz(1.0979872) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4524298) q[2];
sx q[2];
rz(-0.55183691) q[2];
sx q[2];
rz(2.9186644) q[2];
rz(-0.034742268) q[3];
sx q[3];
rz(-1.7545173) q[3];
sx q[3];
rz(-3.070014) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3361622) q[0];
sx q[0];
rz(-0.31328377) q[0];
sx q[0];
rz(-2.0715332) q[0];
rz(-1.7806212) q[1];
sx q[1];
rz(-0.37934163) q[1];
sx q[1];
rz(1.8575352) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59239913) q[0];
sx q[0];
rz(-2.1201029) q[0];
sx q[0];
rz(2.4248289) q[0];
rz(-pi) q[1];
rz(1.1805004) q[2];
sx q[2];
rz(-1.830606) q[2];
sx q[2];
rz(0.63894546) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.2388873) q[1];
sx q[1];
rz(-1.1088015) q[1];
sx q[1];
rz(1.8932896) q[1];
rz(-2.8347557) q[3];
sx q[3];
rz(-2.0862498) q[3];
sx q[3];
rz(-3.0981578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.47508919) q[2];
sx q[2];
rz(-1.44375) q[2];
sx q[2];
rz(2.5040023) q[2];
rz(-0.83667886) q[3];
sx q[3];
rz(-1.1058608) q[3];
sx q[3];
rz(1.7975413) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5752983) q[0];
sx q[0];
rz(-2.7116382) q[0];
sx q[0];
rz(0.56754011) q[0];
rz(0.42770806) q[1];
sx q[1];
rz(-1.5274915) q[1];
sx q[1];
rz(2.2033851) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0465614) q[0];
sx q[0];
rz(-1.9782269) q[0];
sx q[0];
rz(0.33336063) q[0];
x q[1];
rz(-2.5222048) q[2];
sx q[2];
rz(-0.833138) q[2];
sx q[2];
rz(-0.39820652) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.2141014) q[1];
sx q[1];
rz(-0.92988211) q[1];
sx q[1];
rz(-2.6919041) q[1];
rz(-pi) q[2];
rz(-2.4343743) q[3];
sx q[3];
rz(-1.3361317) q[3];
sx q[3];
rz(0.16682391) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.87970916) q[2];
sx q[2];
rz(-1.1738913) q[2];
sx q[2];
rz(1.3860469) q[2];
rz(-1.322768) q[3];
sx q[3];
rz(-1.991792) q[3];
sx q[3];
rz(-3.1125606) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26738527) q[0];
sx q[0];
rz(-0.33339849) q[0];
sx q[0];
rz(-1.4338795) q[0];
rz(-1.2738312) q[1];
sx q[1];
rz(-2.0055983) q[1];
sx q[1];
rz(-2.3103255) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6600663) q[0];
sx q[0];
rz(-2.7114081) q[0];
sx q[0];
rz(0.35920401) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.93584658) q[2];
sx q[2];
rz(-1.4668462) q[2];
sx q[2];
rz(-2.7452591) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.23849328) q[1];
sx q[1];
rz(-1.0817173) q[1];
sx q[1];
rz(-2.3563983) q[1];
rz(-1.3475111) q[3];
sx q[3];
rz(-1.3459473) q[3];
sx q[3];
rz(0.68883483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.5780118) q[2];
sx q[2];
rz(-1.3728023) q[2];
sx q[2];
rz(2.1772299) q[2];
rz(-1.1635121) q[3];
sx q[3];
rz(-2.5301299) q[3];
sx q[3];
rz(-2.4826629) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7678541) q[0];
sx q[0];
rz(-0.73497325) q[0];
sx q[0];
rz(-0.97737616) q[0];
rz(1.3865698) q[1];
sx q[1];
rz(-1.8354548) q[1];
sx q[1];
rz(-1.1057373) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10406636) q[0];
sx q[0];
rz(-2.6797047) q[0];
sx q[0];
rz(0.26432963) q[0];
x q[1];
rz(-2.2374723) q[2];
sx q[2];
rz(-1.6534001) q[2];
sx q[2];
rz(0.77043515) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92377201) q[1];
sx q[1];
rz(-2.8433228) q[1];
sx q[1];
rz(-0.89426269) q[1];
rz(-1.8589784) q[3];
sx q[3];
rz(-1.9087221) q[3];
sx q[3];
rz(-1.6328788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.608312) q[2];
sx q[2];
rz(-2.7567342) q[2];
sx q[2];
rz(2.9917955) q[2];
rz(-1.7685361) q[3];
sx q[3];
rz(-1.7356197) q[3];
sx q[3];
rz(-1.8201374) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6479284) q[0];
sx q[0];
rz(-2.6239008) q[0];
sx q[0];
rz(3.0143484) q[0];
rz(1.6607025) q[1];
sx q[1];
rz(-0.36247411) q[1];
sx q[1];
rz(-0.1677992) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.4807178) q[0];
sx q[0];
rz(-1.6010451) q[0];
sx q[0];
rz(-1.589993) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4160412) q[2];
sx q[2];
rz(-0.96826474) q[2];
sx q[2];
rz(3.0263911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4976616) q[1];
sx q[1];
rz(-1.5776458) q[1];
sx q[1];
rz(2.5531406) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6889257) q[3];
sx q[3];
rz(-1.6932634) q[3];
sx q[3];
rz(0.99692217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.026713513) q[2];
sx q[2];
rz(-0.9388887) q[2];
sx q[2];
rz(2.3804469) q[2];
rz(3.051565) q[3];
sx q[3];
rz(-2.138425) q[3];
sx q[3];
rz(-0.95054039) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5933843) q[0];
sx q[0];
rz(-1.9822639) q[0];
sx q[0];
rz(-0.32250861) q[0];
rz(-2.7535915) q[1];
sx q[1];
rz(-1.7419659) q[1];
sx q[1];
rz(2.3566125) q[1];
rz(2.2446752) q[2];
sx q[2];
rz(-2.2581836) q[2];
sx q[2];
rz(2.7473292) q[2];
rz(2.9705863) q[3];
sx q[3];
rz(-2.1033559) q[3];
sx q[3];
rz(-0.81567473) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];