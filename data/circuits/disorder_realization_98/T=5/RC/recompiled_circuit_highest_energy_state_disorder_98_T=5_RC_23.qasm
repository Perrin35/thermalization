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
rz(-0.67796081) q[0];
sx q[0];
rz(6.0089587) q[0];
sx q[0];
rz(11.275509) q[0];
rz(2.9277247) q[1];
sx q[1];
rz(-0.31615058) q[1];
sx q[1];
rz(1.6419799) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0555178) q[0];
sx q[0];
rz(-2.7788305) q[0];
sx q[0];
rz(-2.4253009) q[0];
rz(-pi) q[1];
rz(-1.2798645) q[2];
sx q[2];
rz(-2.3290344) q[2];
sx q[2];
rz(-1.9951374) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1823481) q[1];
sx q[1];
rz(-2.3492491) q[1];
sx q[1];
rz(-2.6807129) q[1];
rz(-pi) q[2];
rz(-0.47874449) q[3];
sx q[3];
rz(-2.1447721) q[3];
sx q[3];
rz(0.51306242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.9700254) q[2];
sx q[2];
rz(-1.0415404) q[2];
sx q[2];
rz(0.901326) q[2];
rz(-2.6775635) q[3];
sx q[3];
rz(-1.8601067) q[3];
sx q[3];
rz(0.06981167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0932015) q[0];
sx q[0];
rz(-3.1049325) q[0];
sx q[0];
rz(-1.2777591) q[0];
rz(0.094206421) q[1];
sx q[1];
rz(-0.46368805) q[1];
sx q[1];
rz(1.5346079) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8600991) q[0];
sx q[0];
rz(-1.6035542) q[0];
sx q[0];
rz(-2.1707235) q[0];
x q[1];
rz(3.0846527) q[2];
sx q[2];
rz(-1.3888265) q[2];
sx q[2];
rz(1.6873311) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6834641) q[1];
sx q[1];
rz(-0.56741112) q[1];
sx q[1];
rz(2.6415884) q[1];
rz(-pi) q[2];
x q[2];
rz(0.8022763) q[3];
sx q[3];
rz(-0.55219383) q[3];
sx q[3];
rz(1.4331762) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4216807) q[2];
sx q[2];
rz(-1.9399425) q[2];
sx q[2];
rz(0.16033944) q[2];
rz(-2.4028589) q[3];
sx q[3];
rz(-0.84638798) q[3];
sx q[3];
rz(1.7140478) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62683231) q[0];
sx q[0];
rz(-1.4381831) q[0];
sx q[0];
rz(-2.6318188) q[0];
rz(-1.7104507) q[1];
sx q[1];
rz(-1.1809843) q[1];
sx q[1];
rz(0.74657718) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1603161) q[0];
sx q[0];
rz(-1.0183987) q[0];
sx q[0];
rz(0.023560087) q[0];
rz(-pi) q[1];
rz(-0.72991972) q[2];
sx q[2];
rz(-2.0509208) q[2];
sx q[2];
rz(1.4381222) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.67184) q[1];
sx q[1];
rz(-0.79974175) q[1];
sx q[1];
rz(-1.0554764) q[1];
rz(2.0381365) q[3];
sx q[3];
rz(-0.73690542) q[3];
sx q[3];
rz(1.3947427) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.15472445) q[2];
sx q[2];
rz(-0.26686033) q[2];
sx q[2];
rz(-1.9179087) q[2];
rz(2.0679421) q[3];
sx q[3];
rz(-1.7045538) q[3];
sx q[3];
rz(2.2435718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2735485) q[0];
sx q[0];
rz(-0.88669625) q[0];
sx q[0];
rz(2.7622188) q[0];
rz(0.1768449) q[1];
sx q[1];
rz(-1.481448) q[1];
sx q[1];
rz(-2.3462229) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.11957) q[0];
sx q[0];
rz(-1.9384465) q[0];
sx q[0];
rz(3.0421542) q[0];
x q[1];
rz(-1.8537117) q[2];
sx q[2];
rz(-0.77436111) q[2];
sx q[2];
rz(-2.0191569) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6713555) q[1];
sx q[1];
rz(-1.0507601) q[1];
sx q[1];
rz(0.29754559) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5136817) q[3];
sx q[3];
rz(-1.4453672) q[3];
sx q[3];
rz(-0.40298395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.4312326) q[2];
sx q[2];
rz(-1.7962339) q[2];
sx q[2];
rz(1.7809407) q[2];
rz(1.0294754) q[3];
sx q[3];
rz(-1.6607213) q[3];
sx q[3];
rz(1.8195389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65115702) q[0];
sx q[0];
rz(-1.54162) q[0];
sx q[0];
rz(-1.5929476) q[0];
rz(-2.7410638) q[1];
sx q[1];
rz(-1.488204) q[1];
sx q[1];
rz(-1.4097479) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9072394) q[0];
sx q[0];
rz(-0.84915224) q[0];
sx q[0];
rz(-1.3035151) q[0];
x q[1];
rz(-0.39938853) q[2];
sx q[2];
rz(-2.7542973) q[2];
sx q[2];
rz(3.0687817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.21239195) q[1];
sx q[1];
rz(-0.98941313) q[1];
sx q[1];
rz(0.53668944) q[1];
rz(-2.1928113) q[3];
sx q[3];
rz(-1.5153441) q[3];
sx q[3];
rz(0.90622073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3132402) q[2];
sx q[2];
rz(-1.7007549) q[2];
sx q[2];
rz(-2.7102615) q[2];
rz(-3.0865772) q[3];
sx q[3];
rz(-2.8300245) q[3];
sx q[3];
rz(2.5855248) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3738275) q[0];
sx q[0];
rz(-2.8887833) q[0];
sx q[0];
rz(-1.025169) q[0];
rz(-0.45285666) q[1];
sx q[1];
rz(-0.57944524) q[1];
sx q[1];
rz(-2.2377009) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5885279) q[0];
sx q[0];
rz(-0.72719535) q[0];
sx q[0];
rz(-1.0968047) q[0];
x q[1];
rz(-0.63389961) q[2];
sx q[2];
rz(-1.0154775) q[2];
sx q[2];
rz(-0.429053) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.3749668) q[1];
sx q[1];
rz(-1.300525) q[1];
sx q[1];
rz(0.75281669) q[1];
rz(-pi) q[2];
rz(-1.4056095) q[3];
sx q[3];
rz(-1.7398698) q[3];
sx q[3];
rz(2.1220834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.30770939) q[2];
sx q[2];
rz(-1.4557975) q[2];
sx q[2];
rz(-0.21635381) q[2];
rz(0.89573914) q[3];
sx q[3];
rz(-0.68550617) q[3];
sx q[3];
rz(1.640813) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1322587) q[0];
sx q[0];
rz(-2.0348771) q[0];
sx q[0];
rz(-2.4427781) q[0];
rz(1.9578594) q[1];
sx q[1];
rz(-1.7203169) q[1];
sx q[1];
rz(0.85174495) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18139938) q[0];
sx q[0];
rz(-1.532868) q[0];
sx q[0];
rz(2.1071803) q[0];
rz(0.039489517) q[2];
sx q[2];
rz(-2.563884) q[2];
sx q[2];
rz(0.60009391) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.70515436) q[1];
sx q[1];
rz(-2.2346304) q[1];
sx q[1];
rz(-1.1052119) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4999291) q[3];
sx q[3];
rz(-1.9681276) q[3];
sx q[3];
rz(-1.7375377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.304004) q[2];
sx q[2];
rz(-1.3112661) q[2];
sx q[2];
rz(-2.6336929) q[2];
rz(-1.0287644) q[3];
sx q[3];
rz(-2.2977836) q[3];
sx q[3];
rz(-0.60104162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.69865882) q[0];
sx q[0];
rz(-1.826094) q[0];
sx q[0];
rz(0.0096631924) q[0];
rz(1.6161605) q[1];
sx q[1];
rz(-1.3776255) q[1];
sx q[1];
rz(-2.756871) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5670861) q[0];
sx q[0];
rz(-0.28038803) q[0];
sx q[0];
rz(1.0100421) q[0];
x q[1];
rz(1.184395) q[2];
sx q[2];
rz(-1.3829872) q[2];
sx q[2];
rz(1.1593429) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3166805) q[1];
sx q[1];
rz(-0.34873617) q[1];
sx q[1];
rz(0.92479264) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1423903) q[3];
sx q[3];
rz(-0.88547844) q[3];
sx q[3];
rz(0.1911605) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.5821417) q[2];
sx q[2];
rz(-2.1820575) q[2];
sx q[2];
rz(3.1257296) q[2];
rz(-1.9150241) q[3];
sx q[3];
rz(-0.80569402) q[3];
sx q[3];
rz(2.5198643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7568307) q[0];
sx q[0];
rz(-2.4767196) q[0];
sx q[0];
rz(1.0913947) q[0];
rz(0.81820828) q[1];
sx q[1];
rz(-1.3312157) q[1];
sx q[1];
rz(2.4726726) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0398934) q[0];
sx q[0];
rz(-2.153206) q[0];
sx q[0];
rz(-2.3132313) q[0];
x q[1];
rz(1.0878272) q[2];
sx q[2];
rz(-2.3685799) q[2];
sx q[2];
rz(2.9376466) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.131625) q[1];
sx q[1];
rz(-0.56198453) q[1];
sx q[1];
rz(-1.704292) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.091033) q[3];
sx q[3];
rz(-0.63891131) q[3];
sx q[3];
rz(-1.0425488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6341256) q[2];
sx q[2];
rz(-0.52906817) q[2];
sx q[2];
rz(2.1779306) q[2];
rz(-1.4108747) q[3];
sx q[3];
rz(-1.4286634) q[3];
sx q[3];
rz(-1.3655519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1249579) q[0];
sx q[0];
rz(-2.5732915) q[0];
sx q[0];
rz(-1.6868663) q[0];
rz(1.1085054) q[1];
sx q[1];
rz(-1.6051555) q[1];
sx q[1];
rz(-0.32807168) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1301918) q[0];
sx q[0];
rz(-1.5131803) q[0];
sx q[0];
rz(-1.685623) q[0];
x q[1];
rz(0.56888442) q[2];
sx q[2];
rz(-0.21368229) q[2];
sx q[2];
rz(-1.793022) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3076664) q[1];
sx q[1];
rz(-2.1824679) q[1];
sx q[1];
rz(-2.0425379) q[1];
x q[2];
rz(1.0193908) q[3];
sx q[3];
rz(-2.9388722) q[3];
sx q[3];
rz(-0.77714506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.46128094) q[2];
sx q[2];
rz(-1.2993456) q[2];
sx q[2];
rz(2.7247562) q[2];
rz(2.4428115) q[3];
sx q[3];
rz(-2.4650033) q[3];
sx q[3];
rz(2.7509287) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9417435) q[0];
sx q[0];
rz(-1.9244292) q[0];
sx q[0];
rz(1.695965) q[0];
rz(-0.70882123) q[1];
sx q[1];
rz(-2.2292021) q[1];
sx q[1];
rz(3.0880047) q[1];
rz(2.9403953) q[2];
sx q[2];
rz(-0.94379776) q[2];
sx q[2];
rz(-1.1371053) q[2];
rz(2.1772467) q[3];
sx q[3];
rz(-0.83899211) q[3];
sx q[3];
rz(-0.28750026) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
