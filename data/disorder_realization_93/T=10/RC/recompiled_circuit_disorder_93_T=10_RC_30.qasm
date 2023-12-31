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
rz(1.8811037) q[0];
rz(2.1029544) q[1];
sx q[1];
rz(-1.3488052) q[1];
sx q[1];
rz(0.92372149) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0246575) q[0];
sx q[0];
rz(-2.5664461) q[0];
sx q[0];
rz(2.2206109) q[0];
x q[1];
rz(1.853763) q[2];
sx q[2];
rz(-1.4322865) q[2];
sx q[2];
rz(1.9622918) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.87286283) q[1];
sx q[1];
rz(-1.1464719) q[1];
sx q[1];
rz(0.55117589) q[1];
rz(-1.9507017) q[3];
sx q[3];
rz(-2.0348747) q[3];
sx q[3];
rz(-0.78117785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91360056) q[2];
sx q[2];
rz(-1.2486518) q[2];
sx q[2];
rz(-0.16201924) q[2];
rz(0.93531936) q[3];
sx q[3];
rz(-0.98615065) q[3];
sx q[3];
rz(0.71301618) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0733923) q[0];
sx q[0];
rz(-0.22664264) q[0];
sx q[0];
rz(1.9447928) q[0];
rz(-2.4616922) q[1];
sx q[1];
rz(-2.6459243) q[1];
sx q[1];
rz(-1.4555567) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8287795) q[0];
sx q[0];
rz(-1.7377186) q[0];
sx q[0];
rz(2.1524327) q[0];
rz(-0.37462072) q[2];
sx q[2];
rz(-1.4943559) q[2];
sx q[2];
rz(-2.3945216) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1728954) q[1];
sx q[1];
rz(-0.52297938) q[1];
sx q[1];
rz(-1.5437267) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0807642) q[3];
sx q[3];
rz(-2.7292477) q[3];
sx q[3];
rz(-2.6044248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.7130647) q[2];
sx q[2];
rz(-1.45168) q[2];
sx q[2];
rz(1.3519752) q[2];
rz(0.18243608) q[3];
sx q[3];
rz(-2.1648516) q[3];
sx q[3];
rz(0.3119719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75333726) q[0];
sx q[0];
rz(-2.4607846) q[0];
sx q[0];
rz(-2.341111) q[0];
rz(3.1128186) q[1];
sx q[1];
rz(-1.0556227) q[1];
sx q[1];
rz(1.172539) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55484178) q[0];
sx q[0];
rz(-1.2625853) q[0];
sx q[0];
rz(0.59535938) q[0];
rz(-pi) q[1];
x q[1];
rz(2.795479) q[2];
sx q[2];
rz(-0.78595224) q[2];
sx q[2];
rz(-1.2197989) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0769656) q[1];
sx q[1];
rz(-0.58276999) q[1];
sx q[1];
rz(2.0240192) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.058733744) q[3];
sx q[3];
rz(-1.8211094) q[3];
sx q[3];
rz(0.55610031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.0671493) q[2];
sx q[2];
rz(-1.6643486) q[2];
sx q[2];
rz(-2.2303936) q[2];
rz(-2.1905812) q[3];
sx q[3];
rz(-0.8042897) q[3];
sx q[3];
rz(-0.89200154) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7610385) q[0];
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
rz(1.9273705) q[0];
sx q[0];
rz(-0.71338755) q[0];
sx q[0];
rz(-2.5582696) q[0];
x q[1];
rz(-2.8115494) q[2];
sx q[2];
rz(-1.651262) q[2];
sx q[2];
rz(0.45804322) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7586786) q[1];
sx q[1];
rz(-2.3944693) q[1];
sx q[1];
rz(1.1927356) q[1];
rz(-pi) q[2];
rz(-1.773049) q[3];
sx q[3];
rz(-0.60086717) q[3];
sx q[3];
rz(2.3893389) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.52544242) q[2];
sx q[2];
rz(-1.5443065) q[2];
sx q[2];
rz(0.564044) q[2];
rz(-0.28856746) q[3];
sx q[3];
rz(-0.42268649) q[3];
sx q[3];
rz(-2.585876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.48150912) q[0];
sx q[0];
rz(-2.4531589) q[0];
sx q[0];
rz(1.6500641) q[0];
rz(2.2619757) q[1];
sx q[1];
rz(-1.2477701) q[1];
sx q[1];
rz(-2.1496444) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.49682) q[0];
sx q[0];
rz(-0.20475514) q[0];
sx q[0];
rz(-0.55069189) q[0];
rz(0.0089346272) q[2];
sx q[2];
rz(-1.2679456) q[2];
sx q[2];
rz(-2.6645899) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.3494898) q[1];
sx q[1];
rz(-1.8794685) q[1];
sx q[1];
rz(1.5285138) q[1];
rz(-pi) q[2];
rz(1.7080073) q[3];
sx q[3];
rz(-0.76223323) q[3];
sx q[3];
rz(1.2587794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0118959) q[2];
sx q[2];
rz(-2.7719438) q[2];
sx q[2];
rz(-2.8707855) q[2];
rz(-0.21823847) q[3];
sx q[3];
rz(-1.3202347) q[3];
sx q[3];
rz(-2.9158084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4145684) q[0];
sx q[0];
rz(-2.4232061) q[0];
sx q[0];
rz(-1.3487934) q[0];
rz(0.38189608) q[1];
sx q[1];
rz(-2.8254639) q[1];
sx q[1];
rz(1.7165002) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0483635) q[0];
sx q[0];
rz(-1.7943802) q[0];
sx q[0];
rz(-0.50888749) q[0];
rz(-pi) q[1];
rz(-1.9082597) q[2];
sx q[2];
rz(-0.26974264) q[2];
sx q[2];
rz(-0.5459107) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.35364321) q[1];
sx q[1];
rz(-2.0197581) q[1];
sx q[1];
rz(3.06762) q[1];
x q[2];
rz(-0.85961996) q[3];
sx q[3];
rz(-2.4058127) q[3];
sx q[3];
rz(-1.3501292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0075334) q[2];
sx q[2];
rz(-2.7441661) q[2];
sx q[2];
rz(0.56387222) q[2];
rz(0.18051906) q[3];
sx q[3];
rz(-1.518395) q[3];
sx q[3];
rz(-0.40294161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6329704) q[0];
sx q[0];
rz(-0.16212012) q[0];
sx q[0];
rz(-2.7222743) q[0];
rz(1.5527027) q[1];
sx q[1];
rz(-1.2607375) q[1];
sx q[1];
rz(0.82180506) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069698378) q[0];
sx q[0];
rz(-2.0928203) q[0];
sx q[0];
rz(2.4126023) q[0];
rz(-0.49634883) q[2];
sx q[2];
rz(-2.2349572) q[2];
sx q[2];
rz(-1.4585444) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.2591178) q[1];
sx q[1];
rz(-1.7977409) q[1];
sx q[1];
rz(-0.74101733) q[1];
rz(-2.5542198) q[3];
sx q[3];
rz(-1.8837187) q[3];
sx q[3];
rz(1.8250993) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.8043148) q[2];
sx q[2];
rz(-2.3874805) q[2];
sx q[2];
rz(0.24469963) q[2];
rz(-0.129536) q[3];
sx q[3];
rz(-1.1641538) q[3];
sx q[3];
rz(1.5130419) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41641763) q[0];
sx q[0];
rz(-0.019151909) q[0];
sx q[0];
rz(-0.82292557) q[0];
rz(2.8322463) q[1];
sx q[1];
rz(-1.3920709) q[1];
sx q[1];
rz(1.8364505) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1730723) q[0];
sx q[0];
rz(-0.99181306) q[0];
sx q[0];
rz(-1.9645683) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.70456409) q[2];
sx q[2];
rz(-1.8578055) q[2];
sx q[2];
rz(1.9412083) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9040363) q[1];
sx q[1];
rz(-0.95103969) q[1];
sx q[1];
rz(2.0380286) q[1];
rz(-0.023530258) q[3];
sx q[3];
rz(-1.9100034) q[3];
sx q[3];
rz(0.48285218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.7548783) q[0];
sx q[0];
rz(-1.2972378) q[0];
sx q[0];
rz(-0.3219147) q[0];
rz(-1.5362668) q[1];
sx q[1];
rz(-1.9202817) q[1];
sx q[1];
rz(-2.4386491) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2625933) q[0];
sx q[0];
rz(-1.9290553) q[0];
sx q[0];
rz(0.4155638) q[0];
rz(-1.0614971) q[2];
sx q[2];
rz(-1.5361538) q[2];
sx q[2];
rz(-2.408037) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1785537) q[1];
sx q[1];
rz(-2.5759765) q[1];
sx q[1];
rz(1.2600684) q[1];
rz(-pi) q[2];
x q[2];
rz(2.8076257) q[3];
sx q[3];
rz(-2.1595862) q[3];
sx q[3];
rz(2.2759618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.2408509) q[2];
sx q[2];
rz(-0.27975953) q[2];
sx q[2];
rz(-1.8019603) q[2];
rz(2.83589) q[3];
sx q[3];
rz(-1.8140847) q[3];
sx q[3];
rz(1.3302749) q[3];
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
rz(-0.16383485) q[0];
sx q[0];
rz(-2.3840388) q[0];
sx q[0];
rz(1.2257858) q[0];
rz(2.2380791) q[1];
sx q[1];
rz(-0.61360306) q[1];
sx q[1];
rz(2.6729565) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7077431) q[0];
sx q[0];
rz(-1.9540457) q[0];
sx q[0];
rz(1.3462523) q[0];
rz(-2.2655728) q[2];
sx q[2];
rz(-2.6555736) q[2];
sx q[2];
rz(-1.3181869) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.48060265) q[1];
sx q[1];
rz(-2.8669679) q[1];
sx q[1];
rz(-1.3165228) q[1];
rz(2.0300794) q[3];
sx q[3];
rz(-1.1772403) q[3];
sx q[3];
rz(-2.4126088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.6580711) q[2];
sx q[2];
rz(-1.8489685) q[2];
sx q[2];
rz(1.1432077) q[2];
rz(-3.0269567) q[3];
sx q[3];
rz(-2.1879523) q[3];
sx q[3];
rz(1.6121929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
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
rz(-0.86482277) q[2];
sx q[2];
rz(-2.13158) q[2];
sx q[2];
rz(-2.0958015) q[2];
rz(-1.6735531) q[3];
sx q[3];
rz(-2.4874874) q[3];
sx q[3];
rz(-2.7174674) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
