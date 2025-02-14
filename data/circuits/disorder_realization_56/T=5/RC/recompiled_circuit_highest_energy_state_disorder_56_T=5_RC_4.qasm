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
rz(-1.6662958) q[0];
sx q[0];
rz(-1.8721606) q[0];
sx q[0];
rz(2.4694634) q[0];
rz(0.44828662) q[1];
sx q[1];
rz(-1.5223794) q[1];
sx q[1];
rz(-2.3840005) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.282895) q[0];
sx q[0];
rz(-0.13326779) q[0];
sx q[0];
rz(-0.96576502) q[0];
x q[1];
rz(0.50067164) q[2];
sx q[2];
rz(-1.8953875) q[2];
sx q[2];
rz(2.2364834) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.6554101) q[1];
sx q[1];
rz(-1.3453801) q[1];
sx q[1];
rz(0.42014007) q[1];
x q[2];
rz(-0.67528649) q[3];
sx q[3];
rz(-2.3216341) q[3];
sx q[3];
rz(-2.2238071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0932833) q[2];
sx q[2];
rz(-1.8123241) q[2];
sx q[2];
rz(-1.3422356) q[2];
rz(-1.3679158) q[3];
sx q[3];
rz(-1.0932837) q[3];
sx q[3];
rz(-1.5889408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20464483) q[0];
sx q[0];
rz(-2.1277675) q[0];
sx q[0];
rz(2.192705) q[0];
rz(-0.10143796) q[1];
sx q[1];
rz(-1.0600435) q[1];
sx q[1];
rz(-0.92996517) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5407911) q[0];
sx q[0];
rz(-1.7834429) q[0];
sx q[0];
rz(1.4669111) q[0];
x q[1];
rz(-1.5783159) q[2];
sx q[2];
rz(-2.0759656) q[2];
sx q[2];
rz(1.9746321) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.21446642) q[1];
sx q[1];
rz(-8*pi/15) q[1];
sx q[1];
rz(-2.9795477) q[1];
rz(-pi) q[2];
x q[2];
rz(0.92949246) q[3];
sx q[3];
rz(-2.2713823) q[3];
sx q[3];
rz(-1.495958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0019504) q[2];
sx q[2];
rz(-1.4652239) q[2];
sx q[2];
rz(-2.3993717) q[2];
rz(0.46418515) q[3];
sx q[3];
rz(-1.7964541) q[3];
sx q[3];
rz(1.5884885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4699698) q[0];
sx q[0];
rz(-0.18593423) q[0];
sx q[0];
rz(-1.6999014) q[0];
rz(0.16547671) q[1];
sx q[1];
rz(-2.4022357) q[1];
sx q[1];
rz(0.86732078) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.353426) q[0];
sx q[0];
rz(-1.5702308) q[0];
sx q[0];
rz(4.0325469e-05) q[0];
rz(-pi) q[1];
rz(-1.1050842) q[2];
sx q[2];
rz(-0.55464472) q[2];
sx q[2];
rz(1.1797734) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.2056624) q[1];
sx q[1];
rz(-1.8962269) q[1];
sx q[1];
rz(1.2742576) q[1];
rz(0.54384772) q[3];
sx q[3];
rz(-2.2529054) q[3];
sx q[3];
rz(-0.92095514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1778339) q[2];
sx q[2];
rz(-1.9671665) q[2];
sx q[2];
rz(0.7589232) q[2];
rz(-2.3376076) q[3];
sx q[3];
rz(-0.94954973) q[3];
sx q[3];
rz(-2.4132531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3114965) q[0];
sx q[0];
rz(-1.281597) q[0];
sx q[0];
rz(-0.48402825) q[0];
rz(-1.6054035) q[1];
sx q[1];
rz(-1.4151298) q[1];
sx q[1];
rz(-1.7162292) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40544505) q[0];
sx q[0];
rz(-2.510294) q[0];
sx q[0];
rz(0.59595705) q[0];
rz(-pi) q[1];
rz(-0.62230939) q[2];
sx q[2];
rz(-2.7633585) q[2];
sx q[2];
rz(0.44307274) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.64827907) q[1];
sx q[1];
rz(-1.7910622) q[1];
sx q[1];
rz(2.2146398) q[1];
x q[2];
rz(-0.68257733) q[3];
sx q[3];
rz(-0.55844504) q[3];
sx q[3];
rz(2.6598833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.385685) q[2];
sx q[2];
rz(-2.0734831) q[2];
sx q[2];
rz(-2.8483086) q[2];
rz(-3.0934603) q[3];
sx q[3];
rz(-1.3061413) q[3];
sx q[3];
rz(2.8369246) q[3];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83679477) q[0];
sx q[0];
rz(-1.4076819) q[0];
sx q[0];
rz(2.0378713) q[0];
rz(-2.2301105) q[1];
sx q[1];
rz(-1.6171004) q[1];
sx q[1];
rz(-0.10890659) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0287339) q[0];
sx q[0];
rz(-0.92623952) q[0];
sx q[0];
rz(0.95405202) q[0];
rz(0.2997024) q[2];
sx q[2];
rz(-2.5863918) q[2];
sx q[2];
rz(-0.19367684) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8629294) q[1];
sx q[1];
rz(-2.4554774) q[1];
sx q[1];
rz(1.8300232) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.44556983) q[3];
sx q[3];
rz(-2.516149) q[3];
sx q[3];
rz(1.3423962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.64985046) q[2];
sx q[2];
rz(-1.7848585) q[2];
sx q[2];
rz(-0.25807992) q[2];
rz(-2.7598925) q[3];
sx q[3];
rz(-2.2038867) q[3];
sx q[3];
rz(0.55735731) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7447164) q[0];
sx q[0];
rz(-0.98137403) q[0];
sx q[0];
rz(1.1599524) q[0];
rz(-2.5634735) q[1];
sx q[1];
rz(-1.6696397) q[1];
sx q[1];
rz(0.32346183) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.9834546) q[0];
sx q[0];
rz(-1.1295415) q[0];
sx q[0];
rz(2.6652314) q[0];
x q[1];
rz(-2.696229) q[2];
sx q[2];
rz(-0.54071835) q[2];
sx q[2];
rz(0.13810829) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0929133) q[1];
sx q[1];
rz(-1.3447176) q[1];
sx q[1];
rz(2.5552555) q[1];
rz(-pi) q[2];
rz(-1.8927073) q[3];
sx q[3];
rz(-1.1309147) q[3];
sx q[3];
rz(1.3706576) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.20299992) q[2];
sx q[2];
rz(-0.77363571) q[2];
sx q[2];
rz(-2.2933551) q[2];
rz(1.8741459) q[3];
sx q[3];
rz(-1.7402612) q[3];
sx q[3];
rz(-0.69376865) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11442014) q[0];
sx q[0];
rz(-2.4878451) q[0];
sx q[0];
rz(-2.2681336) q[0];
rz(-1.9050441) q[1];
sx q[1];
rz(-2.1658587) q[1];
sx q[1];
rz(0.083560856) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.386451) q[0];
sx q[0];
rz(-1.4487195) q[0];
sx q[0];
rz(0.32213078) q[0];
rz(-pi) q[1];
rz(2.100285) q[2];
sx q[2];
rz(-1.7779967) q[2];
sx q[2];
rz(2.4251314) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.096683) q[1];
sx q[1];
rz(-2.2538141) q[1];
sx q[1];
rz(-2.6234496) q[1];
rz(1.2681581) q[3];
sx q[3];
rz(-1.2250966) q[3];
sx q[3];
rz(-1.3132335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.43625912) q[2];
sx q[2];
rz(-1.5353563) q[2];
sx q[2];
rz(-1.7748888) q[2];
rz(0.27240917) q[3];
sx q[3];
rz(-2.1209769) q[3];
sx q[3];
rz(0.75850707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2347539) q[0];
sx q[0];
rz(-2.315157) q[0];
sx q[0];
rz(-0.93836623) q[0];
rz(0.80288184) q[1];
sx q[1];
rz(-0.96790853) q[1];
sx q[1];
rz(0.82069194) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.11242) q[0];
sx q[0];
rz(-0.39390644) q[0];
sx q[0];
rz(1.9256511) q[0];
rz(-pi) q[1];
rz(0.16751473) q[2];
sx q[2];
rz(-1.2021087) q[2];
sx q[2];
rz(-0.14376727) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-3.0487602) q[1];
sx q[1];
rz(-2.5093749) q[1];
sx q[1];
rz(1.7507751) q[1];
rz(0.16444178) q[3];
sx q[3];
rz(-1.9666082) q[3];
sx q[3];
rz(0.37104097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9390823) q[2];
sx q[2];
rz(-2.9896917) q[2];
sx q[2];
rz(-2.0738156) q[2];
rz(-1.1311401) q[3];
sx q[3];
rz(-2.307297) q[3];
sx q[3];
rz(0.080032674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.15602569) q[0];
sx q[0];
rz(-1.1556867) q[0];
sx q[0];
rz(2.735403) q[0];
rz(1.4093026) q[1];
sx q[1];
rz(-1.0180165) q[1];
sx q[1];
rz(-2.0565775) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70124088) q[0];
sx q[0];
rz(-0.95537649) q[0];
sx q[0];
rz(1.3568272) q[0];
rz(-pi) q[1];
rz(-0.0064956587) q[2];
sx q[2];
rz(-0.61031872) q[2];
sx q[2];
rz(1.4607061) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0250562) q[1];
sx q[1];
rz(-2.3485314) q[1];
sx q[1];
rz(2.039394) q[1];
rz(-pi) q[2];
x q[2];
rz(0.57906753) q[3];
sx q[3];
rz(-2.0430312) q[3];
sx q[3];
rz(3.0621665) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.5558418) q[2];
sx q[2];
rz(-0.61183524) q[2];
sx q[2];
rz(2.2108868) q[2];
rz(-1.3961004) q[3];
sx q[3];
rz(-1.3627108) q[3];
sx q[3];
rz(3.0220368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4569106) q[0];
sx q[0];
rz(-0.97701183) q[0];
sx q[0];
rz(-1.3071625) q[0];
rz(-0.3859418) q[1];
sx q[1];
rz(-1.394505) q[1];
sx q[1];
rz(-1.0345667) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.912303) q[0];
sx q[0];
rz(-1.467839) q[0];
sx q[0];
rz(1.3961755) q[0];
x q[1];
rz(-2.7427086) q[2];
sx q[2];
rz(-0.49408696) q[2];
sx q[2];
rz(-1.8979771) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1423827) q[1];
sx q[1];
rz(-2.2015328) q[1];
sx q[1];
rz(-2.5725065) q[1];
rz(-pi) q[2];
rz(-1.7004556) q[3];
sx q[3];
rz(-0.76414062) q[3];
sx q[3];
rz(2.0779) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4474386) q[2];
sx q[2];
rz(-0.73221451) q[2];
sx q[2];
rz(0.36699692) q[2];
rz(-1.1191818) q[3];
sx q[3];
rz(-2.7435591) q[3];
sx q[3];
rz(-1.6163577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601396) q[0];
sx q[0];
rz(-1.9522788) q[0];
sx q[0];
rz(2.0576394) q[0];
rz(-3.0837334) q[1];
sx q[1];
rz(-1.2670988) q[1];
sx q[1];
rz(1.7115464) q[1];
rz(-2.62769) q[2];
sx q[2];
rz(-2.9436417) q[2];
sx q[2];
rz(-1.8463989) q[2];
rz(1.9290942) q[3];
sx q[3];
rz(-0.90322106) q[3];
sx q[3];
rz(-1.7424104) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
