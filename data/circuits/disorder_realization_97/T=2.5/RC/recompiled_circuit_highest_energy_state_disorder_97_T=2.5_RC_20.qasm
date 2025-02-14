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
rz(2.4801369) q[0];
sx q[0];
rz(-2.8602726) q[0];
sx q[0];
rz(1.9982279) q[0];
rz(0.65302628) q[1];
sx q[1];
rz(4.9622494) q[1];
sx q[1];
rz(9.4402037) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30886762) q[0];
sx q[0];
rz(-2.057353) q[0];
sx q[0];
rz(-0.95054807) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1434606) q[2];
sx q[2];
rz(-2.0734416) q[2];
sx q[2];
rz(-2.5302327) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.28492883) q[1];
sx q[1];
rz(-2.2010815) q[1];
sx q[1];
rz(-0.50215118) q[1];
x q[2];
rz(2.8174345) q[3];
sx q[3];
rz(-2.5796842) q[3];
sx q[3];
rz(-0.066854157) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.36110863) q[2];
sx q[2];
rz(-3.1321654) q[2];
sx q[2];
rz(-2.0375605) q[2];
rz(2.5822254) q[3];
sx q[3];
rz(-2.1909824) q[3];
sx q[3];
rz(2.2573788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1169432) q[0];
sx q[0];
rz(-0.58498061) q[0];
sx q[0];
rz(0.90644932) q[0];
rz(2.8001884) q[1];
sx q[1];
rz(-0.87541348) q[1];
sx q[1];
rz(-0.34814775) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7067862) q[0];
sx q[0];
rz(-0.40074391) q[0];
sx q[0];
rz(-2.4639919) q[0];
rz(-pi) q[1];
rz(-2.1402713) q[2];
sx q[2];
rz(-3.0932326) q[2];
sx q[2];
rz(-1.1011626) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22869884) q[1];
sx q[1];
rz(-1.6457108) q[1];
sx q[1];
rz(0.95864899) q[1];
rz(-2.313638) q[3];
sx q[3];
rz(-1.7418234) q[3];
sx q[3];
rz(-3.1307182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1657095) q[2];
sx q[2];
rz(-2.7996863) q[2];
sx q[2];
rz(0.085414097) q[2];
rz(-0.092770569) q[3];
sx q[3];
rz(-2.2723891) q[3];
sx q[3];
rz(-2.2658277) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3669325) q[0];
sx q[0];
rz(-0.36151883) q[0];
sx q[0];
rz(0.34459484) q[0];
rz(-1.5917646) q[1];
sx q[1];
rz(-1.7235618) q[1];
sx q[1];
rz(-1.4580911) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7994996) q[0];
sx q[0];
rz(-1.5342426) q[0];
sx q[0];
rz(0.7176368) q[0];
rz(-1.6750653) q[2];
sx q[2];
rz(-2.6223256) q[2];
sx q[2];
rz(0.41246513) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.96839206) q[1];
sx q[1];
rz(-1.8138944) q[1];
sx q[1];
rz(0.17420423) q[1];
rz(-pi) q[2];
rz(-0.86784466) q[3];
sx q[3];
rz(-0.72753564) q[3];
sx q[3];
rz(-0.65505469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6452667) q[2];
sx q[2];
rz(-0.93924773) q[2];
sx q[2];
rz(-1.1305031) q[2];
rz(2.1988403) q[3];
sx q[3];
rz(-2.0944984) q[3];
sx q[3];
rz(1.9345136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4261674) q[0];
sx q[0];
rz(-1.3454477) q[0];
sx q[0];
rz(-1.4581534) q[0];
rz(1.0484877) q[1];
sx q[1];
rz(-1.4920934) q[1];
sx q[1];
rz(2.6715211) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94435197) q[0];
sx q[0];
rz(-1.2908123) q[0];
sx q[0];
rz(-0.51843317) q[0];
rz(-pi) q[1];
rz(-2.8747005) q[2];
sx q[2];
rz(-2.0837657) q[2];
sx q[2];
rz(1.5266872) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.64420358) q[1];
sx q[1];
rz(-2.4371689) q[1];
sx q[1];
rz(0.299244) q[1];
rz(-pi) q[2];
rz(0.87155452) q[3];
sx q[3];
rz(-2.1208753) q[3];
sx q[3];
rz(-3.0015903) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5992392) q[2];
sx q[2];
rz(-2.3848644) q[2];
sx q[2];
rz(-2.2139464) q[2];
rz(-1.8574235) q[3];
sx q[3];
rz(-1.410306) q[3];
sx q[3];
rz(-0.94902432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.014378431) q[0];
sx q[0];
rz(-0.40507409) q[0];
sx q[0];
rz(2.6601484) q[0];
rz(2.1638347) q[1];
sx q[1];
rz(-0.6441741) q[1];
sx q[1];
rz(1.0747304) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18431379) q[0];
sx q[0];
rz(-1.9879576) q[0];
sx q[0];
rz(2.8431312) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8406732) q[2];
sx q[2];
rz(-1.3011429) q[2];
sx q[2];
rz(-0.62992879) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.9323651) q[1];
sx q[1];
rz(-1.8098988) q[1];
sx q[1];
rz(-1.1878315) q[1];
rz(-pi) q[2];
rz(-2.7513192) q[3];
sx q[3];
rz(-1.2765795) q[3];
sx q[3];
rz(2.7613784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0398728) q[2];
sx q[2];
rz(-2.0247255) q[2];
sx q[2];
rz(-2.5679585) q[2];
rz(1.6576069) q[3];
sx q[3];
rz(-2.0935757) q[3];
sx q[3];
rz(-0.60605961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(0.023733519) q[0];
sx q[0];
rz(-2.8890299) q[0];
sx q[0];
rz(-2.8299487) q[0];
rz(2.812884) q[1];
sx q[1];
rz(-1.3361822) q[1];
sx q[1];
rz(-2.4005344) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8097654) q[0];
sx q[0];
rz(-0.96967319) q[0];
sx q[0];
rz(-2.580216) q[0];
rz(-3.0799505) q[2];
sx q[2];
rz(-1.7992524) q[2];
sx q[2];
rz(1.018334) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.44081008) q[1];
sx q[1];
rz(-0.68432552) q[1];
sx q[1];
rz(1.6709063) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5696008) q[3];
sx q[3];
rz(-2.5107267) q[3];
sx q[3];
rz(1.6945171) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3662423) q[2];
sx q[2];
rz(-0.58997184) q[2];
sx q[2];
rz(0.74014202) q[2];
rz(-0.38672334) q[3];
sx q[3];
rz(-2.4155278) q[3];
sx q[3];
rz(2.6630785) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1665523) q[0];
sx q[0];
rz(-0.98169011) q[0];
sx q[0];
rz(0.24328406) q[0];
rz(0.89726204) q[1];
sx q[1];
rz(-1.5829395) q[1];
sx q[1];
rz(-0.74403393) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76975856) q[0];
sx q[0];
rz(-0.39880619) q[0];
sx q[0];
rz(-0.21488667) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2076188) q[2];
sx q[2];
rz(-1.2609856) q[2];
sx q[2];
rz(1.8997836) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0571088) q[1];
sx q[1];
rz(-1.2027825) q[1];
sx q[1];
rz(1.012136) q[1];
rz(-2.0970515) q[3];
sx q[3];
rz(-1.730837) q[3];
sx q[3];
rz(2.2621148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9545142) q[2];
sx q[2];
rz(-2.3176471) q[2];
sx q[2];
rz(-0.0079060923) q[2];
rz(1.4963957) q[3];
sx q[3];
rz(-0.073286101) q[3];
sx q[3];
rz(-2.6325398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(0.024648333) q[0];
sx q[0];
rz(-0.032289676) q[0];
sx q[0];
rz(-0.13667983) q[0];
rz(-1.3487123) q[1];
sx q[1];
rz(-1.2929448) q[1];
sx q[1];
rz(-0.50833702) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1404554) q[0];
sx q[0];
rz(-1.9656202) q[0];
sx q[0];
rz(-0.14060256) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2042785) q[2];
sx q[2];
rz(-2.4319639) q[2];
sx q[2];
rz(-0.21096551) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9421025) q[1];
sx q[1];
rz(-2.2245896) q[1];
sx q[1];
rz(2.687017) q[1];
rz(-pi) q[2];
rz(-1.4517205) q[3];
sx q[3];
rz(-1.4919466) q[3];
sx q[3];
rz(0.54311968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.4623744) q[2];
sx q[2];
rz(-2.4269673) q[2];
sx q[2];
rz(-3.083631) q[2];
rz(-0.95311779) q[3];
sx q[3];
rz(-1.0473731) q[3];
sx q[3];
rz(-0.63547772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6111095) q[0];
sx q[0];
rz(-1.7246752) q[0];
sx q[0];
rz(-0.86607754) q[0];
rz(-1.3653612) q[1];
sx q[1];
rz(-1.5581286) q[1];
sx q[1];
rz(-0.80397111) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3063076) q[0];
sx q[0];
rz(-0.87125766) q[0];
sx q[0];
rz(2.7650096) q[0];
rz(-pi) q[1];
rz(-0.78900225) q[2];
sx q[2];
rz(-1.2478634) q[2];
sx q[2];
rz(-0.6168405) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6996346) q[1];
sx q[1];
rz(-0.82773877) q[1];
sx q[1];
rz(-2.1585967) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6220408) q[3];
sx q[3];
rz(-2.7375023) q[3];
sx q[3];
rz(-1.3431637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.08635252) q[2];
sx q[2];
rz(-0.1487727) q[2];
sx q[2];
rz(-1.0521592) q[2];
rz(0.85938984) q[3];
sx q[3];
rz(-0.79901564) q[3];
sx q[3];
rz(2.1990282) q[3];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4683485) q[0];
sx q[0];
rz(-2.4788661) q[0];
sx q[0];
rz(0.41917875) q[0];
rz(-3.1255417) q[1];
sx q[1];
rz(-1.586986) q[1];
sx q[1];
rz(3.0174461) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57174411) q[0];
sx q[0];
rz(-2.1198065) q[0];
sx q[0];
rz(-0.17904539) q[0];
x q[1];
rz(-1.3163465) q[2];
sx q[2];
rz(-1.7702603) q[2];
sx q[2];
rz(-2.2274889) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.1325741) q[1];
sx q[1];
rz(-0.30840519) q[1];
sx q[1];
rz(-2.5243702) q[1];
rz(-pi) q[2];
rz(-2.9591363) q[3];
sx q[3];
rz(-0.92853755) q[3];
sx q[3];
rz(0.69891847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9496256) q[2];
sx q[2];
rz(-2.4032335) q[2];
sx q[2];
rz(2.9730566) q[2];
rz(-1.4847697) q[3];
sx q[3];
rz(-2.6717581) q[3];
sx q[3];
rz(0.43009871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-2.1930502) q[0];
sx q[0];
rz(-0.47946231) q[0];
sx q[0];
rz(-0.23647501) q[0];
rz(2.3169658) q[1];
sx q[1];
rz(-1.5390479) q[1];
sx q[1];
rz(1.8631757) q[1];
rz(-2.3659351) q[2];
sx q[2];
rz(-1.1393329) q[2];
sx q[2];
rz(2.3928497) q[2];
rz(1.0239368) q[3];
sx q[3];
rz(-1.6228068) q[3];
sx q[3];
rz(0.3323298) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
