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
rz(-1.4544571) q[0];
sx q[0];
rz(-1.9782826) q[0];
sx q[0];
rz(-1.8736725) q[0];
rz(1.6549702) q[1];
sx q[1];
rz(-1.854719) q[1];
sx q[1];
rz(0.16965228) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9244512) q[0];
sx q[0];
rz(-2.0471153) q[0];
sx q[0];
rz(-2.7360271) q[0];
rz(-pi) q[1];
rz(-1.0526161) q[2];
sx q[2];
rz(-1.3498326) q[2];
sx q[2];
rz(0.25928869) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.84089609) q[1];
sx q[1];
rz(-1.5920326) q[1];
sx q[1];
rz(1.1752711) q[1];
rz(-pi) q[2];
rz(0.21244399) q[3];
sx q[3];
rz(-1.8768969) q[3];
sx q[3];
rz(-2.3116632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9156645) q[2];
sx q[2];
rz(-1.6253086) q[2];
sx q[2];
rz(0.020326745) q[2];
rz(-0.23975553) q[3];
sx q[3];
rz(-0.25529796) q[3];
sx q[3];
rz(2.3026626) q[3];
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
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13424419) q[0];
sx q[0];
rz(-2.5237995) q[0];
sx q[0];
rz(-2.0452621) q[0];
rz(1.6657375) q[1];
sx q[1];
rz(-1.9225537) q[1];
sx q[1];
rz(2.347167) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5223094) q[0];
sx q[0];
rz(-1.8849775) q[0];
sx q[0];
rz(2.7854325) q[0];
rz(0.75411039) q[2];
sx q[2];
rz(-1.1271141) q[2];
sx q[2];
rz(0.59690969) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.8960171) q[1];
sx q[1];
rz(-1.4652677) q[1];
sx q[1];
rz(0.42550931) q[1];
rz(2.3268324) q[3];
sx q[3];
rz(-1.7069121) q[3];
sx q[3];
rz(-2.5825115) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.085792556) q[2];
sx q[2];
rz(-1.5726568) q[2];
sx q[2];
rz(-0.22573486) q[2];
rz(-0.24383946) q[3];
sx q[3];
rz(-2.1832681) q[3];
sx q[3];
rz(-0.17914151) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8826411) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(-2.7146345) q[0];
rz(0.34307617) q[1];
sx q[1];
rz(-1.1767358) q[1];
sx q[1];
rz(0.25161904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93676567) q[0];
sx q[0];
rz(-0.63049769) q[0];
sx q[0];
rz(2.7555694) q[0];
x q[1];
rz(-1.8015091) q[2];
sx q[2];
rz(-0.78015155) q[2];
sx q[2];
rz(0.095965775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.7921264) q[1];
sx q[1];
rz(-1.9053962) q[1];
sx q[1];
rz(-1.3931331) q[1];
rz(-pi) q[2];
rz(2.2813024) q[3];
sx q[3];
rz(-1.5328457) q[3];
sx q[3];
rz(0.52805985) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.0106657) q[2];
sx q[2];
rz(-1.6341354) q[2];
sx q[2];
rz(2.6962386) q[2];
rz(1.2238067) q[3];
sx q[3];
rz(-1.8755951) q[3];
sx q[3];
rz(-2.7538917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7447516) q[0];
sx q[0];
rz(-2.3362609) q[0];
sx q[0];
rz(-1.1591563) q[0];
rz(-1.2027488) q[1];
sx q[1];
rz(-1.9614599) q[1];
sx q[1];
rz(2.4045827) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.096367717) q[0];
sx q[0];
rz(-2.8857433) q[0];
sx q[0];
rz(-2.8907828) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1221755) q[2];
sx q[2];
rz(-2.4146842) q[2];
sx q[2];
rz(-2.2688933) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1070246) q[1];
sx q[1];
rz(-0.66063297) q[1];
sx q[1];
rz(2.8559277) q[1];
rz(-0.33212157) q[3];
sx q[3];
rz(-2.9654223) q[3];
sx q[3];
rz(-0.15819269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0143339) q[2];
sx q[2];
rz(-0.88982439) q[2];
sx q[2];
rz(-0.8117525) q[2];
rz(-0.74770606) q[3];
sx q[3];
rz(-2.2816608) q[3];
sx q[3];
rz(-3.0809793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49374813) q[0];
sx q[0];
rz(-2.0525377) q[0];
sx q[0];
rz(1.7895948) q[0];
rz(0.79212517) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(1.0101213) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3528004) q[0];
sx q[0];
rz(-1.236434) q[0];
sx q[0];
rz(-1.0092606) q[0];
rz(-pi) q[1];
rz(1.5158122) q[2];
sx q[2];
rz(-1.1786945) q[2];
sx q[2];
rz(-2.7934809) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0466442) q[1];
sx q[1];
rz(-1.6502336) q[1];
sx q[1];
rz(-0.55752505) q[1];
x q[2];
rz(1.406448) q[3];
sx q[3];
rz(-1.6858628) q[3];
sx q[3];
rz(-0.0054520741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1048364) q[2];
sx q[2];
rz(-0.75883055) q[2];
sx q[2];
rz(1.7581615) q[2];
rz(3.1116327) q[3];
sx q[3];
rz(-2.1432917) q[3];
sx q[3];
rz(1.8647319) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5112011) q[0];
sx q[0];
rz(-2.7250405) q[0];
sx q[0];
rz(0.1314441) q[0];
rz(-0.99880544) q[1];
sx q[1];
rz(-1.1591594) q[1];
sx q[1];
rz(-1.3712032) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11267647) q[0];
sx q[0];
rz(-1.925591) q[0];
sx q[0];
rz(-1.2043857) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.68626257) q[2];
sx q[2];
rz(-1.7289844) q[2];
sx q[2];
rz(-0.37294086) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.19325312) q[1];
sx q[1];
rz(-2.8242556) q[1];
sx q[1];
rz(-2.0977519) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5672979) q[3];
sx q[3];
rz(-2.281684) q[3];
sx q[3];
rz(1.9128127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3580631) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(-1.0895458) q[2];
rz(-2.774636) q[3];
sx q[3];
rz(-0.4042545) q[3];
sx q[3];
rz(-2.3914316) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.628196) q[0];
sx q[0];
rz(-0.80444002) q[0];
sx q[0];
rz(-2.2688769) q[0];
rz(1.5645507) q[1];
sx q[1];
rz(-1.6460452) q[1];
sx q[1];
rz(-0.73211342) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66595378) q[0];
sx q[0];
rz(-1.7717965) q[0];
sx q[0];
rz(0.71626407) q[0];
rz(2.9896834) q[2];
sx q[2];
rz(-0.6912125) q[2];
sx q[2];
rz(-0.37746261) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.435675) q[1];
sx q[1];
rz(-2.6443548) q[1];
sx q[1];
rz(-1.9051208) q[1];
x q[2];
rz(-3.0454265) q[3];
sx q[3];
rz(-1.796826) q[3];
sx q[3];
rz(1.9440252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.5319891) q[2];
sx q[2];
rz(-1.0010109) q[2];
sx q[2];
rz(-0.58094376) q[2];
rz(-1.1254492) q[3];
sx q[3];
rz(-1.1939129) q[3];
sx q[3];
rz(1.1552936) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693212) q[0];
sx q[0];
rz(-3.0348365) q[0];
sx q[0];
rz(-2.971055) q[0];
rz(1.2455617) q[1];
sx q[1];
rz(-1.5252557) q[1];
sx q[1];
rz(2.4571498) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67407284) q[0];
sx q[0];
rz(-2.0430123) q[0];
sx q[0];
rz(-1.945482) q[0];
rz(2.9028394) q[2];
sx q[2];
rz(-0.49477067) q[2];
sx q[2];
rz(2.787311) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.0741006) q[1];
sx q[1];
rz(-1.3252186) q[1];
sx q[1];
rz(2.8771993) q[1];
rz(-1.5106279) q[3];
sx q[3];
rz(-0.92272131) q[3];
sx q[3];
rz(-1.5980521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.30274063) q[2];
sx q[2];
rz(-1.799823) q[2];
sx q[2];
rz(2.0212685) q[2];
rz(-2.4773347) q[3];
sx q[3];
rz(-2.7393326) q[3];
sx q[3];
rz(2.6334488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6054194) q[0];
sx q[0];
rz(-2.7689731) q[0];
sx q[0];
rz(-0.63825178) q[0];
rz(-0.0087180184) q[1];
sx q[1];
rz(-2.0254841) q[1];
sx q[1];
rz(-2.7353824) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4753301) q[0];
sx q[0];
rz(-1.8436369) q[0];
sx q[0];
rz(3.0240339) q[0];
rz(1.137144) q[2];
sx q[2];
rz(-1.523345) q[2];
sx q[2];
rz(0.46772568) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.49055028) q[1];
sx q[1];
rz(-0.79572751) q[1];
sx q[1];
rz(2.9814238) q[1];
rz(-pi) q[2];
rz(3.0681455) q[3];
sx q[3];
rz(-1.5961507) q[3];
sx q[3];
rz(-0.0097833477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.093988769) q[2];
sx q[2];
rz(-2.0145907) q[2];
sx q[2];
rz(2.7809533) q[2];
rz(0.7274729) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(2.8235161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3280846) q[0];
sx q[0];
rz(-0.94539517) q[0];
sx q[0];
rz(-0.16950053) q[0];
rz(2.4318579) q[1];
sx q[1];
rz(-1.4163481) q[1];
sx q[1];
rz(-1.08606) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57993488) q[0];
sx q[0];
rz(-2.5626963) q[0];
sx q[0];
rz(0.41335841) q[0];
rz(-pi) q[1];
rz(-1.1342088) q[2];
sx q[2];
rz(-1.2870169) q[2];
sx q[2];
rz(1.6095609) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.12609657) q[1];
sx q[1];
rz(-2.129678) q[1];
sx q[1];
rz(-3.0256767) q[1];
x q[2];
rz(1.754722) q[3];
sx q[3];
rz(-1.1952595) q[3];
sx q[3];
rz(0.25115764) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78709948) q[2];
sx q[2];
rz(-0.084986173) q[2];
sx q[2];
rz(0.3698012) q[2];
rz(1.1981111) q[3];
sx q[3];
rz(-1.5503649) q[3];
sx q[3];
rz(-0.91355356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13854606) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(-2.1239602) q[1];
sx q[1];
rz(-1.3980649) q[1];
sx q[1];
rz(-1.3921888) q[1];
rz(1.5003149) q[2];
sx q[2];
rz(-1.4983836) q[2];
sx q[2];
rz(-0.56624779) q[2];
rz(2.9774278) q[3];
sx q[3];
rz(-1.6995242) q[3];
sx q[3];
rz(-1.8229998) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
