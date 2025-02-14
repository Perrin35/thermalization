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
rz(1.6871356) q[0];
sx q[0];
rz(1.9782826) q[0];
sx q[0];
rz(10.692698) q[0];
rz(-1.4866225) q[1];
sx q[1];
rz(-1.2868737) q[1];
sx q[1];
rz(2.9719404) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.171863) q[0];
sx q[0];
rz(-0.61530441) q[0];
sx q[0];
rz(-2.2236373) q[0];
rz(-pi) q[1];
rz(0.25303475) q[2];
sx q[2];
rz(-1.0664244) q[2];
sx q[2];
rz(1.4358226) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.4625068) q[1];
sx q[1];
rz(-2.7455277) q[1];
sx q[1];
rz(-1.6258662) q[1];
rz(-pi) q[2];
rz(0.21244399) q[3];
sx q[3];
rz(-1.8768969) q[3];
sx q[3];
rz(0.82992947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.22592813) q[2];
sx q[2];
rz(-1.516284) q[2];
sx q[2];
rz(0.020326745) q[2];
rz(2.9018371) q[3];
sx q[3];
rz(-2.8862947) q[3];
sx q[3];
rz(0.83893004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13424419) q[0];
sx q[0];
rz(-2.5237995) q[0];
sx q[0];
rz(1.0963305) q[0];
rz(-1.6657375) q[1];
sx q[1];
rz(-1.9225537) q[1];
sx q[1];
rz(0.79442564) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83704575) q[0];
sx q[0];
rz(-1.908806) q[0];
sx q[0];
rz(-1.90453) q[0];
x q[1];
rz(-2.1486304) q[2];
sx q[2];
rz(-2.2374399) q[2];
sx q[2];
rz(1.3571598) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8960171) q[1];
sx q[1];
rz(-1.4652677) q[1];
sx q[1];
rz(-0.42550931) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3268324) q[3];
sx q[3];
rz(-1.4346806) q[3];
sx q[3];
rz(-0.55908113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0558001) q[2];
sx q[2];
rz(-1.5726568) q[2];
sx q[2];
rz(-2.9158578) q[2];
rz(0.24383946) q[3];
sx q[3];
rz(-2.1832681) q[3];
sx q[3];
rz(-2.9624511) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2589515) q[0];
sx q[0];
rz(-1.3132341) q[0];
sx q[0];
rz(2.7146345) q[0];
rz(0.34307617) q[1];
sx q[1];
rz(-1.9648569) q[1];
sx q[1];
rz(-0.25161904) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7386603) q[0];
sx q[0];
rz(-2.1485747) q[0];
sx q[0];
rz(1.3026122) q[0];
rz(-0.22253665) q[2];
sx q[2];
rz(-2.3250569) q[2];
sx q[2];
rz(0.41513069) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.15089825) q[1];
sx q[1];
rz(-0.37726918) q[1];
sx q[1];
rz(-2.6713085) q[1];
x q[2];
rz(-0.86029025) q[3];
sx q[3];
rz(-1.5328457) q[3];
sx q[3];
rz(-2.6135328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.0106657) q[2];
sx q[2];
rz(-1.5074573) q[2];
sx q[2];
rz(2.6962386) q[2];
rz(1.2238067) q[3];
sx q[3];
rz(-1.2659975) q[3];
sx q[3];
rz(-0.38770097) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39684108) q[0];
sx q[0];
rz(-2.3362609) q[0];
sx q[0];
rz(1.9824363) q[0];
rz(1.9388439) q[1];
sx q[1];
rz(-1.1801327) q[1];
sx q[1];
rz(-2.4045827) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16251462) q[0];
sx q[0];
rz(-1.3231228) q[0];
sx q[0];
rz(1.5059656) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89522712) q[2];
sx q[2];
rz(-1.2784119) q[2];
sx q[2];
rz(1.043373) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4631237) q[1];
sx q[1];
rz(-0.94128525) q[1];
sx q[1];
rz(1.786382) q[1];
rz(2.8094711) q[3];
sx q[3];
rz(-2.9654223) q[3];
sx q[3];
rz(-0.15819269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.1272588) q[2];
sx q[2];
rz(-0.88982439) q[2];
sx q[2];
rz(-2.3298402) q[2];
rz(2.3938866) q[3];
sx q[3];
rz(-0.85993189) q[3];
sx q[3];
rz(-0.060613304) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49374813) q[0];
sx q[0];
rz(-1.0890549) q[0];
sx q[0];
rz(1.3519979) q[0];
rz(2.3494675) q[1];
sx q[1];
rz(-2.242576) q[1];
sx q[1];
rz(2.1314714) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8787694) q[0];
sx q[0];
rz(-2.4973625) q[0];
sx q[0];
rz(-2.1488726) q[0];
rz(-2.7489565) q[2];
sx q[2];
rz(-1.6216039) q[2];
sx q[2];
rz(-1.9399373) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0466442) q[1];
sx q[1];
rz(-1.6502336) q[1];
sx q[1];
rz(-0.55752505) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.185818) q[3];
sx q[3];
rz(-2.9412651) q[3];
sx q[3];
rz(-0.97071394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.036756214) q[2];
sx q[2];
rz(-0.75883055) q[2];
sx q[2];
rz(1.3834312) q[2];
rz(3.1116327) q[3];
sx q[3];
rz(-0.99830097) q[3];
sx q[3];
rz(-1.8647319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6303915) q[0];
sx q[0];
rz(-0.41655219) q[0];
sx q[0];
rz(3.0101486) q[0];
rz(-0.99880544) q[1];
sx q[1];
rz(-1.9824332) q[1];
sx q[1];
rz(-1.7703895) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5509439) q[0];
sx q[0];
rz(-1.9134132) q[0];
sx q[0];
rz(-0.37775535) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4553301) q[2];
sx q[2];
rz(-1.4126083) q[2];
sx q[2];
rz(2.7686518) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.882521) q[1];
sx q[1];
rz(-1.4132199) q[1];
sx q[1];
rz(-1.8474008) q[1];
x q[2];
rz(-0.71089069) q[3];
sx q[3];
rz(-1.5734473) q[3];
sx q[3];
rz(0.34429911) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.78352952) q[2];
sx q[2];
rz(-0.72856599) q[2];
sx q[2];
rz(-1.0895458) q[2];
rz(0.36695668) q[3];
sx q[3];
rz(-2.7373382) q[3];
sx q[3];
rz(-0.75016108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51339665) q[0];
sx q[0];
rz(-0.80444002) q[0];
sx q[0];
rz(-0.87271571) q[0];
rz(1.577042) q[1];
sx q[1];
rz(-1.4955474) q[1];
sx q[1];
rz(-0.73211342) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67938618) q[0];
sx q[0];
rz(-2.4024978) q[0];
sx q[0];
rz(2.8406891) q[0];
rz(-pi) q[1];
x q[1];
rz(0.15190923) q[2];
sx q[2];
rz(-2.4503802) q[2];
sx q[2];
rz(-0.37746261) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.70591761) q[1];
sx q[1];
rz(-2.6443548) q[1];
sx q[1];
rz(1.9051208) q[1];
x q[2];
rz(0.096166178) q[3];
sx q[3];
rz(-1.796826) q[3];
sx q[3];
rz(1.9440252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.6096036) q[2];
sx q[2];
rz(-1.0010109) q[2];
sx q[2];
rz(-0.58094376) q[2];
rz(-1.1254492) q[3];
sx q[3];
rz(-1.9476798) q[3];
sx q[3];
rz(1.9862991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4693212) q[0];
sx q[0];
rz(-3.0348365) q[0];
sx q[0];
rz(0.17053764) q[0];
rz(1.2455617) q[1];
sx q[1];
rz(-1.6163369) q[1];
sx q[1];
rz(0.68444288) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.71971539) q[0];
sx q[0];
rz(-1.2388031) q[0];
sx q[0];
rz(0.50194711) q[0];
rz(-0.23875321) q[2];
sx q[2];
rz(-2.646822) q[2];
sx q[2];
rz(-2.787311) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.43757968) q[1];
sx q[1];
rz(-1.3145169) q[1];
sx q[1];
rz(1.8248454) q[1];
rz(-pi) q[2];
rz(0.6489469) q[3];
sx q[3];
rz(-1.5228378) q[3];
sx q[3];
rz(0.063604442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.30274063) q[2];
sx q[2];
rz(-1.3417696) q[2];
sx q[2];
rz(-2.0212685) q[2];
rz(2.4773347) q[3];
sx q[3];
rz(-2.7393326) q[3];
sx q[3];
rz(0.50814381) q[3];
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
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6054194) q[0];
sx q[0];
rz(-2.7689731) q[0];
sx q[0];
rz(-2.5033409) q[0];
rz(-3.1328746) q[1];
sx q[1];
rz(-2.0254841) q[1];
sx q[1];
rz(-0.4062103) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4753301) q[0];
sx q[0];
rz(-1.8436369) q[0];
sx q[0];
rz(3.0240339) q[0];
x q[1];
rz(1.6833324) q[2];
sx q[2];
rz(-0.43607682) q[2];
sx q[2];
rz(-1.9364408) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.26365852) q[1];
sx q[1];
rz(-0.78805) q[1];
sx q[1];
rz(-1.7321943) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5962192) q[3];
sx q[3];
rz(-1.6442199) q[3];
sx q[3];
rz(1.5591476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0476039) q[2];
sx q[2];
rz(-1.1270019) q[2];
sx q[2];
rz(-2.7809533) q[2];
rz(-2.4141198) q[3];
sx q[3];
rz(-2.45939) q[3];
sx q[3];
rz(-0.31807652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8135081) q[0];
sx q[0];
rz(-2.1961975) q[0];
sx q[0];
rz(-0.16950053) q[0];
rz(-0.70973474) q[1];
sx q[1];
rz(-1.7252445) q[1];
sx q[1];
rz(1.08606) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5616578) q[0];
sx q[0];
rz(-0.57889639) q[0];
sx q[0];
rz(0.41335841) q[0];
x q[1];
rz(-1.1342088) q[2];
sx q[2];
rz(-1.8545758) q[2];
sx q[2];
rz(1.5320317) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0154961) q[1];
sx q[1];
rz(-2.129678) q[1];
sx q[1];
rz(-0.11591594) q[1];
rz(-pi) q[2];
rz(-2.7602155) q[3];
sx q[3];
rz(-1.3998195) q[3];
sx q[3];
rz(1.387763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.78709948) q[2];
sx q[2];
rz(-3.0566065) q[2];
sx q[2];
rz(0.3698012) q[2];
rz(1.9434816) q[3];
sx q[3];
rz(-1.5912278) q[3];
sx q[3];
rz(-0.91355356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0030466) q[0];
sx q[0];
rz(-1.9771165) q[0];
sx q[0];
rz(1.8856915) q[0];
rz(2.1239602) q[1];
sx q[1];
rz(-1.7435278) q[1];
sx q[1];
rz(1.7494038) q[1];
rz(1.5003149) q[2];
sx q[2];
rz(-1.4983836) q[2];
sx q[2];
rz(-0.56624779) q[2];
rz(-1.4403338) q[3];
sx q[3];
rz(-1.7335907) q[3];
sx q[3];
rz(2.8681267) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
