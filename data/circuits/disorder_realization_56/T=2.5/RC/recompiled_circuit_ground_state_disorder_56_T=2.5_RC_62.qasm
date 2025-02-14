OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.174515) q[0];
sx q[0];
rz(-1.6835901) q[0];
sx q[0];
rz(-0.30012497) q[0];
rz(-1.9441654) q[1];
sx q[1];
rz(-1.5265042) q[1];
sx q[1];
rz(1.5010897) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59812058) q[0];
sx q[0];
rz(-1.5859005) q[0];
sx q[0];
rz(1.5669109) q[0];
rz(-pi) q[1];
rz(-2.757171) q[2];
sx q[2];
rz(-0.84058207) q[2];
sx q[2];
rz(-0.79193774) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7607494) q[1];
sx q[1];
rz(-2.2135206) q[1];
sx q[1];
rz(-1.8442783) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80880161) q[3];
sx q[3];
rz(-1.5672657) q[3];
sx q[3];
rz(0.88909066) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.40452051) q[2];
sx q[2];
rz(-1.1417192) q[2];
sx q[2];
rz(0.70297757) q[2];
rz(-2.5718555) q[3];
sx q[3];
rz(-0.8786141) q[3];
sx q[3];
rz(1.5574633) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0074145929) q[0];
sx q[0];
rz(-0.61892048) q[0];
sx q[0];
rz(-1.9084357) q[0];
rz(0.59457072) q[1];
sx q[1];
rz(-2.1023127) q[1];
sx q[1];
rz(2.0436683) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63609517) q[0];
sx q[0];
rz(-1.5042802) q[0];
sx q[0];
rz(1.3696704) q[0];
x q[1];
rz(0.49169831) q[2];
sx q[2];
rz(-1.227081) q[2];
sx q[2];
rz(2.6298486) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6095396) q[1];
sx q[1];
rz(-1.5960448) q[1];
sx q[1];
rz(1.4933426) q[1];
x q[2];
rz(3.0497562) q[3];
sx q[3];
rz(-0.79376924) q[3];
sx q[3];
rz(2.5240999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7132831) q[2];
sx q[2];
rz(-1.8417336) q[2];
sx q[2];
rz(-1.5353047) q[2];
rz(0.71896583) q[3];
sx q[3];
rz(-0.9459559) q[3];
sx q[3];
rz(0.21397056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.83800256) q[0];
sx q[0];
rz(-2.328023) q[0];
sx q[0];
rz(-0.59233061) q[0];
rz(1.8968286) q[1];
sx q[1];
rz(-1.1400305) q[1];
sx q[1];
rz(-0.14561428) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.94731047) q[0];
sx q[0];
rz(-2.6366451) q[0];
sx q[0];
rz(-2.8584252) q[0];
x q[1];
rz(0.8823186) q[2];
sx q[2];
rz(-1.4704508) q[2];
sx q[2];
rz(1.060263) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.023886746) q[1];
sx q[1];
rz(-0.63359208) q[1];
sx q[1];
rz(-1.078245) q[1];
rz(-pi) q[2];
rz(1.0907902) q[3];
sx q[3];
rz(-2.601805) q[3];
sx q[3];
rz(-1.9095498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.40272063) q[2];
sx q[2];
rz(-2.1108184) q[2];
sx q[2];
rz(1.9090451) q[2];
rz(0.23972073) q[3];
sx q[3];
rz(-1.0545878) q[3];
sx q[3];
rz(-2.6565552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0204912) q[0];
sx q[0];
rz(-2.4395269) q[0];
sx q[0];
rz(0.030601587) q[0];
rz(-2.5864511) q[1];
sx q[1];
rz(-1.6629013) q[1];
sx q[1];
rz(2.0487002) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3557778) q[0];
sx q[0];
rz(-0.91453881) q[0];
sx q[0];
rz(0.86571619) q[0];
rz(-pi) q[1];
rz(-2.3762114) q[2];
sx q[2];
rz(-1.5281406) q[2];
sx q[2];
rz(2.9431107) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.79697414) q[1];
sx q[1];
rz(-0.94616468) q[1];
sx q[1];
rz(0.25298869) q[1];
rz(-0.84597702) q[3];
sx q[3];
rz(-2.1095157) q[3];
sx q[3];
rz(-1.2332476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.0838202) q[2];
sx q[2];
rz(-1.7744935) q[2];
sx q[2];
rz(2.8687381) q[2];
rz(-1.3911635) q[3];
sx q[3];
rz(-0.83943668) q[3];
sx q[3];
rz(2.1896037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4680173) q[0];
sx q[0];
rz(-1.8033569) q[0];
sx q[0];
rz(1.2982298) q[0];
rz(2.4507554) q[1];
sx q[1];
rz(-1.1996484) q[1];
sx q[1];
rz(1.615049) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1896601) q[0];
sx q[0];
rz(-0.88610035) q[0];
sx q[0];
rz(0.71345274) q[0];
rz(-pi) q[1];
rz(1.2882907) q[2];
sx q[2];
rz(-1.9108678) q[2];
sx q[2];
rz(-1.1929389) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95350515) q[1];
sx q[1];
rz(-1.5449448) q[1];
sx q[1];
rz(-0.8698747) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4590309) q[3];
sx q[3];
rz(-0.45387156) q[3];
sx q[3];
rz(2.0477432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.7319298) q[2];
sx q[2];
rz(-2.3356428) q[2];
sx q[2];
rz(-0.16732495) q[2];
rz(1.3903728) q[3];
sx q[3];
rz(-2.6087587) q[3];
sx q[3];
rz(2.4840241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.14199) q[0];
sx q[0];
rz(-0.36407343) q[0];
sx q[0];
rz(-1.2784736) q[0];
rz(-1.7059884) q[1];
sx q[1];
rz(-1.6537138) q[1];
sx q[1];
rz(0.37567589) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77462353) q[0];
sx q[0];
rz(-1.3448098) q[0];
sx q[0];
rz(-1.8480728) q[0];
x q[1];
rz(2.9922688) q[2];
sx q[2];
rz(-1.4027601) q[2];
sx q[2];
rz(-2.4150893) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.73114016) q[1];
sx q[1];
rz(-1.2527494) q[1];
sx q[1];
rz(2.9067626) q[1];
rz(-pi) q[2];
rz(-1.814498) q[3];
sx q[3];
rz(-2.3117495) q[3];
sx q[3];
rz(-1.0761225) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0659236) q[2];
sx q[2];
rz(-2.0772987) q[2];
sx q[2];
rz(-0.86090487) q[2];
rz(1.9979477) q[3];
sx q[3];
rz(-2.836561) q[3];
sx q[3];
rz(2.8633269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.025573108) q[0];
sx q[0];
rz(-1.3207859) q[0];
sx q[0];
rz(0.75468165) q[0];
rz(-2.2017551) q[1];
sx q[1];
rz(-0.87516963) q[1];
sx q[1];
rz(-1.2831203) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.59903501) q[0];
sx q[0];
rz(-1.577566) q[0];
sx q[0];
rz(1.2640796) q[0];
x q[1];
rz(0.36410113) q[2];
sx q[2];
rz(-1.1392829) q[2];
sx q[2];
rz(-1.3442775) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.922328) q[1];
sx q[1];
rz(-1.9524876) q[1];
sx q[1];
rz(-1.2569409) q[1];
rz(0.8831034) q[3];
sx q[3];
rz(-1.9839483) q[3];
sx q[3];
rz(-0.93078155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.5963886) q[2];
sx q[2];
rz(-1.117492) q[2];
sx q[2];
rz(-1.0835353) q[2];
rz(-2.3986473) q[3];
sx q[3];
rz(-1.5321621) q[3];
sx q[3];
rz(-1.4325745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.1525986) q[0];
sx q[0];
rz(-0.78482634) q[0];
sx q[0];
rz(-2.3866744) q[0];
rz(-0.49916357) q[1];
sx q[1];
rz(-2.1776336) q[1];
sx q[1];
rz(2.4430433) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7592418) q[0];
sx q[0];
rz(-0.97433358) q[0];
sx q[0];
rz(1.9475219) q[0];
rz(-0.46496181) q[2];
sx q[2];
rz(-1.4413646) q[2];
sx q[2];
rz(1.7120106) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.10741988) q[1];
sx q[1];
rz(-1.0386969) q[1];
sx q[1];
rz(-3.1183395) q[1];
x q[2];
rz(-1.9898765) q[3];
sx q[3];
rz(-1.5718096) q[3];
sx q[3];
rz(0.59194293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.90427202) q[2];
sx q[2];
rz(-0.99232435) q[2];
sx q[2];
rz(-0.80238706) q[2];
rz(-0.55195105) q[3];
sx q[3];
rz(-0.80058432) q[3];
sx q[3];
rz(2.5210181) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17860831) q[0];
sx q[0];
rz(-1.7010138) q[0];
sx q[0];
rz(1.1075903) q[0];
rz(1.3811318) q[1];
sx q[1];
rz(-1.3637204) q[1];
sx q[1];
rz(-1.507087) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3920306) q[0];
sx q[0];
rz(-1.9641341) q[0];
sx q[0];
rz(1.7225207) q[0];
rz(-pi) q[1];
rz(-1.9805895) q[2];
sx q[2];
rz(-1.2772182) q[2];
sx q[2];
rz(0.30337203) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.9642051) q[1];
sx q[1];
rz(-1.2086165) q[1];
sx q[1];
rz(1.7262579) q[1];
x q[2];
rz(-3.0595257) q[3];
sx q[3];
rz(-0.82857271) q[3];
sx q[3];
rz(2.9794793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.0570809) q[2];
sx q[2];
rz(-2.9456186) q[2];
sx q[2];
rz(1.8692807) q[2];
rz(1.48014) q[3];
sx q[3];
rz(-1.1353761) q[3];
sx q[3];
rz(2.541466) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24499527) q[0];
sx q[0];
rz(-2.5801881) q[0];
sx q[0];
rz(1.2835314) q[0];
rz(3.1237579) q[1];
sx q[1];
rz(-2.5051038) q[1];
sx q[1];
rz(-0.12558118) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.784362) q[0];
sx q[0];
rz(-1.0166369) q[0];
sx q[0];
rz(-0.016655075) q[0];
rz(-pi) q[1];
rz(-2.0415885) q[2];
sx q[2];
rz(-0.94816899) q[2];
sx q[2];
rz(2.560844) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.64396137) q[1];
sx q[1];
rz(-0.61986524) q[1];
sx q[1];
rz(-2.3513112) q[1];
rz(-1.3036895) q[3];
sx q[3];
rz(-1.7648762) q[3];
sx q[3];
rz(-1.1052856) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.77832001) q[2];
sx q[2];
rz(-1.2691701) q[2];
sx q[2];
rz(2.3385091) q[2];
rz(-2.965029) q[3];
sx q[3];
rz(-1.8162138) q[3];
sx q[3];
rz(0.77809063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17978996) q[0];
sx q[0];
rz(-0.50146865) q[0];
sx q[0];
rz(1.5487221) q[0];
rz(1.8190307) q[1];
sx q[1];
rz(-1.7772728) q[1];
sx q[1];
rz(-1.5900236) q[1];
rz(2.8990135) q[2];
sx q[2];
rz(-1.5097813) q[2];
sx q[2];
rz(-0.079955242) q[2];
rz(1.6311036) q[3];
sx q[3];
rz(-1.9077488) q[3];
sx q[3];
rz(0.13691402) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
