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
rz(1.5481663) q[0];
sx q[0];
rz(-2.6731773) q[0];
sx q[0];
rz(-0.13844891) q[0];
rz(-1.2568714) q[1];
sx q[1];
rz(-1.0705907) q[1];
sx q[1];
rz(0.021477403) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3119451) q[0];
sx q[0];
rz(-1.2197131) q[0];
sx q[0];
rz(-0.21613516) q[0];
rz(0.98928605) q[2];
sx q[2];
rz(-1.0972736) q[2];
sx q[2];
rz(-0.038528942) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.8565377) q[1];
sx q[1];
rz(-1.3962588) q[1];
sx q[1];
rz(-0.15600295) q[1];
x q[2];
rz(0.50867601) q[3];
sx q[3];
rz(-0.66036105) q[3];
sx q[3];
rz(-2.3213399) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.59637493) q[2];
sx q[2];
rz(-2.6863828) q[2];
sx q[2];
rz(1.5111766) q[2];
rz(3.0916072) q[3];
sx q[3];
rz(-1.2286011) q[3];
sx q[3];
rz(0.25240067) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2216457) q[0];
sx q[0];
rz(-1.9939461) q[0];
sx q[0];
rz(-0.037121437) q[0];
rz(-3.1275753) q[1];
sx q[1];
rz(-0.57828301) q[1];
sx q[1];
rz(2.9055273) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.35821298) q[0];
sx q[0];
rz(-1.5851846) q[0];
sx q[0];
rz(0.074812513) q[0];
rz(-pi) q[1];
rz(2.057197) q[2];
sx q[2];
rz(-0.91569967) q[2];
sx q[2];
rz(0.72593216) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.6746857) q[1];
sx q[1];
rz(-2.094419) q[1];
sx q[1];
rz(0.44498131) q[1];
rz(-pi) q[2];
rz(0.25059741) q[3];
sx q[3];
rz(-1.9204233) q[3];
sx q[3];
rz(0.51167292) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1290805) q[2];
sx q[2];
rz(-1.3690288) q[2];
sx q[2];
rz(0.033626076) q[2];
rz(0.29020894) q[3];
sx q[3];
rz(-0.84607327) q[3];
sx q[3];
rz(0.39053759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6212807) q[0];
sx q[0];
rz(-0.97977591) q[0];
sx q[0];
rz(-0.2990956) q[0];
rz(-1.5215123) q[1];
sx q[1];
rz(-0.027093096) q[1];
sx q[1];
rz(-0.026195899) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55518245) q[0];
sx q[0];
rz(-1.5918709) q[0];
sx q[0];
rz(-3.1350967) q[0];
rz(-pi) q[1];
rz(-0.14339672) q[2];
sx q[2];
rz(-1.2931648) q[2];
sx q[2];
rz(-1.5806701) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.67945665) q[1];
sx q[1];
rz(-2.1894881) q[1];
sx q[1];
rz(-0.87631823) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6947909) q[3];
sx q[3];
rz(-1.7789156) q[3];
sx q[3];
rz(1.758054) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.4574778) q[2];
sx q[2];
rz(-0.44197765) q[2];
sx q[2];
rz(-0.6074062) q[2];
rz(0.36951798) q[3];
sx q[3];
rz(-1.4990467) q[3];
sx q[3];
rz(0.0079689715) q[3];
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
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.58005106) q[0];
sx q[0];
rz(-2.9396368) q[0];
sx q[0];
rz(2.0107021) q[0];
rz(-0.35218969) q[1];
sx q[1];
rz(-0.49030855) q[1];
sx q[1];
rz(2.9255548) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0433885) q[0];
sx q[0];
rz(-1.8134612) q[0];
sx q[0];
rz(-0.46671861) q[0];
x q[1];
rz(1.4540931) q[2];
sx q[2];
rz(-1.1266303) q[2];
sx q[2];
rz(1.5291052) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.8204788) q[1];
sx q[1];
rz(-1.063709) q[1];
sx q[1];
rz(0.45407461) q[1];
rz(-1.400835) q[3];
sx q[3];
rz(-1.483184) q[3];
sx q[3];
rz(-2.5690998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.5242247) q[2];
sx q[2];
rz(-2.2286131) q[2];
sx q[2];
rz(2.6101904) q[2];
rz(1.4517387) q[3];
sx q[3];
rz(-0.51133358) q[3];
sx q[3];
rz(1.0303729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3568929) q[0];
sx q[0];
rz(-1.0837311) q[0];
sx q[0];
rz(-2.8392131) q[0];
rz(2.0284292) q[1];
sx q[1];
rz(-2.1402363) q[1];
sx q[1];
rz(1.0611634) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0264176) q[0];
sx q[0];
rz(-3.0531394) q[0];
sx q[0];
rz(-0.74924107) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.9169754) q[2];
sx q[2];
rz(-1.6664522) q[2];
sx q[2];
rz(-3.0358853) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1969152) q[1];
sx q[1];
rz(-2.3507171) q[1];
sx q[1];
rz(0.18038919) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3187899) q[3];
sx q[3];
rz(-1.2815426) q[3];
sx q[3];
rz(-2.1923329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.82659668) q[2];
sx q[2];
rz(-1.7642517) q[2];
sx q[2];
rz(-2.946741) q[2];
rz(-0.74243122) q[3];
sx q[3];
rz(-2.4104379) q[3];
sx q[3];
rz(-2.8661695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.18405296) q[0];
sx q[0];
rz(-0.66127151) q[0];
sx q[0];
rz(-2.0273965) q[0];
rz(-2.8320352) q[1];
sx q[1];
rz(-0.93300262) q[1];
sx q[1];
rz(1.385744) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0838309) q[0];
sx q[0];
rz(-1.8529735) q[0];
sx q[0];
rz(2.9585005) q[0];
rz(-1.3676332) q[2];
sx q[2];
rz(-2.6250955) q[2];
sx q[2];
rz(-1.1927644) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.42658994) q[1];
sx q[1];
rz(-0.70014435) q[1];
sx q[1];
rz(-0.41058524) q[1];
x q[2];
rz(-2.3666016) q[3];
sx q[3];
rz(-2.0609566) q[3];
sx q[3];
rz(-2.91193) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4625357) q[2];
sx q[2];
rz(-0.44687301) q[2];
sx q[2];
rz(0.40979579) q[2];
rz(3.0637686) q[3];
sx q[3];
rz(-1.9255368) q[3];
sx q[3];
rz(2.3791544) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3536612) q[0];
sx q[0];
rz(-2.9910112) q[0];
sx q[0];
rz(-1.4713564) q[0];
rz(0.81368601) q[1];
sx q[1];
rz(-2.4206471) q[1];
sx q[1];
rz(-2.241316) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0519596) q[0];
sx q[0];
rz(-1.6311627) q[0];
sx q[0];
rz(-3.0230762) q[0];
x q[1];
rz(-1.2834107) q[2];
sx q[2];
rz(-1.6528439) q[2];
sx q[2];
rz(-2.1211185) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1623692) q[1];
sx q[1];
rz(-1.8368533) q[1];
sx q[1];
rz(3.1331269) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6358709) q[3];
sx q[3];
rz(-1.8037919) q[3];
sx q[3];
rz(-2.9386793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.5208931) q[2];
sx q[2];
rz(-0.40105477) q[2];
sx q[2];
rz(0.14191423) q[2];
rz(2.9698931) q[3];
sx q[3];
rz(-1.2407691) q[3];
sx q[3];
rz(-2.5920674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5080344) q[0];
sx q[0];
rz(-1.1548076) q[0];
sx q[0];
rz(-1.5615734) q[0];
rz(-0.70478565) q[1];
sx q[1];
rz(-2.2403658) q[1];
sx q[1];
rz(-0.27552342) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4954105) q[0];
sx q[0];
rz(-2.9998064) q[0];
sx q[0];
rz(-1.1205733) q[0];
x q[1];
rz(-0.24532206) q[2];
sx q[2];
rz(-1.6082014) q[2];
sx q[2];
rz(0.87329162) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6478579) q[1];
sx q[1];
rz(-1.0714515) q[1];
sx q[1];
rz(2.8060421) q[1];
rz(-pi) q[2];
x q[2];
rz(2.4750214) q[3];
sx q[3];
rz(-1.9063933) q[3];
sx q[3];
rz(0.94387142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0602818) q[2];
sx q[2];
rz(-1.5014481) q[2];
sx q[2];
rz(2.8509129) q[2];
rz(-0.79689133) q[3];
sx q[3];
rz(-0.47178888) q[3];
sx q[3];
rz(-2.1577788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31442916) q[0];
sx q[0];
rz(-1.4845347) q[0];
sx q[0];
rz(-0.86781251) q[0];
rz(2.9293291) q[1];
sx q[1];
rz(-0.8808732) q[1];
sx q[1];
rz(2.8689522) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6295687) q[0];
sx q[0];
rz(-1.4958515) q[0];
sx q[0];
rz(-1.7864947) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0402203) q[2];
sx q[2];
rz(-0.82216149) q[2];
sx q[2];
rz(3.0790975) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3301047) q[1];
sx q[1];
rz(-1.6088142) q[1];
sx q[1];
rz(2.2056715) q[1];
rz(2.5288071) q[3];
sx q[3];
rz(-2.6481508) q[3];
sx q[3];
rz(0.32705826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.049204443) q[2];
sx q[2];
rz(-2.7907351) q[2];
sx q[2];
rz(-1.1445047) q[2];
rz(-2.5149964) q[3];
sx q[3];
rz(-0.75855362) q[3];
sx q[3];
rz(0.315061) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.7476244) q[0];
sx q[0];
rz(-0.67248857) q[0];
sx q[0];
rz(2.5482063) q[0];
rz(2.4001832) q[1];
sx q[1];
rz(-0.64677042) q[1];
sx q[1];
rz(-1.5470362) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.55047148) q[0];
sx q[0];
rz(-0.24419366) q[0];
sx q[0];
rz(-1.1039735) q[0];
rz(-2.7266961) q[2];
sx q[2];
rz(-1.3438359) q[2];
sx q[2];
rz(-2.0503873) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.50362557) q[1];
sx q[1];
rz(-0.73439288) q[1];
sx q[1];
rz(-3.0375558) q[1];
rz(-0.029742777) q[3];
sx q[3];
rz(-1.7726726) q[3];
sx q[3];
rz(-2.0824528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0934304) q[2];
sx q[2];
rz(-0.54485816) q[2];
sx q[2];
rz(2.1989934) q[2];
rz(-1.4283098) q[3];
sx q[3];
rz(-1.2497679) q[3];
sx q[3];
rz(-1.1552756) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(2.9485332) q[0];
sx q[0];
rz(-1.5567224) q[0];
sx q[0];
rz(-1.3514883) q[0];
rz(-1.1763186) q[1];
sx q[1];
rz(-0.90570025) q[1];
sx q[1];
rz(-0.66014231) q[1];
rz(-1.0084441) q[2];
sx q[2];
rz(-2.8165419) q[2];
sx q[2];
rz(-1.4680924) q[2];
rz(-0.42849937) q[3];
sx q[3];
rz(-2.5658723) q[3];
sx q[3];
rz(0.65278237) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
