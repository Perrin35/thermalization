OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5105628) q[0];
sx q[0];
rz(-0.50463843) q[0];
sx q[0];
rz(2.7161427) q[0];
rz(0.57234859) q[1];
sx q[1];
rz(4.2387716) q[1];
sx q[1];
rz(7.4833202) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6365189) q[0];
sx q[0];
rz(-2.8905792) q[0];
sx q[0];
rz(-3.1025629) q[0];
rz(2.8699371) q[2];
sx q[2];
rz(-1.2363529) q[2];
sx q[2];
rz(-1.7890499) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.7537743) q[1];
sx q[1];
rz(-1.942831) q[1];
sx q[1];
rz(-1.340999) q[1];
x q[2];
rz(2.9619734) q[3];
sx q[3];
rz(-0.67905513) q[3];
sx q[3];
rz(1.8906901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9878865) q[2];
sx q[2];
rz(-1.1517297) q[2];
sx q[2];
rz(-2.494334) q[2];
rz(2.9636532) q[3];
sx q[3];
rz(-1.2578332) q[3];
sx q[3];
rz(-1.2037163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7385638) q[0];
sx q[0];
rz(-3.0572427) q[0];
sx q[0];
rz(-0.52363288) q[0];
rz(1.2298443) q[1];
sx q[1];
rz(-0.74408999) q[1];
sx q[1];
rz(2.8935166) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1551127) q[0];
sx q[0];
rz(-2.393828) q[0];
sx q[0];
rz(2.98481) q[0];
rz(0.78901864) q[2];
sx q[2];
rz(-0.58525267) q[2];
sx q[2];
rz(-1.326013) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.16399244) q[1];
sx q[1];
rz(-1.6095785) q[1];
sx q[1];
rz(1.7165403) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.66623293) q[3];
sx q[3];
rz(-1.1403475) q[3];
sx q[3];
rz(2.7048049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.5292042) q[2];
sx q[2];
rz(-0.91823053) q[2];
sx q[2];
rz(-2.6017792) q[2];
rz(0.16048935) q[3];
sx q[3];
rz(-1.548998) q[3];
sx q[3];
rz(0.81015712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4983343) q[0];
sx q[0];
rz(-0.67688268) q[0];
sx q[0];
rz(0.13370378) q[0];
rz(-1.5191822) q[1];
sx q[1];
rz(-1.8459873) q[1];
sx q[1];
rz(-0.79230961) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3334219) q[0];
sx q[0];
rz(-1.9229888) q[0];
sx q[0];
rz(-2.9461224) q[0];
rz(-pi) q[1];
rz(1.4390723) q[2];
sx q[2];
rz(-1.4274397) q[2];
sx q[2];
rz(-1.5174587) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.5416273) q[1];
sx q[1];
rz(-2.1481594) q[1];
sx q[1];
rz(-0.75871102) q[1];
rz(-pi) q[2];
rz(0.30012975) q[3];
sx q[3];
rz(-1.5709166) q[3];
sx q[3];
rz(0.67226582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.2866659) q[2];
sx q[2];
rz(-0.67015219) q[2];
sx q[2];
rz(-2.3112042) q[2];
rz(1.3487799) q[3];
sx q[3];
rz(-1.034779) q[3];
sx q[3];
rz(-0.59876284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8778359) q[0];
sx q[0];
rz(-2.1649375) q[0];
sx q[0];
rz(-2.7184955) q[0];
rz(1.1571723) q[1];
sx q[1];
rz(-1.4856228) q[1];
sx q[1];
rz(-2.5194936) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3836648) q[0];
sx q[0];
rz(-0.94449857) q[0];
sx q[0];
rz(-2.6561198) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1272264) q[2];
sx q[2];
rz(-2.6831919) q[2];
sx q[2];
rz(-2.1200402) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9378949) q[1];
sx q[1];
rz(-2.7337498) q[1];
sx q[1];
rz(1.5289343) q[1];
rz(-pi) q[2];
rz(0.027340319) q[3];
sx q[3];
rz(-0.84544824) q[3];
sx q[3];
rz(-0.10408653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.28079924) q[2];
sx q[2];
rz(-1.7223027) q[2];
sx q[2];
rz(-0.61817509) q[2];
rz(1.6587967) q[3];
sx q[3];
rz(-1.3525454) q[3];
sx q[3];
rz(1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(2.9789155) q[0];
sx q[0];
rz(-1.8122883) q[0];
sx q[0];
rz(-0.83258122) q[0];
rz(-2.816448) q[1];
sx q[1];
rz(-2.145576) q[1];
sx q[1];
rz(0.064373374) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58080855) q[0];
sx q[0];
rz(-2.7940895) q[0];
sx q[0];
rz(2.7827713) q[0];
rz(0.60466296) q[2];
sx q[2];
rz(-0.69398601) q[2];
sx q[2];
rz(-2.1460905) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2657703) q[1];
sx q[1];
rz(-1.0604621) q[1];
sx q[1];
rz(-2.2052664) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7835618) q[3];
sx q[3];
rz(-0.97409407) q[3];
sx q[3];
rz(-0.95355129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.23318204) q[2];
sx q[2];
rz(-2.5711214) q[2];
sx q[2];
rz(1.849966) q[2];
rz(-2.6455247) q[3];
sx q[3];
rz(-0.78947624) q[3];
sx q[3];
rz(1.1424278) q[3];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9555776) q[0];
sx q[0];
rz(-2.8514974) q[0];
sx q[0];
rz(-2.7225323) q[0];
rz(-0.31173197) q[1];
sx q[1];
rz(-1.5429976) q[1];
sx q[1];
rz(-2.9980803) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18565047) q[0];
sx q[0];
rz(-2.5505192) q[0];
sx q[0];
rz(-0.67703621) q[0];
rz(1.3355876) q[2];
sx q[2];
rz(-1.303788) q[2];
sx q[2];
rz(-2.5036204) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.12998768) q[1];
sx q[1];
rz(-1.2530348) q[1];
sx q[1];
rz(2.3684273) q[1];
rz(-pi) q[2];
rz(2.8815202) q[3];
sx q[3];
rz(-1.4315769) q[3];
sx q[3];
rz(-1.2633439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.781337) q[2];
sx q[2];
rz(-2.2402253) q[2];
sx q[2];
rz(1.1330913) q[2];
rz(2.0356778) q[3];
sx q[3];
rz(-1.4191041) q[3];
sx q[3];
rz(2.6667986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7853506) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(-1.6957977) q[0];
rz(-0.71714199) q[1];
sx q[1];
rz(-1.278806) q[1];
sx q[1];
rz(-0.76146567) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8728947) q[0];
sx q[0];
rz(-0.94572645) q[0];
sx q[0];
rz(-1.752466) q[0];
rz(0.47777678) q[2];
sx q[2];
rz(-2.226916) q[2];
sx q[2];
rz(0.28766866) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.95132212) q[1];
sx q[1];
rz(-2.1441064) q[1];
sx q[1];
rz(3.0828397) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7588571) q[3];
sx q[3];
rz(-1.4131964) q[3];
sx q[3];
rz(-3.0732875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.7722499) q[2];
sx q[2];
rz(-2.6476634) q[2];
sx q[2];
rz(2.210145) q[2];
rz(-0.61819589) q[3];
sx q[3];
rz(-1.6341011) q[3];
sx q[3];
rz(2.1759694) q[3];
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
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44929993) q[0];
sx q[0];
rz(-1.3626784) q[0];
sx q[0];
rz(1.8071254) q[0];
rz(2.6284699) q[1];
sx q[1];
rz(-0.98315364) q[1];
sx q[1];
rz(-1.9122874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.089565) q[0];
sx q[0];
rz(-2.2917245) q[0];
sx q[0];
rz(-0.57336471) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0487829) q[2];
sx q[2];
rz(-1.3121737) q[2];
sx q[2];
rz(-1.5802204) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.38459435) q[1];
sx q[1];
rz(-1.726686) q[1];
sx q[1];
rz(-1.5916567) q[1];
rz(-pi) q[2];
rz(2.4793998) q[3];
sx q[3];
rz(-1.5569038) q[3];
sx q[3];
rz(0.77869773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.8339771) q[2];
sx q[2];
rz(-0.5270842) q[2];
sx q[2];
rz(0.21044593) q[2];
rz(-3.0757507) q[3];
sx q[3];
rz(-0.47271553) q[3];
sx q[3];
rz(-0.79894799) q[3];
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
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54302067) q[0];
sx q[0];
rz(-0.6544756) q[0];
sx q[0];
rz(-1.2731113) q[0];
rz(1.2741362) q[1];
sx q[1];
rz(-0.68398634) q[1];
sx q[1];
rz(0.88409105) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40128517) q[0];
sx q[0];
rz(-1.5710366) q[0];
sx q[0];
rz(2.1876641) q[0];
rz(-pi) q[1];
x q[1];
rz(0.52427788) q[2];
sx q[2];
rz(-1.345621) q[2];
sx q[2];
rz(0.88722992) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.99779656) q[1];
sx q[1];
rz(-1.9854768) q[1];
sx q[1];
rz(0.30325496) q[1];
rz(-pi) q[2];
x q[2];
rz(0.30852921) q[3];
sx q[3];
rz(-1.8179699) q[3];
sx q[3];
rz(0.17251523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.990443) q[2];
sx q[2];
rz(-2.1572025) q[2];
sx q[2];
rz(-1.3746064) q[2];
rz(2.9511342) q[3];
sx q[3];
rz(-0.74646598) q[3];
sx q[3];
rz(1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16354887) q[0];
sx q[0];
rz(-1.4530285) q[0];
sx q[0];
rz(-1.2646041) q[0];
rz(3.0679852) q[1];
sx q[1];
rz(-2.059748) q[1];
sx q[1];
rz(2.5216865) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0479931) q[0];
sx q[0];
rz(-1.7353829) q[0];
sx q[0];
rz(0.82407804) q[0];
x q[1];
rz(0.61087278) q[2];
sx q[2];
rz(-1.9985285) q[2];
sx q[2];
rz(1.7868702) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.318639) q[1];
sx q[1];
rz(-1.3711437) q[1];
sx q[1];
rz(0.21491383) q[1];
rz(-0.46658559) q[3];
sx q[3];
rz(-1.3562725) q[3];
sx q[3];
rz(1.9746575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.8668883) q[2];
sx q[2];
rz(-2.4394749) q[2];
sx q[2];
rz(-2.9956024) q[2];
rz(-0.22249666) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(2.4116624) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20919007) q[0];
sx q[0];
rz(-1.9575735) q[0];
sx q[0];
rz(0.14458543) q[0];
rz(-1.2296386) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(-1.3521005) q[2];
sx q[2];
rz(-0.98697294) q[2];
sx q[2];
rz(1.0728991) q[2];
rz(-1.4208117) q[3];
sx q[3];
rz(-0.21258988) q[3];
sx q[3];
rz(-1.7439738) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
