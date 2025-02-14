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
rz(-2.1383877) q[0];
sx q[0];
rz(-1.9411074) q[0];
sx q[0];
rz(2.1152273) q[0];
rz(0.29844555) q[1];
sx q[1];
rz(-1.6331853) q[1];
sx q[1];
rz(1.0025947) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5110785) q[0];
sx q[0];
rz(-1.3670232) q[0];
sx q[0];
rz(-3.1183) q[0];
x q[1];
rz(2.4786776) q[2];
sx q[2];
rz(-1.6143245) q[2];
sx q[2];
rz(1.3452665) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.57220578) q[1];
sx q[1];
rz(-1.9574698) q[1];
sx q[1];
rz(1.4264849) q[1];
rz(2.3839124) q[3];
sx q[3];
rz(-2.3509988) q[3];
sx q[3];
rz(-2.9942715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.5528494) q[2];
sx q[2];
rz(-1.0007977) q[2];
sx q[2];
rz(0.0351077) q[2];
rz(1.2037753) q[3];
sx q[3];
rz(-1.3222008) q[3];
sx q[3];
rz(-2.8573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0725919) q[0];
sx q[0];
rz(-0.44206107) q[0];
sx q[0];
rz(-1.0042071) q[0];
rz(-0.12241157) q[1];
sx q[1];
rz(-1.6740084) q[1];
sx q[1];
rz(2.4845128) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9417038) q[0];
sx q[0];
rz(-2.7452069) q[0];
sx q[0];
rz(2.928171) q[0];
rz(-2.9802965) q[2];
sx q[2];
rz(-1.2968924) q[2];
sx q[2];
rz(1.0082001) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2942336) q[1];
sx q[1];
rz(-2.2072133) q[1];
sx q[1];
rz(0.60346236) q[1];
x q[2];
rz(1.3260001) q[3];
sx q[3];
rz(-0.80202937) q[3];
sx q[3];
rz(1.0753461) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.5848026) q[2];
sx q[2];
rz(-0.57791296) q[2];
sx q[2];
rz(2.54971) q[2];
rz(1.6651734) q[3];
sx q[3];
rz(-1.5341325) q[3];
sx q[3];
rz(2.2679451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1378491) q[0];
sx q[0];
rz(-2.9326404) q[0];
sx q[0];
rz(-1.9269706) q[0];
rz(2.1851723) q[1];
sx q[1];
rz(-1.1914057) q[1];
sx q[1];
rz(-1.9699684) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4223683) q[0];
sx q[0];
rz(-2.5107493) q[0];
sx q[0];
rz(-2.2792086) q[0];
x q[1];
rz(-0.20522988) q[2];
sx q[2];
rz(-2.2608888) q[2];
sx q[2];
rz(-3.1179259) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.7293766) q[1];
sx q[1];
rz(-2.4640485) q[1];
sx q[1];
rz(0.55528806) q[1];
rz(-pi) q[2];
rz(-0.089853386) q[3];
sx q[3];
rz(-1.615534) q[3];
sx q[3];
rz(-0.60715946) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.20022915) q[2];
sx q[2];
rz(-1.0503294) q[2];
sx q[2];
rz(-1.7639147) q[2];
rz(-2.7348943) q[3];
sx q[3];
rz(-0.64928693) q[3];
sx q[3];
rz(-2.6810834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3520626) q[0];
sx q[0];
rz(-1.2593513) q[0];
sx q[0];
rz(1.9011185) q[0];
rz(0.71818304) q[1];
sx q[1];
rz(-1.8156464) q[1];
sx q[1];
rz(2.9026418) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0485503) q[0];
sx q[0];
rz(-1.6918381) q[0];
sx q[0];
rz(-2.8516475) q[0];
rz(-pi) q[1];
rz(2.3283655) q[2];
sx q[2];
rz(-1.0095749) q[2];
sx q[2];
rz(-0.49408572) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0480369) q[1];
sx q[1];
rz(-2.844792) q[1];
sx q[1];
rz(2.0778627) q[1];
rz(-pi) q[2];
rz(-0.53455456) q[3];
sx q[3];
rz(-1.6865493) q[3];
sx q[3];
rz(-2.4730269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.51603985) q[2];
sx q[2];
rz(-1.3935139) q[2];
sx q[2];
rz(-1.7986521) q[2];
rz(-2.2602153) q[3];
sx q[3];
rz(-2.9727029) q[3];
sx q[3];
rz(2.3605997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
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
rz(1.7278904) q[0];
sx q[0];
rz(-0.24402937) q[0];
sx q[0];
rz(-1.6572886) q[0];
rz(0.65385747) q[1];
sx q[1];
rz(-1.5212719) q[1];
sx q[1];
rz(-0.13793764) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50820275) q[0];
sx q[0];
rz(-1.7195367) q[0];
sx q[0];
rz(1.2138019) q[0];
x q[1];
rz(1.8191255) q[2];
sx q[2];
rz(-2.2086124) q[2];
sx q[2];
rz(-2.1030104) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.9508507) q[1];
sx q[1];
rz(-0.44176451) q[1];
sx q[1];
rz(-1.8015693) q[1];
rz(-1.3119427) q[3];
sx q[3];
rz(-1.6248684) q[3];
sx q[3];
rz(1.3075706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5518034) q[2];
sx q[2];
rz(-1.9686331) q[2];
sx q[2];
rz(-3.0544082) q[2];
rz(-3.0269571) q[3];
sx q[3];
rz(-0.81511027) q[3];
sx q[3];
rz(-0.4172999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7755985) q[0];
sx q[0];
rz(-1.5394779) q[0];
sx q[0];
rz(-0.48496801) q[0];
rz(-0.94351774) q[1];
sx q[1];
rz(-0.8129932) q[1];
sx q[1];
rz(-1.9224723) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6935558) q[0];
sx q[0];
rz(-2.4139934) q[0];
sx q[0];
rz(1.4240525) q[0];
rz(-1.8928244) q[2];
sx q[2];
rz(-1.5552943) q[2];
sx q[2];
rz(-2.6982234) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.74138) q[1];
sx q[1];
rz(-1.6467491) q[1];
sx q[1];
rz(1.6981359) q[1];
x q[2];
rz(0.7597105) q[3];
sx q[3];
rz(-1.1897161) q[3];
sx q[3];
rz(-2.423513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.24011) q[2];
sx q[2];
rz(-2.0811452) q[2];
sx q[2];
rz(-0.94129747) q[2];
rz(2.8955722) q[3];
sx q[3];
rz(-1.4964208) q[3];
sx q[3];
rz(1.3601607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0435903) q[0];
sx q[0];
rz(-2.0397546) q[0];
sx q[0];
rz(1.164042) q[0];
rz(-1.7835435) q[1];
sx q[1];
rz(-0.86654228) q[1];
sx q[1];
rz(-1.7485626) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1355125) q[0];
sx q[0];
rz(-0.3835668) q[0];
sx q[0];
rz(1.3719029) q[0];
rz(-pi) q[1];
x q[1];
rz(0.4460469) q[2];
sx q[2];
rz(-2.5452633) q[2];
sx q[2];
rz(-1.9976384) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.021780304) q[1];
sx q[1];
rz(-2.0396898) q[1];
sx q[1];
rz(-1.6890668) q[1];
rz(2.249516) q[3];
sx q[3];
rz(-1.3952655) q[3];
sx q[3];
rz(-0.20993671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.4712269) q[2];
sx q[2];
rz(-1.7134075) q[2];
sx q[2];
rz(-2.8426389) q[2];
rz(-2.7361338) q[3];
sx q[3];
rz(-0.36646989) q[3];
sx q[3];
rz(-2.5772742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4295171) q[0];
sx q[0];
rz(-2.8686664) q[0];
sx q[0];
rz(-0.78936973) q[0];
rz(-0.40496597) q[1];
sx q[1];
rz(-1.3507495) q[1];
sx q[1];
rz(0.49496034) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4433561) q[0];
sx q[0];
rz(-2.3670417) q[0];
sx q[0];
rz(-2.2775035) q[0];
rz(2.8414824) q[2];
sx q[2];
rz(-2.1202188) q[2];
sx q[2];
rz(0.97893366) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.61117426) q[1];
sx q[1];
rz(-1.8137167) q[1];
sx q[1];
rz(2.3526885) q[1];
rz(1.0010507) q[3];
sx q[3];
rz(-0.82017878) q[3];
sx q[3];
rz(-2.6626183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.1627545) q[2];
sx q[2];
rz(-1.9163722) q[2];
sx q[2];
rz(2.8863353) q[2];
rz(-0.034959547) q[3];
sx q[3];
rz(-1.5233663) q[3];
sx q[3];
rz(0.24622723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.404945) q[0];
sx q[0];
rz(-1.0628137) q[0];
sx q[0];
rz(-2.431562) q[0];
rz(1.0011477) q[1];
sx q[1];
rz(-0.89718693) q[1];
sx q[1];
rz(-2.5206916) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17711711) q[0];
sx q[0];
rz(-1.510448) q[0];
sx q[0];
rz(-1.9326769) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65695564) q[2];
sx q[2];
rz(-1.1816506) q[2];
sx q[2];
rz(-2.3537113) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.94540751) q[1];
sx q[1];
rz(-1.9176449) q[1];
sx q[1];
rz(-1.4566684) q[1];
x q[2];
rz(-2.8594703) q[3];
sx q[3];
rz(-2.3446313) q[3];
sx q[3];
rz(2.7180501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.8093449) q[2];
sx q[2];
rz(-1.5231909) q[2];
sx q[2];
rz(-2.8524032) q[2];
rz(1.1122164) q[3];
sx q[3];
rz(-2.8465392) q[3];
sx q[3];
rz(-2.1212063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6638829) q[0];
sx q[0];
rz(-1.3417256) q[0];
sx q[0];
rz(0.92673242) q[0];
rz(2.0345188) q[1];
sx q[1];
rz(-1.6278382) q[1];
sx q[1];
rz(0.56799299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3044514) q[0];
sx q[0];
rz(-2.1689609) q[0];
sx q[0];
rz(1.9564081) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6935518) q[2];
sx q[2];
rz(-2.5841641) q[2];
sx q[2];
rz(-0.008226062) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.24577877) q[1];
sx q[1];
rz(-1.2061283) q[1];
sx q[1];
rz(-0.90943894) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0509012) q[3];
sx q[3];
rz(-2.1445159) q[3];
sx q[3];
rz(-0.75999505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.0893112) q[2];
sx q[2];
rz(-2.3312882) q[2];
sx q[2];
rz(-1.7029765) q[2];
rz(-0.29664052) q[3];
sx q[3];
rz(-2.7999925) q[3];
sx q[3];
rz(-0.44704416) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.457837) q[0];
sx q[0];
rz(-2.044027) q[0];
sx q[0];
rz(-2.5720163) q[0];
rz(-0.33869047) q[1];
sx q[1];
rz(-1.6117922) q[1];
sx q[1];
rz(1.502996) q[1];
rz(-1.3647625) q[2];
sx q[2];
rz(-1.3246957) q[2];
sx q[2];
rz(-0.93304721) q[2];
rz(2.972907) q[3];
sx q[3];
rz(-1.8416234) q[3];
sx q[3];
rz(2.4467322) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
