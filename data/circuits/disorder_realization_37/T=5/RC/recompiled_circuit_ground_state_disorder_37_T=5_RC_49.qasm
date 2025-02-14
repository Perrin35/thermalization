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
rz(-2.5692441) q[1];
sx q[1];
rz(-1.0971789) q[1];
sx q[1];
rz(-1.2001349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5050738) q[0];
sx q[0];
rz(-2.8905792) q[0];
sx q[0];
rz(-0.039029718) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2283353) q[2];
sx q[2];
rz(-0.42762556) q[2];
sx q[2];
rz(2.0562125) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.26781005) q[1];
sx q[1];
rz(-1.7846196) q[1];
sx q[1];
rz(-0.38114433) q[1];
rz(-pi) q[2];
rz(2.9619734) q[3];
sx q[3];
rz(-2.4625375) q[3];
sx q[3];
rz(-1.8906901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1537062) q[2];
sx q[2];
rz(-1.989863) q[2];
sx q[2];
rz(2.494334) q[2];
rz(-2.9636532) q[3];
sx q[3];
rz(-1.8837594) q[3];
sx q[3];
rz(-1.2037163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4030289) q[0];
sx q[0];
rz(-0.084349923) q[0];
sx q[0];
rz(2.6179598) q[0];
rz(1.9117484) q[1];
sx q[1];
rz(-0.74408999) q[1];
sx q[1];
rz(-2.8935166) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1988293) q[0];
sx q[0];
rz(-0.83434767) q[0];
sx q[0];
rz(1.7146066) q[0];
x q[1];
rz(2.7047728) q[2];
sx q[2];
rz(-1.9736276) q[2];
sx q[2];
rz(2.6879329) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7404799) q[1];
sx q[1];
rz(-1.4251627) q[1];
sx q[1];
rz(0.039197368) q[1];
rz(-2.5026191) q[3];
sx q[3];
rz(-0.77497122) q[3];
sx q[3];
rz(-0.6460305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5292042) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(2.6017792) q[2];
rz(-0.16048935) q[3];
sx q[3];
rz(-1.5925946) q[3];
sx q[3];
rz(0.81015712) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(-1.6432583) q[0];
sx q[0];
rz(-0.67688268) q[0];
sx q[0];
rz(-3.0078889) q[0];
rz(1.5191822) q[1];
sx q[1];
rz(-1.8459873) q[1];
sx q[1];
rz(-2.349283) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83082047) q[0];
sx q[0];
rz(-1.3874652) q[0];
sx q[0];
rz(1.2123327) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9970005) q[2];
sx q[2];
rz(-1.4404313) q[2];
sx q[2];
rz(3.0693288) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50666675) q[1];
sx q[1];
rz(-2.185195) q[1];
sx q[1];
rz(-2.3022815) q[1];
x q[2];
rz(-2.8414629) q[3];
sx q[3];
rz(-1.570676) q[3];
sx q[3];
rz(-0.67226582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.2866659) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(-2.3112042) q[2];
rz(-1.3487799) q[3];
sx q[3];
rz(-2.1068137) q[3];
sx q[3];
rz(-0.59876284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26375672) q[0];
sx q[0];
rz(-0.97665518) q[0];
sx q[0];
rz(-2.7184955) q[0];
rz(1.1571723) q[1];
sx q[1];
rz(-1.6559699) q[1];
sx q[1];
rz(2.5194936) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.3836648) q[0];
sx q[0];
rz(-0.94449857) q[0];
sx q[0];
rz(-0.48547283) q[0];
x q[1];
rz(-1.1515205) q[2];
sx q[2];
rz(-1.3797233) q[2];
sx q[2];
rz(0.1464459) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8922975) q[1];
sx q[1];
rz(-1.163332) q[1];
sx q[1];
rz(3.1235126) q[1];
rz(-pi) q[2];
rz(-2.29633) q[3];
sx q[3];
rz(-1.5503395) q[3];
sx q[3];
rz(-1.6930228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8607934) q[2];
sx q[2];
rz(-1.7223027) q[2];
sx q[2];
rz(-2.5234176) q[2];
rz(-1.482796) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(-1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(-0.16267714) q[0];
sx q[0];
rz(-1.8122883) q[0];
sx q[0];
rz(-0.83258122) q[0];
rz(0.32514462) q[1];
sx q[1];
rz(-0.99601662) q[1];
sx q[1];
rz(3.0772193) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125809) q[0];
sx q[0];
rz(-1.6906749) q[0];
sx q[0];
rz(-2.8146312) q[0];
rz(-2.5369297) q[2];
sx q[2];
rz(-2.4476066) q[2];
sx q[2];
rz(2.1460905) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.0400914) q[1];
sx q[1];
rz(-1.027193) q[1];
sx q[1];
rz(0.60740791) q[1];
rz(-2.1983002) q[3];
sx q[3];
rz(-1.2767013) q[3];
sx q[3];
rz(-2.7315745) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.9084106) q[2];
sx q[2];
rz(-0.57047129) q[2];
sx q[2];
rz(1.2916267) q[2];
rz(0.49606797) q[3];
sx q[3];
rz(-2.3521164) q[3];
sx q[3];
rz(1.9991649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9555776) q[0];
sx q[0];
rz(-2.8514974) q[0];
sx q[0];
rz(2.7225323) q[0];
rz(0.31173197) q[1];
sx q[1];
rz(-1.5429976) q[1];
sx q[1];
rz(2.9980803) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9559422) q[0];
sx q[0];
rz(-2.5505192) q[0];
sx q[0];
rz(2.4645564) q[0];
rz(0.2741998) q[2];
sx q[2];
rz(-1.7975217) q[2];
sx q[2];
rz(-2.271914) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.12998768) q[1];
sx q[1];
rz(-1.2530348) q[1];
sx q[1];
rz(2.3684273) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.26007248) q[3];
sx q[3];
rz(-1.7100157) q[3];
sx q[3];
rz(1.2633439) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.781337) q[2];
sx q[2];
rz(-2.2402253) q[2];
sx q[2];
rz(2.0085013) q[2];
rz(1.1059149) q[3];
sx q[3];
rz(-1.7224885) q[3];
sx q[3];
rz(-0.47479409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
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
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7853506) q[0];
sx q[0];
rz(-2.3478735) q[0];
sx q[0];
rz(1.6957977) q[0];
rz(-0.71714199) q[1];
sx q[1];
rz(-1.8627867) q[1];
sx q[1];
rz(0.76146567) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035485418) q[0];
sx q[0];
rz(-0.64752827) q[0];
sx q[0];
rz(-2.8962563) q[0];
rz(-0.47777678) q[2];
sx q[2];
rz(-2.226916) q[2];
sx q[2];
rz(2.853924) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.58757979) q[1];
sx q[1];
rz(-1.5214458) q[1];
sx q[1];
rz(0.99669902) q[1];
rz(1.3827356) q[3];
sx q[3];
rz(-1.4131964) q[3];
sx q[3];
rz(0.068305123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.7722499) q[2];
sx q[2];
rz(-0.49392924) q[2];
sx q[2];
rz(0.93144766) q[2];
rz(2.5233968) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(-2.1759694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44929993) q[0];
sx q[0];
rz(-1.7789142) q[0];
sx q[0];
rz(-1.3344673) q[0];
rz(-2.6284699) q[1];
sx q[1];
rz(-2.158439) q[1];
sx q[1];
rz(-1.9122874) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72237678) q[0];
sx q[0];
rz(-2.253719) q[0];
sx q[0];
rz(2.1238668) q[0];
rz(1.3111054) q[2];
sx q[2];
rz(-1.6605111) q[2];
sx q[2];
rz(-3.1083687) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.9521515) q[1];
sx q[1];
rz(-1.5501889) q[1];
sx q[1];
rz(2.9856696) q[1];
rz(-0.66219285) q[3];
sx q[3];
rz(-1.5846888) q[3];
sx q[3];
rz(-0.77869773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.3076155) q[2];
sx q[2];
rz(-0.5270842) q[2];
sx q[2];
rz(0.21044593) q[2];
rz(3.0757507) q[3];
sx q[3];
rz(-0.47271553) q[3];
sx q[3];
rz(0.79894799) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54302067) q[0];
sx q[0];
rz(-0.6544756) q[0];
sx q[0];
rz(1.2731113) q[0];
rz(-1.8674564) q[1];
sx q[1];
rz(-0.68398634) q[1];
sx q[1];
rz(0.88409105) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7403075) q[0];
sx q[0];
rz(-1.570556) q[0];
sx q[0];
rz(-0.9539286) q[0];
rz(-0.52427788) q[2];
sx q[2];
rz(-1.345621) q[2];
sx q[2];
rz(-0.88722992) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.6940004) q[1];
sx q[1];
rz(-1.847637) q[1];
sx q[1];
rz(-2.0029699) q[1];
x q[2];
rz(-1.8296914) q[3];
sx q[3];
rz(-1.2719385) q[3];
sx q[3];
rz(1.3204624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1511496) q[2];
sx q[2];
rz(-0.98439011) q[2];
sx q[2];
rz(1.3746064) q[2];
rz(0.19045842) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16354887) q[0];
sx q[0];
rz(-1.6885641) q[0];
sx q[0];
rz(-1.8769886) q[0];
rz(-3.0679852) q[1];
sx q[1];
rz(-2.059748) q[1];
sx q[1];
rz(-2.5216865) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.093599565) q[0];
sx q[0];
rz(-1.7353829) q[0];
sx q[0];
rz(-2.3175146) q[0];
rz(1.0629547) q[2];
sx q[2];
rz(-1.0217102) q[2];
sx q[2];
rz(2.6428403) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.1312841) q[1];
sx q[1];
rz(-0.29228739) q[1];
sx q[1];
rz(-0.75914219) q[1];
rz(-1.3315174) q[3];
sx q[3];
rz(-1.115723) q[3];
sx q[3];
rz(0.51067715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8668883) q[2];
sx q[2];
rz(-0.70211774) q[2];
sx q[2];
rz(0.14599027) q[2];
rz(2.919096) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(2.4116624) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.20919007) q[0];
sx q[0];
rz(-1.1840191) q[0];
sx q[0];
rz(-2.9970072) q[0];
rz(1.9119541) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(1.7894921) q[2];
sx q[2];
rz(-0.98697294) q[2];
sx q[2];
rz(1.0728991) q[2];
rz(-1.7207809) q[3];
sx q[3];
rz(-2.9290028) q[3];
sx q[3];
rz(1.3976189) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
