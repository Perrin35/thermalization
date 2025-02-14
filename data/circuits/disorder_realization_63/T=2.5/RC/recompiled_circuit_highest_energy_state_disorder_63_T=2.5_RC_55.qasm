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
rz(1.0032049) q[0];
sx q[0];
rz(-1.2004852) q[0];
sx q[0];
rz(1.0263654) q[0];
rz(0.29844555) q[1];
sx q[1];
rz(4.65) q[1];
sx q[1];
rz(10.427373) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5110785) q[0];
sx q[0];
rz(-1.3670232) q[0];
sx q[0];
rz(3.1183) q[0];
rz(-pi) q[1];
x q[1];
rz(0.070656405) q[2];
sx q[2];
rz(-2.4774654) q[2];
sx q[2];
rz(0.16985591) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.5693869) q[1];
sx q[1];
rz(-1.1841229) q[1];
sx q[1];
rz(1.4264849) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.63319541) q[3];
sx q[3];
rz(-2.0811307) q[3];
sx q[3];
rz(1.1309433) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.58874321) q[2];
sx q[2];
rz(-2.1407949) q[2];
sx q[2];
rz(-0.0351077) q[2];
rz(-1.9378174) q[3];
sx q[3];
rz(-1.8193918) q[3];
sx q[3];
rz(2.8573341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.069000706) q[0];
sx q[0];
rz(-0.44206107) q[0];
sx q[0];
rz(2.1373855) q[0];
rz(-0.12241157) q[1];
sx q[1];
rz(-1.4675843) q[1];
sx q[1];
rz(-2.4845128) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1998889) q[0];
sx q[0];
rz(-0.39638576) q[0];
sx q[0];
rz(2.928171) q[0];
x q[1];
rz(-1.8481215) q[2];
sx q[2];
rz(-1.4155626) q[2];
sx q[2];
rz(-0.6065795) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2942336) q[1];
sx q[1];
rz(-0.93437934) q[1];
sx q[1];
rz(-0.60346236) q[1];
rz(0.78389726) q[3];
sx q[3];
rz(-1.3957033) q[3];
sx q[3];
rz(2.8181048) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5848026) q[2];
sx q[2];
rz(-0.57791296) q[2];
sx q[2];
rz(-2.54971) q[2];
rz(-1.4764192) q[3];
sx q[3];
rz(-1.5341325) q[3];
sx q[3];
rz(-0.8736476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0037435) q[0];
sx q[0];
rz(-0.20895222) q[0];
sx q[0];
rz(1.9269706) q[0];
rz(2.1851723) q[1];
sx q[1];
rz(-1.950187) q[1];
sx q[1];
rz(1.9699684) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7192243) q[0];
sx q[0];
rz(-0.63084334) q[0];
sx q[0];
rz(2.2792086) q[0];
x q[1];
rz(-1.3287622) q[2];
sx q[2];
rz(-0.71513745) q[2];
sx q[2];
rz(0.29235833) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.0846935) q[1];
sx q[1];
rz(-2.1325743) q[1];
sx q[1];
rz(-1.9719719) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.615715) q[3];
sx q[3];
rz(-1.6605596) q[3];
sx q[3];
rz(-0.95960754) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9413635) q[2];
sx q[2];
rz(-1.0503294) q[2];
sx q[2];
rz(1.7639147) q[2];
rz(2.7348943) q[3];
sx q[3];
rz(-2.4923057) q[3];
sx q[3];
rz(-2.6810834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7895301) q[0];
sx q[0];
rz(-1.2593513) q[0];
sx q[0];
rz(-1.9011185) q[0];
rz(-0.71818304) q[1];
sx q[1];
rz(-1.3259462) q[1];
sx q[1];
rz(2.9026418) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5137635) q[0];
sx q[0];
rz(-1.2830334) q[0];
sx q[0];
rz(-1.6970558) q[0];
rz(2.3283655) q[2];
sx q[2];
rz(-1.0095749) q[2];
sx q[2];
rz(-0.49408572) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1526665) q[1];
sx q[1];
rz(-1.7133022) q[1];
sx q[1];
rz(1.8320384) q[1];
x q[2];
rz(-1.7051058) q[3];
sx q[3];
rz(-2.1013936) q[3];
sx q[3];
rz(-2.1710896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.51603985) q[2];
sx q[2];
rz(-1.3935139) q[2];
sx q[2];
rz(-1.7986521) q[2];
rz(-0.88137734) q[3];
sx q[3];
rz(-0.16888976) q[3];
sx q[3];
rz(-0.78099293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7278904) q[0];
sx q[0];
rz(-2.8975633) q[0];
sx q[0];
rz(1.6572886) q[0];
rz(0.65385747) q[1];
sx q[1];
rz(-1.6203208) q[1];
sx q[1];
rz(0.13793764) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1178095) q[0];
sx q[0];
rz(-1.9236757) q[0];
sx q[0];
rz(0.15858741) q[0];
rz(-pi) q[1];
rz(1.8191255) q[2];
sx q[2];
rz(-2.2086124) q[2];
sx q[2];
rz(1.0385823) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5522462) q[1];
sx q[1];
rz(-1.668743) q[1];
sx q[1];
rz(2.0022654) q[1];
rz(-pi) q[2];
x q[2];
rz(1.82965) q[3];
sx q[3];
rz(-1.5167243) q[3];
sx q[3];
rz(-1.3075706) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.58978927) q[2];
sx q[2];
rz(-1.9686331) q[2];
sx q[2];
rz(-3.0544082) q[2];
rz(-0.11463556) q[3];
sx q[3];
rz(-2.3264824) q[3];
sx q[3];
rz(-0.4172999) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
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
rz(-0.36599416) q[0];
sx q[0];
rz(-1.5394779) q[0];
sx q[0];
rz(0.48496801) q[0];
rz(-2.1980749) q[1];
sx q[1];
rz(-0.8129932) q[1];
sx q[1];
rz(1.9224723) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6935558) q[0];
sx q[0];
rz(-0.72759923) q[0];
sx q[0];
rz(1.7175402) q[0];
x q[1];
rz(-1.6197422) q[2];
sx q[2];
rz(-2.8192047) q[2];
sx q[2];
rz(-1.9677377) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.3644219) q[1];
sx q[1];
rz(-0.14816532) q[1];
sx q[1];
rz(-1.0309459) q[1];
x q[2];
rz(0.7597105) q[3];
sx q[3];
rz(-1.1897161) q[3];
sx q[3];
rz(0.71807968) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.9014827) q[2];
sx q[2];
rz(-2.0811452) q[2];
sx q[2];
rz(0.94129747) q[2];
rz(2.8955722) q[3];
sx q[3];
rz(-1.6451719) q[3];
sx q[3];
rz(1.7814319) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.098002382) q[0];
sx q[0];
rz(-1.101838) q[0];
sx q[0];
rz(-1.9775506) q[0];
rz(1.3580492) q[1];
sx q[1];
rz(-2.2750504) q[1];
sx q[1];
rz(-1.39303) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9214894) q[0];
sx q[0];
rz(-1.1951726) q[0];
sx q[0];
rz(3.0620248) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.5921408) q[2];
sx q[2];
rz(-1.3260815) q[2];
sx q[2];
rz(-3.0915054) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1198123) q[1];
sx q[1];
rz(-2.0396898) q[1];
sx q[1];
rz(-1.6890668) q[1];
rz(0.89207666) q[3];
sx q[3];
rz(-1.3952655) q[3];
sx q[3];
rz(0.20993671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67036575) q[2];
sx q[2];
rz(-1.4281851) q[2];
sx q[2];
rz(2.8426389) q[2];
rz(2.7361338) q[3];
sx q[3];
rz(-0.36646989) q[3];
sx q[3];
rz(-0.56431842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7120755) q[0];
sx q[0];
rz(-2.8686664) q[0];
sx q[0];
rz(-0.78936973) q[0];
rz(2.7366267) q[1];
sx q[1];
rz(-1.3507495) q[1];
sx q[1];
rz(0.49496034) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5720737) q[0];
sx q[0];
rz(-2.1316307) q[0];
sx q[0];
rz(-2.5755519) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8414824) q[2];
sx q[2];
rz(-2.1202188) q[2];
sx q[2];
rz(-0.97893366) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1972754) q[1];
sx q[1];
rz(-0.8110356) q[1];
sx q[1];
rz(1.2326101) q[1];
rz(-2.6172753) q[3];
sx q[3];
rz(-2.2341407) q[3];
sx q[3];
rz(1.9085705) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.97883812) q[2];
sx q[2];
rz(-1.2252204) q[2];
sx q[2];
rz(-2.8863353) q[2];
rz(-3.1066331) q[3];
sx q[3];
rz(-1.6182263) q[3];
sx q[3];
rz(0.24622723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73664767) q[0];
sx q[0];
rz(-2.078779) q[0];
sx q[0];
rz(0.71003067) q[0];
rz(2.1404449) q[1];
sx q[1];
rz(-0.89718693) q[1];
sx q[1];
rz(2.5206916) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17711711) q[0];
sx q[0];
rz(-1.6311446) q[0];
sx q[0];
rz(1.2089157) q[0];
rz(-pi) q[1];
x q[1];
rz(0.65695564) q[2];
sx q[2];
rz(-1.959942) q[2];
sx q[2];
rz(2.3537113) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8709594) q[1];
sx q[1];
rz(-0.3644202) q[1];
sx q[1];
rz(-2.8363813) q[1];
rz(-1.848351) q[3];
sx q[3];
rz(-0.81349868) q[3];
sx q[3];
rz(-3.1111969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3322477) q[2];
sx q[2];
rz(-1.6184018) q[2];
sx q[2];
rz(0.28918949) q[2];
rz(2.0293763) q[3];
sx q[3];
rz(-2.8465392) q[3];
sx q[3];
rz(2.1212063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4777098) q[0];
sx q[0];
rz(-1.7998671) q[0];
sx q[0];
rz(2.2148602) q[0];
rz(2.0345188) q[1];
sx q[1];
rz(-1.5137545) q[1];
sx q[1];
rz(-0.56799299) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0999789) q[0];
sx q[0];
rz(-1.2547412) q[0];
sx q[0];
rz(0.63412447) q[0];
rz(2.6935518) q[2];
sx q[2];
rz(-0.5574286) q[2];
sx q[2];
rz(0.008226062) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.5957678) q[1];
sx q[1];
rz(-2.1819356) q[1];
sx q[1];
rz(2.6910438) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0509012) q[3];
sx q[3];
rz(-0.9970768) q[3];
sx q[3];
rz(-2.3815976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.052281436) q[2];
sx q[2];
rz(-0.81030446) q[2];
sx q[2];
rz(-1.4386162) q[2];
rz(0.29664052) q[3];
sx q[3];
rz(-2.7999925) q[3];
sx q[3];
rz(-2.6945485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.457837) q[0];
sx q[0];
rz(-1.0975657) q[0];
sx q[0];
rz(0.56957635) q[0];
rz(2.8029022) q[1];
sx q[1];
rz(-1.6117922) q[1];
sx q[1];
rz(1.502996) q[1];
rz(-2.8903932) q[2];
sx q[2];
rz(-1.770537) q[2];
sx q[2];
rz(0.58687576) q[2];
rz(-2.972907) q[3];
sx q[3];
rz(-1.2999693) q[3];
sx q[3];
rz(-0.6948605) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
