OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.2405038) q[0];
sx q[0];
rz(3.3189964) q[0];
sx q[0];
rz(11.431974) q[0];
rz(1.1881243) q[1];
sx q[1];
rz(-2.1048574) q[1];
sx q[1];
rz(-0.66361767) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21270574) q[0];
sx q[0];
rz(-1.1097254) q[0];
sx q[0];
rz(1.5443718) q[0];
rz(2.3637949) q[2];
sx q[2];
rz(-1.3133089) q[2];
sx q[2];
rz(-0.68629005) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9259778) q[1];
sx q[1];
rz(-1.3091334) q[1];
sx q[1];
rz(1.2222626) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8144242) q[3];
sx q[3];
rz(-1.4654136) q[3];
sx q[3];
rz(-1.6107314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2628281) q[2];
sx q[2];
rz(-2.6972013) q[2];
sx q[2];
rz(3.0905241) q[2];
rz(-0.55705327) q[3];
sx q[3];
rz(-2.3414108) q[3];
sx q[3];
rz(1.5548271) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.58650815) q[0];
sx q[0];
rz(-2.3095135) q[0];
sx q[0];
rz(2.5449975) q[0];
rz(-2.3157628) q[1];
sx q[1];
rz(-1.700371) q[1];
sx q[1];
rz(-1.2260431) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7694089) q[0];
sx q[0];
rz(-2.6499977) q[0];
sx q[0];
rz(0.43222897) q[0];
x q[1];
rz(-2.366757) q[2];
sx q[2];
rz(-0.92445395) q[2];
sx q[2];
rz(1.6006084) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.023149816) q[1];
sx q[1];
rz(-2.9217302) q[1];
sx q[1];
rz(1.5727732) q[1];
x q[2];
rz(2.1434104) q[3];
sx q[3];
rz(-1.5787573) q[3];
sx q[3];
rz(-1.5496467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.79919672) q[2];
sx q[2];
rz(-1.1725972) q[2];
sx q[2];
rz(2.8105695) q[2];
rz(-0.80667574) q[3];
sx q[3];
rz(-1.2253864) q[3];
sx q[3];
rz(1.7139009) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1598635) q[0];
sx q[0];
rz(-1.9112497) q[0];
sx q[0];
rz(1.8925517) q[0];
rz(-0.088009134) q[1];
sx q[1];
rz(-2.0188589) q[1];
sx q[1];
rz(-1.0294611) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35917657) q[0];
sx q[0];
rz(-2.4826907) q[0];
sx q[0];
rz(-1.2078148) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2842032) q[2];
sx q[2];
rz(-0.9409875) q[2];
sx q[2];
rz(1.0857925) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.50476915) q[1];
sx q[1];
rz(-1.4154589) q[1];
sx q[1];
rz(1.7117281) q[1];
rz(1.449552) q[3];
sx q[3];
rz(-0.83449927) q[3];
sx q[3];
rz(0.27437011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.133698) q[2];
sx q[2];
rz(-1.7241314) q[2];
sx q[2];
rz(-2.6339445) q[2];
rz(1.7525904) q[3];
sx q[3];
rz(-0.32998431) q[3];
sx q[3];
rz(1.1631789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4841109) q[0];
sx q[0];
rz(-0.27814516) q[0];
sx q[0];
rz(-1.5456276) q[0];
rz(2.0987299) q[1];
sx q[1];
rz(-1.9680126) q[1];
sx q[1];
rz(-1.5159336) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31947485) q[0];
sx q[0];
rz(-1.428077) q[0];
sx q[0];
rz(2.4583754) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.707294) q[2];
sx q[2];
rz(-2.578306) q[2];
sx q[2];
rz(1.0922645) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.2559291) q[1];
sx q[1];
rz(-0.71223488) q[1];
sx q[1];
rz(2.7182012) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7901778) q[3];
sx q[3];
rz(-0.99321584) q[3];
sx q[3];
rz(-2.1882309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.76379124) q[2];
sx q[2];
rz(-2.2968473) q[2];
sx q[2];
rz(-1.2949004) q[2];
rz(-0.14136782) q[3];
sx q[3];
rz(-0.54261345) q[3];
sx q[3];
rz(0.98658371) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35448733) q[0];
sx q[0];
rz(-1.961345) q[0];
sx q[0];
rz(1.6217344) q[0];
rz(-0.63201085) q[1];
sx q[1];
rz(-2.4193587) q[1];
sx q[1];
rz(-2.246726) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8653523) q[0];
sx q[0];
rz(-1.9812752) q[0];
sx q[0];
rz(-0.68666896) q[0];
rz(-0.61770265) q[2];
sx q[2];
rz(-0.82759826) q[2];
sx q[2];
rz(-2.655381) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.27069651) q[1];
sx q[1];
rz(-2.0082698) q[1];
sx q[1];
rz(-0.84769627) q[1];
rz(-2.8551293) q[3];
sx q[3];
rz(-1.5781919) q[3];
sx q[3];
rz(1.8252107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.616509) q[2];
sx q[2];
rz(-2.775165) q[2];
sx q[2];
rz(-1.1425225) q[2];
rz(2.3948495) q[3];
sx q[3];
rz(-1.2544422) q[3];
sx q[3];
rz(-1.07553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.557945) q[0];
sx q[0];
rz(-2.8201411) q[0];
sx q[0];
rz(1.7640132) q[0];
rz(2.6691943) q[1];
sx q[1];
rz(-2.6230085) q[1];
sx q[1];
rz(-0.46498743) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0156292) q[0];
sx q[0];
rz(-1.1383801) q[0];
sx q[0];
rz(0.36884357) q[0];
x q[1];
rz(2.9372413) q[2];
sx q[2];
rz(-0.34826476) q[2];
sx q[2];
rz(-1.6153) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.28593674) q[1];
sx q[1];
rz(-2.2054513) q[1];
sx q[1];
rz(0.059307701) q[1];
rz(-pi) q[2];
rz(1.3607849) q[3];
sx q[3];
rz(-0.87267733) q[3];
sx q[3];
rz(2.7426646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.49089367) q[2];
sx q[2];
rz(-1.8785672) q[2];
sx q[2];
rz(-2.6965551) q[2];
rz(-2.2079091) q[3];
sx q[3];
rz(-1.4327587) q[3];
sx q[3];
rz(2.8745108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0141107) q[0];
sx q[0];
rz(-1.5662136) q[0];
sx q[0];
rz(-1.4200462) q[0];
rz(0.02380112) q[1];
sx q[1];
rz(-2.5271466) q[1];
sx q[1];
rz(-2.9856317) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10931817) q[0];
sx q[0];
rz(-0.30800691) q[0];
sx q[0];
rz(0.59662915) q[0];
x q[1];
rz(0.51257001) q[2];
sx q[2];
rz(-1.1935496) q[2];
sx q[2];
rz(-0.23114983) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.8487726) q[1];
sx q[1];
rz(-0.18310586) q[1];
sx q[1];
rz(-1.8715026) q[1];
x q[2];
rz(2.6050623) q[3];
sx q[3];
rz(-2.0091669) q[3];
sx q[3];
rz(2.8560672) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.069313958) q[2];
sx q[2];
rz(-2.5645655) q[2];
sx q[2];
rz(2.0689266) q[2];
rz(-0.32564751) q[3];
sx q[3];
rz(-1.9653392) q[3];
sx q[3];
rz(-1.6252888) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1919365) q[0];
sx q[0];
rz(-2.7392116) q[0];
sx q[0];
rz(2.8038213) q[0];
rz(2.0514964) q[1];
sx q[1];
rz(-0.97507674) q[1];
sx q[1];
rz(2.8930194) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0524806) q[0];
sx q[0];
rz(-2.5626474) q[0];
sx q[0];
rz(-2.6738033) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48034251) q[2];
sx q[2];
rz(-2.0311653) q[2];
sx q[2];
rz(0.72887052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.92330248) q[1];
sx q[1];
rz(-2.3308838) q[1];
sx q[1];
rz(-2.2680125) q[1];
rz(1.9931273) q[3];
sx q[3];
rz(-0.18581192) q[3];
sx q[3];
rz(1.2687792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.2934072) q[2];
sx q[2];
rz(-0.53331393) q[2];
sx q[2];
rz(3.0548813) q[2];
rz(2.6596206) q[3];
sx q[3];
rz(-0.50589365) q[3];
sx q[3];
rz(-1.988525) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50865737) q[0];
sx q[0];
rz(-1.0077227) q[0];
sx q[0];
rz(-0.28276643) q[0];
rz(0.70156082) q[1];
sx q[1];
rz(-0.82156721) q[1];
sx q[1];
rz(-1.3185906) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0193664) q[0];
sx q[0];
rz(-1.1302233) q[0];
sx q[0];
rz(0.088935436) q[0];
rz(-pi) q[1];
rz(2.4852072) q[2];
sx q[2];
rz(-1.0644541) q[2];
sx q[2];
rz(0.93751794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9642155) q[1];
sx q[1];
rz(-2.4417158) q[1];
sx q[1];
rz(1.7978976) q[1];
rz(-pi) q[2];
x q[2];
rz(0.41140822) q[3];
sx q[3];
rz(-0.79380006) q[3];
sx q[3];
rz(2.0421162) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.41137722) q[2];
sx q[2];
rz(-0.29835478) q[2];
sx q[2];
rz(0.39917699) q[2];
rz(-2.2579851) q[3];
sx q[3];
rz(-1.6688321) q[3];
sx q[3];
rz(-1.2214899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3783962) q[0];
sx q[0];
rz(-0.77532399) q[0];
sx q[0];
rz(1.5989074) q[0];
rz(2.0762766) q[1];
sx q[1];
rz(-0.97384802) q[1];
sx q[1];
rz(1.8803966) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6529918) q[0];
sx q[0];
rz(-1.4654667) q[0];
sx q[0];
rz(2.7306042) q[0];
rz(-pi) q[1];
rz(1.5930428) q[2];
sx q[2];
rz(-0.97264475) q[2];
sx q[2];
rz(-1.0215789) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.65834261) q[1];
sx q[1];
rz(-1.634942) q[1];
sx q[1];
rz(-0.49883962) q[1];
rz(-pi) q[2];
rz(0.68393771) q[3];
sx q[3];
rz(-1.4468907) q[3];
sx q[3];
rz(0.77139664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5836872) q[2];
sx q[2];
rz(-2.5186899) q[2];
sx q[2];
rz(2.5718001) q[2];
rz(1.2184881) q[3];
sx q[3];
rz(-2.4945538) q[3];
sx q[3];
rz(2.5789554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31310836) q[0];
sx q[0];
rz(-2.2611571) q[0];
sx q[0];
rz(-1.8631998) q[0];
rz(2.5333511) q[1];
sx q[1];
rz(-0.47641644) q[1];
sx q[1];
rz(0.48412916) q[1];
rz(3.1153395) q[2];
sx q[2];
rz(-0.99925169) q[2];
sx q[2];
rz(3.0796438) q[2];
rz(-2.9363587) q[3];
sx q[3];
rz(-2.5406465) q[3];
sx q[3];
rz(3.061486) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
