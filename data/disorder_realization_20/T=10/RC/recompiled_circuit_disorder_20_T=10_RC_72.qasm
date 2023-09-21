OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.5469172) q[0];
sx q[0];
rz(4.1424799) q[0];
sx q[0];
rz(9.6371798) q[0];
rz(-2.4266333) q[1];
sx q[1];
rz(-0.7874878) q[1];
sx q[1];
rz(1.8600872) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40385383) q[0];
sx q[0];
rz(-1.7716584) q[0];
sx q[0];
rz(-1.9181197) q[0];
rz(1.6002866) q[2];
sx q[2];
rz(-0.87848308) q[2];
sx q[2];
rz(-0.32683795) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.9077397) q[1];
sx q[1];
rz(-2.6820161) q[1];
sx q[1];
rz(1.185226) q[1];
x q[2];
rz(-3.078457) q[3];
sx q[3];
rz(-1.6992484) q[3];
sx q[3];
rz(-2.0383143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.32221258) q[2];
sx q[2];
rz(-0.11495049) q[2];
sx q[2];
rz(1.6248576) q[2];
rz(2.8939698) q[3];
sx q[3];
rz(-1.6678436) q[3];
sx q[3];
rz(-0.65555278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4441929) q[0];
sx q[0];
rz(-1.6635165) q[0];
sx q[0];
rz(0.15629388) q[0];
rz(2.5022751) q[1];
sx q[1];
rz(-1.0126065) q[1];
sx q[1];
rz(-1.4888391) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1600906) q[0];
sx q[0];
rz(-1.0233876) q[0];
sx q[0];
rz(1.7574969) q[0];
rz(-2.4376051) q[2];
sx q[2];
rz(-2.1805602) q[2];
sx q[2];
rz(-3.0733382) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.21178791) q[1];
sx q[1];
rz(-2.1165407) q[1];
sx q[1];
rz(0.8916698) q[1];
rz(-2.923521) q[3];
sx q[3];
rz(-3.0017188) q[3];
sx q[3];
rz(-1.1034031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2911825) q[2];
sx q[2];
rz(-3.0427142) q[2];
sx q[2];
rz(-0.29671159) q[2];
rz(-0.3324278) q[3];
sx q[3];
rz(-1.5184831) q[3];
sx q[3];
rz(1.9077574) q[3];
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
rz(1.6333703) q[0];
sx q[0];
rz(-1.8828266) q[0];
sx q[0];
rz(0.53031522) q[0];
rz(-0.5439533) q[1];
sx q[1];
rz(-2.0201611) q[1];
sx q[1];
rz(1.189032) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8231359) q[0];
sx q[0];
rz(-1.728115) q[0];
sx q[0];
rz(2.8419028) q[0];
x q[1];
rz(-2.8934115) q[2];
sx q[2];
rz(-1.5319053) q[2];
sx q[2];
rz(1.789202) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.5570453) q[1];
sx q[1];
rz(-1.0439596) q[1];
sx q[1];
rz(0.79750632) q[1];
rz(-pi) q[2];
rz(0.53718209) q[3];
sx q[3];
rz(-0.77386412) q[3];
sx q[3];
rz(-2.8230132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.236078) q[2];
sx q[2];
rz(-1.4730075) q[2];
sx q[2];
rz(-1.2515602) q[2];
rz(-2.6873612) q[3];
sx q[3];
rz(-2.7082704) q[3];
sx q[3];
rz(0.35513487) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.77804756) q[0];
sx q[0];
rz(-0.30968928) q[0];
sx q[0];
rz(-1.3409412) q[0];
rz(-0.60607934) q[1];
sx q[1];
rz(-0.93901912) q[1];
sx q[1];
rz(1.9569424) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.143648) q[0];
sx q[0];
rz(-2.8638683) q[0];
sx q[0];
rz(2.4585637) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0834288) q[2];
sx q[2];
rz(-1.9789654) q[2];
sx q[2];
rz(-2.4525016) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9222316) q[1];
sx q[1];
rz(-2.1582099) q[1];
sx q[1];
rz(-3.0012793) q[1];
rz(-2.437856) q[3];
sx q[3];
rz(-2.6061213) q[3];
sx q[3];
rz(0.4243917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.23345315) q[2];
sx q[2];
rz(-0.88580695) q[2];
sx q[2];
rz(-2.502029) q[2];
rz(-0.8574287) q[3];
sx q[3];
rz(-1.9968417) q[3];
sx q[3];
rz(2.6517984) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.63067591) q[0];
sx q[0];
rz(-0.54172051) q[0];
sx q[0];
rz(2.6211991) q[0];
rz(-0.10294542) q[1];
sx q[1];
rz(-2.091145) q[1];
sx q[1];
rz(-2.4127814) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2136912) q[0];
sx q[0];
rz(-1.3662845) q[0];
sx q[0];
rz(-1.6708899) q[0];
x q[1];
rz(-1.7545627) q[2];
sx q[2];
rz(-0.73256058) q[2];
sx q[2];
rz(-1.5866605) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0923857) q[1];
sx q[1];
rz(-1.2485463) q[1];
sx q[1];
rz(-2.9258123) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.8647996) q[3];
sx q[3];
rz(-2.3949361) q[3];
sx q[3];
rz(-1.696123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0985428) q[2];
sx q[2];
rz(-2.7816732) q[2];
sx q[2];
rz(-2.4318802) q[2];
rz(0.63306159) q[3];
sx q[3];
rz(-1.148162) q[3];
sx q[3];
rz(0.13171296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1028041) q[0];
sx q[0];
rz(-3.0392201) q[0];
sx q[0];
rz(0.88371712) q[0];
rz(-0.9206413) q[1];
sx q[1];
rz(-1.0620774) q[1];
sx q[1];
rz(-0.19764915) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5699428) q[0];
sx q[0];
rz(-2.2075704) q[0];
sx q[0];
rz(2.3417579) q[0];
rz(0.10651393) q[2];
sx q[2];
rz(-1.649756) q[2];
sx q[2];
rz(-0.67895141) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.16285322) q[1];
sx q[1];
rz(-1.1161242) q[1];
sx q[1];
rz(-1.5233558) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2570409) q[3];
sx q[3];
rz(-2.2042847) q[3];
sx q[3];
rz(1.3214878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(3.1137696) q[2];
sx q[2];
rz(-1.7980857) q[2];
sx q[2];
rz(-3.0659884) q[2];
rz(-1.6522225) q[3];
sx q[3];
rz(-2.7421156) q[3];
sx q[3];
rz(2.3125125) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.41912115) q[0];
sx q[0];
rz(-0.88547456) q[0];
sx q[0];
rz(-0.042073123) q[0];
rz(-1.178859) q[1];
sx q[1];
rz(-1.983164) q[1];
sx q[1];
rz(1.1694318) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.021335) q[0];
sx q[0];
rz(-1.7354256) q[0];
sx q[0];
rz(3.0740601) q[0];
rz(-pi) q[1];
rz(-0.39016907) q[2];
sx q[2];
rz(-1.3069659) q[2];
sx q[2];
rz(2.2389776) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.52160727) q[1];
sx q[1];
rz(-2.1072525) q[1];
sx q[1];
rz(-1.1779838) q[1];
rz(0.38798214) q[3];
sx q[3];
rz(-1.531732) q[3];
sx q[3];
rz(-2.8475873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.823267) q[2];
sx q[2];
rz(-1.1151168) q[2];
sx q[2];
rz(-0.014483359) q[2];
rz(-1.5197586) q[3];
sx q[3];
rz(-1.8448011) q[3];
sx q[3];
rz(2.5884132) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13977215) q[0];
sx q[0];
rz(-2.791239) q[0];
sx q[0];
rz(-0.56234223) q[0];
rz(-1.7100122) q[1];
sx q[1];
rz(-1.7446012) q[1];
sx q[1];
rz(0.45773488) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2248174) q[0];
sx q[0];
rz(-0.48368925) q[0];
sx q[0];
rz(-1.7513357) q[0];
rz(-pi) q[1];
rz(-0.87403239) q[2];
sx q[2];
rz(-1.1296141) q[2];
sx q[2];
rz(0.87768427) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8843958) q[1];
sx q[1];
rz(-2.307502) q[1];
sx q[1];
rz(-0.019330545) q[1];
rz(-1.0192972) q[3];
sx q[3];
rz(-0.5987474) q[3];
sx q[3];
rz(-2.7703354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.83595014) q[2];
sx q[2];
rz(-0.93985158) q[2];
sx q[2];
rz(2.8590554) q[2];
rz(0.95747581) q[3];
sx q[3];
rz(-1.7493533) q[3];
sx q[3];
rz(0.47874513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1829421) q[0];
sx q[0];
rz(-0.68529737) q[0];
sx q[0];
rz(-1.4022934) q[0];
rz(0.69933403) q[1];
sx q[1];
rz(-1.3032841) q[1];
sx q[1];
rz(-1.8797849) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5926873) q[0];
sx q[0];
rz(-1.3015916) q[0];
sx q[0];
rz(2.6972428) q[0];
rz(-pi) q[1];
rz(-2.8724573) q[2];
sx q[2];
rz(-2.4705774) q[2];
sx q[2];
rz(-0.035949635) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0314434) q[1];
sx q[1];
rz(-1.3182536) q[1];
sx q[1];
rz(-2.4367711) q[1];
rz(-0.34187596) q[3];
sx q[3];
rz(-2.5856527) q[3];
sx q[3];
rz(-0.13560175) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4385779) q[2];
sx q[2];
rz(-1.5038749) q[2];
sx q[2];
rz(1.4578488) q[2];
rz(2.4750989) q[3];
sx q[3];
rz(-1.4278744) q[3];
sx q[3];
rz(0.75171793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2465729) q[0];
sx q[0];
rz(-0.81037766) q[0];
sx q[0];
rz(-1.9816459) q[0];
rz(0.054140422) q[1];
sx q[1];
rz(-1.663798) q[1];
sx q[1];
rz(1.0704401) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0265855) q[0];
sx q[0];
rz(-0.75782776) q[0];
sx q[0];
rz(1.3269781) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.530933) q[2];
sx q[2];
rz(-0.76518744) q[2];
sx q[2];
rz(2.6305692) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4813303) q[1];
sx q[1];
rz(-2.3412625) q[1];
sx q[1];
rz(-0.34923133) q[1];
x q[2];
rz(0.6108547) q[3];
sx q[3];
rz(-0.29963747) q[3];
sx q[3];
rz(2.6328997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.76876172) q[2];
sx q[2];
rz(-0.54905811) q[2];
sx q[2];
rz(1.5268415) q[2];
rz(-2.5620143) q[3];
sx q[3];
rz(-1.7947581) q[3];
sx q[3];
rz(0.65210623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.51331818) q[0];
sx q[0];
rz(-1.5341298) q[0];
sx q[0];
rz(-2.4232724) q[0];
rz(-2.1075481) q[1];
sx q[1];
rz(-2.4517192) q[1];
sx q[1];
rz(2.413961) q[1];
rz(0.037739567) q[2];
sx q[2];
rz(-0.83913091) q[2];
sx q[2];
rz(-0.057374949) q[2];
rz(-0.53229971) q[3];
sx q[3];
rz(-1.1333864) q[3];
sx q[3];
rz(-0.65896853) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
