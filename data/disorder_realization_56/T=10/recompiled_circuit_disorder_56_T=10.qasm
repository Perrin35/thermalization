OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg meas[4];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.38874414563179) q[0];
sx q[0];
rz(3.67787704070146) q[0];
sx q[0];
rz(10.3725467681806) q[0];
rz(-1.32879686355591) q[1];
sx q[1];
rz(4.40899446805055) q[1];
sx q[1];
rz(10.4525273799817) q[1];
cx q[1],q[0];
rz(2.00048136711121) q[0];
sx q[0];
rz(1.94654646714265) q[0];
sx q[0];
rz(12.541833138458) q[0];
rz(-0.489333122968674) q[2];
sx q[2];
rz(4.51641193230683) q[2];
sx q[2];
rz(9.08422148822948) q[2];
cx q[2],q[1];
rz(-0.487840503454208) q[1];
sx q[1];
rz(2.23347053130204) q[1];
sx q[1];
rz(9.91878641246959) q[1];
rz(-2.15604376792908) q[3];
sx q[3];
rz(2.42607763608033) q[3];
sx q[3];
rz(10.5948987960736) q[3];
cx q[3],q[2];
rz(-1.12969219684601) q[2];
sx q[2];
rz(7.99017110665376) q[2];
sx q[2];
rz(8.37449333666965) q[2];
rz(2.02832794189453) q[3];
sx q[3];
rz(5.3914751132303) q[3];
sx q[3];
rz(9.35667095928594) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-0.0724090784788132) q[0];
sx q[0];
rz(1.24531117280061) q[0];
sx q[0];
rz(12.2686178445737) q[0];
rz(0.619666635990143) q[1];
sx q[1];
rz(4.14877346356446) q[1];
sx q[1];
rz(14.5998682737271) q[1];
cx q[1],q[0];
rz(0.715826869010925) q[0];
sx q[0];
rz(5.32154146035249) q[0];
sx q[0];
rz(11.2513293981473) q[0];
rz(-2.53331875801086) q[2];
sx q[2];
rz(4.35246685345704) q[2];
sx q[2];
rz(8.69665620326206) q[2];
cx q[2],q[1];
rz(3.12167096138) q[1];
sx q[1];
rz(2.60560754139955) q[1];
sx q[1];
rz(10.2517261266629) q[1];
rz(-0.868222653865814) q[3];
sx q[3];
rz(1.36300233204896) q[3];
sx q[3];
rz(9.03214181064769) q[3];
cx q[3],q[2];
rz(0.462532818317413) q[2];
sx q[2];
rz(2.16475370724732) q[2];
sx q[2];
rz(11.5910725355069) q[2];
rz(2.22359275817871) q[3];
sx q[3];
rz(1.8564329465204) q[3];
sx q[3];
rz(6.57937452792331) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.17920303344727) q[0];
sx q[0];
rz(2.29008916218812) q[0];
sx q[0];
rz(8.88186774253055) q[0];
rz(-2.25935983657837) q[1];
sx q[1];
rz(4.27692464192445) q[1];
sx q[1];
rz(11.6015224218289) q[1];
cx q[1],q[0];
rz(-0.0977926254272461) q[0];
sx q[0];
rz(1.85051968892152) q[0];
sx q[0];
rz(13.8022479772489) q[0];
rz(1.27386081218719) q[2];
sx q[2];
rz(1.61162165005738) q[2];
sx q[2];
rz(9.64097709058925) q[2];
cx q[2],q[1];
rz(-4.25449562072754) q[1];
sx q[1];
rz(-1.635105458898) q[1];
sx q[1];
rz(15.5865673780362) q[1];
rz(-0.425905674695969) q[3];
sx q[3];
rz(0.452554615336009) q[3];
sx q[3];
rz(11.3602392434995) q[3];
cx q[3],q[2];
rz(-3.04700112342834) q[2];
sx q[2];
rz(3.75244483550126) q[2];
sx q[2];
rz(8.29163739680454) q[2];
rz(2.90996932983398) q[3];
sx q[3];
rz(4.41461315949494) q[3];
sx q[3];
rz(8.66750732659503) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.421067237854) q[0];
sx q[0];
rz(3.13114875846589) q[0];
sx q[0];
rz(11.1898550748746) q[0];
rz(0.518494129180908) q[1];
sx q[1];
rz(5.01876750786836) q[1];
sx q[1];
rz(9.18265501259967) q[1];
cx q[1],q[0];
rz(0.748944163322449) q[0];
sx q[0];
rz(0.90396514733369) q[0];
sx q[0];
rz(12.3588008642118) q[0];
rz(3.15466904640198) q[2];
sx q[2];
rz(7.42438283761079) q[2];
sx q[2];
rz(10.4440825939099) q[2];
cx q[2],q[1];
rz(-5.94887733459473) q[1];
sx q[1];
rz(1.1242927630716) q[1];
sx q[1];
rz(9.61307098566695) q[1];
rz(0.141080677509308) q[3];
sx q[3];
rz(5.12670055230195) q[3];
sx q[3];
rz(12.2679598092954) q[3];
cx q[3],q[2];
rz(-3.1597855091095) q[2];
sx q[2];
rz(5.31454864342744) q[2];
sx q[2];
rz(9.51963113843604) q[2];
rz(-1.79939603805542) q[3];
sx q[3];
rz(4.53884509404237) q[3];
sx q[3];
rz(13.0054845571439) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.78385329246521) q[0];
sx q[0];
rz(1.8308323939615) q[0];
sx q[0];
rz(12.2055775880735) q[0];
rz(1.75332534313202) q[1];
sx q[1];
rz(4.95239082177217) q[1];
sx q[1];
rz(7.41777608393832) q[1];
cx q[1],q[0];
rz(-3.87730526924133) q[0];
sx q[0];
rz(7.23050466378266) q[0];
sx q[0];
rz(8.80316183566257) q[0];
rz(-5.85343837738037) q[2];
sx q[2];
rz(1.89987126191194) q[2];
sx q[2];
rz(10.2491805314939) q[2];
cx q[2],q[1];
rz(1.51350092887878) q[1];
sx q[1];
rz(2.66736319859559) q[1];
sx q[1];
rz(13.7908558607022) q[1];
rz(4.27559852600098) q[3];
sx q[3];
rz(4.63632133801515) q[3];
sx q[3];
rz(14.4738449811856) q[3];
cx q[3],q[2];
rz(5.27789688110352) q[2];
sx q[2];
rz(4.56239250500733) q[2];
sx q[2];
rz(2.5689382314603) q[2];
rz(4.07034921646118) q[3];
sx q[3];
rz(3.66321882803971) q[3];
sx q[3];
rz(8.25722453593417) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.81444334983826) q[0];
sx q[0];
rz(2.04540184338624) q[0];
sx q[0];
rz(8.17543814181491) q[0];
rz(1.91847395896912) q[1];
sx q[1];
rz(1.52531162102754) q[1];
sx q[1];
rz(11.4141380548398) q[1];
cx q[1],q[0];
rz(1.45706796646118) q[0];
sx q[0];
rz(-3.24764188925689) q[0];
sx q[0];
rz(9.54486064463063) q[0];
rz(0.41209352016449) q[2];
sx q[2];
rz(5.1788069327646) q[2];
sx q[2];
rz(8.99594775437518) q[2];
cx q[2],q[1];
rz(-1.8730947971344) q[1];
sx q[1];
rz(0.226462038355418) q[1];
sx q[1];
rz(10.5116475582044) q[1];
rz(0.860299825668335) q[3];
sx q[3];
rz(2.32642546494538) q[3];
sx q[3];
rz(10.6239551067273) q[3];
cx q[3],q[2];
rz(2.67270088195801) q[2];
sx q[2];
rz(4.35585442383821) q[2];
sx q[2];
rz(5.29222247599765) q[2];
rz(0.647836029529572) q[3];
sx q[3];
rz(7.24872270424897) q[3];
sx q[3];
rz(6.63072392939731) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(3.39954113960266) q[0];
sx q[0];
rz(6.51028195221955) q[0];
sx q[0];
rz(9.48706343620225) q[0];
rz(-6.46899080276489) q[1];
sx q[1];
rz(7.96807638009126) q[1];
sx q[1];
rz(6.6779618024747) q[1];
cx q[1],q[0];
rz(-0.617907226085663) q[0];
sx q[0];
rz(3.60804453690583) q[0];
sx q[0];
rz(11.5430979490201) q[0];
rz(1.19455564022064) q[2];
sx q[2];
rz(4.41240123112733) q[2];
sx q[2];
rz(12.2645506620328) q[2];
cx q[2],q[1];
rz(3.31932306289673) q[1];
sx q[1];
rz(4.28675976594026) q[1];
sx q[1];
rz(12.8383729219358) q[1];
rz(-1.41733705997467) q[3];
sx q[3];
rz(4.74891916115815) q[3];
sx q[3];
rz(12.5342261552732) q[3];
cx q[3],q[2];
rz(-5.63468933105469) q[2];
sx q[2];
rz(2.81153699954087) q[2];
sx q[2];
rz(15.4349398374478) q[2];
rz(1.83888828754425) q[3];
sx q[3];
rz(4.96989193757112) q[3];
sx q[3];
rz(12.2543291807096) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(0.397723138332367) q[0];
sx q[0];
rz(1.13451543648774) q[0];
sx q[0];
rz(10.7098887920301) q[0];
rz(-4.64311170578003) q[1];
sx q[1];
rz(8.03382793267305) q[1];
sx q[1];
rz(4.482294535629) q[1];
cx q[1],q[0];
rz(-3.33638715744019) q[0];
sx q[0];
rz(1.18009820778901) q[0];
sx q[0];
rz(8.20636901854678) q[0];
rz(-1.63166832923889) q[2];
sx q[2];
rz(6.97991576989228) q[2];
sx q[2];
rz(11.8074941396634) q[2];
cx q[2],q[1];
rz(5.0816011428833) q[1];
sx q[1];
rz(1.92189470131929) q[1];
sx q[1];
rz(9.81007648109599) q[1];
rz(-0.442052125930786) q[3];
sx q[3];
rz(7.56732621987397) q[3];
sx q[3];
rz(9.52118246852561) q[3];
cx q[3],q[2];
rz(1.03547072410583) q[2];
sx q[2];
rz(3.94806251128251) q[2];
sx q[2];
rz(11.5427653551023) q[2];
rz(0.184939548373222) q[3];
sx q[3];
rz(0.390263231592723) q[3];
sx q[3];
rz(12.324483370773) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-2.09693646430969) q[0];
sx q[0];
rz(2.14502951701219) q[0];
sx q[0];
rz(11.0460623264234) q[0];
rz(-2.81144905090332) q[1];
sx q[1];
rz(5.07548538048799) q[1];
sx q[1];
rz(14.8708295583646) q[1];
cx q[1],q[0];
rz(0.225682079792023) q[0];
sx q[0];
rz(4.51056209405) q[0];
sx q[0];
rz(8.54135987757846) q[0];
rz(0.296600937843323) q[2];
sx q[2];
rz(5.26857748826081) q[2];
sx q[2];
rz(7.38450214862033) q[2];
cx q[2],q[1];
rz(-2.37714624404907) q[1];
sx q[1];
rz(6.46810022194917) q[1];
sx q[1];
rz(8.82090500592395) q[1];
rz(4.42134046554565) q[3];
sx q[3];
rz(1.0927199443155) q[3];
sx q[3];
rz(7.21109817027255) q[3];
cx q[3],q[2];
rz(3.74283671379089) q[2];
sx q[2];
rz(4.22656813462312) q[2];
sx q[2];
rz(9.60440552829906) q[2];
rz(-2.14586973190308) q[3];
sx q[3];
rz(8.17600265343721) q[3];
sx q[3];
rz(10.7357370614926) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(2.37784647941589) q[0];
sx q[0];
rz(2.79599309166009) q[0];
sx q[0];
rz(8.36755571364566) q[0];
rz(-0.107488475739956) q[1];
sx q[1];
rz(5.02974596818025) q[1];
sx q[1];
rz(8.44479916094943) q[1];
cx q[1],q[0];
rz(1.55325925350189) q[0];
sx q[0];
rz(6.14757409890229) q[0];
sx q[0];
rz(9.03934792279407) q[0];
rz(2.08851480484009) q[2];
sx q[2];
rz(5.15483346779878) q[2];
sx q[2];
rz(18.595395064346) q[2];
cx q[2],q[1];
rz(0.926147818565369) q[1];
sx q[1];
rz(4.283151896792) q[1];
sx q[1];
rz(8.11703834532901) q[1];
rz(3.43411803245544) q[3];
sx q[3];
rz(5.94620529015596) q[3];
sx q[3];
rz(4.98042628764316) q[3];
cx q[3],q[2];
rz(-1.30486273765564) q[2];
sx q[2];
rz(7.48799339135224) q[2];
sx q[2];
rz(7.31098172663852) q[2];
rz(-4.52843141555786) q[3];
sx q[3];
rz(8.11567178566987) q[3];
sx q[3];
rz(5.99957058428928) q[3];
cx q[3],q[2];
cx q[2],q[1];
cx q[1],q[0];
rz(-1.63593125343323) q[0];
sx q[0];
rz(5.5951179583841) q[0];
sx q[0];
rz(7.21329972743198) q[0];
rz(pi/2) q[0];
sx q[0];
rz(3*pi/2) q[0];
sx q[0];
rz(5*pi/2) q[0];
rz(-1.34941327571869) q[1];
sx q[1];
rz(4.43259445031221) q[1];
sx q[1];
rz(13.7613439321439) q[1];
rz(pi/2) q[1];
sx q[1];
rz(3*pi/2) q[1];
sx q[1];
rz(5*pi/2) q[1];
rz(-3.28445053100586) q[2];
sx q[2];
rz(4.23937812645967) q[2];
sx q[2];
rz(10.0544916152875) q[2];
rz(pi/2) q[2];
sx q[2];
rz(3*pi/2) q[2];
sx q[2];
rz(5*pi/2) q[2];
rz(-2.1909499168396) q[3];
sx q[3];
rz(7.64515176613862) q[3];
sx q[3];
rz(16.7159328222196) q[3];
rz(pi/2) q[3];
sx q[3];
rz(3*pi/2) q[3];
sx q[3];
rz(5*pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> meas[0];
measure q[1] -> meas[1];
measure q[2] -> meas[2];
measure q[3] -> meas[3];
