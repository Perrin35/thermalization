OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.8024017) q[0];
sx q[0];
rz(-1.985476) q[0];
sx q[0];
rz(2.6240786) q[0];
rz(1.327688) q[1];
sx q[1];
rz(7.1751243) q[1];
sx q[1];
rz(8.8827477) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9487614) q[0];
sx q[0];
rz(-1.7418712) q[0];
sx q[0];
rz(-0.77617742) q[0];
x q[1];
rz(2.3752604) q[2];
sx q[2];
rz(-1.7558985) q[2];
sx q[2];
rz(1.9232149) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.4229976) q[1];
sx q[1];
rz(-1.7684775) q[1];
sx q[1];
rz(2.343178) q[1];
rz(0.23302257) q[3];
sx q[3];
rz(-1.4015504) q[3];
sx q[3];
rz(-0.47692933) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.0077670495) q[2];
sx q[2];
rz(-2.0962891) q[2];
sx q[2];
rz(-2.2110151) q[2];
rz(-3.0104356) q[3];
sx q[3];
rz(-2.7635062) q[3];
sx q[3];
rz(-1.650943) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5432878) q[0];
sx q[0];
rz(-0.85619339) q[0];
sx q[0];
rz(-1.244586) q[0];
rz(2.2159131) q[1];
sx q[1];
rz(-1.8857748) q[1];
sx q[1];
rz(1.6474887) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5410575) q[0];
sx q[0];
rz(-0.78406683) q[0];
sx q[0];
rz(1.1009786) q[0];
x q[1];
rz(0.41183128) q[2];
sx q[2];
rz(-0.71626012) q[2];
sx q[2];
rz(1.0541944) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.0665613) q[1];
sx q[1];
rz(-1.9157456) q[1];
sx q[1];
rz(-0.48081545) q[1];
rz(-pi) q[2];
rz(-3.1140987) q[3];
sx q[3];
rz(-2.8373233) q[3];
sx q[3];
rz(-2.9472443) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1458448) q[2];
sx q[2];
rz(-2.6971942) q[2];
sx q[2];
rz(2.2399529) q[2];
rz(-1.3580648) q[3];
sx q[3];
rz(-1.8243022) q[3];
sx q[3];
rz(0.54734126) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3667592) q[0];
sx q[0];
rz(-1.9623373) q[0];
sx q[0];
rz(-1.9955848) q[0];
rz(2.215812) q[1];
sx q[1];
rz(-2.130276) q[1];
sx q[1];
rz(-0.19736966) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1903052) q[0];
sx q[0];
rz(-1.4337375) q[0];
sx q[0];
rz(0.97676116) q[0];
rz(-1.8756798) q[2];
sx q[2];
rz(-1.8499057) q[2];
sx q[2];
rz(-0.32372083) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0377522) q[1];
sx q[1];
rz(-1.3430129) q[1];
sx q[1];
rz(-0.040018602) q[1];
x q[2];
rz(0.54346187) q[3];
sx q[3];
rz(-1.3980143) q[3];
sx q[3];
rz(1.9256189) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.4027412) q[2];
sx q[2];
rz(-0.98029843) q[2];
sx q[2];
rz(-2.951238) q[2];
rz(-0.061577408) q[3];
sx q[3];
rz(-0.96934167) q[3];
sx q[3];
rz(-1.0524582) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15678081) q[0];
sx q[0];
rz(-1.4182014) q[0];
sx q[0];
rz(-2.5132827) q[0];
rz(-1.5208987) q[1];
sx q[1];
rz(-1.2077121) q[1];
sx q[1];
rz(2.3627088) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.223004) q[0];
sx q[0];
rz(-2.1445334) q[0];
sx q[0];
rz(2.6911435) q[0];
x q[1];
rz(2.7812296) q[2];
sx q[2];
rz(-1.833263) q[2];
sx q[2];
rz(-1.2059463) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4877704) q[1];
sx q[1];
rz(-1.5842591) q[1];
sx q[1];
rz(-3.0938992) q[1];
x q[2];
rz(-0.67653894) q[3];
sx q[3];
rz(-2.153844) q[3];
sx q[3];
rz(-2.9522715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6598307) q[2];
sx q[2];
rz(-2.063282) q[2];
sx q[2];
rz(-4.3241186e-05) q[2];
rz(1.0846042) q[3];
sx q[3];
rz(-1.1559887) q[3];
sx q[3];
rz(0.46561766) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4701009) q[0];
sx q[0];
rz(-1.4441613) q[0];
sx q[0];
rz(1.8430365) q[0];
rz(2.5647054) q[1];
sx q[1];
rz(-1.2737609) q[1];
sx q[1];
rz(-2.6947122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90335411) q[0];
sx q[0];
rz(-2.2132259) q[0];
sx q[0];
rz(-1.7967671) q[0];
rz(-1.8700897) q[2];
sx q[2];
rz(-1.6422049) q[2];
sx q[2];
rz(-1.4324607) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.943534) q[1];
sx q[1];
rz(-1.3683649) q[1];
sx q[1];
rz(-1.979503) q[1];
rz(-pi) q[2];
rz(0.3077363) q[3];
sx q[3];
rz(-2.947181) q[3];
sx q[3];
rz(2.3178063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.4345066) q[2];
sx q[2];
rz(-0.27927566) q[2];
sx q[2];
rz(0.2198098) q[2];
rz(1.6541727) q[3];
sx q[3];
rz(-1.6916964) q[3];
sx q[3];
rz(0.69948227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
rz(2.952482) q[0];
sx q[0];
rz(-2.7087152) q[0];
sx q[0];
rz(0.89749807) q[0];
rz(1.3709566) q[1];
sx q[1];
rz(-1.0400583) q[1];
sx q[1];
rz(-0.8955566) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33226704) q[0];
sx q[0];
rz(-2.4054962) q[0];
sx q[0];
rz(-2.9311137) q[0];
x q[1];
rz(2.1826971) q[2];
sx q[2];
rz(-2.09225) q[2];
sx q[2];
rz(0.34470169) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8297146) q[1];
sx q[1];
rz(-2.6304881) q[1];
sx q[1];
rz(2.4097777) q[1];
x q[2];
rz(-0.47116739) q[3];
sx q[3];
rz(-1.1352603) q[3];
sx q[3];
rz(-0.35126951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7728277) q[2];
sx q[2];
rz(-2.0389098) q[2];
sx q[2];
rz(-0.15711288) q[2];
rz(-2.469192) q[3];
sx q[3];
rz(-1.3998569) q[3];
sx q[3];
rz(3.0290643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2812578) q[0];
sx q[0];
rz(-0.66547886) q[0];
sx q[0];
rz(1.459664) q[0];
rz(0.36059391) q[1];
sx q[1];
rz(-1.4692042) q[1];
sx q[1];
rz(-1.2999387) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0361745) q[0];
sx q[0];
rz(-1.4486533) q[0];
sx q[0];
rz(2.7226177) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5333161) q[2];
sx q[2];
rz(-0.71675863) q[2];
sx q[2];
rz(0.83502095) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.056562076) q[1];
sx q[1];
rz(-1.6675648) q[1];
sx q[1];
rz(0.29778778) q[1];
rz(2.5266493) q[3];
sx q[3];
rz(-1.9087026) q[3];
sx q[3];
rz(-2.8470521) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2099057) q[2];
sx q[2];
rz(-1.0301215) q[2];
sx q[2];
rz(-2.9295861) q[2];
rz(2.0906406) q[3];
sx q[3];
rz(-1.7651599) q[3];
sx q[3];
rz(-3.0534548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6050922) q[0];
sx q[0];
rz(-0.78992805) q[0];
sx q[0];
rz(-2.9686046) q[0];
rz(2.0269003) q[1];
sx q[1];
rz(-2.098691) q[1];
sx q[1];
rz(3.0317422) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3152471) q[0];
sx q[0];
rz(-1.7638806) q[0];
sx q[0];
rz(-0.56435926) q[0];
rz(-pi) q[1];
rz(-2.8809261) q[2];
sx q[2];
rz(-0.83999604) q[2];
sx q[2];
rz(0.8651917) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.50625769) q[1];
sx q[1];
rz(-1.3999363) q[1];
sx q[1];
rz(-1.8299163) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9494826) q[3];
sx q[3];
rz(-0.56280901) q[3];
sx q[3];
rz(0.61432381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4032119) q[2];
sx q[2];
rz(-1.9387559) q[2];
sx q[2];
rz(-2.476725) q[2];
rz(-2.6730149) q[3];
sx q[3];
rz(-1.5237619) q[3];
sx q[3];
rz(-0.81336462) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.68798962) q[0];
sx q[0];
rz(-0.60992321) q[0];
sx q[0];
rz(-0.74158057) q[0];
rz(-1.954151) q[1];
sx q[1];
rz(-1.6148184) q[1];
sx q[1];
rz(0.56914079) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81108196) q[0];
sx q[0];
rz(-1.8966604) q[0];
sx q[0];
rz(1.261607) q[0];
x q[1];
rz(-2.9958565) q[2];
sx q[2];
rz(-2.6837641) q[2];
sx q[2];
rz(2.482058) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7267159) q[1];
sx q[1];
rz(-1.1099713) q[1];
sx q[1];
rz(-2.0445092) q[1];
rz(-2.3027116) q[3];
sx q[3];
rz(-1.0195707) q[3];
sx q[3];
rz(-1.6987926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.49447122) q[2];
sx q[2];
rz(-0.98860604) q[2];
sx q[2];
rz(-0.76623255) q[2];
rz(1.1726441) q[3];
sx q[3];
rz(-1.6021043) q[3];
sx q[3];
rz(0.22181454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
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
rz(0.93429339) q[0];
sx q[0];
rz(-0.3485637) q[0];
sx q[0];
rz(-0.19185129) q[0];
rz(2.3244997) q[1];
sx q[1];
rz(-2.8849738) q[1];
sx q[1];
rz(-2.4339035) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1908042) q[0];
sx q[0];
rz(-0.30256264) q[0];
sx q[0];
rz(-2.5214419) q[0];
rz(-pi) q[1];
x q[1];
rz(0.97918503) q[2];
sx q[2];
rz(-1.8509097) q[2];
sx q[2];
rz(1.4794939) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.62218828) q[1];
sx q[1];
rz(-0.71336245) q[1];
sx q[1];
rz(0.98825561) q[1];
x q[2];
rz(0.58587425) q[3];
sx q[3];
rz(-2.088024) q[3];
sx q[3];
rz(0.45211238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.91844687) q[2];
sx q[2];
rz(-0.94474363) q[2];
sx q[2];
rz(0.72500149) q[2];
rz(-2.9265192) q[3];
sx q[3];
rz(-2.8408065) q[3];
sx q[3];
rz(0.71896368) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9216777) q[0];
sx q[0];
rz(-2.324993) q[0];
sx q[0];
rz(-2.8391229) q[0];
rz(1.1014145) q[1];
sx q[1];
rz(-1.9093724) q[1];
sx q[1];
rz(-2.2473635) q[1];
rz(3.0240622) q[2];
sx q[2];
rz(-1.0385658) q[2];
sx q[2];
rz(-0.081296878) q[2];
rz(0.054417944) q[3];
sx q[3];
rz(-0.46811737) q[3];
sx q[3];
rz(0.56437592) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
