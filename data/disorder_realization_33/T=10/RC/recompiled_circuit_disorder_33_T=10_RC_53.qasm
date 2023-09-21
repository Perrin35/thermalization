OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5053951) q[0];
sx q[0];
rz(-2.8656821) q[0];
sx q[0];
rz(1.8338058) q[0];
rz(1.1360599) q[1];
sx q[1];
rz(-0.93568957) q[1];
sx q[1];
rz(-1.5712665) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9309064) q[0];
sx q[0];
rz(-1.2331729) q[0];
sx q[0];
rz(2.7772285) q[0];
rz(-pi) q[1];
rz(0.76531305) q[2];
sx q[2];
rz(-1.0368477) q[2];
sx q[2];
rz(3.090976) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(3.1085514) q[1];
sx q[1];
rz(-1.812495) q[1];
sx q[1];
rz(-2.7944195) q[1];
rz(-pi) q[2];
rz(1.2075495) q[3];
sx q[3];
rz(-1.7540635) q[3];
sx q[3];
rz(-2.7024384) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.87542614) q[2];
sx q[2];
rz(-0.29310075) q[2];
sx q[2];
rz(-1.1323294) q[2];
rz(-1.6752361) q[3];
sx q[3];
rz(-1.8050067) q[3];
sx q[3];
rz(1.0124538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19673008) q[0];
sx q[0];
rz(-0.20962993) q[0];
sx q[0];
rz(2.9557513) q[0];
rz(-0.56022412) q[1];
sx q[1];
rz(-1.2954243) q[1];
sx q[1];
rz(-2.9247608) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29138716) q[0];
sx q[0];
rz(-2.4225525) q[0];
sx q[0];
rz(-1.1262116) q[0];
x q[1];
rz(-2.6685643) q[2];
sx q[2];
rz(-1.066726) q[2];
sx q[2];
rz(2.8216528) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.195897) q[1];
sx q[1];
rz(-2.2356114) q[1];
sx q[1];
rz(2.871454) q[1];
rz(-2.5012245) q[3];
sx q[3];
rz(-2.159517) q[3];
sx q[3];
rz(2.0224188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.310114) q[2];
sx q[2];
rz(-2.3159413) q[2];
sx q[2];
rz(-1.8537834) q[2];
rz(-0.76256049) q[3];
sx q[3];
rz(-1.1688787) q[3];
sx q[3];
rz(-0.30502239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4644311) q[0];
sx q[0];
rz(-2.7966249) q[0];
sx q[0];
rz(-0.60423869) q[0];
rz(1.3263946) q[1];
sx q[1];
rz(-1.3605958) q[1];
sx q[1];
rz(0.93260971) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7639097) q[0];
sx q[0];
rz(-1.253486) q[0];
sx q[0];
rz(0.056563932) q[0];
x q[1];
rz(2.9595397) q[2];
sx q[2];
rz(-1.3472054) q[2];
sx q[2];
rz(1.6857266) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0886791) q[1];
sx q[1];
rz(-2.1267849) q[1];
sx q[1];
rz(-1.4312137) q[1];
rz(-1.014939) q[3];
sx q[3];
rz(-1.1851289) q[3];
sx q[3];
rz(2.8876497) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.147826) q[2];
sx q[2];
rz(-2.0596762) q[2];
sx q[2];
rz(-1.0926584) q[2];
rz(2.5993733) q[3];
sx q[3];
rz(-1.0850302) q[3];
sx q[3];
rz(2.1742163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-1.3595235) q[0];
sx q[0];
rz(-0.096465915) q[0];
sx q[0];
rz(-2.6413667) q[0];
rz(-2.3362828) q[1];
sx q[1];
rz(-1.1601245) q[1];
sx q[1];
rz(-1.4979699) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.362975) q[0];
sx q[0];
rz(-2.551429) q[0];
sx q[0];
rz(2.1182563) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.285032) q[2];
sx q[2];
rz(-0.19837241) q[2];
sx q[2];
rz(0.49516585) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5768891) q[1];
sx q[1];
rz(-0.35208811) q[1];
sx q[1];
rz(0.7301773) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.3718932) q[3];
sx q[3];
rz(-0.61004988) q[3];
sx q[3];
rz(-1.9609914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3952289) q[2];
sx q[2];
rz(-0.56240288) q[2];
sx q[2];
rz(0.70181075) q[2];
rz(2.3102405) q[3];
sx q[3];
rz(-0.96389198) q[3];
sx q[3];
rz(2.519616) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9005301) q[0];
sx q[0];
rz(-0.59589544) q[0];
sx q[0];
rz(-2.3262614) q[0];
rz(-1.5218081) q[1];
sx q[1];
rz(-2.3074469) q[1];
sx q[1];
rz(-2.0933847) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5883023) q[0];
sx q[0];
rz(-0.36670812) q[0];
sx q[0];
rz(0.81272965) q[0];
rz(-pi) q[1];
rz(-1.5585209) q[2];
sx q[2];
rz(-2.2222812) q[2];
sx q[2];
rz(-2.2001681) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.293562) q[1];
sx q[1];
rz(-1.3333496) q[1];
sx q[1];
rz(-2.8192987) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7311677) q[3];
sx q[3];
rz(-0.99567185) q[3];
sx q[3];
rz(0.23469532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.52577019) q[2];
sx q[2];
rz(-0.56695357) q[2];
sx q[2];
rz(-1.099951) q[2];
rz(0.82529092) q[3];
sx q[3];
rz(-1.0422948) q[3];
sx q[3];
rz(0.88551372) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0681756) q[0];
sx q[0];
rz(-2.5475579) q[0];
sx q[0];
rz(-0.90240479) q[0];
rz(1.0166608) q[1];
sx q[1];
rz(-2.0817751) q[1];
sx q[1];
rz(0.12983233) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31045612) q[0];
sx q[0];
rz(-1.9400915) q[0];
sx q[0];
rz(-0.62311689) q[0];
rz(0.36826276) q[2];
sx q[2];
rz(-1.0266745) q[2];
sx q[2];
rz(2.1899109) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8772014) q[1];
sx q[1];
rz(-1.0069205) q[1];
sx q[1];
rz(1.1370204) q[1];
x q[2];
rz(-2.8270709) q[3];
sx q[3];
rz(-0.57146996) q[3];
sx q[3];
rz(2.275327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8292024) q[2];
sx q[2];
rz(-0.94909334) q[2];
sx q[2];
rz(0.20425805) q[2];
rz(-1.2060818) q[3];
sx q[3];
rz(-1.6198502) q[3];
sx q[3];
rz(2.9061785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7234574) q[0];
sx q[0];
rz(-1.3292987) q[0];
sx q[0];
rz(1.4468505) q[0];
rz(1.8824668) q[1];
sx q[1];
rz(-0.99021688) q[1];
sx q[1];
rz(-0.68626219) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9335564) q[0];
sx q[0];
rz(-2.4718923) q[0];
sx q[0];
rz(0.74525381) q[0];
x q[1];
rz(-0.94481988) q[2];
sx q[2];
rz(-0.84569028) q[2];
sx q[2];
rz(-3.0996029) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.4916617) q[1];
sx q[1];
rz(-1.4964536) q[1];
sx q[1];
rz(-1.7564303) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8920184) q[3];
sx q[3];
rz(-1.1676844) q[3];
sx q[3];
rz(1.5461127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.69616047) q[2];
sx q[2];
rz(-1.7636718) q[2];
sx q[2];
rz(-0.0017722842) q[2];
rz(-2.5799675) q[3];
sx q[3];
rz(-2.2300945) q[3];
sx q[3];
rz(1.5047489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5381662) q[0];
sx q[0];
rz(-2.4551233) q[0];
sx q[0];
rz(1.4461393) q[0];
rz(-2.360545) q[1];
sx q[1];
rz(-1.3054409) q[1];
sx q[1];
rz(-1.6400281) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7088889) q[0];
sx q[0];
rz(-2.6484657) q[0];
sx q[0];
rz(1.6119484) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.976622) q[2];
sx q[2];
rz(-0.067194447) q[2];
sx q[2];
rz(-2.7213328) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.20128076) q[1];
sx q[1];
rz(-0.33946013) q[1];
sx q[1];
rz(0.27075726) q[1];
x q[2];
rz(2.4092259) q[3];
sx q[3];
rz(-2.3673956) q[3];
sx q[3];
rz(0.68294169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.35187307) q[2];
sx q[2];
rz(-1.3616273) q[2];
sx q[2];
rz(-1.3191351) q[2];
rz(1.2119279) q[3];
sx q[3];
rz(-1.8550248) q[3];
sx q[3];
rz(2.8222728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33655745) q[0];
sx q[0];
rz(-0.55258495) q[0];
sx q[0];
rz(1.9375027) q[0];
rz(-0.38326344) q[1];
sx q[1];
rz(-2.6158694) q[1];
sx q[1];
rz(0.35167545) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4708913) q[0];
sx q[0];
rz(-2.0235217) q[0];
sx q[0];
rz(-0.41505138) q[0];
rz(-0.36231626) q[2];
sx q[2];
rz(-3*pi/13) q[2];
sx q[2];
rz(-2.97646) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5896776) q[1];
sx q[1];
rz(-0.81112408) q[1];
sx q[1];
rz(3.0685436) q[1];
rz(-pi) q[2];
rz(0.10961253) q[3];
sx q[3];
rz(-1.622756) q[3];
sx q[3];
rz(-2.6250641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7982771) q[2];
sx q[2];
rz(-2.0337992) q[2];
sx q[2];
rz(1.2822255) q[2];
rz(1.4964237) q[3];
sx q[3];
rz(-1.5346425) q[3];
sx q[3];
rz(-1.055868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6431817) q[0];
sx q[0];
rz(-1.2675985) q[0];
sx q[0];
rz(2.9472651) q[0];
rz(-2.1037897) q[1];
sx q[1];
rz(-0.56832814) q[1];
sx q[1];
rz(2.1077572) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72458306) q[0];
sx q[0];
rz(-1.4405182) q[0];
sx q[0];
rz(0.91086046) q[0];
rz(-pi) q[1];
x q[1];
rz(0.50844426) q[2];
sx q[2];
rz(-0.93548453) q[2];
sx q[2];
rz(2.9490162) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.6453213) q[1];
sx q[1];
rz(-1.4693345) q[1];
sx q[1];
rz(-0.99781499) q[1];
rz(-pi) q[2];
rz(-0.71318993) q[3];
sx q[3];
rz(-1.4048368) q[3];
sx q[3];
rz(-0.99456577) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0795435) q[2];
sx q[2];
rz(-2.1958308) q[2];
sx q[2];
rz(-2.5058084) q[2];
rz(-2.87129) q[3];
sx q[3];
rz(-0.79939866) q[3];
sx q[3];
rz(-1.5283782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4476267) q[0];
sx q[0];
rz(-1.8287369) q[0];
sx q[0];
rz(1.0736314) q[0];
rz(1.4355961) q[1];
sx q[1];
rz(-1.5789079) q[1];
sx q[1];
rz(0.78067738) q[1];
rz(-1.6171261) q[2];
sx q[2];
rz(-0.6033069) q[2];
sx q[2];
rz(2.6848007) q[2];
rz(-1.2407606) q[3];
sx q[3];
rz(-1.5041372) q[3];
sx q[3];
rz(2.1048673) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];