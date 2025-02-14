OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.7527591) q[0];
sx q[0];
rz(-1.4491117) q[0];
sx q[0];
rz(1.4504855) q[0];
rz(-2.1679572) q[1];
sx q[1];
rz(-1.4373625) q[1];
sx q[1];
rz(0.91926423) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0231173) q[0];
sx q[0];
rz(-3.1065515) q[0];
sx q[0];
rz(-2.6316597) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4953142) q[2];
sx q[2];
rz(-2.9913104) q[2];
sx q[2];
rz(0.26963797) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.011292) q[1];
sx q[1];
rz(-0.71349547) q[1];
sx q[1];
rz(-1.43047) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1000865) q[3];
sx q[3];
rz(-1.9506427) q[3];
sx q[3];
rz(-2.0947411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.3322525) q[2];
sx q[2];
rz(-1.6925749) q[2];
sx q[2];
rz(-1.7285041) q[2];
rz(0.20279065) q[3];
sx q[3];
rz(-1.7594124) q[3];
sx q[3];
rz(0.10281674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24421144) q[0];
sx q[0];
rz(-1.0845217) q[0];
sx q[0];
rz(2.4308423) q[0];
rz(-2.3731025) q[1];
sx q[1];
rz(-2.0748506) q[1];
sx q[1];
rz(2.1220727) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1508038) q[0];
sx q[0];
rz(-1.5312563) q[0];
sx q[0];
rz(0.44724748) q[0];
rz(-pi) q[1];
rz(1.1395733) q[2];
sx q[2];
rz(-1.6399472) q[2];
sx q[2];
rz(-0.56996843) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.72538917) q[1];
sx q[1];
rz(-2.7692911) q[1];
sx q[1];
rz(1.6561693) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7048265) q[3];
sx q[3];
rz(-2.5833327) q[3];
sx q[3];
rz(-0.2414862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3780313) q[2];
sx q[2];
rz(-1.8969994) q[2];
sx q[2];
rz(1.9192609) q[2];
rz(-1.9289121) q[3];
sx q[3];
rz(-1.0168394) q[3];
sx q[3];
rz(-0.71162629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0843622) q[0];
sx q[0];
rz(-0.98452345) q[0];
sx q[0];
rz(1.9236176) q[0];
rz(-0.51586622) q[1];
sx q[1];
rz(-2.5936544) q[1];
sx q[1];
rz(-0.9224433) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3271752) q[0];
sx q[0];
rz(-1.6285768) q[0];
sx q[0];
rz(-1.5301276) q[0];
x q[1];
rz(2.5848021) q[2];
sx q[2];
rz(-0.80711952) q[2];
sx q[2];
rz(-0.082376235) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6906185) q[1];
sx q[1];
rz(-0.96267525) q[1];
sx q[1];
rz(-2.4324424) q[1];
x q[2];
rz(1.4018784) q[3];
sx q[3];
rz(-1.2460105) q[3];
sx q[3];
rz(-1.6095773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.018365232) q[2];
sx q[2];
rz(-1.0845228) q[2];
sx q[2];
rz(1.8017192) q[2];
rz(-0.54723048) q[3];
sx q[3];
rz(-2.7180143) q[3];
sx q[3];
rz(2.4303998) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0860586) q[0];
sx q[0];
rz(-2.9965897) q[0];
sx q[0];
rz(0.15922971) q[0];
rz(0.010146443) q[1];
sx q[1];
rz(-1.0228913) q[1];
sx q[1];
rz(-0.25746447) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2672294) q[0];
sx q[0];
rz(-0.82634514) q[0];
sx q[0];
rz(2.6494725) q[0];
rz(-pi) q[1];
rz(3.0765216) q[2];
sx q[2];
rz(-0.79539585) q[2];
sx q[2];
rz(-2.3415023) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.5415693) q[1];
sx q[1];
rz(-0.63508717) q[1];
sx q[1];
rz(2.9734008) q[1];
rz(-pi) q[2];
rz(0.80735029) q[3];
sx q[3];
rz(-1.0991384) q[3];
sx q[3];
rz(-2.014267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3782392) q[2];
sx q[2];
rz(-1.8469609) q[2];
sx q[2];
rz(0.47284687) q[2];
rz(-0.7044479) q[3];
sx q[3];
rz(-1.7770146) q[3];
sx q[3];
rz(2.4956467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6351629) q[0];
sx q[0];
rz(-1.1744873) q[0];
sx q[0];
rz(-0.065486431) q[0];
rz(2.7323515) q[1];
sx q[1];
rz(-2.0114653) q[1];
sx q[1];
rz(1.7154891) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5201841) q[0];
sx q[0];
rz(-1.4778293) q[0];
sx q[0];
rz(-0.21958242) q[0];
x q[1];
rz(-0.3887811) q[2];
sx q[2];
rz(-1.6531303) q[2];
sx q[2];
rz(-1.8052693) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.1115426) q[1];
sx q[1];
rz(-0.36007127) q[1];
sx q[1];
rz(0.67782016) q[1];
rz(-pi) q[2];
rz(0.17284837) q[3];
sx q[3];
rz(-1.2842872) q[3];
sx q[3];
rz(0.89927538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.9466729) q[2];
sx q[2];
rz(-1.0871004) q[2];
sx q[2];
rz(1.8348414) q[2];
rz(-1.170018) q[3];
sx q[3];
rz(-0.40478671) q[3];
sx q[3];
rz(-0.10725966) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4829247) q[0];
sx q[0];
rz(-1.985745) q[0];
sx q[0];
rz(2.7401127) q[0];
rz(-0.67277706) q[1];
sx q[1];
rz(-2.3428226) q[1];
sx q[1];
rz(-2.2875517) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6333165) q[0];
sx q[0];
rz(-0.22249732) q[0];
sx q[0];
rz(-1.359324) q[0];
rz(1.4208904) q[2];
sx q[2];
rz(-1.3873563) q[2];
sx q[2];
rz(0.73018398) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.8484162) q[1];
sx q[1];
rz(-0.59976116) q[1];
sx q[1];
rz(-0.17677115) q[1];
x q[2];
rz(2.5307114) q[3];
sx q[3];
rz(-1.1958836) q[3];
sx q[3];
rz(1.5883816) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.65471571) q[2];
sx q[2];
rz(-0.87149039) q[2];
sx q[2];
rz(1.1775449) q[2];
rz(2.5901637) q[3];
sx q[3];
rz(-1.5827554) q[3];
sx q[3];
rz(2.3059755) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7198782) q[0];
sx q[0];
rz(-2.8604909) q[0];
sx q[0];
rz(1.9792492) q[0];
rz(1.1516736) q[1];
sx q[1];
rz(-1.3538227) q[1];
sx q[1];
rz(0.91845671) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3716482) q[0];
sx q[0];
rz(-1.838883) q[0];
sx q[0];
rz(1.2275342) q[0];
rz(-pi) q[1];
rz(-0.26084857) q[2];
sx q[2];
rz(-2.2493304) q[2];
sx q[2];
rz(-1.4229753) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.7734997) q[1];
sx q[1];
rz(-1.5775023) q[1];
sx q[1];
rz(-3.1342602) q[1];
rz(-pi) q[2];
rz(0.17042589) q[3];
sx q[3];
rz(-1.9099351) q[3];
sx q[3];
rz(-0.91811524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.1071757) q[2];
sx q[2];
rz(-1.3849266) q[2];
sx q[2];
rz(0.50298634) q[2];
rz(-2.3668187) q[3];
sx q[3];
rz(-0.46190244) q[3];
sx q[3];
rz(1.3431965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4485432) q[0];
sx q[0];
rz(-0.22668426) q[0];
sx q[0];
rz(-1.5013129) q[0];
rz(2.8003108) q[1];
sx q[1];
rz(-1.4309859) q[1];
sx q[1];
rz(-1.1192082) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5401181) q[0];
sx q[0];
rz(-1.0154503) q[0];
sx q[0];
rz(-2.0106273) q[0];
x q[1];
rz(-1.0985435) q[2];
sx q[2];
rz(-0.53016463) q[2];
sx q[2];
rz(1.7165754) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.71257617) q[1];
sx q[1];
rz(-1.6916995) q[1];
sx q[1];
rz(2.4790133) q[1];
x q[2];
rz(-1.4358892) q[3];
sx q[3];
rz(-1.9622318) q[3];
sx q[3];
rz(-1.5086482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1559653) q[2];
sx q[2];
rz(-2.1817744) q[2];
sx q[2];
rz(-2.7723374) q[2];
rz(0.30820942) q[3];
sx q[3];
rz(-0.9674558) q[3];
sx q[3];
rz(-1.6109899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8592598) q[0];
sx q[0];
rz(-1.8349324) q[0];
sx q[0];
rz(-2.7985213) q[0];
rz(-1.9725017) q[1];
sx q[1];
rz(-1.7770551) q[1];
sx q[1];
rz(1.4287359) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0853303) q[0];
sx q[0];
rz(-1.3199727) q[0];
sx q[0];
rz(-2.73231) q[0];
rz(-pi) q[1];
rz(0.83447225) q[2];
sx q[2];
rz(-1.021046) q[2];
sx q[2];
rz(-2.2817734) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.0020121) q[1];
sx q[1];
rz(-1.5095469) q[1];
sx q[1];
rz(-1.9467926) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2637706) q[3];
sx q[3];
rz(-1.1980499) q[3];
sx q[3];
rz(-0.28723785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.94065654) q[2];
sx q[2];
rz(-0.20874615) q[2];
sx q[2];
rz(-1.7500056) q[2];
rz(-0.40677795) q[3];
sx q[3];
rz(-1.4429561) q[3];
sx q[3];
rz(0.87882915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.10130356) q[0];
sx q[0];
rz(-0.60281301) q[0];
sx q[0];
rz(-1.783675) q[0];
rz(-1.2376002) q[1];
sx q[1];
rz(-2.1289181) q[1];
sx q[1];
rz(-0.79992574) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023066085) q[0];
sx q[0];
rz(-0.79341054) q[0];
sx q[0];
rz(1.2454459) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2679891) q[2];
sx q[2];
rz(-2.1005512) q[2];
sx q[2];
rz(2.985266) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8854196) q[1];
sx q[1];
rz(-0.90682632) q[1];
sx q[1];
rz(-1.3732984) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0217729) q[3];
sx q[3];
rz(-0.52996892) q[3];
sx q[3];
rz(-1.3254904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7633729) q[2];
sx q[2];
rz(-1.7174481) q[2];
sx q[2];
rz(3.0685032) q[2];
rz(2.8857005) q[3];
sx q[3];
rz(-2.3292694) q[3];
sx q[3];
rz(-1.6528486) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.273461) q[0];
sx q[0];
rz(-1.68597) q[0];
sx q[0];
rz(-1.2706533) q[0];
rz(1.801626) q[1];
sx q[1];
rz(-1.7908962) q[1];
sx q[1];
rz(-2.4790196) q[1];
rz(0.70941464) q[2];
sx q[2];
rz(-1.5723036) q[2];
sx q[2];
rz(0.2607762) q[2];
rz(0.70684915) q[3];
sx q[3];
rz(-2.3656188) q[3];
sx q[3];
rz(-1.594365) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
