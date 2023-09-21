OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.049790073) q[0];
sx q[0];
rz(3.0135305) q[0];
sx q[0];
rz(11.749) q[0];
rz(-2.1583537) q[1];
sx q[1];
rz(-2.6020738) q[1];
sx q[1];
rz(-1.9411545) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.565633) q[0];
sx q[0];
rz(-2.6607249) q[0];
sx q[0];
rz(2.5984882) q[0];
x q[1];
rz(3.009216) q[2];
sx q[2];
rz(-1.0356324) q[2];
sx q[2];
rz(-1.7420499) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5719205) q[1];
sx q[1];
rz(-1.8351646) q[1];
sx q[1];
rz(-0.7260679) q[1];
rz(-pi) q[2];
rz(-2.3463763) q[3];
sx q[3];
rz(-1.8024369) q[3];
sx q[3];
rz(2.8455646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.0212705) q[2];
sx q[2];
rz(-2.5734084) q[2];
sx q[2];
rz(1.5585287) q[2];
rz(2.1448686) q[3];
sx q[3];
rz(-0.45209) q[3];
sx q[3];
rz(2.7157917) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5834171) q[0];
sx q[0];
rz(-2.3893864) q[0];
sx q[0];
rz(0.054071991) q[0];
rz(1.1955098) q[1];
sx q[1];
rz(-2.1046488) q[1];
sx q[1];
rz(-0.53584677) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5862522) q[0];
sx q[0];
rz(-2.147701) q[0];
sx q[0];
rz(2.9931086) q[0];
rz(-pi) q[1];
rz(0.29157721) q[2];
sx q[2];
rz(-0.73086408) q[2];
sx q[2];
rz(2.4183194) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(3.0332196) q[1];
sx q[1];
rz(-2.6544826) q[1];
sx q[1];
rz(-2.5750722) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4278533) q[3];
sx q[3];
rz(-1.5385475) q[3];
sx q[3];
rz(1.3168207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1198931) q[2];
sx q[2];
rz(-2.0844441) q[2];
sx q[2];
rz(2.1832441) q[2];
rz(0.066453233) q[3];
sx q[3];
rz(-1.5840014) q[3];
sx q[3];
rz(0.44979969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7217343) q[0];
sx q[0];
rz(-1.2092713) q[0];
sx q[0];
rz(0.15047519) q[0];
rz(2.6843605) q[1];
sx q[1];
rz(-0.91232863) q[1];
sx q[1];
rz(-0.025807468) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8919825) q[0];
sx q[0];
rz(-1.5907856) q[0];
sx q[0];
rz(1.32094) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8823207) q[2];
sx q[2];
rz(-1.6561964) q[2];
sx q[2];
rz(-0.66750079) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.64297134) q[1];
sx q[1];
rz(-1.9830623) q[1];
sx q[1];
rz(-0.69795124) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2529536) q[3];
sx q[3];
rz(-2.7405973) q[3];
sx q[3];
rz(0.81438118) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0187443) q[2];
sx q[2];
rz(-0.3652502) q[2];
sx q[2];
rz(-2.5562111) q[2];
rz(-2.9600926) q[3];
sx q[3];
rz(-1.3320965) q[3];
sx q[3];
rz(-1.5649149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.240775) q[0];
sx q[0];
rz(-2.5188991) q[0];
sx q[0];
rz(-2.9649819) q[0];
rz(2.2606842) q[1];
sx q[1];
rz(-1.0662339) q[1];
sx q[1];
rz(-0.53612971) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49318424) q[0];
sx q[0];
rz(-0.83183653) q[0];
sx q[0];
rz(2.1302845) q[0];
rz(-pi) q[1];
x q[1];
rz(1.5136112) q[2];
sx q[2];
rz(-1.8071037) q[2];
sx q[2];
rz(1.1255217) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0587479) q[1];
sx q[1];
rz(-1.1943598) q[1];
sx q[1];
rz(2.46545) q[1];
rz(2.4651277) q[3];
sx q[3];
rz(-1.8026661) q[3];
sx q[3];
rz(-1.7959309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.6716016) q[2];
sx q[2];
rz(-1.4208379) q[2];
sx q[2];
rz(1.0446576) q[2];
rz(2.4345543) q[3];
sx q[3];
rz(-2.1576594) q[3];
sx q[3];
rz(-0.35693359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9064643) q[0];
sx q[0];
rz(-1.3236073) q[0];
sx q[0];
rz(-2.0902324) q[0];
rz(1.6479187) q[1];
sx q[1];
rz(-0.58902478) q[1];
sx q[1];
rz(-3.0984745) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.608425) q[0];
sx q[0];
rz(-2.5072917) q[0];
sx q[0];
rz(-1.8932635) q[0];
rz(-pi) q[1];
rz(2.0897323) q[2];
sx q[2];
rz(-1.3216002) q[2];
sx q[2];
rz(3.0254226) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.67152126) q[1];
sx q[1];
rz(-0.38362353) q[1];
sx q[1];
rz(2.3684711) q[1];
rz(-pi) q[2];
rz(2.352036) q[3];
sx q[3];
rz(-1.6766747) q[3];
sx q[3];
rz(-1.7136128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.0126426) q[2];
sx q[2];
rz(-0.49762112) q[2];
sx q[2];
rz(2.7894003) q[2];
rz(0.59018618) q[3];
sx q[3];
rz(-2.6679109) q[3];
sx q[3];
rz(-0.56110704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6234289) q[0];
sx q[0];
rz(-1.2151027) q[0];
sx q[0];
rz(1.9859001) q[0];
rz(-0.75025264) q[1];
sx q[1];
rz(-0.9393839) q[1];
sx q[1];
rz(2.0828784) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4689894) q[0];
sx q[0];
rz(-0.4840584) q[0];
sx q[0];
rz(2.6812535) q[0];
rz(-pi) q[1];
rz(2.6201453) q[2];
sx q[2];
rz(-1.2505184) q[2];
sx q[2];
rz(-0.59200586) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.9204273) q[1];
sx q[1];
rz(-1.7943008) q[1];
sx q[1];
rz(-1.8799053) q[1];
rz(-pi) q[2];
rz(-2.433957) q[3];
sx q[3];
rz(-1.71873) q[3];
sx q[3];
rz(-1.9572452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.47026149) q[2];
sx q[2];
rz(-1.417421) q[2];
sx q[2];
rz(1.9227825) q[2];
rz(1.9865215) q[3];
sx q[3];
rz(-0.23854908) q[3];
sx q[3];
rz(-1.4412122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0999775) q[0];
sx q[0];
rz(-1.3168553) q[0];
sx q[0];
rz(-1.6301427) q[0];
rz(1.7639683) q[1];
sx q[1];
rz(-2.8306077) q[1];
sx q[1];
rz(0.84164936) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8182939) q[0];
sx q[0];
rz(-0.57384402) q[0];
sx q[0];
rz(0.42563514) q[0];
x q[1];
rz(-1.7633348) q[2];
sx q[2];
rz(-0.95270573) q[2];
sx q[2];
rz(2.3877909) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6517666) q[1];
sx q[1];
rz(-3.0511599) q[1];
sx q[1];
rz(2.138278) q[1];
x q[2];
rz(1.8089201) q[3];
sx q[3];
rz(-0.75546414) q[3];
sx q[3];
rz(1.8734224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(3.1402309) q[2];
sx q[2];
rz(-1.7732239) q[2];
sx q[2];
rz(0.049953071) q[2];
rz(0.66155457) q[3];
sx q[3];
rz(-0.52246061) q[3];
sx q[3];
rz(2.9522827) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.89896232) q[0];
sx q[0];
rz(-1.0268509) q[0];
sx q[0];
rz(1.4021953) q[0];
rz(0.095480355) q[1];
sx q[1];
rz(-1.1663368) q[1];
sx q[1];
rz(2.7239674) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4293489) q[0];
sx q[0];
rz(-2.0040383) q[0];
sx q[0];
rz(-2.594125) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8629486) q[2];
sx q[2];
rz(-0.38049618) q[2];
sx q[2];
rz(-0.7064864) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.046125267) q[1];
sx q[1];
rz(-1.0091262) q[1];
sx q[1];
rz(2.5820877) q[1];
rz(1.9440218) q[3];
sx q[3];
rz(-1.2959359) q[3];
sx q[3];
rz(-1.1653792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.3461356) q[2];
sx q[2];
rz(-1.4435578) q[2];
sx q[2];
rz(-1.0428838) q[2];
rz(2.4677094) q[3];
sx q[3];
rz(-1.6582812) q[3];
sx q[3];
rz(2.2414482) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3089356) q[0];
sx q[0];
rz(-1.7333663) q[0];
sx q[0];
rz(-2.5119264) q[0];
rz(0.57485238) q[1];
sx q[1];
rz(-1.8300627) q[1];
sx q[1];
rz(-0.94690698) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1614721) q[0];
sx q[0];
rz(-1.8959909) q[0];
sx q[0];
rz(1.8878216) q[0];
rz(-pi) q[1];
rz(1.5309179) q[2];
sx q[2];
rz(-1.4738184) q[2];
sx q[2];
rz(3.1181042) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.8967594) q[1];
sx q[1];
rz(-1.6549126) q[1];
sx q[1];
rz(0.33478488) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0128653) q[3];
sx q[3];
rz(-1.3978492) q[3];
sx q[3];
rz(0.63463075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5809014) q[2];
sx q[2];
rz(-1.1297444) q[2];
sx q[2];
rz(-1.2488731) q[2];
rz(2.4272264) q[3];
sx q[3];
rz(-1.8656105) q[3];
sx q[3];
rz(-2.8760288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0666075) q[0];
sx q[0];
rz(-2.5191436) q[0];
sx q[0];
rz(-2.4865436) q[0];
rz(-0.89637268) q[1];
sx q[1];
rz(-2.2348149) q[1];
sx q[1];
rz(2.4972829) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5376741) q[0];
sx q[0];
rz(-0.18225056) q[0];
sx q[0];
rz(-0.3084348) q[0];
x q[1];
rz(2.9115453) q[2];
sx q[2];
rz(-2.2736079) q[2];
sx q[2];
rz(2.7237797) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.58562467) q[1];
sx q[1];
rz(-1.2871736) q[1];
sx q[1];
rz(0.3507627) q[1];
rz(1.0993016) q[3];
sx q[3];
rz(-1.9715371) q[3];
sx q[3];
rz(2.5221962) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.3506938) q[2];
sx q[2];
rz(-1.472298) q[2];
sx q[2];
rz(-2.541686) q[2];
rz(-2.24263) q[3];
sx q[3];
rz(-2.9581684) q[3];
sx q[3];
rz(1.3658587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8469289) q[0];
sx q[0];
rz(-1.8142721) q[0];
sx q[0];
rz(-0.66692764) q[0];
rz(-0.22944336) q[1];
sx q[1];
rz(-2.2506917) q[1];
sx q[1];
rz(-3.0058203) q[1];
rz(-1.7462126) q[2];
sx q[2];
rz(-0.63175628) q[2];
sx q[2];
rz(0.14887688) q[2];
rz(2.7502144) q[3];
sx q[3];
rz(-2.2073675) q[3];
sx q[3];
rz(-1.7195306) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];